#include <iostream>
#include <map>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "external/cxxopts.hpp"

#include "include/config.h"
#include "include/util.h"
#include "include/event/ore/raw_trigger_event.h"
#include "include/smelt/sifting.h"

using namespace forgerib;

int ReadXiaTimes(
	const char *path,
	std::vector<double> &times,
	std::vector<double> &external_times,
	std::vector<long long> &entries,
	bool report = false
) {
	TFile ipf(path, "read");
	TTree *ipt = static_cast<TTree*>(ipf.Get("tree"));
	if (!ipt) {
		std::cerr << "Error: Get tree from " << path << " failed.\n";
		return -1;
	}
	forgerib::RawTriggerEvent trigger_event;
	forgerib::SetupInput(ipt, trigger_event);

	times.clear();
	external_times.clear();
	entries.clear();
	if (report) {
		printf("Reading XIA times   0%%");
		fflush(stdout);
	}
	long long total = ipt->GetEntries();
	long long last_percentage = 0;
	for (long long entry = 0; entry < total; ++entry) {
		if (report && entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		if (trigger_event.type == 1) {
			times.push_back(trigger_event.time);              // ns
			external_times.push_back(trigger_event.external_time);
			entries.push_back(entry);
		}
	}
	if (report) {
		printf("\b\b\b\b100%%\n");
	}
	return 0;
}


int ReadVmeTimes(
	const char *path,
	std::vector<double> &times,          // 输出：ns
	std::vector<long long> &entries,
	double time_threshold_ns,            // 输入：ns
	long long &total_valid,              // 输出：有效且通过阈值的总数（用于打印）
	bool report = false
) {
	TFile ipf(path, "read");
	if (!ipf.IsOpen()) {
		std::cerr << "Error: Cannot open " << path << ".\n";
		return -1;
	}
	TTree *ipt = static_cast<TTree*>(ipf.Get("tree"));
	if (!ipt) {
		std::cerr << "Error: get tree from " << path << " failed.\n";
		return -1;
	}

	bool valid = false;
	long long time = 0;                  // us
	ipt->SetBranchAddress("valid", &valid);
	ipt->SetBranchAddress("time", &time);

	times.clear();
	entries.clear();
	total_valid = 0;
	long long total = ipt->GetEntries();
	long long last_percentage = 0;
	if (report) {
		printf("Reading VME trigger times   0%%");
		fflush(stdout);
	}
	for (long long entry = 0; entry < total; ++entry) {
		if (report && entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		if (!valid) continue;
		++total_valid;

		// time 单位 us -> ns
		double time_ns = static_cast<double>(time) * 1000.0;
		if (time_ns <= time_threshold_ns) continue;

		times.push_back(time_ns);
		entries.push_back(entry);
	}
	if (report) printf("\b\b\b\b100%%\n");
	ipf.Close();
	return 0;
}


int main(int argc, char **argv) {
	cxxopts::Options options("sift", "sift trigger alignment");
	options.add_options()
		("h,help", "Print usage")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"path"
		)
		(
			"x,xia_run",
			"XIA run number, required.",
			cxxopts::value<int>(),
			"run"
		)
		(
			"v,vme_run",
			"VME run number, required.",
			cxxopts::value<int>(),
			"run"
		)
		(
			"e,external",
			"Use XIA's external timestamp.",
			cxxopts::value<bool>()->default_value("false")
		)
		(
			"time",
			"VME absolute time threshold in seconds. VME times (us) are converted to ns and only those > threshold are kept. Default 0 s.",
			cxxopts::value<double>()->default_value("0.0"),
			"seconds"
		)
		(
			"dt",
			"Match window = dt * 1000 (ns). Default 15.",
			cxxopts::value<double>()->default_value("15.0"),
			"value"
		);
	auto parse_result = options.parse(argc, argv);

	if (
		parse_result.count("help")
		|| !parse_result.count("xia_run")
		|| !parse_result.count("vme_run")
	) {
		std::cout << options.help() << std::endl;
		return 0;
	}
	int xia_run = parse_result["xia_run"].as<int>();
	int vme_run = parse_result["vme_run"].as<int>();
	bool external = parse_result["external"].as<bool>();
	double time_threshold_s = parse_result["time"].as<double>();
	double dt_value         = parse_result["dt"].as<double>();
	const std::string config_path = parse_result["config"].as<std::string>();

	// 统一到 ns
	double vme_time_threshold_ns = time_threshold_s * 1e9;
	double window                = dt_value * 1000.0;

	AppConfig config;
	if (LoadConfig(config_path, config)) {
		return -1;
	}
	const std::string grit_dir  = JoinPath(config.workspace, config.paths.grit);
	const std::string grain_dir = JoinPath(config.workspace, config.paths.grain);

	// ---- XIA ----
	std::vector<double> xia_times;
	std::vector<double> xia_external_times;
	std::vector<long long> xia_entries;
	TString xia_filename = TString::Format(
		"%s/trigger_%04d.root",
		grit_dir.c_str(),
		xia_run
	);
	if (ReadXiaTimes(
		xia_filename.Data(),
		xia_times,
		xia_external_times,
		xia_entries,
		true
	)) {
		std::cerr << "Error: Read XIA times failed.\n";
		return -2;
	}

	// ---- VME（读取时按阈值过滤）----
	std::vector<double> vme_times;       // ns
	std::vector<long long> vme_entries;
	long long vme_total_valid = 0;
	TString vme_filename = TString::Format(
		"%s/trigger_vme_%04d.root",
		grit_dir.c_str(),
		vme_run
	);
	if (ReadVmeTimes(
		vme_filename.Data(),
		vme_times,
		vme_entries,
		vme_time_threshold_ns,
		vme_total_valid,
		true
	)) {
		std::cerr << "Error: Read VME trigger times failed.\n";
		return -2;
	}

	TString output_file = TString::Format(
		"%s/grain_%04d.root",
		grain_dir.c_str(),
		xia_run
	);

	int result = Sift(
		(external ? xia_external_times : xia_times),
		xia_entries,
		vme_times,
		vme_entries,
		output_file,
		window,
		10,    // lookhead
		1000,
		1000,
		true
	);

	// ---- 结尾统一打印 ----
	std::cout << "\n===== sift summary =====\n";
	std::cout << "XIA run            : " << xia_run << "\n";
	std::cout << "VME run            : " << vme_run << "\n";
	std::cout << "Use external XIA   : " << (external ? "yes" : "no") << "\n";
	std::cout << "XIA entries        : " << xia_times.size() << "\n";
	std::cout << "VME valid entries  : " << vme_total_valid << "\n";
	std::cout << "VME threshold      : " << time_threshold_s << " s ("
	          << vme_time_threshold_ns << " ns)\n";
	std::cout << "VME kept entries   : " << vme_times.size() << "\n";
	std::cout << "Match window       : " << window << " ns (dt="
	          << dt_value << " * 1000)\n";
	std::cout << "Output file        : " << output_file.Data() << "\n";
	if (result) {
		std::cerr << "Error: Sift failed.\n";
		return -1;
	}
	std::cout << "Status             : OK\n";
	return 0;
}