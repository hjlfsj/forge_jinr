#include <algorithm>
#include <array>
#include <cctype>
#include <dirent.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include "external/cxxopts.hpp"
#include "external/toml.hpp"
#include "include/config.h"
#include "include/util.h"

namespace {

using namespace forgerib;

// 从目录列出 *.root，供交互选择（POSIX，无 ROOT 依赖）
std::vector<std::string> ListRootFiles(const std::string &dir) {
	std::vector<std::string> files;
	DIR *d = opendir(dir.c_str());
	if (!d) return files;
	struct dirent *ent = nullptr;
	while ((ent = readdir(d)) != nullptr) {
		std::string s(ent->d_name);
		if (s.size() > 5 && s.compare(s.size() - 5, 5, ".root") == 0) {
			files.push_back(dir + "/" + s);
		}
	}
	closedir(d);
	std::sort(files.begin(), files.end());
	return files;
}

// 终端交互：询问输入文件。返回空串表示取消。
std::string PromptInputFile(const std::string &raw_vme_dir) {
	std::vector<std::string> files = ListRootFiles(raw_vme_dir);
	std::cout << "\n=== Input file selection ===\n";
	if (!files.empty()) {
		std::cout << "Files found in " << raw_vme_dir << ":\n";
		for (size_t i = 0; i < files.size(); ++i) {
			std::cout << "  [" << i << "] " << files[i] << "\n";
		}
	} else {
		std::cout << "No *.root found in " << raw_vme_dir
			<< ". Please input a full path.\n";
	}
	std::cout << "Enter index or path (empty to abort): ";
	std::string line;
	if (!std::getline(std::cin, line)) return {};
	// trim
	auto l = line.find_first_not_of(" \t\r\n");
	auto r = line.find_last_not_of(" \t\r\n");
	if (l == std::string::npos) return {};
	line = line.substr(l, r - l + 1);
	if (line.empty()) return {};

	// 纯数字 -> 索引
	bool all_digit = std::all_of(line.begin(), line.end(),
		[](unsigned char c) { return std::isdigit(c) != 0; });
	if (all_digit) {
		size_t idx = static_cast<size_t>(std::stoul(line));
		if (idx < files.size()) return files[idx];
		std::cerr << "Error: index " << idx << " out of range.\n";
		return {};
	}
	return line;
}

// 读取 Go4 AnalysisxTree 的 NeEvent.scalerXIA，填充 times。
int ReadScalerXia(const std::string &path, std::vector<long long> &times) {
	TFile input_file(path.c_str(), "read");
	if (!input_file.IsOpen()) {
		std::cerr << "Error: Cannot open " << path << ".\n";
		return -1;
	}
	TTree *tree = static_cast<TTree *>(input_file.Get("AnalysisxTree"));
	if (!tree) {
		std::cerr << "Error: Get AnalysisxTree from " << path << " failed.\n";
		return -1;
	}

	TTreeReader reader(tree);
	TTreeReaderValue<ULong64_t> scaler_xia(reader, "NeEvent.scalerXIA");

	const Long64_t total = tree->GetEntries();
	times.clear();
	times.reserve(static_cast<size_t>(total));

	long long last_percentage = 0;
	printf("Reading scalerXIA   0%%");
	fflush(stdout);
	for (Long64_t entry = 0; entry < total; ++entry) {
		if (entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		reader.Next();
		if (reader.GetEntryStatus() != TTreeReader::kEntryValid) {
			std::cerr << "\nError: TTreeReader failed at entry " << entry << ".\n";
			return -1;
		}
		times.push_back(static_cast<long long>(*scaler_xia));
	}
	printf("\b\b\b\b100%%\n");

	input_file.Close();
	return 0;
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options(
		"crush_vme",
		"convert Go4 scalerXIA to trigger_vme_XXXX.root"
	);
	options.add_options()
		("h,help", "Print usage")
		(
			"run",
			"VME run number (used for output file name).",
			cxxopts::value<int>(),
			"run"
		)
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"path"
		)
		// 保留 threshold 选项以保持接口不变；本程序不再使用
		(
			"threshold",
			"(unused) Threshold file path, kept for interface compatibility.",
			cxxopts::value<std::string>()->default_value("vme_threshold.toml"),
			"path"
		)
		;
	options.parse_positional({"run"});
	auto parse_result = options.parse(argc, argv);
	if (parse_result.count("help") || !parse_result.count("run")) {
		std::cout << options.help() << std::endl;
		return 0;
	}

	const int vme_run = parse_result["run"].as<int>();
	const std::string config_path = parse_result["config"].as<std::string>();

	AppConfig config;
	if (LoadConfig(config_path, config)) {
		return -1;
	}
	const std::string output_dir = JoinPath(config.workspace, config.paths.grit);

	// 终端交互：询问输入文件
	const std::string input_path = PromptInputFile(config.paths.raw_vme);
	if (input_path.empty()) {
		std::cerr << "Aborted: no input file selected.\n";
		return -1;
	}
	std::cout << "Input : " << input_path << "\n";

	// 读取 scalerXIA
	std::vector<long long> trigger_times;
	if (ReadScalerXia(input_path, trigger_times)) {
		return -1;
	}
	std::cout << "Entries: " << trigger_times.size() << "\n";

	// 输出
	TString trigger_path = TString::Format(
		"%s/trigger_vme_%04d.root",
		output_dir.c_str(),
		vme_run
	);
	std::cout << "Output: " << trigger_path << "\n";

	TFile trigger_file(trigger_path, "recreate");
	if (!trigger_file.IsOpen()) {
		std::cerr << "Error: Cannot create " << trigger_path << ".\n";
		return -1;
	}
	TTree trigger_tree("tree", "trigger vme");
	bool trigger_valid = false;
	long long trigger_time = 0;
	trigger_tree.Branch("valid", &trigger_valid, "v/O");
	trigger_tree.Branch("time", &trigger_time, "time/L");

	for (long long t : trigger_times) {
		trigger_valid = (t != 0);
		trigger_time = t;
		trigger_tree.Fill();
	}

	trigger_file.cd();
	trigger_tree.Write();
	trigger_file.Close();
	return 0;
}
