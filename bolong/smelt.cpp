#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>

#include "external/cxxopts.hpp"

#include "include/config.h"
#include "include/util.h"
#include "include/event/ingot/trigger_event.h"
#include "include/event/ore/raw_trigger_event.h"
#include "include/event/ingot/csi_event.h"
#include "include/event/ingot/dssd_event.h"
#include "include/event/ingot/tafd_event.h"
#include "include/event/ingot/ppac_event.h"
#include "include/smelt/smelt_beam.h"
#include "include/smelt/smelt_detector.h"
#include "include/smelt/smelt_ppac.h"
#include "include/smelt/smelt_ssd.h"
#include "include/smelt/smelt_t0csi.h"
#include "include/smelt/smelt_t0csi_bloom.h"
#include "include/smelt/smelt_gagg.h"
#include "include/smelt/smelt_t0d1.h"
#include "include/smelt/smelt_t0d2.h"
#include "include/smelt/smelt_t0d3.h"
#include "include/smelt/smelt_t0d4.h"

namespace {

using forgerib::AlignEvent;

std::string TriggerStem(const std::string &trigger_type) {
	return trigger_type == "main" ? "" : trigger_type + "_";
}

TString GritInputPath(
	const std::string &grit_dir,
	const char *detector,
	int run
) {
	return TString::Format("%s/%s_%04d.root", grit_dir.c_str(), detector, run);
}

TString GritVmeInputPath(
	const std::string &grit_dir,
	const char *detector,
	int run
) {
	return TString::Format("%s/%s_vme_%04d.root", grit_dir.c_str(), detector, run);
}

TString IngotOutputPath(
	const std::string &ingot_dir,
	const char *detector,
	const std::string &trigger_type,
	int run
) {
	return TString::Format(
		"%s/%s_%s%04d.root",
		ingot_dir.c_str(),
		detector,
		TriggerStem(trigger_type).c_str(),
		run
	);
}

TString IngotTriggeredPath(
	const std::string &ingot_dir,
	const char *detector,
	const std::string &trigger_type,
	int run
) {
	return TString::Format(
		"%s/%s_%s%04d.root",
		ingot_dir.c_str(),
		detector,
		TriggerStem(trigger_type).c_str(),
		run
	);
}

// ---------------------------------------------------------------------------
// 把刻度分支（cal）补进 smelt 的 ingot 输出。
// forgerib 的 TafdEvent / DssdEvent / CsiEvent 里没有 cal 成员，所以
// SmeltDetector<T> 不会写它们；这里按 vme_entries 的对应关系，从
// <det>_vme_XXXX.root 逐条复制，顺序与输出事例一一对应。
// ---------------------------------------------------------------------------

int CopyTafdCal(
	const char *out_path,
	const char *vme_path,
	const std::vector<long long> &vme_entries
) {
	TFile of(out_path, "update");
	TTree *ot = static_cast<TTree*>(of.Get("tree"));
	TFile vf(vme_path, "read");
	TTree *vt = static_cast<TTree*>(vf.Get("tree"));
	if (!of.IsOpen() || !ot || !vf.IsOpen() || !vt) {
		std::cerr << "Error: open cal source/output failed: "
			<< out_path << " / " << vme_path << "\n";
		return -1;
	}

	double fc[6] = {};
	double bc[6] = {};
	vt->SetBranchAddress("front_cal", fc);
	vt->SetBranchAddress("back_cal", bc);

	TBranch *b_fc = ot->Branch("front_cal", fc, "fc[6]/D");
	TBranch *b_bc = ot->Branch("back_cal", bc, "bc[6]/D");
	if (!b_fc || !b_bc) {
		std::cerr << "Error: create cal branches in " << out_path << "\n";
		return -2;
	}

	const Long64_t n = ot->GetEntries();
	if (static_cast<size_t>(n) != vme_entries.size()) {
		std::cerr << "Warning: " << out_path << " has " << n
			<< " entries, trigger has " << vme_entries.size() << "\n";
	}
	for (Long64_t i = 0; i < n; ++i) {
		ot->GetEntry(i);
		const long long ve = static_cast<size_t>(i) < vme_entries.size()
			? vme_entries[i] : -1;
		if (ve >= 0 && ve < vt->GetEntries()) {
			vt->GetEntry(ve);
		} else {
			for (int k = 0; k < 6; ++k) {
				fc[k] = 0;
				bc[k] = 0;
			}
		}
		b_fc->Fill();
		b_bc->Fill();
	}
	// TTree::Write 写到当前目录，必须先切回输出文件（否则会写进只读的输入文件）
	of.cd();
	ot->Write("", TObject::kOverwrite);
	of.Close();
	return 0;
}

int CopyDssdCal(
	const char *out_path,
	const char *vme_path,
	const std::vector<long long> &vme_entries
) {
	TFile of(out_path, "update");
	TTree *ot = static_cast<TTree*>(of.Get("tree"));
	TFile vf(vme_path, "read");
	TTree *vt = static_cast<TTree*>(vf.Get("tree"));
	if (!of.IsOpen() || !ot || !vf.IsOpen() || !vt) {
		std::cerr << "Error: open cal source/output failed: "
			<< out_path << " / " << vme_path << "\n";
		return -1;
	}

	double fc[32] = {};
	double bc[32] = {};
	vt->SetBranchAddress("front_cal", fc);
	vt->SetBranchAddress("back_cal", bc);

	TBranch *b_fc = ot->Branch("front_cal", fc, "fc[fn]/D");
	TBranch *b_bc = ot->Branch("back_cal", bc, "bc[bn]/D");
	if (!b_fc || !b_bc) {
		std::cerr << "Error: create cal branches in " << out_path << "\n";
		return -2;
	}

	const Long64_t n = ot->GetEntries();
	if (static_cast<size_t>(n) != vme_entries.size()) {
		std::cerr << "Warning: " << out_path << " has " << n
			<< " entries, trigger has " << vme_entries.size() << "\n";
	}
	for (Long64_t i = 0; i < n; ++i) {
		ot->GetEntry(i);
		const long long ve = static_cast<size_t>(i) < vme_entries.size()
			? vme_entries[i] : -1;
		if (ve >= 0 && ve < vt->GetEntries()) {
			vt->GetEntry(ve);
		} else {
			for (int k = 0; k < 32; ++k) {
				fc[k] = 0;
				bc[k] = 0;
			}
		}
		b_fc->Fill();
		b_bc->Fill();
	}
	// TTree::Write 写到当前目录，必须先切回输出文件（否则会写进只读的输入文件）
	of.cd();
	ot->Write("", TObject::kOverwrite);
	of.Close();
	return 0;
}

int CopyCsiCal(
	const char *out_path,
	const char *vme_path,
	const std::vector<long long> &vme_entries,
	int n_csi
) {
	TFile of(out_path, "update");
	TTree *ot = static_cast<TTree*>(of.Get("tree"));
	TFile vf(vme_path, "read");
	TTree *vt = static_cast<TTree*>(vf.Get("tree"));
	if (!of.IsOpen() || !ot || !vf.IsOpen() || !vt) {
		std::cerr << "Error: open cal source/output failed: "
			<< out_path << " / " << vme_path << "\n";
		return -1;
	}

	double cal[16] = {};
	vt->SetBranchAddress("cal", cal);

	TString leaf = TString::Format("c[%d]/D", n_csi);
	TBranch *b_cal = ot->Branch("cal", cal, leaf);
	if (!b_cal) {
		std::cerr << "Error: create cal branch in " << out_path << "\n";
		return -2;
	}

	const Long64_t n = ot->GetEntries();
	if (static_cast<size_t>(n) != vme_entries.size()) {
		std::cerr << "Warning: " << out_path << " has " << n
			<< " entries, trigger has " << vme_entries.size() << "\n";
	}
	for (Long64_t i = 0; i < n; ++i) {
		ot->GetEntry(i);
		const long long ve = static_cast<size_t>(i) < vme_entries.size()
			? vme_entries[i] : -1;
		if (ve >= 0 && ve < vt->GetEntries()) {
			vt->GetEntry(ve);
		} else {
			for (int k = 0; k < 16; ++k) {
				cal[k] = 0;
			}
		}
		b_cal->Fill();
	}
	// TTree::Write 写到当前目录，必须先切回输出文件（否则会写进只读的输入文件）
	of.cd();
	ot->Write("", TObject::kOverwrite);
	of.Close();
	return 0;
}

// TString GrainPath(
// 	const std::string &grain_dir,
// 	const std::string &trigger_type,
// 	int run
// ) {
// 	return TString::Format(
// 		"%s/grain_%s%04d.root",
// 		grain_dir.c_str(),
// 		TriggerStem(trigger_type).c_str(),
// 		run
// 	);
// }

int ReadTrigger(
	const char *path,
	std::vector<double> &trigger,
	std::vector<long long> &vme_entries,
	bool report = false
) {
	TFile ipf(path, "read");
	TTree *ipt = static_cast<TTree*>(ipf.Get("tree"));
	if (!ipt) {
		std::cerr << "Error: Get tree from " << path << " failed.\n";
		return -1;
	}

	forgerib::TriggerEvent trigger_event;
	forgerib::SetupInput(ipt, trigger_event);

	trigger.clear();
	vme_entries.clear();
	if (report) {
		printf("Reading trigger   0%%");
		fflush(stdout);
	}
	long long total = ipt->GetEntries();
	long long last_percentage = 0;
	for (long long entry = 0; entry < total; ++entry) {
		if (report && total > 0 && entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		if (!trigger_event.valid[forgerib::kTriggerMain]) continue;     
		trigger.push_back(trigger_event.time[forgerib::kTriggerMain]);
		vme_entries.push_back(trigger_event.vme_entry);
	}
	if (report) {
		printf("\b\b\b\b100%%\n");
	}
	ipf.Close();
	return 0;
}

// long long GetXiaTotal(const std::string &ingot_dir, int run) {
// 	TString filename = TString::Format("%s/trigger_%04d.root", ingot_dir.c_str(), run);
// 	TFile ipf(filename, "read");
// 	TTree *ipt = static_cast<TTree*>(ipf.Get("tree"));
// 	if (!ipt) {
// 		std::cerr << "Error: Get tree from " << filename << " failed.\n";
// 		return -1;
// 	}
// 	const long long total = ipt->GetEntries();
// 	ipf.Close();
// 	return total;
// }

bool NeedsVmeRun(const std::string &detector) {
	return detector == "t0csi"
		|| detector == "t1du"
		|| detector == "t1ud"
		|| detector == "t1dd"
		|| detector == "t1csiu"
		|| detector == "t1csid"
		|| detector == "tafd"
		|| detector == "tafcsi"
		|| detector == "vppac";
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("smelt", "smelt detectors into ingot output");
	options.add_options()
		("h,help", "Print usage")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"path"
		)
		(
			"t,trigger",
			"Trigger type: main, taf, or t1.",
			cxxopts::value<std::string>(),
			"trigger"
		)
		(
			"r,run",
			"XIA and VME run number.",
			cxxopts::value<int>(),
			"run"
		)
		(
			"x,xia-run",
			"XIA run number.",
			cxxopts::value<int>(),
			"run"
		)
		(
			"v,vme-run",
			"VME run number for fused detectors.",
			cxxopts::value<int>(),
			"run"
		)
		(
			"detectors",
			"Detector list: ppac, vppac, t0d1, t0d2, t0d3, t0d4, t0s, t0s_trace, t1su, t1sd, t0csi, t0csi_trace, beam, gagg, t1du, t1dd, t1csiu, t1csid, tafd, tafcsi",
			cxxopts::value<std::vector<std::string>>(),
			"detector"
		);
	options.parse_positional({"detectors"});
	auto parse_result = options.parse(argc, argv);

	if (parse_result.count("help")) {
		std::cout << options.help() << "\n";
		return 0;
	}

	const bool has_run = parse_result.count("run") > 0;
	const bool has_xia_run = parse_result.count("xia-run") > 0;
	const bool has_vme_run = parse_result.count("vme-run") > 0;
	if (!has_run && !has_xia_run) {
		std::cout << options.help() << "\n";
		return 0;
	}

	int xia_run = has_run
		? parse_result["run"].as<int>()
		: parse_result["xia-run"].as<int>();
	if (has_run && has_xia_run) {
		const int alias_run = parse_result["xia-run"].as<int>();
		if (xia_run != alias_run) {
			std::cerr << "Error: --run and --xia-run must match.\n";
			return -1;
		}
	}

	const std::vector<std::string> detectors = parse_result.count("detectors")
		? parse_result["detectors"].as<std::vector<std::string>>()
		: std::vector<std::string>();
	if (detectors.empty()) {
		std::cout << options.help() << "\n";
		return 0;
	}

	const std::string config_path = parse_result["config"].as<std::string>();
	forgerib::AppConfig config;
	if (forgerib::LoadConfig(config_path, config)) {
		return -1;
	}

	const std::string trigger_type = parse_result.count("trigger")
		? parse_result["trigger"].as<std::string>()
		: config.trigger;
	if (trigger_type != "main" && trigger_type != "t1") {
		std::cerr << "Error: Unsupported trigger " << trigger_type
			<< ", use main or t1.\n";
		return -1;
	}

	const bool require_vme = std::any_of(
		detectors.begin(),
		detectors.end(),
		NeedsVmeRun
	);
	if (require_vme && !has_vme_run && !has_run) {
		std::cerr << "Error: --vme-run is required for VME detectors.\n";
		return -1;
	}
	const int vme_run = parse_result.count("vme-run")
		? parse_result["vme-run"].as<int>()
		: parse_result["run"].as<int>();

	const std::string grit_dir = forgerib::JoinPath(config.workspace, config.paths.grit);
	const std::string grain_dir = forgerib::JoinPath(config.workspace, config.paths.grain);
	const std::string ingot_dir = forgerib::JoinPath(config.workspace, config.paths.ingot);

	std::vector<double> triggers;
	std::vector<long long> vme_entries;
	TString trigger_path = TString::Format(
		"%s/trigger_%s%04d.root",
		ingot_dir.c_str(),
		TriggerStem(trigger_type).c_str(),
		xia_run
	);
	if (ReadTrigger(trigger_path.Data(), triggers, vme_entries, true)) {
		std::cerr << "Error: Read trigger failed.\n";
		return -1;
	}

	for (std::string detector : detectors) {
		if (detector == "t1ud") detector = "t1du";

		int result = 0;
		if (detector == "ppac") {
			result = forgerib::SmeltPpac(
				triggers,
				GritInputPath(grit_dir, "ppac", xia_run).Data(),
				IngotOutputPath(ingot_dir, "ppac", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t0d1") {
			result = forgerib::SmeltT0d1(
				triggers,
				GritInputPath(grit_dir, "t0d1", xia_run).Data(),
				IngotOutputPath(ingot_dir, "t0d1", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t0d2") {
			result = forgerib::SmeltT0d2(
				triggers,
				GritInputPath(grit_dir, "t0d2", xia_run).Data(),
				IngotOutputPath(ingot_dir, "t0d2", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t0d3") {
			result = forgerib::SmeltT0d3(
				triggers,
				GritInputPath(grit_dir, "t0d3", xia_run).Data(),
				IngotOutputPath(ingot_dir, "t0d3", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t0d4") {
			result = forgerib::SmeltT0d4(
				triggers,
				GritInputPath(grit_dir, "t0d4", xia_run).Data(),
				IngotOutputPath(ingot_dir, "t0d4", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t0s") {
			result = forgerib::SmeltT0s(
				triggers,
				GritInputPath(grit_dir, "t0s", xia_run).Data(),
				IngotOutputPath(ingot_dir, "t0s", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t0s_trace") {
			result = forgerib::SmeltT0s_trace(
				triggers,
				GritInputPath(grit_dir, "t0s_trace", xia_run).Data(),
				IngotOutputPath(ingot_dir, "t0s_trace", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t1su") {
			result = forgerib::SmeltT1su(
				triggers,
				GritInputPath(grit_dir, "t1su", xia_run).Data(),
				IngotOutputPath(ingot_dir, "t1su", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t1sd") {
			result = forgerib::SmeltT1sd(
				triggers,
				GritInputPath(grit_dir, "t1sd", xia_run).Data(),
				IngotOutputPath(ingot_dir, "t1sd", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "beam") {
			result = forgerib::SmeltBeam(
				triggers,
				GritInputPath(grit_dir, "beam", xia_run).Data(),
				IngotOutputPath(ingot_dir, "beam", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "gagg") {
			result = forgerib::SmeltGagg(
				triggers,
				GritInputPath(grit_dir, "gagg", xia_run).Data(),
				IngotOutputPath(ingot_dir, "gagg", trigger_type, xia_run).Data(),
				1000.0,
				true
			);
		} else if (detector == "t0csi") {
			const TString bloom_path = IngotTriggeredPath(
				ingot_dir,
				"t0csi_bloom",
				trigger_type,
				xia_run
			);
			result = forgerib::SmeltT0CsiBloom(
				triggers,
				GritInputPath(grit_dir, "t0csi", xia_run).Data(),
				bloom_path.Data(),
				1000.0,
				true
			);
			if (!result) {
				const TString output_path = IngotTriggeredPath(
					ingot_dir,
					"t0csi",
					trigger_type,
					xia_run
				);
				result = forgerib::SmeltT0Csi(
					vme_entries,
					bloom_path.Data(),
					GritVmeInputPath(grit_dir, "t0csi", vme_run).Data(),
					output_path.Data(),
					true
				);
			}
		} else if (detector == "t0csi_trace") {
			std::cout << "T0csi trace.\n";
			const TString bloom_path = IngotTriggeredPath(
				ingot_dir,
				"t0csi_trace",
				trigger_type,
				xia_run
			);
			// std::cout << bloom_path.Data() << "\n";
			// std::cout << GritInputPath(grit_dir, "t0csi_trace", xia_run).Data() << "\n";
			// std::cout << __LINE__ << "\n";
			result = forgerib::SmeltT0Csi_trace(
				triggers,
				GritInputPath(grit_dir, "t0csi_trace", xia_run).Data(),
				bloom_path.Data(),
				1000.0,
				true
			);
			// std::cout << __LINE__ << "\n";

		} else if (detector == "t1du") {
			result = forgerib::SmeltDetector<forgerib::DssdEvent>(
				vme_entries,
				GritVmeInputPath(grit_dir, "t1du", vme_run).Data(),
				IngotOutputPath(ingot_dir, "t1du", trigger_type, xia_run).Data(),
				"t1du",
				true
			);
			if (!result) {
				result = CopyDssdCal(
					IngotOutputPath(ingot_dir, "t1du", trigger_type, xia_run).Data(),
					GritVmeInputPath(grit_dir, "t1du", vme_run).Data(),
					vme_entries
				);
			}
		} else if (detector == "t1dd") {
			result = forgerib::SmeltDetector<forgerib::DssdEvent>(
				vme_entries,
				GritVmeInputPath(grit_dir, "t1dd", vme_run).Data(),
				IngotOutputPath(ingot_dir, "t1dd", trigger_type, xia_run).Data(),
				"t1dd",
				true
			);
		} else if (detector == "t1csiu") {
			result = forgerib::SmeltDetector<forgerib::CsiEvent<4>>(
				vme_entries,
				GritVmeInputPath(grit_dir, "t1csiu", vme_run).Data(),
				IngotOutputPath(ingot_dir, "t1csiu", trigger_type, xia_run).Data(),
				"t1csiu",
				true
			);
			if (!result) {
				result = CopyCsiCal(
					IngotOutputPath(ingot_dir, "t1csiu", trigger_type, xia_run).Data(),
					GritVmeInputPath(grit_dir, "t1csiu", vme_run).Data(),
					vme_entries,
					4
				);
			}
		} else if (detector == "t1csid") {
			result = forgerib::SmeltDetector<forgerib::CsiEvent<13>>(
				vme_entries,
				GritVmeInputPath(grit_dir, "t1csid", vme_run).Data(),
				IngotOutputPath(ingot_dir, "t1csid", trigger_type, xia_run).Data(),
				"t1csid",
				true
			);
		} else if (detector == "tafd") {
			result = forgerib::SmeltDetector<forgerib::TafdEvent>(
				vme_entries,
				GritVmeInputPath(grit_dir, "tafd", vme_run).Data(),
				IngotOutputPath(ingot_dir, "tafd", trigger_type, xia_run).Data(),
				"tafd",
				true
			);
			if (!result) {
				result = CopyTafdCal(
					IngotOutputPath(ingot_dir, "tafd", trigger_type, xia_run).Data(),
					GritVmeInputPath(grit_dir, "tafd", vme_run).Data(),
					vme_entries
				);
			}
		} else if (detector == "tafcsi") {
			result = forgerib::SmeltDetector<forgerib::CsiEvent<12>>(
				vme_entries,
				GritVmeInputPath(grit_dir, "tafcsi", vme_run).Data(),
				IngotOutputPath(ingot_dir, "tafcsi", trigger_type, xia_run).Data(),
				"tafcsi",
				true
			);
			if (!result) {
				result = CopyCsiCal(
					IngotOutputPath(ingot_dir, "tafcsi", trigger_type, xia_run).Data(),
					GritVmeInputPath(grit_dir, "tafcsi", vme_run).Data(),
					vme_entries,
					12
				);
			}
		} else if (detector == "vppac") {
			result = forgerib::SmeltDetector<forgerib::PpacEvent>(
				vme_entries,
				GritVmeInputPath(grit_dir, "ppac", vme_run).Data(),
				IngotOutputPath(ingot_dir, "vppac", trigger_type, xia_run).Data(),
				"vppac",
				true
			);
		} else {
			std::cerr << "Error: Unsupported detector " << detector << ".\n";
			return -1;
		}

		if (result) {
			std::cerr << "Error: Smelt " << detector << " failed.\n";
			return result;
		}
	}

	return 0;
}