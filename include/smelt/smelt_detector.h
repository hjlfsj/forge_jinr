#ifndef __FUSE_DETECTOR_H__
#define __FUSE_DETECTOR_H__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <TFile.h>
#include <TTree.h>

#include "include/event/ingot/align_event.h"

namespace forgerib {

// 在 trigger_vme_XXXX.root 里按时间匹配 vme_times，得到对应的 VME entry 号
// vme_times : 来自 grain 的 vme_time，已从小到大排序
// 匹配不到的位置填 -1
// 返回 0 成功，-1 失败
inline int MatchVmeEntries(
	const std::vector<double> &vme_times,
	const char *trigger_vme_path,
	std::vector<long long> &vme_entries,
	double tolerance = 1.0,
	bool report = false
) {
	TFile ipf(trigger_vme_path, "read");
	if (!ipf.IsOpen()) {
		std::cerr << "Error: Cannot open " << trigger_vme_path << ".\n";
		return -1;
	}
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from " << trigger_vme_path << " failed.\n";
		ipf.Close();
		return -1;
	}

	bool valid = false;
	long long time = 0;
	if (ipt->SetBranchAddress("valid", &valid) < 0
		|| ipt->SetBranchAddress("time", &time) < 0) {
		std::cerr << "Error: cannot bind valid/time in "
			<< trigger_vme_path << ".\n";
		ipf.Close();
		return -1;
	}

	const long long total = (long long)vme_times.size();
	const long long vme_total = ipt->GetEntries();
	vme_entries.assign(total, -1);

	long long last_percentage = 0;
	if (report) {
		printf("Matching VME entries   0%%");
		fflush(stdout);
	}

	long long vme_cursor = 0;
	for (long long entry = 0; entry < total; ++entry) {
		if (report && total > 0
			&& entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}

		const double target = vme_times[entry];

		// 跳过明显小于 target - tolerance 的 entry（只考虑 valid）
		while (vme_cursor < vme_total) {
			ipt->GetEntry(vme_cursor);
			if (valid && (double)time >= target - tolerance) break;
			++vme_cursor;
		}
		if (vme_cursor >= vme_total) continue;

		// 比较 vme_cursor 和 vme_cursor-1 哪个更接近 target（只考虑 valid）
		long long best_entry = -1;
		double best_diff = 1e18;

		ipt->GetEntry(vme_cursor);
		if (valid) {
			const double d = std::fabs((double)time - target);
			if (d < best_diff) { best_diff = d; best_entry = vme_cursor; }
		}
		if (vme_cursor > 0) {
			ipt->GetEntry(vme_cursor - 1);
			if (valid) {
				const double d = std::fabs((double)time - target);
				if (d < best_diff) { best_diff = d; best_entry = vme_cursor - 1; }
			}
		}

		if (best_entry >= 0 && best_diff < tolerance) {
			vme_entries[entry] = best_entry;
		}
	}

	if (report) printf("\b\b\b\b100%%\n");
	ipf.Close();

	long long matched = 0;
	for (auto e : vme_entries) if (e >= 0) ++matched;
	std::cout << "VME entry matched: " << matched << " / " << total << "\n";
	return 0;
}

// 用 vme_entries 从探测器 VME 文件里取数
template<typename Event>
int SmeltDetector(
	const std::vector<long long> &vme_entries,
	const char *vme_path,
	const char *output_path,
	const char *detector_name,
	bool report
) {
	TFile vme_file(vme_path, "read");
	TTree *vme_tree = (TTree*)vme_file.Get("tree");
	if (!vme_tree) {
		std::cerr << "Error: Get tree from " << vme_path << " failed.\n";
		return -2;
	}
	Event vme_event;
	SetupInput(vme_tree, vme_event);

	TFile opf(output_path, "recreate");
	TTree opt("tree", TString::Format("fused %s", detector_name));
	Event event;
	SetupOutput(&opt, event);

	if (report) {
		printf("Fusing %s   0%%", detector_name);
		fflush(stdout);
	}
	long long total = (long long)vme_entries.size();
	long long last_percentage = 0;
	for (long long entry = 0; entry < total; ++entry) {
		if (report && total > 0
			&& entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		if (vme_entries[entry] >= 0) {
			vme_tree->GetEntry(vme_entries[entry]);
			event = vme_event;
		} else {
			Reset(event);
		}
		opt.Fill();
	}
	if (report) printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	opf.Close();
	vme_file.Close();
	return 0;
}

}

#endif