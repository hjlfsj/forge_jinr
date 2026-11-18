#include "include/smelt/smelt_gagg.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

#include "include/event/ingot/csi_event.h"
#include "include/event/ore/raw_csi_event.h"

namespace forgerib {

namespace {

constexpr size_t kSlotNum = 1000;
constexpr int kGaggNum = 14;

void ResetGaggEvent(CsiEvent<kGaggNum> &csi) {
	csi.flag = 0;
	for (int i = 0; i < kGaggNum; ++i) {
		csi.valid[i] = false;
		csi.time[i] = 0.0;
		csi.energy[i] = 0;
	}
}

void UpdateGaggEvent(const RawCsiEvent &raw, CsiEvent<kGaggNum> &csi) {
	if (raw.index < 0 || raw.index >= kGaggNum) return;
	if (!csi.valid[raw.index] || raw.energy > csi.energy[raw.index]) {
		csi.valid[raw.index] = true;
		csi.flag |= 1ull << raw.index;
		csi.time[raw.index] = raw.time;
		csi.energy[raw.index] = raw.energy;
	}
}

int SmeltWithTrigger(
	const std::vector<double> &trigger_time,
	const char *path,
	const char *output_path,
	const double window,
	bool report
) {
	TFile opf(output_path, "recreate");
	TH1F forge_window("fw", "forge window", 200, -window, window);
	TTree opt("tree", "forged gagg");
	CsiEvent<kGaggNum> csi;
	SetupOutput(&opt, csi);

	TFile ipf(path, "read");
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from " << path << " failed.\n";
		return -1;
	}
	RawCsiEvent raw;
	SetupInput(ipt, raw);

	CsiEvent<kGaggNum> slots[kSlotNum];
	for (size_t i = 0; i < kSlotNum; ++i) {
		ResetGaggEvent(slots[i]);
	}
	size_t tofill_entry = 0;

	long long total = ipt->GetEntriesFast();
	long long last_percentage = 0;
	if (report) {
		printf("Forging gagg with trigger   0%%");
		fflush(stdout);
	}
	for (long long entry = 0; entry < total; ++entry) {
		if (report && entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		double min_time = 2*window;
		size_t min_time_entry = 0;
		for (
			auto iter = std::lower_bound(
				trigger_time.begin()+tofill_entry,
				trigger_time.end(),
				raw.time - window
			);
			iter != trigger_time.end();
			++iter
		) {
			if (*iter > raw.time + window) break;
			forge_window.Fill(raw.time - *iter);
			if (std::fabs(*iter - raw.time) < min_time) {
				min_time = std::fabs(*iter - raw.time);
				min_time_entry = iter - trigger_time.begin();
			}
		}
		if (min_time > window) continue;
		if (min_time_entry < tofill_entry) continue;
		if (min_time_entry - tofill_entry >= kSlotNum) {
			for (size_t fill = tofill_entry; fill <= min_time_entry-kSlotNum; ++fill) {
				csi = slots[fill%kSlotNum];
				opt.Fill();
				ResetGaggEvent(slots[fill%kSlotNum]);
			}
			tofill_entry = min_time_entry - kSlotNum + 1;
		}
		UpdateGaggEvent(raw, slots[min_time_entry%kSlotNum]);
	}
	for (size_t fill = tofill_entry; fill < trigger_time.size(); ++fill) {
		csi = slots[fill%kSlotNum];
		opt.Fill();
	}
	if (report) printf("\b\b\b\b100%%\n");

	opf.cd();
	forge_window.Write();
	opt.Write();
	ipf.Close();
	opf.Close();
	return 0;
}

int SmeltWithoutTrigger(
	const char *path,
	const char *output_path,
	const double window,
	bool report
) {
	TFile opf(output_path, "recreate");
	TH1F forge_window("fw", "forge window", 200, -window, window);
	TTree opt("tree", "forged gagg");
	CsiEvent<kGaggNum> csi;
	SetupOutput(&opt, csi);

	TFile ipf(path, "read");
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from " << path << " failed.\n";
		return -1;
	}
	RawCsiEvent raw;
	SetupInput(ipt, raw);

	ResetGaggEvent(csi);
	double ref_time = -1.0;
	long long total_entries = ipt->GetEntriesFast();
	long long last_percentage = 0;
	if (report) {
		printf("Forging gagg without trigger   0%%");
		fflush(stdout);
	}
	for (long long entry = 0; entry < total_entries; ++entry) {
		if (report && entry * 100ll / total_entries > last_percentage) {
			last_percentage = entry * 100ll / total_entries;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		if (csi.flag == 0) {
			UpdateGaggEvent(raw, csi);
			ref_time = raw.time;
		} else if (raw.time - ref_time < window) {
			forge_window.Fill(raw.time - ref_time);
			UpdateGaggEvent(raw, csi);
		} else {
			opt.Fill();
			ResetGaggEvent(csi);
			UpdateGaggEvent(raw, csi);
			ref_time = raw.time;
		}
	}
	if (csi.flag != 0) opt.Fill();
	if (report) printf("\b\b\b\b100%%\n");

	opf.cd();
	forge_window.Write();
	opt.Write();
	ipf.Close();
	opf.Close();
	return 0;
}

} // namespace

int SmeltGagg(
	const std::vector<double> &trigger_time,
	const char *path,
	const char *output_path,
	const double window,
	bool report
) {
	if (trigger_time.empty()) {
		return SmeltWithoutTrigger(path, output_path, window, report);
	}
	return SmeltWithTrigger(trigger_time, path, output_path, window, report);
}

}