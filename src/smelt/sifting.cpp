#include "include/smelt/sifting.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include <TGraph.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include "include/event/ingot/align_event.h"

namespace forgerib {

int Sift(
	const std::vector<double> &xia_times,          // ns
	const std::vector<long long> &xia_entries,
	const std::vector<double> &vme_times,          // ns（已在读取时过滤）
	const std::vector<long long> &vme_entries,
	const char *output_path,
	const double window,                           // ns
	const size_t range,
	const size_t vme_window,
	const size_t xia_window,
	bool report
) {
	TFile opf(output_path, "recreate");
	if (!opf.IsOpen()) {
		std::cerr << "Error: Cannot create output file " << output_path
			<< ". Check permissions/disk.\n";
		return -1;
	}
	TGraph offset_vs_group[3];
	TGraph window_vs_group[3];
	TTree opt("tree", "alignment result");
	AlignEvent event;
	SetupOutput(&opt, event);

	if (xia_times.size() < 4 || vme_times.size() < 4) {
		std::cerr << "Error: Not enough events. XIA=" << xia_times.size()
			<< ", VME=" << vme_times.size() << ".\n";
		return -1;
	}
	if (xia_entries.size() != xia_times.size()
		|| vme_entries.size() != vme_times.size()) {
		std::cerr << "Error: times/entries size mismatch.\n";
		return -1;
	}

	std::vector<double> vme_dt;
	for (size_t i = 1; i < vme_times.size(); ++i) {
		vme_dt.push_back(vme_times[i] - vme_times[i-1]);
	}

	size_t anchor_xia = static_cast<size_t>(-1);
	size_t anchor_vme = static_cast<size_t>(-1);
	bool found = false;

	const size_t xia_search_end   = std::min(xia_times.size() / 2, xia_times.size());
	const size_t vme_search_begin = 3;
	const size_t vme_search_end   = vme_times.size();

	if (vme_search_end <= vme_search_begin + 3) {
		std::cerr << "Error: VME search range too small: ["
			<< vme_search_begin << ", " << vme_search_end << ").\n";
		return -1;
	}

	for (anchor_xia = 3; anchor_xia < xia_search_end; ++anchor_xia) {
		double xia_dt0 = xia_times[anchor_xia-2] - xia_times[anchor_xia-3];
		double xia_dt1 = xia_times[anchor_xia-1] - xia_times[anchor_xia-2];
		double xia_dt2 = xia_times[anchor_xia]   - xia_times[anchor_xia-1];
		for (anchor_vme = vme_search_begin;
		     anchor_vme + 2 < vme_search_end; ++anchor_vme) {
			double vme_dt0 = vme_dt[anchor_vme-3];
			double vme_dt1 = vme_dt[anchor_vme-2];
			double vme_dt2 = vme_dt[anchor_vme-1];
			if (fabs(vme_dt0 - xia_dt0) < window
				&& fabs(vme_dt1 - xia_dt1) < window
				&& fabs(vme_dt2 - xia_dt2) < window) {
				found = true;
				break;
			}
		}
		if (found) break;
	}
	if (!found) {
		std::cerr << "Error: Could not find initial anchor (4 consecutive events) "
			<< "in XIA window [3," << xia_search_end << "), "
			<< "VME window [" << vme_search_begin << ","
			<< vme_search_end << ").\n";
		return -1;
	}

	if (report) {
		std::cout << "Anchor at XIA entry " << anchor_xia
			<< ", VME entry " << anchor_vme << " (4-event anchor)\n";
	}

	size_t align_events = 0;
	double last_matched_vme_time_ns = -1.0;

	for (int i = 3; i >= 0; --i) {
		event.vme_entry = vme_entries[anchor_vme-i];
		event.vme_time  = vme_times[anchor_vme-i] / 1000.0;   // ns -> us
		event.xia_entry = xia_entries[anchor_xia-i];
		event.xia_time  = xia_times[anchor_xia-i];            // ns
		opt.Fill();
		++align_events;
	}
	last_matched_vme_time_ns = vme_times[anchor_vme];

	if (report) {
		printf("Sifting   0%%");
		fflush(stdout);
	}
	size_t last_percentage = 0;
	bool is_lookhead = true;
	while (anchor_xia < xia_times.size()-1 && anchor_vme < vme_times.size()-1 && is_lookhead) {
		if (report && anchor_xia * 100 / xia_times.size() > last_percentage) {
			last_percentage = anchor_xia * 100 / xia_times.size();
			printf("\b\b\b\b%3lu%%", last_percentage);
			fflush(stdout);
		}
		bool match = false;
		const size_t xia_lim = std::min(anchor_xia + range, xia_times.size() - 1);
		const size_t vme_lim = std::min(anchor_vme + range, vme_times.size() - 1);
		for (size_t i = anchor_xia+1; i <= xia_lim; ++i) {
			double xia_dt = xia_times[i] - xia_times[anchor_xia];
			for (size_t j = anchor_vme+1; j <= vme_lim; ++j) {
				double vme_dt_val = vme_times[j] - vme_times[anchor_vme];
				if (fabs(xia_dt - vme_dt_val) < window) {
					event.vme_entry = vme_entries[j];
					event.vme_time  = vme_times[j] / 1000.0;   // ns -> us
					event.xia_entry = xia_entries[i];
					event.xia_time  = xia_times[i];            // ns
					opt.Fill();
					++align_events;
					last_matched_vme_time_ns = vme_times[j];
					anchor_vme = j;
					anchor_xia = i;
					match = true;
					break;
				}
			}
			if (match) break;
		}
		if (!match && !(anchor_xia+vme_window < xia_times.size()-1 && anchor_vme+xia_window < vme_times.size()-1))
		{
			is_lookhead = false;
		}
		found = false;
		while (!match && !found && anchor_xia+vme_window < xia_times.size()-1 && anchor_vme+xia_window < vme_times.size()-1)
		{
			std::cout << "Matching stalled at XIA idx " << anchor_xia << ". Re-anchoring..." << std::endl;
			size_t ref_xia = anchor_xia + 1;
			size_t ref_vme = anchor_vme + 1;
			found = false;
			const size_t xia_re_end = std::min(xia_window + ref_xia, xia_times.size());
			const size_t vme_re_end = std::min(vme_window + ref_vme, vme_times.size());
			for (anchor_xia = ref_xia + 2; anchor_xia < xia_re_end; ++anchor_xia)
			{
				double xia_dt0 = xia_times[anchor_xia-1] - xia_times[anchor_xia-2];
				double xia_dt1 = xia_times[anchor_xia] - xia_times[anchor_xia-1];
				for (anchor_vme = ref_vme; anchor_vme + 1 < vme_re_end; ++anchor_vme)
				{
					double vme_dt0 = vme_times[anchor_vme+1] - vme_times[anchor_vme];
					double vme_dt1 = vme_times[anchor_vme+2] - vme_times[anchor_vme+1];
					if (fabs(vme_dt0-xia_dt0) < window && fabs(vme_dt1-xia_dt1) < window)
					{
						found = true;
						anchor_vme += 2;
						break;
					}
				}
				if (found)
				{
					for (int i = 2; i >= 0; --i)
					{
						event.vme_entry = vme_entries[anchor_vme-i];
						event.vme_time  = vme_times[anchor_vme-i] / 1000.0;  // ns -> us
						event.xia_entry = xia_entries[anchor_xia-i];
						event.xia_time  = xia_times[anchor_xia-i];           // ns
						opt.Fill();
						++align_events;
					}
					last_matched_vme_time_ns = vme_times[anchor_vme];
					std::cout << "Anchor at XIA entry " << anchor_xia
						<< ", VME entry " << anchor_vme << "\n";
					break;
				}
			}

			if (!found)
			{
				std::cerr << "Could not find anchor in XIA window "
					<< xia_window << ", VME window " << vme_window << ".\n";
			}
		}
	}
	if (report) printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	opf.Close();

	if (report) {
		std::cout
			<< "XIA align rate "
			<< align_events << " / " << xia_times.size() << " "
			<< double(align_events) / double(xia_times.size()) << "\n"
			<< "VME align rate "
			<< align_events << " / " << vme_times.size() << " "
			<< double(align_events) / double(vme_times.size()) << "\n";
		if (last_matched_vme_time_ns >= 0.0) {
			std::cout << "Last matched VME time = "
				<< last_matched_vme_time_ns / 1e9
				<< " s\n";
		}
	}
	return 0;
}

} // namespace forgerib