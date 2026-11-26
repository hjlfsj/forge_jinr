#ifndef __SMELT_SIFTING_H__
#define __SMELT_SIFTING_H__

#include <vector>
#include <cstddef>

namespace forgerib {

int Sift(
	const std::vector<double> &xia_times,
	const std::vector<long long> &xia_entries,
	const std::vector<double> &vme_times,
	const std::vector<long long> &vme_entries,
	const char *output_path,
	const double window,
	const size_t range,
	const size_t vme_window,
	const size_t xia_window,
	bool report
);

}

#endif