#include "include/event/ingot/gagg_event.h"

namespace forgerib {

void SetupInput(TTree *tree, GaggEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "num").c_str(), &event.num);
	tree->SetBranchAddress((prefix + "energy").c_str(), event.energy);
	tree->SetBranchAddress((prefix + "time").c_str(), event.time);
	tree->SetBranchAddress((prefix + "index").c_str(), event.index);
}

void SetupOutput(TTree *tree, GaggEvent &event) {
	tree->Branch("num", &event.num, "num/I");
	tree->Branch("energy", event.energy, "energy[num]/I");
	tree->Branch("time", event.time, "time[num]/D");
	tree->Branch("index", event.index, "index[num]/I");
}

void Reset(GaggEvent &event) {
	event.num = 0;
}

void Update(GaggEvent &event, const int idx, const int e, const double t) {
	// 若 index 已存在：能量更大则替换；否则追加
	for (int i = 0; i < event.num; ++i) {
		if (event.index[i] == idx) {
			if (e > event.energy[i]) {
				event.energy[i] = e;
				event.time[i] = t;
			}
			return;
		}
	}
	if (event.num >= GaggEvent::kMaxNum) return;
	event.index[event.num] = idx;
	event.energy[event.num] = e;
	event.time[event.num] = t;
	++event.num;
}

} // namespace forgerib