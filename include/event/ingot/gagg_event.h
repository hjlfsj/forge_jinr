#ifndef __GAGG_EVENT_H__
#define __GAGG_EVENT_H__

#include <string>

#include <TTree.h>

namespace forgerib {

// GAGG 变长数组事件
// num      : 命中通道数
// energy[] : 每个命中的能量
// time[]   : 每个命中的时间
// index[]  : 每个命中的通道号
struct GaggEvent {
	static constexpr int kMaxNum = 256;   // 单事件缓存上限，按需调大

	int num;
	int energy[kMaxNum];
	double time[kMaxNum];
	int index[kMaxNum];
};

// 声明（实现放在 gagg_event.cpp）
void SetupInput(TTree *tree, GaggEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, GaggEvent &event);
void Reset(GaggEvent &event);
void Update(GaggEvent &event, const int index, const int energy, const double time);

} // namespace forgerib

#endif