//#pragma once
#ifndef __RAW_BEAM_EVENT__
#define __RAW_BEAM_EVENT__
#include <string>

#include <TTree.h>

namespace forgerib {

constexpr int kBeamTof3_1 = 0;
constexpr int kBeamTof3_2 = 1;
constexpr int kBeamTof3_3 = 2;
constexpr int kBeamTof3_4 = 3;
constexpr int kBeamTof5_1 = 4;
constexpr int kBeamTof5_2 = 5;
constexpr int kBeamTof5_3 = 6;
constexpr int kBeamTof5_4 = 7;

struct RawBeamEvent {
	int type;
	int energy;
	double time;
	bool cv;
};

void SetupInput(TTree *tree, RawBeamEvent &event, const std::string &prefix = "");

void SetupOutput(TTree *tree, RawBeamEvent &event);


}

#endif