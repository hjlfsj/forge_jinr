#pragma once

#include <string>

#include <TTree.h>

namespace forgerib {

struct BeamEvent {
	unsigned long long flag;
	bool valid[8];
	double time[8];
	double tof;
};

void SetupInput(TTree *tree, BeamEvent &event, const std::string &prefix = "");

void SetupOutput(TTree *tree, BeamEvent &event);

}