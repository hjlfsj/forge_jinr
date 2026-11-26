#include <algorithm>
#include <array>
#include <cctype>
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
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

// 列出目录下所有 *.root（POSIX，无 ROOT 依赖）
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

// 在文件名（不含目录）中查找 run 名，格式如 run08_01-04.root。
// 判定规则：文件名中出现 "_<run_tag>_" 或 "-<run_tag>-" 或 以 "<run_tag>_" 开头段
// 这里直接用 "run<run_tag>_" 更贴合实际命名，避免误匹配。
bool FileMatchesRun(const std::string &basename, const std::string &run_tag) {
    // 去掉 .root 后缀
    std::string stem = basename;
    if (stem.size() > 5 && stem.compare(stem.size() - 5, 5, ".root") == 0) {
        stem = stem.substr(0, stem.size() - 5);
    }

    // 形如 run08_01-04
    const std::string prefix = "run" + run_tag + "_";
    if (stem.compare(0, prefix.size(), prefix) == 0) return true;

    // 兜底：文件名中任意位置出现 "_08_" / "-08-"
    const std::string n1 = "_" + run_tag + "_";
    const std::string n2 = "-" + run_tag + "-";
    if (stem.find(n1) != std::string::npos) return true;
    if (stem.find(n2) != std::string::npos) return true;
    return false;
}

// 从目录中筛选出属于指定 run 的所有 *.root，按文件名排序返回
std::vector<std::string> CollectRunFiles(const std::string &dir, int vme_run) {
    char run_tag[8];
    std::snprintf(run_tag, sizeof(run_tag), "%02d", vme_run);  // "08"

    std::vector<std::string> all = ListRootFiles(dir);
    std::vector<std::string> matched;
    matched.reserve(all.size());
    for (auto &f : all) {
        std::string base = f.substr(f.find_last_of('/') + 1);
        if (FileMatchesRun(base, run_tag)) {
            matched.push_back(f);
        }
    }
    return matched;
}

// ---------------------------------------------------------------------------
// Go4 AnalysisxTree 的 NeEvent 内容（输入缓冲，分支名 = "NeEvent." + 成员对应的原名）
// 输出时压缩成 NeOutEvent，只保留超过阈值的道
// ---------------------------------------------------------------------------
struct NeEventData {
	// S1 - S6：每层 8 条（a: ADC，t: TDC，c: 刻度后能量）
	UShort_t s1a[8];   UShort_t s1t[8];
	UShort_t r1a[16];  UShort_t r1t[16];
	UShort_t s2a[8];   UShort_t s2t[8];
	UShort_t r2a[16];  UShort_t r2t[16];
	UShort_t s3a[8];   UShort_t s3t[8];
	UShort_t r3a[16];  UShort_t r3t[16];
	UShort_t s4a[8];   UShort_t s4t[8];
	UShort_t r4a[16];  UShort_t r4t[16];
	UShort_t s5a[8];   UShort_t s5t[8];
	UShort_t r5a[16];  UShort_t r5t[16];
	UShort_t s6a[8];   UShort_t s6t[8];
	UShort_t r6a[16];  UShort_t r6t[16];
	// C / X / Y：16 条
	UShort_t csa[16];  UShort_t cst[16];
	UShort_t xa[16];   UShort_t xt[16];
	UShort_t ya[16];   UShort_t yt[16];
	// 闪烁体 / MWPC
	UShort_t f3[4];    UShort_t tf3[4];
	UShort_t f5[4];    UShort_t tf5[4];
	UShort_t tmwpc[4];
	// 标量
	UInt_t   scaler[16];
	ULong64_t scaler_xia;
	// 全局量
	Double_t af3;
	Double_t af5;
	Double_t tof;
	Double_t xbt;
	Double_t ybt;
	Double_t xbd;
	Double_t ybd;
	// 刻度后能量
	Double_t s1c[8];   Double_t r1c[16];
	Double_t s2c[8];   Double_t r2c[16];
	Double_t s3c[8];   Double_t r3c[16];
	Double_t s4c[8];   Double_t r4c[16];
	Double_t s5c[8];   Double_t r5c[16];
	Double_t s6c[8];   Double_t r6c[16];
	Double_t csc[16];  Double_t xc[16];  Double_t yc[16];
	// x1/y1（32 条）与命中数
	UShort_t nx1;      UShort_t ny1;
	UShort_t nx2;      UShort_t ny2;
	UShort_t x1[32];   UShort_t y1[32];
	UShort_t x2[32];   UShort_t y2[32];
	UInt_t   rx1;      UInt_t   rx2;
	UInt_t   ry1;      UInt_t   ry2;
	Int_t    nevent;
	Int_t    trigger;
	Int_t    subevents;
	Int_t    evsize;
	UInt_t   mtime;
	UInt_t   utime;
	UInt_t   utime_old;
};

// 输出分支的叶子类型代码（与 TTree::Branch 的 leaf list 一致）
// 加 Ne 前缀，避免与本程序其它头文件中的同名模板冲突
template <typename T>
const char *NeLeafCode() {
	if constexpr (std::is_same_v<T, UShort_t>) return "s";
	else if constexpr (std::is_same_v<T, UInt_t>) return "i";
	else if constexpr (std::is_same_v<T, ULong64_t>) return "l";
	else if constexpr (std::is_same_v<T, Int_t>) return "I";
	else if constexpr (std::is_same_v<T, Double_t>) return "D";
	else return "";
}

template <typename T>
void AddNeBranch(TTree *tree, const char *name, T &value) {
	TBranch *br = nullptr;
	if constexpr (std::is_array_v<T>) {
		constexpr size_t N = std::extent_v<T>;
		using Elem = std::remove_extent_t<T>;
		br = tree->Branch(
			name,
			static_cast<void *>(value),
			TString::Format("%s[%zu]/%s", name, N, NeLeafCode<Elem>()).Data()
		);
	} else {
		br = tree->Branch(
			name,
			&value,
			TString::Format("%s/%s", name, NeLeafCode<T>()).Data()
		);
	}
	if (!br) std::cerr << "Warning: cannot create branch " << name << "\n";
}

// 变长数组分支：长度由同树中的计数分支 count 决定（只写前 num 个元素）
template <typename T, size_t N>
void AddNeBranchV(TTree *tree, const char *name, const char *count, T (&arr)[N]) {
	TBranch *br = tree->Branch(
		name,
		static_cast<void *>(arr),
		TString::Format("%s[%s]/%s", name, count, NeLeafCode<T>()).Data()
	);
	if (!br) std::cerr << "Warning: cannot create branch " << name << "\n";
}

// 字段表：X(成员名, 原分支名)，用于绑定输入（分支名 = "NeEvent." + 名字）
#define NE_EVENT_FIELDS(X) \
	X(s1a, "S1a") X(s1t, "S1t") X(r1a, "R1a") X(r1t, "R1t") \
	X(s2a, "S2a") X(s2t, "S2t") X(r2a, "R2a") X(r2t, "R2t") \
	X(s3a, "S3a") X(s3t, "S3t") X(r3a, "R3a") X(r3t, "R3t") \
	X(s4a, "S4a") X(s4t, "S4t") X(r4a, "R4a") X(r4t, "R4t") \
	X(s5a, "S5a") X(s5t, "S5t") X(r5a, "R5a") X(r5t, "R5t") \
	X(s6a, "S6a") X(s6t, "S6t") X(r6a, "R6a") X(r6t, "R6t") \
	X(csa, "Csa") X(cst, "Cst") X(xa, "Xa") X(xt, "Xt") \
	X(ya, "Ya") X(yt, "Yt") \
	X(f3, "F3") X(tf3, "tF3") X(f5, "F5") X(tf5, "tF5") \
	X(tmwpc, "tMWPC") \
	X(scaler, "scaler") X(scaler_xia, "scalerXIA") \
	X(af3, "af3") X(af5, "af5") X(tof, "ToF") \
	X(xbt, "xbt") X(ybt, "ybt") X(xbd, "xbd") X(ybd, "ybd") \
	X(s1c, "S1c") X(r1c, "R1c") X(s2c, "S2c") X(r2c, "R2c") \
	X(s3c, "S3c") X(r3c, "R3c") X(s4c, "S4c") X(r4c, "R4c") \
	X(s5c, "S5c") X(r5c, "R5c") X(s6c, "S6c") X(r6c, "R6c") \
	X(csc, "Csc") X(xc, "Xc") X(yc, "Yc") \
	X(nx1, "nx1") X(ny1, "ny1") X(nx2, "nx2") X(ny2, "ny2") \
	X(x1, "x1") X(y1, "y1") X(x2, "x2") X(y2, "y2") \
	X(rx1, "rx1") X(rx2, "rx2") X(ry1, "ry1") X(ry2, "ry2") \
	X(nevent, "nevent") X(trigger, "trigger") \
	X(subevents, "subevents") X(evsize, "evsize") \
	X(mtime, "mtime") X(utime, "utime") X(utime_old, "utime_old")

// ---------------------------------------------------------------------------
// 压缩输出结构：像 crush_vme_old 一样只保留有事件（超过阈值）的道
//   每个探测器的命中数组都按能量从大到小排列，num 为命中数
//   输出分支为变长数组 [num]，strip/energy/time/cal 与排好序的 energy 一一对应
//   6 个 ADSSD 分开存储（a1…a6）
// 阈值直接写在程序里，见 DefaultNeThreshold()
// ---------------------------------------------------------------------------
struct NeThreshold {
	std::array<double, 16> adssd_front{};
	std::array<double, 8>  adssd_back{};
	std::array<double, 16> csi{};
	std::array<double, 16> dssd_front{};
	std::array<double, 16> dssd_back{};
};

// 6 个 ADSSD 分开存储，每个的命中按能量从大到小排序
constexpr int kAdssdCount = 6;
constexpr int kAdssdFrontStrips = 16;
constexpr int kAdssdBackStrips = 8;

struct AdssdHits {
	Int_t    front_num = 0;
	UShort_t front_strip[kAdssdFrontStrips] = {};
	UShort_t front_energy[kAdssdFrontStrips] = {};
	UShort_t front_time[kAdssdFrontStrips] = {};
	Double_t front_cal[kAdssdFrontStrips] = {};
	Int_t    back_num = 0;
	UShort_t back_strip[kAdssdBackStrips] = {};
	UShort_t back_energy[kAdssdBackStrips] = {};
	UShort_t back_time[kAdssdBackStrips] = {};
	Double_t back_cal[kAdssdBackStrips] = {};
};

struct CsiHits {
	Int_t    num = 0;
	UShort_t strip[16] = {};
	UShort_t energy[16] = {};
	UShort_t time[16] = {};
	Double_t cal[16] = {};
};

struct DssdHits {
	Int_t    front_num = 0;
	UShort_t front_strip[16] = {};
	UShort_t front_energy[16] = {};
	UShort_t front_time[16] = {};
	Double_t front_cal[16] = {};
	Int_t    back_num = 0;
	UShort_t back_strip[16] = {};
	UShort_t back_energy[16] = {};
	UShort_t back_time[16] = {};
	Double_t back_cal[16] = {};
};

struct NeOutEvent {
	AdssdHits adssd[kAdssdCount];
	CsiHits   csi;
	DssdHits  dssd;
	// 以下分支按原样保存
	UShort_t  f3[4] = {};    UShort_t tf3[4] = {};
	UShort_t  f5[4] = {};    UShort_t tf5[4] = {};
	UShort_t  tmwpc[4] = {};
	UInt_t    scaler[16] = {};
	ULong64_t scaler_xia = 0;
	Double_t  af3 = 0;
	Double_t  af5 = 0;
	Double_t  tof = 0;
	Double_t  xbt = 0;
	Double_t  ybt = 0;
	Double_t  xbd = 0;
	Double_t  ybd = 0;
	UShort_t  nx1 = 0;
	UShort_t  ny1 = 0;
	UShort_t  nx2 = 0;
	UShort_t  ny2 = 0;
	UShort_t  x1[32] = {};
	UShort_t  y1[32] = {};
	UShort_t  x2[32] = {};
	UShort_t  y2[32] = {};
	UInt_t    rx1 = 0;
	UInt_t    rx2 = 0;
	UInt_t    ry1 = 0;
	UInt_t    ry2 = 0;
	Int_t     nevent = 0;
	Int_t     trigger = 0;
	Int_t     subevents = 0;
	Int_t     evsize = 0;
	UInt_t    mtime = 0;
	UInt_t    utime = 0;
	UInt_t    utime_old = 0;
};

// 原样保存的字段：X(成员名, 分支名)
#define NE_REST_FIELDS(X) \
	X(f3, "F3") X(tf3, "tF3") X(f5, "F5") X(tf5, "tF5") \
	X(tmwpc, "tMWPC") \
	X(scaler, "scaler") X(scaler_xia, "scalerXIA") \
	X(af3, "af3") X(af5, "af5") X(tof, "ToF") \
	X(xbt, "xbt") X(ybt, "ybt") X(xbd, "xbd") X(ybd, "ybd") \
	X(nx1, "nx1") X(ny1, "ny1") X(nx2, "nx2") X(ny2, "ny2") \
	X(x1, "x1") X(y1, "y1") X(x2, "x2") X(y2, "y2") \
	X(rx1, "rx1") X(rx2, "rx2") X(ry1, "ry1") X(ry2, "ry2") \
	X(nevent, "nevent") X(trigger, "trigger") \
	X(subevents, "subevents") X(evsize, "evsize") \
	X(mtime, "mtime") X(utime, "utime") X(utime_old, "utime_old")

// 在树中查找分支：兼容 "NeEvent.X"、"X" 以及带其它前缀的写法
// 先按 "NeEvent.X" 和 "X" 直接查，再遍历所有分支（含 branch element 的子分支）
// 比较名字的最后一段，避免前缀/层级差异导致整份数据读不到
TBranch *FindNeBranch(TTree *tree, const char *suffix) {
	TString full = TString::Format("NeEvent.%s", suffix);
	TBranch *br = tree->GetBranch(full);
	if (!br) br = tree->FindBranch(full);
	if (!br) br = tree->GetBranch(suffix);
	if (!br) br = tree->FindBranch(suffix);
	if (br) return br;

	const TString tail = TString::Format(".%s", suffix);
	TObjArray *top = tree->GetListOfBranches();
	const Int_t ntop = top ? top->GetEntries() : 0;
	for (Int_t i = 0; i < ntop; ++i) {
		TBranch *b = static_cast<TBranch *>(top->At(i));
		if (!b) continue;
		const TString nm = b->GetName();
		if (nm == suffix || nm.EndsWith(tail)) return b;
		TObjArray *sub = b->GetListOfBranches();
		const Int_t nsub = sub ? sub->GetEntries() : 0;
		for (Int_t k = 0; k < nsub; ++k) {
			TBranch *sb = static_cast<TBranch *>(sub->At(k));
			if (!sb) continue;
			const TString snm = sb->GetName();
			if (snm == suffix || snm.EndsWith(tail)) return sb;
		}
	}
	return nullptr;
}

// 绑定输入分支：直接指向输入结构体成员，读入后无需再拷贝
// 返回绑定失败的个数（分支缺失/类型不匹配时会失败）
// 注意：AnalysisxTree 里的 NeEvent 是 split 对象分支，必须先 SetMakeClass(1)
//       把它拆成纯叶子分支，否则对子分支 SetAddress 不会生效（读到的是未初始化内存）
int BindNeEventInput(TTree *tree, NeEventData &ev) {
	tree->SetMakeClass(1);
	tree->SetBranchStatus("*", 1);

	int failed = 0;
#define NE_BIND_IN(field, name) \
	{ \
		TBranch *ne_br = FindNeBranch(tree, name); \
		if (!ne_br) { \
			++failed; \
			std::cerr << "Warning: cannot find branch " name << "\n"; \
		} else if (tree->SetBranchAddress(ne_br->GetName(), \
		           static_cast<void *>(&ev.field)) < 0) { \
			++failed; \
			std::cerr << "Warning: cannot bind " << ne_br->GetName() << "\n"; \
		} \
	}
	NE_EVENT_FIELDS(NE_BIND_IN)
#undef NE_BIND_IN
	return failed;
}

// 建立压缩后的输出分支（数组为变长数组，长度由对应的 *_num 分支决定）
// 命名缩写：a1…a6=ADSSD，d=DSSD，c=CsI；f/b=front/back；n/s/e/t/c=num/strip/energy/time/cal
void SetupNeOutBranches(TTree *tree, NeOutEvent &ev) {
	for (int d = 0; d < kAdssdCount; ++d) {
		AdssdHits &h = ev.adssd[d];
		TString p = TString::Format("a%d", d + 1);
		// 计数分支必须先建，变长数组才有长度可引用
		AddNeBranch(tree, (p + "fn").Data(), h.front_num);
		AddNeBranchV(tree, (p + "fs").Data(), (p + "fn").Data(), h.front_strip);
		AddNeBranchV(tree, (p + "fe").Data(), (p + "fn").Data(), h.front_energy);
		AddNeBranchV(tree, (p + "ft").Data(), (p + "fn").Data(), h.front_time);
		AddNeBranchV(tree, (p + "fc").Data(), (p + "fn").Data(), h.front_cal);
		AddNeBranch(tree, (p + "bn").Data(), h.back_num);
		AddNeBranchV(tree, (p + "bs").Data(), (p + "bn").Data(), h.back_strip);
		AddNeBranchV(tree, (p + "be").Data(), (p + "bn").Data(), h.back_energy);
		AddNeBranchV(tree, (p + "bt").Data(), (p + "bn").Data(), h.back_time);
		AddNeBranchV(tree, (p + "bc").Data(), (p + "bn").Data(), h.back_cal);
	}

	AddNeBranch(tree, "cn", ev.csi.num);
	AddNeBranchV(tree, "cs", "cn", ev.csi.strip);
	AddNeBranchV(tree, "ce", "cn", ev.csi.energy);
	AddNeBranchV(tree, "ct", "cn", ev.csi.time);
	AddNeBranchV(tree, "cc", "cn", ev.csi.cal);

	AddNeBranch(tree, "dfn", ev.dssd.front_num);
	AddNeBranchV(tree, "dfs", "dfn", ev.dssd.front_strip);
	AddNeBranchV(tree, "dfe", "dfn", ev.dssd.front_energy);
	AddNeBranchV(tree, "dft", "dfn", ev.dssd.front_time);
	AddNeBranchV(tree, "dfc", "dfn", ev.dssd.front_cal);
	AddNeBranch(tree, "dbn", ev.dssd.back_num);
	AddNeBranchV(tree, "dbs", "dbn", ev.dssd.back_strip);
	AddNeBranchV(tree, "dbe", "dbn", ev.dssd.back_energy);
	AddNeBranchV(tree, "dbt", "dbn", ev.dssd.back_time);
	AddNeBranchV(tree, "dbc", "dbn", ev.dssd.back_cal);

#define NE_ADD_REST(field, name) AddNeBranch(tree, name, ev.field);
	NE_REST_FIELDS(NE_ADD_REST)
#undef NE_ADD_REST
}

// 命中能量的有效范围（ADC 道数），超出范围的道不记录
constexpr double kNeEnergyMin = 100.0;
constexpr double kNeEnergyMax = 8000.0;

inline bool NeEnergyInRange(double energy) {
	return energy >= kNeEnergyMin && energy <= kNeEnergyMax;
}

// 按能量从大到小插入一个命中（数组容量 N）
// strip / energy / time / cal 一起搬移，保证与排序后的 energy 一一对应
template <size_t N>
void InsertNeHit(
	int &num,
	UShort_t *strip, UShort_t *energy, UShort_t *time, Double_t *cal,
	int s, UShort_t e, UShort_t t, Double_t c
) {
	int pos = 0;
	while (pos < num && energy[pos] >= e) ++pos;
	if (pos >= static_cast<int>(N)) return;
	int limit = std::min<int>(num, static_cast<int>(N) - 1);
	for (int i = limit; i > pos; --i) {
		strip[i] = strip[i-1];
		energy[i] = energy[i-1];
		time[i] = time[i-1];
		cal[i] = cal[i-1];
	}
	strip[pos] = static_cast<UShort_t>(s);
	energy[pos] = e;
	time[pos] = t;
	cal[pos] = c;
	if (num < static_cast<int>(N)) ++num;
}

// 把输入的 NeEvent 压缩成只含命中的结构，其余字段原样拷贝
void FillNeOut(const NeEventData &in, const NeThreshold &thr, NeOutEvent &out) {
	// ADSSD：r1-r6 为 front（16 条），s1-s6 为 back（8 条），6 个分开存储
	const UShort_t *fa[6] = { in.r1a, in.r2a, in.r3a, in.r4a, in.r5a, in.r6a };
	const UShort_t *ft[6] = { in.r1t, in.r2t, in.r3t, in.r4t, in.r5t, in.r6t };
	const Double_t *fc[6] = { in.r1c, in.r2c, in.r3c, in.r4c, in.r5c, in.r6c };
	const UShort_t *ba[6] = { in.s1a, in.s2a, in.s3a, in.s4a, in.s5a, in.s6a };
	const UShort_t *bt[6] = { in.s1t, in.s2t, in.s3t, in.s4t, in.s5t, in.s6t };
	const Double_t *bc[6] = { in.s1c, in.s2c, in.s3c, in.s4c, in.s5c, in.s6c };

	for (int d = 0; d < kAdssdCount; ++d) {
		AdssdHits &h = out.adssd[d];
		h.front_num = 0;
		h.back_num = 0;
		for (int s = 0; s < kAdssdFrontStrips; ++s) {
			if (fa[d][s] > thr.adssd_front[s] && NeEnergyInRange(fa[d][s])) {
				InsertNeHit<kAdssdFrontStrips>(
					h.front_num, h.front_strip, h.front_energy,
					h.front_time, h.front_cal,
					s, fa[d][s], ft[d][s], fc[d][s]
				);
			}
		}
		for (int s = 0; s < kAdssdBackStrips; ++s) {
			if (ba[d][s] > thr.adssd_back[s] && NeEnergyInRange(ba[d][s])) {
				InsertNeHit<kAdssdBackStrips>(
					h.back_num, h.back_strip, h.back_energy,
					h.back_time, h.back_cal,
					s, ba[d][s], bt[d][s], bc[d][s]
				);
			}
		}
	}

	// CsI
	out.csi.num = 0;
	for (int s = 0; s < 16; ++s) {
		if (in.csa[s] > thr.csi[s] && NeEnergyInRange(in.csa[s])) {
			InsertNeHit<16>(
				out.csi.num, out.csi.strip, out.csi.energy,
				out.csi.time, out.csi.cal,
				s, in.csa[s], in.cst[s], in.csc[s]
			);
		}
	}

	// DSSD：x 为 front，y 为 back
	out.dssd.front_num = 0;
	out.dssd.back_num = 0;
	for (int s = 0; s < 16; ++s) {
		if (in.xa[s] > thr.dssd_front[s] && NeEnergyInRange(in.xa[s])) {
			InsertNeHit<16>(
				out.dssd.front_num, out.dssd.front_strip,
				out.dssd.front_energy, out.dssd.front_time,
				out.dssd.front_cal, s, in.xa[s], in.xt[s], in.xc[s]
			);
		}
		if (in.ya[s] > thr.dssd_back[s] && NeEnergyInRange(in.ya[s])) {
			InsertNeHit<16>(
				out.dssd.back_num, out.dssd.back_strip,
				out.dssd.back_energy, out.dssd.back_time,
				out.dssd.back_cal, s, in.ya[s], in.yt[s], in.yc[s]
			);
		}
	}

	// 其余分支原样保存
#define NE_COPY_REST(field, name) std::memcpy(&out.field, &in.field, sizeof(out.field));
	NE_REST_FIELDS(NE_COPY_REST)
#undef NE_COPY_REST
}

// VME 触发有效性：NeEvent.trigger 记录该事例的触发字
// 若 t1 只对应其中某一位，把这里改成 (in.trigger & (1 << kT1Bit)) != 0
inline bool NeTriggerValid(const NeEventData &in) {
	return in.trigger != 0;
}

// 命中阈值（ADC 道数，严格大于才算命中），直接写在程序里
// 需要逐条阈值时，把下面的 fill 替换成逐元素赋值，例如：
//   t.adssd_front = { 50, 50, 60, ... };   // 16 个
constexpr double kNeHitThreshold = 50.0;
// DSSD 单独使用更高阈值
constexpr double kNeDssdThreshold = 250.0;

NeThreshold DefaultNeThreshold() {
	NeThreshold t;
	t.adssd_front.fill(kNeHitThreshold);
	t.adssd_back.fill(kNeHitThreshold);
	t.csi.fill(kNeHitThreshold);
	t.dssd_front.fill(kNeDssdThreshold);
	t.dssd_back.fill(kNeDssdThreshold);
	return t;
}

// ---------------------------------------------------------------------------
// ingot 兼容输出：分支名/类型与 include/event/ingot/*.h 的事件结构一一对应，
//   供 smelt 的 SmeltDetector<T>（内部按成员名绑定分支）直接读取
//   tafd_vme_XXXX.root   : a1-a6 六个 ADSSD   -> TafdEvent
//   t1du_vme_XXXX.root   : dssd 的 x(front)/y(back) -> DssdEvent
//   tafcsi_vme_XXXX.root : CsI 第 0-11 条     -> CsiEvent<12>
//   t1csiu_vme_XXXX.root : CsI 第 12-15 条    -> CsiEvent<4>
// 同时保留新加的刻度分支（front_cal / back_cal / cal）
// ---------------------------------------------------------------------------
constexpr int kDssdOutHits = 32;    // 与 forgerib::kDssdMaxHits 一致
constexpr int kTafdCount = 6;
constexpr int kCsiHighCount = 12;
constexpr int kCsiLowCount = 4;

struct IngotTafdEvent {
	int    flag = 0;
	bool   valid[kTafdCount] = {};
	double front_energy[kTafdCount] = {};
	double back_energy[kTafdCount] = {};
	int    front_strip[kTafdCount] = {};
	int    back_strip[kTafdCount] = {};
	double front_time[kTafdCount] = {};
	double front_cal[kTafdCount] = {};
	double back_cal[kTafdCount] = {};
};

struct IngotDssdEvent {
	int    front_num = 0;
	int    front_strip[kDssdOutHits] = {};
	double front_energy[kDssdOutHits] = {};
	double front_time[kDssdOutHits] = {};
	double front_cal[kDssdOutHits] = {};
	int    back_num = 0;
	int    back_strip[kDssdOutHits] = {};
	double back_energy[kDssdOutHits] = {};
	double back_time[kDssdOutHits] = {};
	double back_cal[kDssdOutHits] = {};
};

template <int N>
struct IngotCsiEvent {
	ULong64_t flag = 0;
	bool   valid[N] = {};
	double time[N] = {};
	int    energy[N] = {};
	double cal[N] = {};
};

void SetupTafdOut(TTree *tree, IngotTafdEvent &ev) {
	tree->Branch("flag", &ev.flag, "flag/I");
	tree->Branch("valid", ev.valid, "v[6]/O");
	tree->Branch("front_energy", ev.front_energy, "fe[6]/D");
	tree->Branch("back_energy", ev.back_energy, "be[6]/D");
	tree->Branch("front_strip", ev.front_strip, "fs[6]/I");
	tree->Branch("back_strip", ev.back_strip, "bs[6]/I");
	tree->Branch("front_time", ev.front_time, "ft[6]/D");
	tree->Branch("front_cal", ev.front_cal, "fc[6]/D");
	tree->Branch("back_cal", ev.back_cal, "bc[6]/D");
}

// 与 forgerib::DssdEvent::SetupOutput 完全一致：前/背均为变长数组
void SetupDssdOut(TTree *tree, IngotDssdEvent &ev) {
	tree->Branch("front_num", &ev.front_num, "fn/I");
	tree->Branch("front_strip", ev.front_strip, "fs[fn]/I");
	tree->Branch("front_energy", ev.front_energy, "fe[fn]/D");
	tree->Branch("front_time", ev.front_time, "ft[fn]/D");
	tree->Branch("front_cal", ev.front_cal, "fc[fn]/D");
	tree->Branch("back_num", &ev.back_num, "bn/I");
	tree->Branch("back_strip", ev.back_strip, "bs[bn]/I");
	tree->Branch("back_energy", ev.back_energy, "be[bn]/D");
	tree->Branch("back_time", ev.back_time, "bt[bn]/D");
	tree->Branch("back_cal", ev.back_cal, "bc[bn]/D");
}

template <int N>
void SetupCsiOut(TTree *tree, IngotCsiEvent<N> &ev) {
	tree->Branch("flag", &ev.flag, "flag/l");
	tree->Branch("valid", ev.valid, TString::Format("v[%d]/O", N).Data());
	tree->Branch("time", ev.time, TString::Format("t[%d]/D", N).Data());
	tree->Branch("energy", ev.energy, TString::Format("e[%d]/I", N).Data());
	tree->Branch("cal", ev.cal, TString::Format("c[%d]/D", N).Data());
}

// 每个 ADSSD 只取能量最大的命中（命中已按能量降序），对应 TafdEvent 的单命中语义
void FillTafdOut(const NeOutEvent &in, IngotTafdEvent &ev) {
	ev.flag = 0;
	for (int d = 0; d < kTafdCount; ++d) {
		ev.valid[d] = false;
		ev.front_energy[d] = 0;
		ev.back_energy[d] = 0;
		ev.front_strip[d] = -1;
		ev.back_strip[d] = -1;
		ev.front_time[d] = 0;
		ev.front_cal[d] = 0;
		ev.back_cal[d] = 0;
		const AdssdHits &h = in.adssd[d];
		if (h.front_num > 0) {
			ev.front_energy[d] = h.front_energy[0];
			ev.front_strip[d] = h.front_strip[0];
			ev.front_time[d] = h.front_time[0];
			ev.front_cal[d] = h.front_cal[0];
			ev.valid[d] = true;
			ev.flag |= (1 << d);
		}
		if (h.back_num > 0) {
			ev.back_energy[d] = h.back_energy[0];
			ev.back_strip[d] = h.back_strip[0];
			ev.back_cal[d] = h.back_cal[0];
			ev.valid[d] = true;
			ev.flag |= (1 << d);
		}
	}
}

void FillDssdOut(const NeOutEvent &in, IngotDssdEvent &ev) {
	const DssdHits &h = in.dssd;
	ev.front_num = std::min(h.front_num, kDssdOutHits);
	for (int i = 0; i < ev.front_num; ++i) {
		ev.front_strip[i] = h.front_strip[i];
		ev.front_energy[i] = h.front_energy[i];
		ev.front_time[i] = h.front_time[i];
		ev.front_cal[i] = h.front_cal[i];
	}
	ev.back_num = std::min(h.back_num, kDssdOutHits);
	for (int i = 0; i < ev.back_num; ++i) {
		ev.back_strip[i] = h.back_strip[i];
		ev.back_energy[i] = h.back_energy[i];
		ev.back_time[i] = h.back_time[i];
		ev.back_cal[i] = h.back_cal[i];
	}
}

// first 为这段 CsI 的起始条号（0-11 或 12-15）
template <int N>
void FillCsiOut(const NeOutEvent &in, int first, IngotCsiEvent<N> &ev) {
	ev.flag = 0;
	for (int i = 0; i < N; ++i) {
		ev.valid[i] = false;
		ev.time[i] = 0;
		ev.energy[i] = 0;
		ev.cal[i] = 0;
	}
	for (int i = 0; i < in.csi.num; ++i) {
		const int s = in.csi.strip[i];
		if (s < first || s >= first + N) continue;
		const int k = s - first;
		ev.valid[k] = true;
		ev.time[k] = in.csi.time[i];
		ev.energy[k] = in.csi.energy[i];
		ev.cal[k] = in.csi.cal[i];
		ev.flag |= (1ull << k);
	}
}

// 四个 ingot 兼容文件的集合，与 vme_XXXX.root 同步逐事例填充
struct IngotOutputs {
	TFile *tafd_file = nullptr;
	TFile *t1du_file = nullptr;
	TFile *tafcsi_file = nullptr;
	TFile *t1csiu_file = nullptr;
	TTree *tafd_tree = nullptr;
	TTree *t1du_tree = nullptr;
	TTree *tafcsi_tree = nullptr;
	TTree *t1csiu_tree = nullptr;
	IngotTafdEvent tafd;
	IngotDssdEvent t1du;
	IngotCsiEvent<kCsiHighCount> tafcsi;
	IngotCsiEvent<kCsiLowCount> t1csiu;

	bool Open(const std::string &dir, int run) {
		tafd_file = new TFile(
			TString::Format("%s/tafd_vme_%04d.root", dir.c_str(), run),
			"recreate"
		);
		t1du_file = new TFile(
			TString::Format("%s/t1du_vme_%04d.root", dir.c_str(), run),
			"recreate"
		);
		tafcsi_file = new TFile(
			TString::Format("%s/tafcsi_vme_%04d.root", dir.c_str(), run),
			"recreate"
		);
		t1csiu_file = new TFile(
			TString::Format("%s/t1csiu_vme_%04d.root", dir.c_str(), run),
			"recreate"
		);
		if (!tafd_file->IsOpen() || !t1du_file->IsOpen()
			|| !tafcsi_file->IsOpen() || !t1csiu_file->IsOpen()) {
			std::cerr << "Error: cannot create ingot vme files in "
			          << dir << "\n";
			return false;
		}
		tafd_file->cd();
		tafd_tree = new TTree("tree", "tafd (a1-a6) from vme");
		SetupTafdOut(tafd_tree, tafd);
		t1du_file->cd();
		t1du_tree = new TTree("tree", "t1du dssd from vme");
		SetupDssdOut(t1du_tree, t1du);
		tafcsi_file->cd();
		tafcsi_tree = new TTree("tree", "tafcsi (csi 0-11) from vme");
		SetupCsiOut(tafcsi_tree, tafcsi);
		t1csiu_file->cd();
		t1csiu_tree = new TTree("tree", "t1csiu (csi 12-15) from vme");
		SetupCsiOut(t1csiu_tree, t1csiu);
		std::cout << "Output: " << tafd_file->GetName() << "\n"
		          << "Output: " << t1du_file->GetName() << "\n"
		          << "Output: " << tafcsi_file->GetName() << "\n"
		          << "Output: " << t1csiu_file->GetName() << "\n";
		return true;
	}

	void Fill(const NeOutEvent &in) {
		FillTafdOut(in, tafd);
		FillDssdOut(in, t1du);
		FillCsiOut(in, 0, tafcsi);
		FillCsiOut(in, kCsiHighCount, t1csiu);
		tafd_tree->Fill();
		t1du_tree->Fill();
		tafcsi_tree->Fill();
		t1csiu_tree->Fill();
	}

	void Write() {
		tafd_file->cd();
		tafd_tree->Write();
		t1du_file->cd();
		t1du_tree->Write();
		tafcsi_file->cd();
		tafcsi_tree->Write();
		t1csiu_file->cd();
		t1csiu_tree->Write();
	}

	void Close() {
		tafd_file->Close();
		t1du_file->Close();
		tafcsi_file->Close();
		t1csiu_file->Close();
	}
};

// 把 runXX_*.root 中的 AnalysisxTree 转换压缩结构后写入 vme_XXXX.root
int WriteVmeTree(
	const std::vector<std::string> &input_files,
	const std::string &output_dir,
	int vme_run
) {
	TString vme_path = TString::Format(
		"%s/vme_%04d.root",
		output_dir.c_str(),
		vme_run
	);
	std::cout << "Output: " << vme_path << "\n";

	TFile vme_file(vme_path, "recreate");
	if (!vme_file.IsOpen()) {
		std::cerr << "Error: Cannot create " << vme_path << ".\n";
		return -1;
	}
	TTree vme_tree("tree", "vme converted from Go4 AnalysisxTree");
	NeEventData ne_event;
	NeOutEvent ne_out;
	NeThreshold ne_threshold = DefaultNeThreshold();
	SetupNeOutBranches(&vme_tree, ne_out);
	vme_file.cd();

	// VME t1 触发文件：每条事例一个标量 valid，顺序与 vme_XXXX.root 完全一致
	// coke 用它决定哪些 VME entry 写入 trigger_t1_XXXX.root
	TString t1_path = TString::Format(
		"%s/trigger_vme_t1_%04d.root",
		output_dir.c_str(),
		vme_run
	);
	TFile t1_file(t1_path, "recreate");
	if (!t1_file.IsOpen()) {
		std::cerr << "Error: Cannot create " << t1_path << ".\n";
		vme_file.Close();
		return -1;
	}
	bool t1_valid = false;
	t1_file.cd();
	TTree t1_tree("tree", "vme t1 trigger");
	t1_tree.Branch("valid", &t1_valid, "valid/O");
	vme_file.cd();
	std::cout << "Output: " << t1_path << "\n";

	// 与 vme_XXXX.root 同步输出的 ingot 兼容文件（tafd / t1du / tafcsi / t1csiu）
	IngotOutputs ingot;
	if (!ingot.Open(output_dir, vme_run)) {
		vme_file.Close();
		return -1;
	}

	long long converted = 0;
	long long t1_valid_count = 0;
	for (size_t i = 0; i < input_files.size(); ++i) {
		TFile input_file(input_files[i].c_str(), "read");
		TTree *input_tree =
			static_cast<TTree *>(input_file.Get("AnalysisxTree"));
		if (!input_tree) {
			std::cerr << "Error: Get AnalysisxTree from "
			          << input_files[i] << " failed.\n";
			return -1;
		}
		const int nfail = BindNeEventInput(input_tree, ne_event);
		if (nfail > 0) {
			std::cerr << "Warning: " << nfail
			          << " input branch(es) failed to bind in "
			          << input_files[i] << "\n";
		}
		if (i == 0) {
			TObjArray *all = input_tree->GetListOfBranches();
			const Int_t nall = all ? all->GetEntries() : 0;
			std::cout << "Input branches (" << nall
			          << " top-level, MakeClass=1), first names:";
			for (Int_t k = 0; k < nall && k < 8; ++k) {
				TBranch *b = static_cast<TBranch *>(all->At(k));
				if (b) std::cout << " " << b->GetName();
			}
			std::cout << "\n";
		}

		const Long64_t total = input_tree->GetEntries();
		// 自检：打印第一个 entry 的原始值，确认输入确实读到位
		if (total > 0) {
			input_tree->GetEntry(0);
			std::cout << "Check " << input_files[i]
			          << " entry0: nevent=" << ne_event.nevent
			          << " R1a[0]=" << ne_event.r1a[0]
			          << " R1a[1]=" << ne_event.r1a[1]
			          << " S1a[0]=" << ne_event.s1a[0]
			          << " Xa[0]=" << ne_event.xa[0]
			          << " scalerXIA=" << ne_event.scaler_xia << "\n";
		}
		printf("Converting %s   0%%", input_files[i].c_str());
		fflush(stdout);
		long long last_percentage = 0;
		for (Long64_t entry = 0; entry < total; ++entry) {
			if (total > 0
				&& entry * 100ll / total > last_percentage) {
				last_percentage = entry * 100ll / total;
				printf("\b\b\b\b%3lld%%", last_percentage);
				fflush(stdout);
			}
			input_tree->GetEntry(entry);
			FillNeOut(ne_event, ne_threshold, ne_out);
			vme_tree.Fill();
			ingot.Fill(ne_out);
			// 与 vme_XXXX.root 同序写一条 valid，coke 按 entry 号索引
			t1_valid = NeTriggerValid(ne_event);
			if (t1_valid) ++t1_valid_count;
			t1_tree.Fill();
		}
		printf("\b\b\b\b100%%\n");
		converted += total;

		input_file.Close();
		vme_file.cd();
	}
	std::cout << "Converted entries: " << converted << "\n";
	std::cout << "VME t1 valid entries: " << t1_valid_count
	          << " / " << converted << "\n";
	if (converted > 0 && t1_valid_count == 0) {
		std::cerr << "Warning: all VME events have trigger == 0, "
		          << "trigger_vme_t1 is empty (t1 selection drops everything).\n";
	}

	ingot.Write();
	ingot.Close();

	t1_file.cd();
	t1_tree.Write();
	t1_file.Close();

	vme_file.cd();
	vme_tree.Write();
	vme_file.Close();

	// 自检：重新打开输出文件，打印前几个 entry，确认写入的数据随事例变化
	{
		TFile check_file(vme_path, "read");
		TTree *check_tree = static_cast<TTree *>(check_file.Get("tree"));
		if (!check_tree) {
			std::cerr << "Warning: cannot reopen " << vme_path << " for check.\n";
			return 0;
		}
		Int_t a1fn = 0, dfn = 0, nevent = 0;
		UShort_t a1fe[kAdssdFrontStrips] = {};
		UShort_t a1fs[kAdssdFrontStrips] = {};
		check_tree->SetBranchAddress("a1fn", &a1fn);
		check_tree->SetBranchAddress("a1fe", a1fe);
		check_tree->SetBranchAddress("a1fs", a1fs);
		check_tree->SetBranchAddress("dfn", &dfn);
		check_tree->SetBranchAddress("nevent", &nevent);
		const Long64_t nshow =
			std::min<Long64_t>(check_tree->GetEntries(), 3);
		for (Long64_t k = 0; k < nshow; ++k) {
			check_tree->GetEntry(k);
			std::cout << "Check output entry " << k
			          << ": nevent=" << nevent
			          << " a1fn=" << a1fn
			          << " a1fs[0]=" << (a1fn > 0 ? a1fs[0] : 0)
			          << " a1fe[0]=" << (a1fn > 0 ? a1fe[0] : 0)
			          << " dfn=" << dfn << "\n";
		}
	}
	return 0;
}

// 字段表宏只在本文件使用，用完取消定义，避免污染其它编译单元
#undef NE_EVENT_FIELDS
#undef NE_REST_FIELDS

} // namespace

int main(int argc, char **argv) {
    cxxopts::Options options(
        "crush_vme",
        "convert Go4 AnalysisxTree of all runXX_*.root into vme_XXXX.root "
        "(plus tafd/t1du/tafcsi/t1csiu_vme_XXXX.root for smelt)"
    );
    options.add_options()
        ("h,help", "Print usage")
        (
            "run",
            "VME run number (matches runXX_*.root, used in output file name).",
            cxxopts::value<int>(),
            "run"
        )
        (
            "c,config",
            "Config file path.",
            cxxopts::value<std::string>()->default_value("config.toml"),
            "path"
        )
        (
            "threshold",
            "(unused) Hit thresholds are defined in DefaultNeThreshold().",
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

    // ---- 1. 收集所有 runXX_*.root ----
    char run_tag[8];
    std::snprintf(run_tag, sizeof(run_tag), "%02d", vme_run);
    std::cout << "Scanning VME raw dir: " << config.paths.raw_vme
              << " for run" << run_tag << "_*.root ...\n";

    std::vector<std::string> input_files =
        CollectRunFiles(config.paths.raw_vme, vme_run);
    if (input_files.empty()) {
        std::cerr << "Error: No *.root matching run" << run_tag
                  << "_* in " << config.paths.raw_vme << ".\n";
        return -1;
    }
    std::cout << "Matched " << input_files.size() << " file(s):\n";
    for (auto &f : input_files) std::cout << "  " << f << "\n";

    // ---- 2. 将 AnalysisxTree 转换结构后写入 vme_XXXX.root ----
    //   trigger 文件由单独的匹配程序生成（trigger_XXXX.root），这里不写
    if (WriteVmeTree(input_files, output_dir, vme_run) != 0) {
        return -1;
    }

    return 0;
}