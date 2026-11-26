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

// ---------------- 文件收集 ----------------
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

bool FileMatchesRun(const std::string &basename, const std::string &run_tag) {
    std::string stem = basename;
    if (stem.size() > 5 && stem.compare(stem.size() - 5, 5, ".root") == 0) {
        stem = stem.substr(0, stem.size() - 5);
    }
    const std::string prefix = "run" + run_tag + "_";
    if (stem.compare(0, prefix.size(), prefix) == 0) return true;
    const std::string n1 = "_" + run_tag + "_";
    const std::string n2 = "-" + run_tag + "-";
    if (stem.find(n1) != std::string::npos) return true;
    if (stem.find(n2) != std::string::npos) return true;
    return false;
}

std::vector<std::string> CollectRunFiles(const std::string &dir, int vme_run) {
    char run_tag[8];
    std::snprintf(run_tag, sizeof(run_tag), "%02d", vme_run);
    std::vector<std::string> all = ListRootFiles(dir);
    std::vector<std::string> matched;
    matched.reserve(all.size());
    for (auto &f : all) {
        std::string base = f.substr(f.find_last_of('/') + 1);
        if (FileMatchesRun(base, run_tag)) matched.push_back(f);
    }
    return matched;
}

// ---------------- 读 scalerXIA（原逻辑） ----------------
int ReadScalerXia(const std::string &path,
                  std::vector<long long> &times,
                  bool report = false) {
    TFile input_file(path.c_str(), "read");
    if (!input_file.IsOpen()) {
        std::cerr << "Error: Cannot open " << path << ".\n";
        return -1;
    }
    TTree *tree = static_cast<TTree *>(input_file.Get("AnalysisxTree"));
    if (!tree) {
        std::cerr << "Error: Get AnalysisxTree from " << path << " failed.\n";
        return -1;
    }
    TTreeReader reader(tree);
    TTreeReaderValue<ULong64_t> scaler_xia(reader, "NeEvent.scalerXIA");
    const Long64_t total = tree->GetEntries();
    if (report) {
        printf("Reading %s   0%%", path.c_str());
        fflush(stdout);
    }
    long long last_percentage = 0;
    for (Long64_t entry = 0; entry < total; ++entry) {
        if (report && total > 0
            && entry * 100ll / total > last_percentage) {
            last_percentage = entry * 100ll / total;
            printf("\b\b\b\b%3lld%%", last_percentage);
            fflush(stdout);
        }
        reader.Next();
        if (reader.GetEntryStatus() != TTreeReader::kEntryValid) {
            std::cerr << "\nError: TTreeReader failed at entry "
                      << entry << " in " << path << ".\n";
            return -1;
        }
        times.push_back(static_cast<long long>(*scaler_xia));
    }
    if (report) printf("\b\b\b\b100%%\n");
    input_file.Close();
    return 0;
}

// ---------------- 新增：NeEvent 输入缓冲 ----------------
struct NeEventData {
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
    UShort_t csa[16];  UShort_t cst[16];
    UShort_t xa[16];   UShort_t xt[16];
    UShort_t ya[16];   UShort_t yt[16];
    Double_t s1c[8];   Double_t r1c[16];
    Double_t s2c[8];   Double_t r2c[16];
    Double_t s3c[8];   Double_t r3c[16];
    Double_t s4c[8];   Double_t r4c[16];
    Double_t s5c[8];   Double_t r5c[16];
    Double_t s6c[8];   Double_t r6c[16];
    Double_t csc[16];  Double_t xc[16];  Double_t yc[16];
    Int_t    trigger;
};

#define NE_EVENT_FIELDS(X) \
    X(s1a, "S1a") X(s1t, "S1t") X(r1a, "R1a") X(r1t, "R1t") \
    X(s2a, "S2a") X(s2t, "S2t") X(r2a, "R2a") X(r2t, "R2t") \
    X(s3a, "S3a") X(s3t, "S3t") X(r3a, "R3a") X(r3t, "R3t") \
    X(s4a, "S4a") X(s4t, "S4t") X(r4a, "R4a") X(r4t, "R4t") \
    X(s5a, "S5a") X(s5t, "S5t") X(r5a, "R5a") X(r5t, "R5t") \
    X(s6a, "S6a") X(s6t, "S6t") X(r6a, "R6a") X(r6t, "R6t") \
    X(csa, "Csa") X(cst, "Cst") X(xa, "Xa") X(xt, "Xt") \
    X(ya, "Ya") X(yt, "Yt") \
    X(s1c, "S1c") X(r1c, "R1c") X(s2c, "S2c") X(r2c, "R2c") \
    X(s3c, "S3c") X(r3c, "R3c") X(s4c, "S4c") X(r4c, "R4c") \
    X(s5c, "S5c") X(r5c, "R5c") X(s6c, "S6c") X(r6c, "R6c") \
    X(csc, "Csc") X(xc, "Xc") X(yc, "Yc") \
    X(trigger, "trigger")

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

inline bool NeTriggerValid(const NeEventData &in) {
    return in.trigger != 0;
}

// ---------------- 新增：按探测器拆分的 _t1.root 输出 ----------------
constexpr int kAdssdCount = 6;
constexpr int kAdssdFrontStrips = 16;
constexpr int kAdssdBackStrips = 8;
constexpr int kDssdOutHits = 32;
constexpr int kCsiHighCount = 12;
constexpr int kCsiLowCount = 4;

struct TafdT1Event {
    int    flag = 0;
    bool   valid[kAdssdCount] = {};
    double front_energy[kAdssdCount] = {};
    double back_energy[kAdssdCount] = {};
    int    front_strip[kAdssdCount] = {};
    int    back_strip[kAdssdCount] = {};
    double front_time[kAdssdCount] = {};
    double front_cal[kAdssdCount] = {};
    double back_cal[kAdssdCount] = {};
};

struct DssdT1Event {
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
struct CsiT1Event {
    ULong64_t flag = 0;
    bool   valid[N] = {};
    double time[N] = {};
    int    energy[N] = {};
    double cal[N] = {};
};

void SetupTafdT1(TTree *tree, TafdT1Event &ev) {
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

void SetupDssdT1(TTree *tree, DssdT1Event &ev) {
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
void SetupCsiT1(TTree *tree, CsiT1Event<N> &ev) {
    tree->Branch("flag", &ev.flag, "flag/l");
    tree->Branch("valid", ev.valid, TString::Format("v[%d]/O", N).Data());
    tree->Branch("time", ev.time, TString::Format("t[%d]/D", N).Data());
    tree->Branch("energy", ev.energy, TString::Format("e[%d]/I", N).Data());
    tree->Branch("cal", ev.cal, TString::Format("c[%d]/D", N).Data());
}

// 每个 ADSSD 取能量最大的 front/back 命中（按 strip 直接扫描，取最大）
void FillTafdT1(const NeEventData &in, TafdT1Event &ev) {
    const UShort_t *fa[6] = { in.r1a, in.r2a, in.r3a, in.r4a, in.r5a, in.r6a };
    const UShort_t *ft[6] = { in.r1t, in.r2t, in.r3t, in.r4t, in.r5t, in.r6t };
    const Double_t *fc[6] = { in.r1c, in.r2c, in.r3c, in.r4c, in.r5c, in.r6c };
    const UShort_t *ba[6] = { in.s1a, in.s2a, in.s3a, in.s4a, in.s5a, in.s6a };
    const Double_t *bc[6] = { in.s1c, in.s2c, in.s3c, in.s4c, in.s5c, in.s6c };

    ev.flag = 0;
    for (int d = 0; d < kAdssdCount; ++d) {
        ev.valid[d] = false;
        ev.front_energy[d] = 0;
        ev.back_energy[d] = 0;
        ev.front_strip[d] = -1;
        ev.back_strip[d] = -1;
        ev.front_time[d] = 0;
        ev.front_cal[d] = 0;
        ev.back_cal[d] = 0;

        int best_s = -1;
        UShort_t best_e = 0;
        for (int s = 0; s < kAdssdFrontStrips; ++s) {
            if (fa[d][s] > best_e) { best_e = fa[d][s]; best_s = s; }
        }
        if (best_s >= 0) {
            ev.front_energy[d] = best_e;
            ev.front_strip[d] = best_s;
            ev.front_time[d] = ft[d][best_s];
            ev.front_cal[d] = fc[d][best_s];
            ev.valid[d] = true;
            ev.flag |= (1 << d);
        }

        best_s = -1;
        best_e = 0;
        for (int s = 0; s < kAdssdBackStrips; ++s) {
            if (ba[d][s] > best_e) { best_e = ba[d][s]; best_s = s; }
        }
        if (best_s >= 0) {
            ev.back_energy[d] = best_e;
            ev.back_strip[d] = best_s;
            ev.back_cal[d] = bc[d][best_s];
            ev.valid[d] = true;
            ev.flag |= (1 << d);
        }
    }
}

// DSSD：x 为 front，y 为 back，保留全部命中
void FillDssdT1(const NeEventData &in, DssdT1Event &ev) {
    ev.front_num = 0;
    for (int s = 0; s < 16; ++s) {
        if (in.xa[s] == 0) continue;
        if (ev.front_num >= kDssdOutHits) break;
        ev.front_strip[ev.front_num] = s;
        ev.front_energy[ev.front_num] = in.xa[s];
        ev.front_time[ev.front_num] = in.xt[s];
        ev.front_cal[ev.front_num] = in.xc[s];
        ++ev.front_num;
    }
    ev.back_num = 0;
    for (int s = 0; s < 16; ++s) {
        if (in.ya[s] == 0) continue;
        if (ev.back_num >= kDssdOutHits) break;
        ev.back_strip[ev.back_num] = s;
        ev.back_energy[ev.back_num] = in.ya[s];
        ev.back_time[ev.back_num] = in.yt[s];
        ev.back_cal[ev.back_num] = in.yc[s];
        ++ev.back_num;
    }
}

// CsI：first 为起始条号（0-11 或 12-15）
template <int N>
void FillCsiT1(const NeEventData &in, int first, CsiT1Event<N> &ev) {
    ev.flag = 0;
    for (int i = 0; i < N; ++i) {
        ev.valid[i] = false;
        ev.time[i] = 0;
        ev.energy[i] = 0;
        ev.cal[i] = 0;
    }
    for (int i = 0; i < N; ++i) {
        const int s = first + i;
        if (s >= 16) break;
        if (in.csa[s] == 0) continue;
        ev.valid[i] = true;
        ev.time[i] = in.cst[s];
        ev.energy[i] = in.csa[s];
        ev.cal[i] = in.csc[s];
        ev.flag |= (1ull << i);
    }
}

// 四个 _t1.root 的集合
struct T1Outputs {
    TFile *tafd_file = nullptr;
    TFile *t1du_file = nullptr;
    TFile *tafcsi_file = nullptr;
    TFile *t1csiu_file = nullptr;
    TTree *tafd_tree = nullptr;
    TTree *t1du_tree = nullptr;
    TTree *tafcsi_tree = nullptr;
    TTree *t1csiu_tree = nullptr;
    TafdT1Event tafd;
    DssdT1Event t1du;
    CsiT1Event<kCsiHighCount> tafcsi;
    CsiT1Event<kCsiLowCount> t1csiu;

    bool Open(const std::string &dir, int run) {
        tafd_file = new TFile(TString::Format(
            "%s/tafd_t1_%04d.root", dir.c_str(), run), "recreate");
        t1du_file = new TFile(TString::Format(
            "%s/t1du_t1_%04d.root", dir.c_str(), run), "recreate");
        tafcsi_file = new TFile(TString::Format(
            "%s/tafcsi_t1_%04d.root", dir.c_str(), run), "recreate");
        t1csiu_file = new TFile(TString::Format(
            "%s/t1csiu_t1_%04d.root", dir.c_str(), run), "recreate");
        if (!tafd_file->IsOpen() || !t1du_file->IsOpen()
            || !tafcsi_file->IsOpen() || !t1csiu_file->IsOpen()) {
            std::cerr << "Error: cannot create t1 files in " << dir << "\n";
            return false;
        }
        tafd_file->cd();
        tafd_tree = new TTree("tree", "tafd (a1-a6) from vme t1");
        SetupTafdT1(tafd_tree, tafd);
        t1du_file->cd();
        t1du_tree = new TTree("tree", "t1du dssd from vme t1");
        SetupDssdT1(t1du_tree, t1du);
        tafcsi_file->cd();
        tafcsi_tree = new TTree("tree", "tafcsi (csi 0-11) from vme t1");
        SetupCsiT1(tafcsi_tree, tafcsi);
        t1csiu_file->cd();
        t1csiu_tree = new TTree("tree", "t1csiu (csi 12-15) from vme t1");
        SetupCsiT1(t1csiu_tree, t1csiu);
        std::cout << "Output: " << tafd_file->GetName() << "\n"
                  << "Output: " << t1du_file->GetName() << "\n"
                  << "Output: " << tafcsi_file->GetName() << "\n"
                  << "Output: " << t1csiu_file->GetName() << "\n";
        return true;
    }

    void Fill(const NeEventData &in) {
        FillTafdT1(in, tafd);
        FillDssdT1(in, t1du);
        FillCsiT1(in, 0, tafcsi);
        FillCsiT1(in, kCsiHighCount, t1csiu);
        tafd_tree->Fill();
        t1du_tree->Fill();
        tafcsi_tree->Fill();
        t1csiu_tree->Fill();
    }

    void Write() {
        tafd_file->cd();   tafd_tree->Write();
        t1du_file->cd();   t1du_tree->Write();
        tafcsi_file->cd(); tafcsi_tree->Write();
        t1csiu_file->cd(); t1csiu_tree->Write();
    }

    void Close() {
        tafd_file->Close();
        t1du_file->Close();
        tafcsi_file->Close();
        t1csiu_file->Close();
    }
};

} // namespace

int main(int argc, char **argv) {
    cxxopts::Options options(
        "crush_vme",
        "merge Go4 scalerXIA of all runXX_*.root into trigger_vme_XXXX.root, "
        "and write tafd/t1du/tafcsi/t1csiu_t1_XXXX.root for t1 events"
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
            "(unused) Threshold file path, kept for interface compatibility.",
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

    // ---- 2. 写 trigger_vme_XXXX.root（原逻辑，保持不变）----
    std::vector<long long> trigger_times;
    for (size_t i = 0; i < input_files.size(); ++i) {
        std::cout << "[" << i+1 << "/" << input_files.size() << "] ";
        if (ReadScalerXia(input_files[i], trigger_times, true)) {
            return -1;
        }
    }
    std::cout << "Total entries: " << trigger_times.size() << "\n";

    TString trigger_path = TString::Format(
        "%s/trigger_vme_%04d.root",
        output_dir.c_str(),
        vme_run
    );
    std::cout << "Output: " << trigger_path << "\n";
    TFile trigger_file(trigger_path, "recreate");
    if (!trigger_file.IsOpen()) {
        std::cerr << "Error: Cannot create " << trigger_path << ".\n";
        return -1;
    }
    TTree trigger_tree("tree", "trigger vme");
    bool trigger_valid = false;
    long long trigger_time = 0;
    trigger_tree.Branch("valid", &trigger_valid, "v/O");
    trigger_tree.Branch("time", &trigger_time, "time/L");
    for (long long t : trigger_times) {
        trigger_valid = (t != 0);
        trigger_time = t;
        trigger_tree.Fill();
    }
    trigger_file.cd();
    trigger_tree.Write();
    trigger_file.Close();

    // ---- 3. 新增：写 tafd/t1du/tafcsi/t1csiu_t1_XXXX.root ----
    T1Outputs t1;
    if (!t1.Open(output_dir, vme_run)) {
        return -1;
    }

    long long total_entries = 0;
    long long t1_valid_count = 0;
    NeEventData ne_event;

    for (size_t i = 0; i < input_files.size(); ++i) {
        TFile input_file(input_files[i].c_str(), "read");
        TTree *input_tree =
            static_cast<TTree *>(input_file.Get("AnalysisxTree"));
        if (!input_tree) {
            std::cerr << "Error: Get AnalysisxTree from "
                      << input_files[i] << " failed.\n";
            t1.Close();
            return -1;
        }
        const int nfail = BindNeEventInput(input_tree, ne_event);
        if (nfail > 0) {
            std::cerr << "Warning: " << nfail
                      << " input branch(es) failed to bind in "
                      << input_files[i] << "\n";
        }

        const Long64_t total = input_tree->GetEntries();
        printf("Writing t1 from %s   0%%", input_files[i].c_str());
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
            ++total_entries;
            if (!NeTriggerValid(ne_event)) continue;
            ++t1_valid_count;
            t1.Fill(ne_event);
        }
        printf("\b\b\b\b100%%\n");
        input_file.Close();
    }

    t1.Write();
    t1.Close();

    std::cout << "Total VME entries: " << total_entries << "\n";
    std::cout << "t1 valid entries: " << t1_valid_count
              << " / " << total_entries << "\n";
    if (total_entries > 0 && t1_valid_count == 0) {
        std::cerr << "Warning: all VME events have trigger == 0, "
                  << "t1 output is empty.\n";
    }
    return 0;
}