#include <algorithm>
#include <array>
#include <cctype>
#include <dirent.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
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

// 读取一个 Go4 AnalysisxTree 的 NeEvent.scalerXIA，追加到 times。
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

} // namespace

int main(int argc, char **argv) {
    cxxopts::Options options(
        "crush_vme",
        "merge Go4 scalerXIA of all runXX_*.root into trigger_vme_XXXX.root"
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

    // ---- 2. 逐个读取 scalerXIA 并合并 ----
    std::vector<long long> trigger_times;
    for (size_t i = 0; i < input_files.size(); ++i) {
        std::cout << "[" << i+1 << "/" << input_files.size() << "] ";
        if (ReadScalerXia(input_files[i], trigger_times, true)) {
            return -1;
        }
    }
    std::cout << "Total entries: " << trigger_times.size() << "\n";

    // ---- 3. 写 trigger_vme_XXXX.root ----
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
    return 0;
}