#include "newcode/pipeline_runner.hpp"
#include "newcode/params.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace newcode;

static void print_usage(const char* prog) {
  std::cout
    << "Usage: " << prog << " [--ebn0 2,4,6,8]\n"
    << "  Run the oFEC pipeline without interleaving (plain).\n"
    << "  --ebn0   Comma-separated Eb/N0 list in dB (default: 2,4,6,8)\n";
}

static std::vector<float> parse_ebn0_list(const std::string& s) {
  std::vector<float> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (!tok.empty()) out.push_back(std::strtof(tok.c_str(), nullptr));
  }
  return out;
}

int main(int argc, char** argv) {
  // 默认的 Eb/N0 列表
  std::vector<float> ebn0_list = {2.f, 4.f, 6.f, 8.f};

  // 解析命令行
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      print_usage(argv[0]);
      return 0;
    } else if (arg == "--ebn0" && i + 1 < argc) {
      ebn0_list = parse_ebn0_list(argv[++i]);
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      print_usage(argv[0]);
      return 1;
    }
  }

  // 运行流水线：不交织（用标签区分场景）
  Params p; // 使用默认参数（块尺寸、调制、帧长等都走默认）
  for (float eb : ebn0_list) {
    std::cout << "\n=== [plain] Eb/N0 = " << eb << " dB ===\n";
    // run_pipeline 内部会完成编码/调制/加噪/LLR/译码并打印 BER（若你在 pipeline 中实现了打印）
    // 第2个参数是场景标签，便于统计/输出区分
    (void)run_pipeline(p, /*label=*/"plain", /*ebn0_db=*/eb);
  }
  return 0;
}
