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
    << "  Run the oFEC pipeline WITH interleaving (ofec interleaver).\n"
    << "  --ebn0       Comma-separated Eb/N0 list in dB (default: 2,4,6,8)\n"
    << "  --no-norm    Disable extrinsic normalization (default: on)\n";
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
  std::vector<float> ebn0_list = {2.f, 4.f, 6.f, 8.f};
  bool normalize_extrinsic = true;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      print_usage(argv[0]);
      return 0;
    } else if (arg == "--ebn0" && i + 1 < argc) {
      ebn0_list = parse_ebn0_list(argv[++i]);
    } else if (arg == "--no-norm") {
      normalize_extrinsic = false;
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      print_usage(argv[0]);
      return 1;
    }
  }

  Params p;
  PipelineConfig cfg;
  cfg.interleaver_name = "ofec";
  cfg.decoder_name = "ebchPF";
  cfg.normalize_extrinsic = normalize_extrinsic;
  for (float eb : ebn0_list) {
    std::cout << "\n=== [interleaved] Eb/N0 = " << eb << " dB ===\n";
    (void)run_pipeline(p, cfg, /*label=*/"interleaved", /*ebn0_db=*/eb);
  }
  return 0;
}
