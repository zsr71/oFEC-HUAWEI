#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "newcode/pipeline_runner.hpp"

namespace {

// oFEC 编码器要求输入比特数是 32*111 的整数倍。
// 这里固定取 2 个矩形块，尽量把回归测试做得足够小、足够快。
constexpr std::size_t kRectBits = 32u * 111u;

newcode::Params make_base_params() {
  newcode::Params p;
  // 用一套很小但完整的参数跑真实 pipeline，
  // 目的不是看 BER，而是让 tile/window 调度逻辑完整走一遍。
  p.NUM_INFO_BITS = 2u * kRectBits;
  p.BITGEN_SEED = 20260506;
  p.CHANNEL_SEED = 20260507;
  p.BITGEN_RANDOM_BITS = true;
  p.LLR_BITS = 6;
  p.LLR_CLIP_RATIO = 0.5f;

  p.ENABLE_EARLY_STOP = true;
  p.EARLY_STOP_ENABLE_LIST = {1, 1, 1, 1, 1, 1};
  p.EARLY_STOP_CONDITION_MODE = 1;
  p.EARLY_STOP_ACTION_MODE = 1;
  p.EARLY_STOP_BIND_GROUP_SIZE = 1;
  p.EARLY_STOP_COND_V1_REQUIRE_BCH = true;
  p.EARLY_STOP_COND_V1_REQUIRE_OVERALL = true;

  p.SISO_ACTIVE_LIST = {32, 24, 16, 12, 8, 8};
  p.MUX_GROUP_G = 1;
  p.MUX_SCHEDULING_MODE = 0;
  p.MUX_EARLY_STOP_PRIORITY_RULE = 0;
  p.MUX_ENABLE_RECONFIG = false;

  p.ALPHA_LIST = {0.428571f, 0.447738f, 0.482782f, 0.528162f, 0.581902f, 0.642857f};
  p.beta_list = {2.857143f, 6.179301f, 12.253626f, 20.119585f, 29.434408f, 40.000000f};
  p.ALPHA = p.ALPHA_LIST.front();
  p.beta = p.beta_list.front();
  p.EARLY_STOP_ACTION_SIGN_BETA_LIST = p.beta_list;
  p.EARLY_STOP_ACTION_SIGN_BETA = p.EARLY_STOP_ACTION_SIGN_BETA_LIST.front();

  p.HYBRID_ENABLE = false;
  p.HYBRID_USE_FAST_CLASSIFIER = false;
  p.HYBRID_NORMALIZE_SOFT_ONLY = false;
  // 打开这个开关后，soft tile 会在新路径结束后额外再跑一次旧路径做逐 tile 对比。
  p.HYBRID_VERIFY_DISABLED_MATCH_LEGACY = true;
  return p;
}

newcode::PipelineConfig make_pipeline_config() {
  newcode::PipelineConfig cfg;
  cfg.decoder_name = "chase_baseline";
  cfg.interleaver_name = "identity";
  cfg.normalize_extrinsic = false;
  cfg.bits_per_symbol = 1;
  cfg.quiet = true;
  return cfg;
}

newcode::PipelineResult run_case(const std::string& label,
                                 newcode::Params params,
                                 float ebn0_db) {
  const auto cfg = make_pipeline_config();
  // 若新旧路径存在任何差异，run_pipeline 内部会直接抛异常，
  // 这里能走到打印 PASS，就说明整条验证链通过了。
  const auto result = newcode::run_pipeline(params, cfg, label, ebn0_db);
  std::cout << "[PASS] " << label
            << " pre_fec=" << result.pre_fec.ber
            << " post_fec=" << result.post_fec.ber;
  if (!result.tile_hard_finish_count.empty()) {
    std::cout << " hard_finish_counts=[";
    for (std::size_t i = 0; i < result.tile_hard_finish_count.size(); ++i) {
      if (i) std::cout << ", ";
      std::cout << result.tile_hard_finish_count[i];
    }
    std::cout << "]";
  }
  std::cout << "\n";
  return result;
}

}  // namespace

int main() {
  try {
    // case 1: 只验证 “MUX 路径” 一致性，不让 early-stop 参与。
    auto mux_only = make_base_params();
    mux_only.ENABLE_EARLY_STOP = false;
    mux_only.EARLY_STOP_ENABLE_LIST = {0, 0, 0, 0, 0, 0};
    run_case("hybrid_disabled_regression_mux_only", mux_only, 3.05f);

    // case 2: 让 early-stop 和 MUX 一起参与，覆盖更接近真实配置的路径。
    auto early_stop_and_mux = make_base_params();
    run_case("hybrid_disabled_regression_early_stop", early_stop_and_mux, 3.05f);

    // case 3: 打开方案三和 S0/S1/S3 快速分类器。
    // 这个 case 不和旧路径做等价比较，只确认 fast-classifier 分支能完整跑通。
    auto fast_classifier = make_base_params();
    fast_classifier.HYBRID_ENABLE = true;
    fast_classifier.HYBRID_USE_FAST_CLASSIFIER = true;
    fast_classifier.HYBRID_VERIFY_DISABLED_MATCH_LEGACY = false;
    run_case("hybrid_fast_classifier_smoke", fast_classifier, 3.05f);

    // case 4: 只让后几个 tile 启用方案三，前几个 tile 必须保持旧 soft 路径。
    auto hybrid_tail_only = make_base_params();
    hybrid_tail_only.HYBRID_ENABLE = false;
    hybrid_tail_only.HYBRID_ENABLE_LIST = {0, 0, 0, 1, 1, 1};
    hybrid_tail_only.HYBRID_USE_FAST_CLASSIFIER = true;
    hybrid_tail_only.HYBRID_VERIFY_DISABLED_MATCH_LEGACY = false;
    const auto tail_result =
        run_case("hybrid_tail_only_smoke", hybrid_tail_only, 3.05f);
    for (std::size_t i = 0; i < 3 && i < tail_result.tile_hard_finish_count.size(); ++i) {
      if (tail_result.tile_hard_finish_count[i] != 0) {
        throw std::runtime_error(
            "hybrid_tail_only_smoke: front tiles must not produce hard-finish rows");
      }
    }
  } catch (const std::exception& ex) {
    std::cerr << "[FAIL] hybrid regression check: " << ex.what() << "\n";
    return 1;
  }
  return 0;
}
