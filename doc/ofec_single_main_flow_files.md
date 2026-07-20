可以，按 `ofec_single` 一次运行的主链路，主要代码文件是这些：

1. 入口与参数组装  
- [apps/ofec_single.cpp](/home/zsr71/projects/newcode/apps/ofec_single.cpp): `main`，你平时改的运行参数都在这里。  
- [include/newcode/ofec_single_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_single_runner.hpp): `Config` 定义。  
- [src/ofec_single/ofec_single_runner.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_runner.cpp): `run_ofec_single` 主入口。  
- [src/ofec_single/ofec_single_params.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_params.cpp): 把 `Config` 转成 `Params`，并校验 `SISO_ACTIVE_LIST/ALPHA_LIST/beta_list`。  

2. pipeline 总流程（发端到收端）  
- [include/newcode/pipeline_runner.hpp](/home/zsr71/projects/newcode/include/newcode/pipeline_runner.hpp): pipeline 接口。  
- [src/common/pipeline/pipeline_runner.cpp](/home/zsr71/projects/newcode/src/common/pipeline/pipeline_runner.cpp): 完整流程：bitgen -> ofec encode -> 交织 -> 调制 -> AWGN -> 解调LLR -> 反交织 -> decoder。  

3. decoder 工厂与分发  
- [include/newcode/decoder_api.hpp](/home/zsr71/projects/newcode/include/newcode/decoder_api.hpp): 解码器抽象接口。  
- [src/rx/decoder/factories/plain_decoder_factory.cpp](/home/zsr71/projects/newcode/src/rx/decoder/factories/plain_decoder_factory.cpp): `plain` 解码器入口。  
- [src/rx/decoder/factories/ebchPF_decoder_factory.cpp](/home/zsr71/projects/newcode/src/rx/decoder/factories/ebchPF_decoder_factory.cpp): `ebchPF` 解码器入口。  

4. oFEC frame/window/tile 主解码链  
- [src/rx/ofec/ofec_frame_decode.cpp](/home/zsr71/projects/newcode/src/rx/ofec/ofec_frame_decode.cpp): `ofec_decode_llr_plain/ebchPF`。  
- [src/rx/ofec/detail/ofec_decode_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_decode_impl.ipp): frame 级循环窗口，调用 `process_window_impl`。  
- [src/rx/ofec/detail/ofec_window_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_window_impl.ipp): window 内逐 tile 处理，按 tile 取 `ALPHA/beta/SISO_ACTIVE_LIST`。  
- [src/rx/ofec/detail/ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp): tile 主流程（prepare -> earlystop/mux -> decode -> writeback）。  
- [src/rx/ofec/detail/ofec_tile_input.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_input.ipp): 构造每行 256 输入。  
- [src/rx/ofec/detail/ofec_tile_decode.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_decode.ipp): 调 core decoder、归一化/量化。  
- [src/rx/ofec/detail/ofec_tile_writeback.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_writeback.ipp): 写回到 tile/work_llr。  

5. early-stop 与 mux（你最近在改的核心）  
- [src/rx/ofec/earlystop/tile_early_stop_stats.ipp](/home/zsr71/projects/newcode/src/rx/ofec/earlystop/tile_early_stop_stats.ipp): 计算 `TileEarlyStopResult` 和 `row_passed_flags`。  
- [src/rx/ofec/mux/mux_state_builder.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_state_builder.cpp): `row_passed_flags -> state(0/1/2)`。  
- [src/rx/ofec/mux/mux_siso_budget.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_siso_budget.cpp): 按 tile 的 SISO 配额把超出的 `0` 改成 `2`。  
- [src/rx/ofec/mux/mux_config_validate.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_config_validate.cpp): 校验 `SISO_ACTIVE_LIST`。  

6. 行级 core decoder（最终执行 Chase/早停分支）  
- [src/rx/ofec/ofec_row_decoder_core.cpp](/home/zsr71/projects/newcode/src/rx/ofec/ofec_row_decoder_core.cpp): 根据 `mux_state` 选择 `skip / early-stop / chase`。  
- [src/rx/ofec/earlystop/row_early_stop_process_2.ipp](/home/zsr71/projects/newcode/src/rx/ofec/earlystop/row_early_stop_process_2.ipp): 当前早停输出路径。  
- [src/rx/ofec/plain/chase256_plain.cpp](/home/zsr71/projects/newcode/src/rx/ofec/plain/chase256_plain.cpp): plain Chase。  
- [src/rx/ofec/ebchPF/chase256_ebchPF.cpp](/home/zsr71/projects/newcode/src/rx/ofec/ebchPF/chase256_ebchPF.cpp): ebchPF Chase。  

如果你愿意，我可以下一步给你画一版“只看 mux/early-stop 的最短阅读路径（5个文件）”。
