#pragma once

namespace new_float_only {
namespace detail {

static void dump_chase_csv(const new_float_only::Params::DebugTraceConfig& trace,
                           const float* y,
                           const uint8_t* hard_ch,
                           const uint8_t* ML,
                           const float* omega)
{
    if (!trace.enable || !trace.dump_chase_csv) return;
    if (trace.chase_tile_index < 0 || trace.chase_invocation < 0) return;

    std::vector<new_float_only::Params::DebugTraceConfig::ChaseTraceEntry> entries =
        trace.active_chase_entries;
    if (entries.empty()) {
        if (trace.chase_decoder_col < 0) return;
        new_float_only::Params::DebugTraceConfig::ChaseTraceEntry fallback;
        fallback.row_index = trace.chase_decoder_row;
        fallback.k = trace.chase_decoder_col;
        fallback.global_row = trace.row;
        fallback.global_col = trace.col;
        fallback.bit_index = -1;
        entries.push_back(fallback);
    }

    const int tile_idx = trace.chase_tile_index;
    const int inv_id   = trace.chase_invocation;

    namespace fs = std::filesystem;
    const fs::path dir = trace.chase_csv_dir.empty()
                           ? fs::path("data/chase_csv")
                           : fs::path(trace.chase_csv_dir);
    std::error_code ec;
    fs::create_directories(dir, ec);

    for (const auto& entry : entries) {
        if (entry.row_index < 0 || entry.k < 0) continue;
        std::string label = entry.label;
        if (label.empty()) {
            if (entry.bit_index >= 0) {
                label = "bit" + std::to_string(entry.bit_index);
            } else {
                label = "row" + std::to_string(entry.global_row) +
                        "_col" + std::to_string(entry.global_col);
            }
        }
        std::string safe_label = label;
        for (char& ch : safe_label) {
            if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '-') {
                ch = '_';
            }
        }
        const fs::path file = dir / ("tile" + std::to_string(tile_idx)
                                     + "_call" + std::to_string(inv_id)
                                     + "_" + safe_label + ".csv");
        std::ofstream out(file);
        if (!out) continue;
        out << "target_label," << label << '\n';
        out << "target_bit_index," << entry.bit_index << '\n';
        out << "target_global_row," << entry.global_row << '\n';
        out << "target_global_col," << entry.global_col << '\n';
        out << "target_expected_bit," << entry.expected_bit << '\n';
        out << "lin_matrix_row_index," << entry.row_index << '\n';
        out << "lin_matrix_k," << entry.k << '\n';
        const std::vector<int8_t>* expected_bits_row = trace.chase_expected_bits_row;
        const bool has_competing =
            entry.cplus_bits.size() == detail::BCH_N_TOTAL &&
            entry.cminus_bits.size() == detail::BCH_N_TOTAL;
        out << "Index,LLR,Omega,ML,HardCh,ExpectedBit";
        if (has_competing) {
            out << ",CplusBit,CminusBit";
        }
        out << '\n';
        for (int k = 0; k < BCH_N_TOTAL; ++k) {
            int expected_bit = -1;
            if (expected_bits_row && static_cast<size_t>(k) < expected_bits_row->size()) {
                expected_bit = (*expected_bits_row)[static_cast<size_t>(k)];
            }
            out << k << ',' << y[k] << ',' << omega[k] << ','
                << int(ML[k]) << ',' << int(hard_ch[k]) << ','
                << expected_bit;
            if (has_competing) {
                out << ',' << int(entry.cplus_bits[static_cast<size_t>(k)])
                    << ',' << int(entry.cminus_bits[static_cast<size_t>(k)]);
            }
            out << '\n';
        }
    }
}

} // namespace detail
} // namespace new_float_only

