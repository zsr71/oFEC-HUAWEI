function stats = compute_and_print_ber(ref_bits, rx_bits, label, p)
%COMPUTE_AND_PRINT_BER Wrapper around compute_ber with logging.
stats = compute_ber(ref_bits, rx_bits, p);
L = min(numel(ref_bits), numel(rx_bits));
row_bits = Params_BCH_K() - p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM;
win_bits = min(L, win_height_rows(p) * row_bits);
cut_total = min(L, win_bits) + min(max(L - win_bits, 0), win_bits);
fprintf('[RESULT] %s BER=%.6g  (errs=%u / %u compared, cut=%u of %u)\n', ...
        label, stats.ber, stats.errors, stats.total, cut_total, L);
end
