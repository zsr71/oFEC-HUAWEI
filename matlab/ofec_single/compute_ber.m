function stats = compute_ber(ref_bits, rx_bits, p)
%COMPUTE_BER Compare bit sequences with edge-window trimming.
L = min(numel(ref_bits), numel(rx_bits));
row_bits = p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM;
win_bits = min(L, win_height_rows(p) * row_bits);
skip_prefix = min(L, 0 * win_bits);
skip_suffix = min(L - skip_prefix, 0*win_bits);
start_idx = skip_prefix + 1;
stop_idx = L - skip_suffix;
if stop_idx < start_idx
    total = 0;
    errors = 0;
else
    cmp_ref = ref_bits(start_idx:stop_idx);
    cmp_rx  = rx_bits(start_idx:stop_idx);
    total = numel(cmp_ref);
    errors = sum(bitxor(uint8(cmp_ref), uint8(cmp_rx)));
end
stats.errors = errors;
stats.total = total;
if total == 0
    stats.ber = 0;
else
    stats.ber = double(errors) / double(total);
end
end
tpcdec