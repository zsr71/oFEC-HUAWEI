function info = rx_info_from_bit_llr(bit_llr_mat, p)
%RX_INFO_FROM_BIT_LLR Extract flattened info bits from LLR matrix.
B = p.BITS_PER_SUBBLOCK_DIM;
N = p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM;
K = Params_BCH_K();
TAKE_BITS = K - N;
if mod(N, B) ~= 0
    error('rx_info_from_bit_llr:InvalidDimensions', 'N must be multiple of B.');
end
[rows, cols] = size(bit_llr_mat);
if cols ~= N
    error('rx_info_from_bit_llr:InvalidCols', 'Unexpected column count.');
end
warmup_rows = min(win_height_rows(p), rows);
useful_rows = rows - warmup_rows;
info = zeros(1, useful_rows * TAKE_BITS, 'uint8');
idx = 1;
for r = warmup_rows + 1:rows
    for c = 1:TAKE_BITS
        info(idx) = uint8(bit_llr_mat(r, c) < 0);
        idx = idx + 1;
    end
end
end
