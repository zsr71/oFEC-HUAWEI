function mat = ofec_encode(bits, p)
%OFEC_ENCODE Mirror of newcode::ofec_encode -> matrix of uint8 bits.
B = p.BITS_PER_SUBBLOCK_DIM;
N = p.NUM_SUBBLOCK_COLS * B;
G = p.NUM_GUARD_SUBROWS;
K = Params_BCH_K();
PAR_LEN = Params_BCH_PARITY_BITS();
TAKE_BITS = K - N;

warmup_rows = win_height_rows(p);
mat = zeros(warmup_rows, N, 'uint8');

bit_pos = 1;
global_row = warmup_rows;
bits_len = numel(bits);

while bit_pos <= bits_len
    R = floor(global_row / B);
    r = mod(global_row, B);

    row = zeros(1, N, 'uint8');
    info239 = zeros(1, K, 'uint8');
    idx = 1;

    for k = 0:(N - 1)
        br = bitxor(R, 1) - 2 * G - 2 * (N / B) + 2 * floor(k / B);
        bc = floor(k / B);
        bit_row_in_block = bitxor(mod(k, B), r);
        bit_col_in_block = r;

        rr = br * B + bit_row_in_block;
        cc = bc * B + bit_col_in_block;

        info239(idx) = mat(rr + 1, cc + 1);
        idx = idx + 1;
    end

    for i = 0:(TAKE_BITS - 1)
        if bit_pos <= bits_len
            v = bits(bit_pos);
            bit_pos = bit_pos + 1;
        else
            v = uint8(0);
        end
        info239(idx) = v;
        idx = idx + 1;
        row(i + 1) = v;
    end

    code_bits = bch_255_239_encode(info239);
    parity16 = code_bits(K + 1:end);
    row(TAKE_BITS + (1:PAR_LEN)) = parity16;

    overall = uint8(mod(sum(code_bits), 2));
    row(TAKE_BITS + PAR_LEN + 1) = overall;

    global_row = global_row + 1;
    mat(global_row, :) = row;
end
end
