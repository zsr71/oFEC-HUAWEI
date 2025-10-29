function [decoded, tile_stats] = ofec_decode_llr(llr_mat, p)
%OFEC_DECODE_LLR Top-level LLR decoding mirroring newcode::ofec_decode_llr.
N = p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM;
[RROWS, CCOLS] = size(llr_mat);
if CCOLS ~= N
    error('ofec_decode_llr:BadCols', 'llr_mat columns must equal N (%d).', N);
end

tile_stats = repmat(struct('triggered', 0, 'total', 0), 1, p.TILES_PER_WIN);
TILE_HEIGHT_ROWS = tile_height_rows(p);
TILE_STRIDE_ROWS = tile_stride_rows(p);
WIN_HEIGHT_ROWS  = win_height_rows(p);
POP_PUSH_ROWS    = pop_push_rows(p);
TILES_PER_WIN    = p.TILES_PER_WIN;

channel_llr = llr_mat;
work_llr = zeros(RROWS, N, 'single');

if RROWS < WIN_HEIGHT_ROWS
    decoded = channel_llr;
    return;
end

win_start = initial_win_start_rows(p);
last_ws = RROWS - WIN_HEIGHT_ROWS;

while win_start <= last_ws
    win_end = win_start + WIN_HEIGHT_ROWS - 1;
    [work_llr, tile_stats] = process_window(work_llr, channel_llr, win_start, win_end, p, ...
                                            TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN, tile_stats);
    win_start = win_start + POP_PUSH_ROWS;
end

decoded = channel_llr + work_llr;
end

function [work_llr, tile_stats] = process_window(work_llr, channel_llr, win_start, win_end, p, ...
                                                 tile_height_rows_val, ~, TILES_PER_WIN, tile_stats)
N = p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM;
for t = 0:(TILES_PER_WIN - 1)
    tile_bottom_row = win_end - t * tile_height_rows_val;
    tile_top_row = tile_bottom_row - tile_height_rows_val + 1;
    rows = tile_bottom_row - tile_top_row + 1;
    idx_rows = tile_top_row + 1:tile_bottom_row + 1;
    tile_in = work_llr(idx_rows, :);
    ch_tile = channel_llr(idx_rows, :);

    tile_params = p;
    tile_params.beta = pick_or_default(p.beta_list, t + 1, p.beta);
    tile_params.ALPHA = pick_or_default(p.ALPHA_LIST, t + 1, p.ALPHA);
    use_hard = pick_or_default(p.HARD_TILE_LIST, t + 1, double(p.HARD_DECODE_DEFAULT)) ~= 0;

    [tile_result, early_stop] = process_tile(tile_in, ch_tile, tile_params, tile_top_row, use_hard);

    tile_stats(t + 1).total = tile_stats(t + 1).total + 1;
    if early_stop
        tile_stats(t + 1).triggered = tile_stats(t + 1).triggered + 1;
    end

    work_llr(idx_rows, :) = tile_result(1:rows, 1:N);
end
end

function value = pick_or_default(list, idx, fallback)
if idx >= 1 && idx <= numel(list)
    value = list(idx);
else
    value = fallback;
end
end

function [tile_out, early_stop_triggered] = process_tile(tile_in, ch_tile, p, tile_top_row_global, use_hard)
B = p.BITS_PER_SUBBLOCK_DIM;
N = p.NUM_SUBBLOCK_COLS * B;
K = Params_BCH_K();
TAKE_BITS = K - N;
BCH_PAR = Params_BCH_PARITY_BITS();
OVR_IDX = Params_BCH_N();

[H, W] = size(tile_in);
if W ~= N
    error('process_tile:BadWidth', 'Unexpected tile width.');
end

tile_out = tile_in;
SBR = p.CHASE_SBR;
if SBR ~= 1 && SBR ~= 2
    error('process_tile:InvalidSBR', 'CHASE_SBR must be 1 or 2.');
end

rows_to_decode = SBR * B;
decoder_cols = 2 * N;
lin_matrix = zeros(rows_to_decode, decoder_cols, 'single');
lch_matrix = zeros(rows_to_decode, decoder_cols, 'single');
row_local_lookup = zeros(rows_to_decode, 1);
row_global_lookup = zeros(rows_to_decode, 1);

row_idx = 0;
for s = 0:(SBR - 1)
    sbr_row0_local = H - (SBR - s) * B;
    for r_off = 0:(B - 1)
        row_idx = row_idx + 1;
        row_local = sbr_row0_local + r_off;
        row_global = tile_top_row_global + row_local;

        row_local_lookup(row_idx) = row_local;
        row_global_lookup(row_idx) = row_global;

        R = floor(row_global / B);
        r = mod(row_global, B);

        for k = 0:(N - 1)
            br = bitxor(R, 1) - 2 * p.NUM_GUARD_SUBROWS - 2 * (N / B) + 2 * floor(k / B);
            bc = floor(k / B);
            bit_row_in_block = bitxor(mod(k, B), r);
            bit_col_in_block = r;

            rr_global = br * B + bit_row_in_block;
            cc_global = bc * B + bit_col_in_block;

            rr_local2 = rr_global - tile_top_row_global;
            cc_local2 = cc_global;

            Lch = ch_tile(rr_local2 + 1, cc_local2 + 1);
            La  = tile_in(rr_local2 + 1, cc_local2 + 1);
            lin_matrix(row_idx, k + 1) = single(Lch + La);
            lch_matrix(row_idx, k + 1) = single(Lch);
        end

        for i = 0:(TAKE_BITS - 1)
            col = N + i;
            Lch = ch_tile(row_local + 1, i + 1);
            La  = tile_in(row_local + 1, i + 1);
            lin_matrix(row_idx, col + 1) = single(Lch + La);
            lch_matrix(row_idx, col + 1) = single(Lch);
        end

        for j = 0:(BCH_PAR - 1)
            col = K + j;
            src_col = TAKE_BITS + j;
            Lch = ch_tile(row_local + 1, src_col + 1);
            La  = tile_in(row_local + 1, src_col + 1);
            lin_matrix(row_idx, col + 1) = single(Lch + La);
            lch_matrix(row_idx, col + 1) = single(Lch);
        end

        Lch = ch_tile(row_local + 1, TAKE_BITS + BCH_PAR + 1);
        La  = tile_in(row_local + 1, TAKE_BITS + BCH_PAR + 1);
        lin_matrix(row_idx, OVR_IDX) = single(Lch + La);
        lch_matrix(row_idx, OVR_IDX) = single(Lch);
    end
end

early_stop_triggered = tile_should_early_stop(lin_matrix);

decoder_res = Decoder_Core(lin_matrix, lch_matrix, use_hard, p);

for s = 0:(SBR - 1)
    for r_off = 0:(B - 1)
        row_idx = s * B + r_off + 1;
        if row_idx > size(decoder_res.lout, 1)
            continue;
        end
        if ~decoder_res.produced_rows(row_idx)
            continue;
        end
        row_local = row_local_lookup(row_idx);
        row_global = row_global_lookup(row_idx);
        lout_row = decoder_res.lout(row_idx, :);

        for i = 0:(TAKE_BITS - 1)
            tile_out(row_local + 1, i + 1) = lout_row(N + i + 1);
        end
        for j = 0:(BCH_PAR - 1)
            tile_out(row_local + 1, TAKE_BITS + j + 1) = lout_row(K + j + 1);
        end
        tile_out(row_local + 1, TAKE_BITS + BCH_PAR + 1) = lout_row(OVR_IDX);

        R = floor(row_global / B);
        r = mod(row_global, B);

        for k = 0:(N - 1)
            br = bitxor(R, 1) - 2 * p.NUM_GUARD_SUBROWS - 2 * (N / B) + 2 * floor(k / B);
            bc = floor(k / B);
            bit_row_in_block = bitxor(mod(k, B), r);
            bit_col_in_block = r;

            rr_global = br * B + bit_row_in_block;
            cc_global = bc * B + bit_col_in_block;

            rr_local2 = rr_global - tile_top_row_global;
            cc_local2 = cc_global;
            tile_out(rr_local2 + 1, cc_local2 + 1) = lout_row(k + 1);
        end
    end
end
end

function res = Decoder_Core(lin_matrix, lch_matrix, use_hard_decode, p)
expected_cols = 2 * p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM;
[rows, cols] = size(lin_matrix);
if cols ~= expected_cols
    error('Decoder_Core:UnexpectedCols', 'Unexpected column count.');
end
res.lout = zeros(rows, cols, 'single');
res.produced_rows = false(rows, 1);
for row = 1:rows
    LinVec = lin_matrix(row, :);
    LchVec = lch_matrix(row, :);
    if use_hard_decode
        error('Decoder_Core:HardDecodeNotSupported', 'Hard decode path not implemented.');
    else
        Y2 = chase_decode_256(LinVec, LchVec, p);
        produced = true;
    end
    if produced
        res.produced_rows(row) = true;
        res.lout(row, :) = Y2;
    end
end
end

function Y2 = chase_decode_256(Lin256, ~, p)
beta = p.beta;
alpha = p.ALPHA;
L = max(1, p.CHASE_L);
NTEST = max(1, p.CHASE_NTEST);

y = double(Lin256(:).');
hard_ch = y < 0;

abs_vals = abs(y(1:255));
[~, order] = sort(abs_vals, 'ascend');
L_eff = min(L, numel(order));
lrp_pos = order(1:L_eff) - 1;

patt = generate_test_patterns(L_eff, NTEST);

CW_all = zeros(NTEST, 256, 'uint8');
ok_mask = false(1, NTEST);
for c = 1:NTEST
    tmp = hard_ch;
    for j = 1:L_eff
        if patt(c, j)
            pos = lrp_pos(j) + 1;
            tmp(pos) = ~tmp(pos);
        end
    end
    [decoded, ok] = bch_255_239_decode_candidate(tmp);
    CW_all(c, :) = decoded;
    ok_mask(c) = ok;
end

codeset = build_unique_codeset_pm1(CW_all, ok_mask);
softin = y;
mean_soft = mean(abs(softin));
if mean_soft > 0
    softin = softin / mean_soft;
end

if isempty(codeset)
    fallback = ones(1, 256);
    fallback(1:255) = (~hard_ch(1:255)) * 2 - 1;
    fallback(256) = (parity256_from255(uint8(hard_ch(1:255))) == 0) * 2 - 1;
    codeset = fallback;
end

M = size(codeset, 1);
dists = zeros(M, 1);
for m = 1:M
    dists(m) = sum((softin - codeset(m, :)).^2);
end
[distD, d_idx] = min(dists);
d = codeset(d_idx, :);

omega = zeros(1, 256);
if M == 1
    omega = beta * d;
else
    for j = 1:256
        mask = codeset(:, j) ~= d(j);
        if any(mask)
            bestComp = min(dists(mask));
            delta = bestComp - distD;
            omega(j) = 0.25 * delta * d(j) - softin(j);
        else
            omega(j) = beta * d(j);
        end
    end
end

mask_not_beta = abs(omega) ~= beta;
if any(mask_not_beta)
    g = mean(abs(omega(mask_not_beta)));
    if g > 0
        omega = omega / g;
    end
end
Y2 = single(alpha * omega);
end

function patt = generate_test_patterns(L, n_test)
if n_test < 1
    n_test = 1;
end
patt = false(n_test, L);
if L >= 0 && L < 31
    full = 2^L;
    if n_test <= full
        for c = 0:(n_test - 1)
            bits = bitget(uint32(c), 1:L);
            patt(c + 1, :) = bits;
        end
        return;
    end
end
c = 1;
if c <= n_test
    c = c + 1;
end
for layer = 1:L
    for j = layer:L
        if c > n_test
            return;
        end
        patt(c, 1:L) = false;
        patt(c, 1:layer-1) = true;
        patt(c, j) = true;
        c = c + 1;
    end
end
end

function [decoded, ok] = bch_255_239_decode_candidate(bits256)
decoded = zeros(1, 256, 'uint8');
decoded(1:255) = uint8(bits256(1:255));
[ok, out255] = bch_255_239_decode_hiho_cw_255(decoded(1:255));
decoded(1:255) = out255;
decoded(256) = parity256_from255(out255);
end

function parity = parity256_from255(cw255)
parity = mod(sum(cw255), 2);
end

function codeset = build_unique_codeset_pm1(CW_all, ok_mask)
mask_idx = find(ok_mask);
if isempty(mask_idx)
    codeset = [];
    return;
end
candidates = (1 - 2 * double(CW_all(mask_idx, :)));
[~, ia] = unique(candidates, 'rows');
codeset = candidates(ia, :);
end

function early = tile_should_early_stop(lin_matrix)
rows = size(lin_matrix, 1);
hard255 = zeros(1, 255, 'uint8');
decoded255 = zeros(1, 255, 'uint8');
for r = 1:rows
    for j = 1:255
        hard255(j) = uint8(lin_matrix(r, j) < 0);
    end
    [ok, decoded255] = bch_255_239_decode_hiho_cw_255(hard255);
    if ~ok
        early = false;
        return;
    end
    parity255 = mod(sum(hard255), 2);
    overall = uint8(lin_matrix(r, 256) < 0);
    if bitxor(parity255, overall) ~= 0
        early = false;
        return;
    end
end
early = true;
end
