function result = run_ofec_single(varargin)
%RUN_OFEC_SINGLE MATLAB replica of apps/ofec_single.cpp pipeline.
%   RESULT = RUN_OFEC_SINGLE() executes the default scenario. Optional
%   name/value pairs mirror the C++ constants, e.g.:
%       run_ofec_single('Label','debug_L6','EbN0dB',3.2,'ChaseL',6);


if mod(numel(varargin), 2) ~= 0
    error('run_ofec_single:InvalidArgs', 'Name/value arguments must come in pairs.');
end

opts = struct('Label', 'debug_L6', ...
              'EbN0dB', 3.67, ...
              'ChaseL', 6, ...
              'AlphaFill', 1.0, ...
              'BetaFill', 0.80, ...
              'AlphaExplicit', [0.3,0.4,0.5,0.6], ...
              'BetaExplicit', [0.6,0.7,0.9,1.0]);

for k = 1:2:numel(varargin)
    name = varargin{k};
    value = varargin{k + 1};
    if ~(ischar(name) || isstring(name))
        error('run_ofec_single:InvalidOptionName', 'Option names must be strings.');
    end
    switch lower(string(name))
        case "label"
            opts.Label = char(value);
        case {"ebn0", "ebn0db"}
            opts.EbN0dB = double(value);
        case {"chasel", "chasel_override"}
            opts.ChaseL = double(value);
        case "alphafill"
            opts.AlphaFill = double(value);
        case "betafill"
            opts.BetaFill = double(value);
        case {"alphaexplicit", "alphalist"}
            opts.AlphaExplicit = double(value(:)');
        case {"betaexplicit", "betalist"}
            opts.BetaExplicit = double(value(:)');
        otherwise
            error('run_ofec_single:UnsupportedOption', ...
                  'Unsupported option "%s".', name);
    end
end

fprintf('--- run_ofec_single ---\n');
fprintf('[SETUP] label="%s", Eb/N0=%.2f dB\n', opts.Label, opts.EbN0dB);

p = default_params();
if opts.ChaseL >= 0
    p.CHASE_L = opts.ChaseL;
    p.CHASE_NTEST = 2^p.CHASE_L;
end

T = p.TILES_PER_WIN;

if ~isempty(opts.AlphaExplicit)
    if numel(opts.AlphaExplicit) ~= T
        error('run_ofec_single:AlphaExplicitLength', ...
              'AlphaExplicit length must equal TILES_PER_WIN (%d).', T);
    end
    p.ALPHA_LIST = opts.AlphaExplicit;
else
    p.ALPHA_LIST = repmat(opts.AlphaFill, 1, T);
end
if ~isempty(p.ALPHA_LIST), p.ALPHA = p.ALPHA_LIST(1); end

if ~isempty(opts.BetaExplicit)
    if numel(opts.BetaExplicit) ~= T
        error('run_ofec_single:BetaExplicitLength', ...
              'BetaExplicit length must equal TILES_PER_WIN (%d).', T);
    end
    p.beta_list = opts.BetaExplicit;
else
    p.beta_list = repmat(opts.BetaFill, 1, T);
end
if ~isempty(p.beta_list), p.beta = p.beta_list(1); end

fprintf('[RUN] Scenario: %s\n', opts.Label);

info_bits = generate_bits(p);
fprintf('[INFO] Generated bits: %u\n', numel(info_bits));

code_matrix = ofec_encode(info_bits, p);
[n_rows, n_cols] = size(code_matrix);
fprintf('[INFO] oFEC matrix: %u x %u\n', n_rows, n_cols);

coded_bits = reshape(code_matrix.', 1, []);
fprintf('[INFO] Coded bits (flattened): %u\n', numel(coded_bits));

n_bps = 2; % QPSK
tx_syms = qam_modulate(coded_bits, n_bps);
fprintf('[INFO] Modulated symbols: %u (Es≈1)\n', numel(tx_syms));

N = p.NUM_SUBBLOCK_COLS * p.BITS_PER_SUBBLOCK_DIM;
K = Params_BCH_K();
TAKEBITS = K - N;
code_rate = double(TAKEBITS) / double(N);

awgn_seed = uint32(p.BITGEN_SEED + 100);
rx_syms = add_awgn(tx_syms, opts.EbN0dB, n_bps, awgn_seed, code_rate);

fprintf('[INFO] Eb/N0 set to %.2f dB\n', opts.EbN0dB);
print_example_symbols(tx_syms, rx_syms);

llr = qam_llr_from_ebn0(rx_syms, n_bps, opts.EbN0dB, code_rate);
fprintf('[INFO] LLR count: %u\n', numel(llr));
print_first_llrs(llr);

llr_mat = reshape(llr, n_cols, n_rows).'; % row-major reshape
llr_mat = apply_known_zero_prefix(llr_mat, p);

[decoded_llr, tile_stats] = ofec_decode_llr(llr_mat, p);

rx_info_pre  = rx_info_from_bit_llr(llr_mat, p);
rx_info_post = rx_info_from_bit_llr(decoded_llr, p);
fprintf('[INFO] rx_info_bits: %u (flattened, warmup skipped)\n', numel(rx_info_pre));

pre_stats  = compute_and_print_ber(info_bits, rx_info_pre,  sprintf('%s Pre-FEC',  opts.Label), p);
post_stats = compute_and_print_ber(info_bits, rx_info_post, sprintf('%s Post-FEC', opts.Label), p);

tile_pct = compute_early_stop_percentages(tile_stats);

result = struct();
result.label = opts.Label;
result.ebn0_db = opts.EbN0dB;
result.pre_fec = pre_stats;
result.post_fec = post_stats;
result.tile_early_stop_pct = tile_pct;

fprintf('[DONE] Pipeline complete. Tile early-stop (%%): %s\n', sprintf('%.1f ', tile_pct));
end
