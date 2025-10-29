function syms = qam_modulate(bits, n_bps)
%QAM_MODULATE Square QAM (Gray mapping) with Es=1 normalisation.
if n_bps <= 0 || mod(n_bps, 2) ~= 0
    error('qam_modulate:InvalidNBPS', 'n_bps must be a positive even integer.');
end
bits = bitand(uint8(bits), 1);
n_sym = ceil(numel(bits) / n_bps);
syms = complex(zeros(1, n_sym, 'single'));
buf = zeros(1, n_bps, 'uint8');
sqrt_es = sqrt(2 * (2^n_bps - 1) / 3);
for s = 1:n_sym
    start_idx = (s - 1) * n_bps + 1;
    end_idx = min(start_idx + n_bps - 1, numel(bits));
    buf(:) = 0;
    buf(1:(end_idx - start_idx + 1)) = bits(start_idx:end_idx);
    m = n_bps / 2;
    I = 1 - 2 * double(buf(1));
    Q = 1 - 2 * double(buf(m + 1));
    for j = 2:m
        I = (1 - 2 * double(buf(j))) * (2^(j - 1) - I);
        Q = (1 - 2 * double(buf(m + j))) * (2^(j - 1) - Q);
    end
    syms(s) = complex(I / sqrt_es, Q / sqrt_es);
end
end
