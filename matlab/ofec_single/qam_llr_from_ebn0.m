function llr = qam_llr_from_ebn0(y, n_bps, ebn0_dB, code_rate)
%QAM_LLR_FROM_EBN0 Compute bit LLRs via exact log-sum-exp demod.
sigma = ebn0_to_sigma(ebn0_dB, n_bps, code_rate);
llr = qam_llr_logsumexp(y, n_bps, sigma);
end

function llr = qam_llr_logsumexp(y, n_bps, sigma)
if mod(n_bps, 2) ~= 0
    error('qam_llr_logsumexp:InvalidNBPS', 'n_bps must be a positive even integer.');
end
if sigma <= 0
    error('qam_llr_logsumexp:InvalidSigma', 'sigma must be > 0.');
end
M = 2^n_bps;
constellation = build_constellation(n_bps);
inv_sigma2 = 1 / (2 * sigma^2);
N_bits = numel(y) * n_bps;
llr = zeros(1, N_bits, 'single');
for n = 1:N_bits
    b = mod(n - 1, n_bps);
    k = floor((n - 1) / n_bps) + 1;
    yk = y(k);
    m0 = -Inf;
    m1 = -Inf;
    metrics = zeros(1, M);
    for j = 0:(M - 1)
        d2 = abs(yk - constellation(j + 1))^2;
        met = -d2 * inv_sigma2;
        metrics(j + 1) = met;
        if bitand(j, bitshift(1, b)) == 0
            m0 = max(m0, met);
        else
            m1 = max(m1, met);
        end
    end
    sum0 = 0;
    sum1 = 0;
    for j = 0:(M - 1)
        met = metrics(j + 1);
        if bitand(j, bitshift(1, b)) == 0
            sum0 = sum0 + exp(double(met - m0));
        else
            sum1 = sum1 + exp(double(met - m1));
        end
    end
    llr(n) = single((m0 + log(sum0)) - (m1 + log(sum1)));
end
end

function constellation = build_constellation(n_bps)
M = 2^n_bps;
bits = zeros(1, M * n_bps, 'uint8');
for j = 0:(M - 1)
    for b = 0:(n_bps - 1)
        bits(j * n_bps + b + 1) = bitand(bitshift(j, -b), 1);
    end
end
constellation = qam_modulate(bits, n_bps);
end
