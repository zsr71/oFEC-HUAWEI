function y = add_awgn(x, ebn0_dB, bits_per_symbol, seed, code_rate)
%ADD_AWGN Add AWGN with variance computed from Eb/N0 and code rate.
if nargin < 5 || isempty(code_rate)
    code_rate = 1.0;
end
sigma = ebn0_to_sigma(ebn0_dB, bits_per_symbol, code_rate);
rng_state = rng;
cleanup = onCleanup(@() rng(rng_state));
rng(double(seed), 'twister');
noise = sigma * (randn(size(x), 'single') + 1i * randn(size(x), 'single'));
y = x + noise;
end
