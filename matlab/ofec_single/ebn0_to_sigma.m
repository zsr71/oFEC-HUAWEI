function sigma = ebn0_to_sigma(ebn0_dB, bits_per_symbol, code_rate)
%EBN0_TO_SIGMA Convert Eb/N0 to Gaussian std dev (per dimension).
if bits_per_symbol <= 0
    error('ebn0_to_sigma:InvalidM', 'bits_per_symbol must be > 0.');
end
if code_rate <= 0 || code_rate > 1
    error('ebn0_to_sigma:InvalidRate', 'code_rate must lie in (0,1].');
end
ebn0_lin = 10^(ebn0_dB / 10);
esn0_lin = ebn0_lin * bits_per_symbol * code_rate;
N0 = 1 / esn0_lin;
sigma = sqrt(N0 / 2);
end
