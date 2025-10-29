function val = Params_BCH_PARITY_BITS()
%PARAMS_BCH_PARITY_BITS Number of BCH parity bits (excluding overall parity).
val = Params_BCH_N() - Params_BCH_K() - 1;
end
