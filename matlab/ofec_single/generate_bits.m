function info_bits = generate_bits(p)
%GENERATE_BITS Draw random info bits using Params.BITGEN_SEED.
num_bits = p.NUM_INFO_BITS;
rng_state = rng;
cleanup = onCleanup(@() rng(rng_state));
rng(p.BITGEN_SEED, 'twister');
info_bits = uint8(randi([0, 1], 1, num_bits));
end
