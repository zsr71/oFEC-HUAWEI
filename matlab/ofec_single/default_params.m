function p = default_params()
%DEFAULT_PARAMS Populate struct mirroring newcode::Params defaults.
p.NUM_SUBBLOCK_COLS     = 8;
p.INFO_SUBROWS_PER_CODE = 16;
p.BITS_PER_SUBBLOCK_DIM = 16;

p.NUM_INFO_BITS     = (1 * 88+2) * 16 * 111;
p.BITGEN_SEED       = 43;
p.NUM_GUARD_SUBROWS = 2;

p.TILES_PER_WIN   = 4;
p.TILE_OVERLAP_BR = 0;
p.TILE_HEIGHT_BR  = 22;
p.WINDOW_POP_PUSH = 2;

p.LLR_BITS = 16;
p.LLR_CLIP = 8.0;

p.CHASE_L     = 6;
p.CHASE_NTEST = 64;
p.CHASE_SBR   = 2;
p.CHASE_TP    = 1;

p.beta = 0.35;
p.ALPHA = 1.0;
p.ALPHA_LIST = [0.3, 0.4, 0.5, 0.6];
p.beta_list  = [0.9, 1.0, 1.1, 1.2];

p.HARD_DECODE_DEFAULT = false;
p.HARD_TILE_LIST = [0, 0, 0, 0, 0];
p.HARD_LLR_MAG = 12.0;
end
