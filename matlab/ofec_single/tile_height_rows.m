function val = tile_height_rows(p)
%TILE_HEIGHT_ROWS Height of one tile in bit-rows.
val = p.TILE_HEIGHT_BR * p.BITS_PER_SUBBLOCK_DIM;
end
