function val = tile_stride_rows(p)
%TILE_STRIDE_ROWS Tile stride in bit-rows.
val = (p.TILE_HEIGHT_BR - p.TILE_OVERLAP_BR) * p.BITS_PER_SUBBLOCK_DIM;
end
