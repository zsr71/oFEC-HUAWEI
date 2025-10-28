function val = win_height_rows(p)
%WIN_HEIGHT_ROWS Convert Params window height to bit-rows.
val = p.TILES_PER_WIN * p.TILE_HEIGHT_BR * p.BITS_PER_SUBBLOCK_DIM ...
      - (p.TILES_PER_WIN - 1) * p.TILE_OVERLAP_BR * p.BITS_PER_SUBBLOCK_DIM;
end
