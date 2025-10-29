function val = pop_push_rows(p)
%POP_PUSH_ROWS Pop/push stride in bit-rows for sliding window.
val = p.WINDOW_POP_PUSH * p.BITS_PER_SUBBLOCK_DIM;
end
