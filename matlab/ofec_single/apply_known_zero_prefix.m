function mat = apply_known_zero_prefix(mat, p)
%APPLY_KNOWN_ZERO_PREFIX Clamp known prefix rows to bit=0 and normalise tail.
R = size(mat, 1);
if R == 0
    return;
end
known_rows = min(win_height_rows(p), R);
if known_rows > 0
    mat(1:known_rows, :) = cast(1.0, 'like', mat);
end
if known_rows < R
    tail = mat(known_rows + 1:end, :);
    if ~isempty(tail)
        llrMean = mean(abs(tail), 'all');
        if llrMean > 0
            mat(known_rows + 1:end, :) = tail / llrMean;
        end
    end
end
end
