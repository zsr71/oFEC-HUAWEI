function print_first_llrs(llr)
%PRINT_FIRST_LLRS Dump first few LLR values for inspection.
count = min(8, numel(llr));
fprintf('[INFO] First few LLRs: ');
for i = 1:count
    if i < count
        fprintf('%.4f, ', llr(i));
    else
        fprintf('%.4f\n', llr(i));
    end
end
end
