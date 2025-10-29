function print_example_symbols(tx_syms, rx_syms)
%PRINT_EXAMPLE_SYMBOLS Display first few TX/RX symbol pairs.
count = min(3, numel(tx_syms));
fprintf('[INFO] Example symbols (TX -> RX):\n');
for i = 1:count
    fprintf('  %d: (%.4f, %.4f) -> (%.4f, %.4f)\n', i - 1, real(tx_syms(i)), imag(tx_syms(i)), ...
            real(rx_syms(i)), imag(rx_syms(i)));
end
end
