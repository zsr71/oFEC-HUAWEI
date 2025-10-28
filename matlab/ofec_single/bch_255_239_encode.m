function code255 = bch_255_239_encode(info239)
%BCH_255_239_ENCODE Encode 239-bit message to 255-bit BCH codeword.

n_core = Params_BCH_N() - 1;
k = Params_BCH_K();

persistent enc gp_coeffs
if isempty(enc)
    gp = bchgenpoly(n_core, k);
    gp_coeffs = double(gp.x);
    enc = comm.BCHEncoder('CodewordLength', n_core, ...
                          'MessageLength', k, ...
                          'GeneratorPolynomialSource', 'Property', ...
                          'GeneratorPolynomial', gp_coeffs);
end

msg = logical(info239(:));
code_col = enc(msg);
code255 = uint8(code_col(:).');
end
