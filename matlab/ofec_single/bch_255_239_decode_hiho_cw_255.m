function [success, out255] = bch_255_239_decode_hiho_cw_255(in255)
%BCH_255_239_DECODE_HIHO_CW_255 Decode 255-bit BCH codeword via System object.

n_core = Params_BCH_N() - 1;
k = Params_BCH_K();

persistent dec gp_coeffs
if isempty(dec)
    gp = bchgenpoly(n_core, k);
    gp_coeffs = double(gp.x);
    dec = comm.BCHDecoder('CodewordLength', n_core, ...
                          'MessageLength', k, ...
                          'GeneratorPolynomialSource', 'Property', ...
                          'GeneratorPolynomial', gp_coeffs, ...
                          'NumCorrectedErrorsOutputPort', true);
end

code_vec = logical(uint8(bitand(in255(:), 1)));
[decoded_msg, num_err] = dec(code_vec);
success = all(~isnan(num_err));
if success
    out255 = bch_255_239_encode(decoded_msg.');
else
    out255 = uint8(code_vec(:).');
end
end
