#include "newcode/bch_255_239.hpp"

#include <array>
#include <cstdio>
#include <cstdint>

int main()
{
  std::array<uint8_t, 255> tmp_in{};
  std::array<uint8_t, 255> decoded{};
  int corrected_errors = 0;

  const bool ok =
      newcode::bch_255_239_decode_hiho_cw_255(tmp_in.data(),
                                              decoded.data(),
                                              &corrected_errors);

  std::printf("Testing bch_255_239_decode_hiho_cw_255 with all-zero input\n");
  std::printf("ok = %s, corrected_errors = %d\n", ok ? "true" : "false",
              corrected_errors);
  std::printf("Decoded codeword (first 32 bits): ");
  for (std::size_t i = 0; i < 32 && i < decoded.size(); ++i) {
    std::printf("%u", decoded[i] & 1u);
  }
  std::printf(" ...\n");

  return ok ? 0 : 1;
}

