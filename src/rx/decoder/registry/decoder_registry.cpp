#include "newcode/decoder_api.hpp"

#include <mutex>
#include <vector>

namespace newcode {
namespace {

std::vector<DecoderFactory>& factories() {
  static auto* storage = new std::vector<DecoderFactory>();
  return *storage;
}

std::mutex& factory_mutex() {
  static auto* m = new std::mutex();
  return *m;
}

} // namespace

void register_decoder_factory(DecoderFactory factory) {
  if (!factory) return;
  std::lock_guard<std::mutex> lock(factory_mutex());
  factories().push_back(std::move(factory));
}

std::unique_ptr<IDecoder> make_decoder(const std::string& name) {
  std::lock_guard<std::mutex> lock(factory_mutex());
  auto& list = factories();
  for (auto it = list.rbegin(); it != list.rend(); ++it) {
    if (!*it) continue;
    if (auto decoder = (*it)(name)) {
      return decoder;
    }
  }
  return nullptr;
}

} // namespace newcode
