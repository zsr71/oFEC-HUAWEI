#pragma once

namespace interleaver {

template <typename T>
std::vector<T> Interleaver::Handle::interleave(const std::vector<T>& x) const {
  if (idx_in.empty()) return x;
  if (x.size() != idx_in.size())
    throw std::invalid_argument("interleave: input size mismatch");
  std::vector<T> y(idx_in.size());
  for (std::size_t i = 0; i < idx_in.size(); ++i)
    y[i] = x[static_cast<std::size_t>(idx_in[i])];
  return y;
}

template <typename T>
std::vector<T> Interleaver::Handle::deinterleave(const std::vector<T>& y) const {
  if (idx_out.empty()) return y;
  if (y.size() != idx_out.size())
    throw std::invalid_argument("deinterleave: input size mismatch");
  std::vector<T> x(idx_out.size());
  for (std::size_t i = 0; i < idx_out.size(); ++i)
    x[i] = y[static_cast<std::size_t>(idx_out[i])];
  return x;
}

template <typename T>
std::vector<T> Interleaver::Handle::interleave_chunks(const std::vector<T>& v) const {
  return interleave(v);
}

template <typename T>
std::vector<T> Interleaver::Handle::deinterleave_chunks(const std::vector<T>& v) const {
  return deinterleave(v);
}

} // namespace interleaver
