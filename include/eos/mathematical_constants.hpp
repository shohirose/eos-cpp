#pragma once

#include <cmath>
#include <type_traits>

namespace eos {

template <typename T>
inline constexpr auto sqrtTwo() noexcept
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
#ifndef M_SQRT2
#error M_SQRT2 is not defined!
#endif
  return static_cast<T>(M_SQRT2);
}

template <typename T>
inline constexpr auto pi() noexcept
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
#ifndef M_PI
#error M_PI is not defined!
#endif
  return static_cast<T>(M_PI);
}

}  // namespace eos
