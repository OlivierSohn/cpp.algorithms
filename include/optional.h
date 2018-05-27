
namespace imajuscule
{
  template<typename T>
  using Optional =
#if __has_include(<optional>)
      std::optional<T>;
#elif __has_include(<experimental/optional>)
      std::experimental::optional<T>;
#else
#   error Must have an optional type, either from <optional> or if not supported from <experimental/optional>.
#endif

  // inspired from https://stackoverflow.com/questions/44217316/how-do-i-use-stdoptional-in-c
  template<typename T>
  T const& get_value(Optional<T>const &opt)
  {
      Assert(opt);
      return *opt;
  }
}
