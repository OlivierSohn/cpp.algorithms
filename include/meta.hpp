/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#define NON_COPYABLE_NOR_MOVABLE(T) \
T(T const &) = delete; \
T(T &&) = delete; \
T& operator=(T const &) = delete; \
T& operator=(T &&) = delete;

namespace imajuscule {

  //////////////////////////////////////////////////////////////////////////////
  // fold
  //////////////////////////////////////////////////////////////////////////////

  template <typename Reducer,
            template <typename> typename Accessor>
  constexpr auto fold(typename Reducer::value_type m) {
    return m;
  }

  template <typename Reducer,
            template <typename> typename Accessor,
            typename T1>
  constexpr auto fold(typename Reducer::value_type m) {
    return Reducer{}(m, Accessor<T1>::value());
  }

  template<typename Reducer,
           template <typename> typename Accessor,
           typename T1,
           typename T2,
           typename...TNs>
  constexpr auto fold(typename Reducer::value_type m) {
    return fold<Reducer, Accessor, T2, TNs...>(fold<Reducer, Accessor, T1>(m));
  }

  template <template <typename> typename Accessor, typename T1, typename ...TNs>
  constexpr auto minValue() {
    using V = decltype(Accessor<T1>::value());
    struct Reducer {
      using value_type = V;
      constexpr V operator () (V a, V b) {
        return std::min(a, b);
      }
    };
    return fold<Reducer, Accessor, T1, TNs...>(std::numeric_limits<V>::max());
  }

  template <template <typename> typename Accessor, typename T1, typename ...TNs>
  constexpr auto maxValue() {
    using V = decltype(Accessor<T1>::value());
    struct Reducer {
      using value_type = V;
      constexpr V operator () (V a, V b) {
        return std::max(a, b);
      }
    };
    return fold<Reducer, Accessor, T1, TNs...>(std::numeric_limits<V>::min());
  }

  template <template <typename> typename Accessor, typename T1, typename ...TNs>
  constexpr auto sumValues() {
    using V = decltype(Accessor<T1>::value());
    struct Reducer {
      using value_type = V;
      constexpr V operator () (V a, V b) {
        return a + b;
      }
    };
    return fold<Reducer, Accessor, T1, TNs...>(0);
  }

  template <template <typename> typename Accessor, typename T1, typename ...TNs>
  constexpr auto multiplyValues() {
    using V = decltype(Accessor<T1>::value());
    struct Reducer {
      using value_type = V;
      constexpr V operator () (V a, V b) {
        return a * b;
      }
    };
    return fold<Reducer, Accessor, T1, TNs...>(1);
  }

  //////////////////////////////////////////////////////////////////////////////


  template <class... Formats, class T, size_t N, size_t... Is>
  std::tuple<Formats...> as_tuple_int(std::array<T, N> & arr,
                                      std::index_sequence<Is...>) {
    return std::make_tuple(Formats{arr[Is]}...);
  }

  template <class... Formats, class T, size_t N,
  class = std::enable_if_t<(N == sizeof...(Formats))>>
  std::tuple<Formats...> as_tuple(std::array<T, N> & arr) {
    return as_tuple_int<Formats...>(arr, std::make_index_sequence<N>{});
  }


    template<class F, class TUPLE, std::size_t...Is>
    void for_each(TUPLE const & tuple, F func, std::index_sequence<Is...>){
        using expander = int[];
        (void)expander { 0, ((void)func(std::get<Is>(tuple)), 0)... };
    }

    template<class F, class TUPLE>
    void for_each(TUPLE const & tuple, F func){
        for_each(tuple, func, std::make_index_sequence<std::tuple_size<TUPLE>::value>());
    }

    template<class F, class TUPLE, std::size_t...Is>
    void for_each_i(TUPLE const & tuple, F func, std::index_sequence<Is...>){
        using expander = int[];
        (void)expander { 0, ((void)func(Is,std::get<Is>(tuple)), 0)... };
    }

    template<class F, class TUPLE>
    void for_each_i(TUPLE const  & tuple, F func){
        for_each_i(tuple, func, std::make_index_sequence<std::tuple_size<TUPLE>::value>());
    }

    template<int N, typename... Ts>
    using NthTypeOf = typename std::tuple_element<N, std::tuple<Ts...>>::type;

    // http://www.gotw.ca/gotw/071.htm
    template<class D, class B>
    class IsDerivedFrom
    {
    private:
        class Yes { char a[1]; };
        class No { char a[10]; };

        static Yes Test( B* ); // undefined
        static No Test( ... ); // undefined

    public:
        static constexpr bool value = sizeof(Test(safe_cast<D*>(nullptr))) == sizeof(Yes) ? 1 : 0;
    };


} // NS imajuscule
