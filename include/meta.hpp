/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#define NON_COPYABLE_NOR_MOVABLE(T) \
T(T const &) = delete; \
void operator=(T const &t) = delete; \
T(T &&) = delete;

namespace imajuscule {
    
    template<class F, class TUPLE, std::size_t...Is>
    void for_each(TUPLE & tuple, F func, std::index_sequence<Is...>){
        using expander = int[];
        (void)expander { 0, ((void)func(std::get<Is>(tuple)), 0)... };
    }
    
    template<class F, class TUPLE>
    void for_each(TUPLE & tuple, F func){
        for_each(tuple, func, std::make_index_sequence<std::tuple_size<TUPLE>::value>());
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
        static constexpr bool Is = sizeof(Test(safe_cast<D*>(nullptr))) == sizeof(Yes) ? 1 : 0;
    };
    
    template<typename F>
    struct Periodically {
        
        NON_COPYABLE_NOR_MOVABLE(Periodically);
        
        Periodically(int period, F f) : f(f), period(period) {}
        ~Periodically() {
            if(i) {
                // the last period is shorter
                f();
            }
        }
        
        void step() {
            ++i;
            if(i >= period) {
                i = 0;
                f();
            }
        }
    private:
        int i=0;
        int period;
        F f;
    };
    

} // NS imajuscule

