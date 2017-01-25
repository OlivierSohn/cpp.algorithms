
namespace imajuscule {
    
    template<class F, class...Ts, std::size_t...Is>
    void for_each(std::tuple<Ts...> & tuple, F func, std::index_sequence<Is...>){
        using expander = int[];
        (void)expander { 0, ((void)func(std::get<Is>(tuple)), 0)... };
    }
    
    template<class F, class...Ts>
    void for_each(std::tuple<Ts...> & tuple, F func){
        for_each(tuple, func, std::make_index_sequence<sizeof...(Ts)>());
    }
    
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
    
} // NS imajuscule

