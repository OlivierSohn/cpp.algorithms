
namespace imajuscule {

template<typename T>
struct pointed {
    using type = T;
    static auto const & get(T const & t) {
        return t;
    }
};
template<typename T>
struct pointed<std::unique_ptr<T>> {
    using type = T;
    static auto const & get(std::unique_ptr<T> const & t) {
        return *(t.get());
    }
};

template<typename C>
double epsilonOfNaiveSummation(C const & cont) {
    using Pointed = pointed<typename C::value_type>;
    using FPT = typename std::remove_reference_t<typename Pointed::type>::FPT;
    double err = {};
    // inner epsilons:
    for(auto const & c : cont) {
        err += Pointed::get(c).getEpsilon();
    }
    // and there are cont.size() additions :
    err += cont.size() * std::numeric_limits<FPT>::epsilon();
    return err;
}

}
