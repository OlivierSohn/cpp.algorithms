
namespace imajuscule::audio {

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

template<typename T, typename U>
struct pointed<std::pair<T, U>> {
    using type = U;
    static auto const & get(std::pair<T, U> const & t) {
        return t.second;
    }
};

template<typename C, typename ...Args>
double epsilonOfNaiveSummation(C const & cont, Args&&... args) {
    using Pointed = pointed<typename C::value_type>;
    using FPT = typename std::remove_reference_t<typename Pointed::type>::FPT;
    double err = {};
    // inner epsilons:
    for(auto const & c : cont) {
        err += Pointed::get(c).getEpsilon(std::forward<Args>(args)...);
    }
    // and there are cont.size() additions :
    err += cont.size() * std::numeric_limits<FPT>::epsilon();
    return err;
}

}
