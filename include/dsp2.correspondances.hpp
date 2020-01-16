namespace imajuscule {

template<typename C>
struct corresponding_legacy_dsp {
    using type = C;
};

template<typename C>
using corresponding_legacy_dsp_t = typename corresponding_legacy_dsp<C>::type;

}
