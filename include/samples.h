namespace imajuscule {
    namespace audio {
        
        
        template<int n_bytes>
        struct Uint8array_to_int32;
        
        template<>
        struct Uint8array_to_int32<1> {
            static int32_t run(uint8_t d[1]) {
                return (d[0] << 24) >> 24;
            }
        };
        
        template<>
        struct Uint8array_to_int32<2> {
            static int32_t run(uint8_t d[2]) {
                return ((d[1] << 24) | (d[0] << 16)) >> 16;
            }
        };
        
        template<>
        struct Uint8array_to_int32<3> {
            static int32_t run(uint8_t d[3]) {
                return ((d[2] << 24) | (d[1] << 16) | (d[0] << 8)) >> 8;
            }
        };
        
        template<>
        struct Uint8array_to_int32<4> {
            static int32_t run(uint8_t d[4]) {
                return ((d[3] << 24) | (d[2] << 16) | (d[1] << 8) | d[0]);
            }
        };
        
        template<int n_bytes>
        constexpr int32_t uint8array_to_int32(uint8_t array[n_bytes]) {
            return Uint8array_to_int32<n_bytes>::run(array);
        }
        
        struct int24_t {
            static constexpr auto n_bytes = 3;
            static constexpr int64_t n_different_values = pow2<8 * n_bytes>();
            static constexpr int32_t m = -n_different_values/2;
            static constexpr int32_t M = n_different_values/2 - 1;

            int24_t() = default;
            
            void fromInt32(int32_t i) {
                a[0] = i;
                a[1] = i >> 8;
                a[2] = i >> 16;
                assert(uint8array_to_int32<n_bytes>(a.data()) == i);
            }

            int24_t(float flt) {
                fromInt32(flt);
            }
            
            std::array<uint8_t, 3> a;
            operator int32_t() const {
                int32_t v = ((a[2] << 24) | (a[1] << 16) | (a[0] << 8)) >> 8;
                return v;
            }
        };
        
        template<typename T>
        struct NumericLimits {
            static constexpr auto min() { return std::numeric_limits<T>::min(); }
            static constexpr auto max() { return std::numeric_limits<T>::max(); }
        };
        
        template<>
        struct NumericLimits<int24_t> {

            static constexpr int32_t min() { return int24_t::m; }
            static constexpr int32_t max() { return int24_t::M; }
        };
        
        template<
        typename SIGNED,
        typename FLT,
        decltype(NumericLimits<SIGNED>::max()) M = NumericLimits<SIGNED>::max(),
        decltype(NumericLimits<SIGNED>::min()) m = NumericLimits<SIGNED>::min()
        >
        static SIGNED float_to_signed(FLT flt) {
            static_assert(std::is_floating_point<FLT>());

            if(flt > NumTraits<FLT>::zero()) {
                if(flt > NumTraits<FLT>::one()) {
                    return static_cast<SIGNED>(M);
                }
                constexpr auto mult = NumTraits<FLT>::half() + M;
                return static_cast<SIGNED>(flt * mult);
            }
            else {
                if(flt < -NumTraits<FLT>::one()) {
                    return static_cast<SIGNED>(m);
                }
                constexpr auto mult = NumTraits<FLT>::half() - m;
                return static_cast<SIGNED>(flt * mult);
            }
        }
        
        template<
        typename FLT,
        typename SIGNED,
        decltype(NumericLimits<SIGNED>::max()) M = NumericLimits<SIGNED>::max(),
        decltype(NumericLimits<SIGNED>::min()) m = NumericLimits<SIGNED>::min()
        >
        static auto signed_to_float(SIGNED s) -> FLT {
            static_assert(std::is_floating_point<FLT>());

            if(s > 0) {
                return s / static_cast<FLT>(M);
            }
            else {
                return -s / static_cast<FLT>(m);
            }
        }

        /*
         * Denormals can appear in reverb algorithm, when signal becomes close to 0.
         */
        void disableDenormals();

    } // NS audio
} // NS imajuscule
