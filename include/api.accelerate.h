
namespace imajuscule {
    namespace accelerate {

        template<typename T>
        struct API;
        
        template<>
        struct API<float> {
            using T_DSPComplex      = DSPComplex;
            using T_DSPSplitComplex = DSPSplitComplex;
            
            using T_FFTSetup        = FFTSetup;
            
            constexpr static auto f_create_fftsetup  = vDSP_create_fftsetup;
            constexpr static auto f_destroy_fftsetup = vDSP_destroy_fftsetup;

            constexpr static auto f_zvfill   = vDSP_zvfill;
            constexpr static auto f_vmul     = vDSP_vmul;
            constexpr static auto f_vasm     = vDSP_vasm;
            constexpr static auto f_mmov     = vDSP_mmov;
            constexpr static auto f_vcpy     = cblas_scopy;
            constexpr static auto f_zvmul    = vDSP_zvmul;
            constexpr static auto f_zvma     = vDSP_zvma;
            constexpr static auto f_ctoz     = vDSP_ctoz;
            constexpr static auto f_ztoc     = vDSP_ztoc;
            constexpr static auto f_fft_zrip = vDSP_fft_zrip;
        };
        
        template<>
        struct API<double> {
            using T_DSPComplex      = DSPDoubleComplex;
            using T_DSPSplitComplex = DSPDoubleSplitComplex;

            using T_FFTSetup        = FFTSetupD;
            
            constexpr static auto f_create_fftsetup  = vDSP_create_fftsetupD;
            constexpr static auto f_destroy_fftsetup = vDSP_destroy_fftsetupD;

            constexpr static auto f_zvfill   = vDSP_zvfillD;
            constexpr static auto f_vmul     = vDSP_vmulD;
            constexpr static auto f_vasm     = vDSP_vasmD;
            constexpr static auto f_mmov     = vDSP_mmovD;
            constexpr static auto f_vcpy     = cblas_dcopy;
            constexpr static auto f_zvmul    = vDSP_zvmulD;
            constexpr static auto f_zvma     = vDSP_zvmaD;
            constexpr static auto f_ctoz     = vDSP_ctozD;
            constexpr static auto f_ztoc     = vDSP_ztocD;
            constexpr static auto f_fft_zrip = vDSP_fft_zripD;
        };

        template<typename T>
        using SplitComplex = typename API<T>::T_DSPSplitComplex;
        
        template<typename T>
        using Complex = typename API<T>::T_DSPComplex;

        template<typename T>
        using FFTSetup_ = typename API<T>::T_FFTSetup;
        
        template<typename T, typename ...Args>
        auto create_fftsetup(Args&&... args) {
            return API<T>::f_create_fftsetup(std::forward<Args>(args)...);
        }

        template<typename T, typename ...Args>
        auto destroy_fftsetup(Args&&... args) {
            return API<T>::f_destroy_fftsetup(std::forward<Args>(args)...);
        }

        template<typename T, typename ...Args>
        auto zvmul(Args&&... args) {
            return API<T>::f_zvmul(std::forward<Args>(args)...);
        }

        template<typename T, typename ...Args>
        auto vmul(Args&&... args) {
            return API<T>::f_vmul(std::forward<Args>(args)...);
        }
        
        template<typename T, typename ...Args>
        auto ctoz(Args&&... args) {
            return API<T>::f_ctoz(std::forward<Args>(args)...);
        }
        
        template<typename T, typename ...Args>
        auto ztoc(Args&&... args) {
            return API<T>::f_ztoc(std::forward<Args>(args)...);
        }
        
        template<typename T, typename ...Args>
        auto fft_zrip(Args&&... args) {
            return API<T>::f_fft_zrip(std::forward<Args>(args)...);
        }
    }
}