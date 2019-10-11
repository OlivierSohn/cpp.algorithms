
namespace imajuscule {
    namespace audio {

        /*
        * Denormals can appear in reverb algorithm, when signal becomes close to 0.
        */
        void disableDenormals() {


#if __has_include(<fenv.h>)
#  if defined FE_DFL_DISABLE_SSE_DENORMS_ENV
          fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);
#  endif
#  if TARGET_OS_IOS
#    if defined FE_DFL_DISABLE_DENORMS_ENV
          fesetenv(FE_DFL_DISABLE_DENORMS_ENV);
#    endif
#  endif
#endif


#if __has_include(<xmmintrin.h>)
#  if defined _MM_SET_FLUSH_ZERO_MODE && defined _MM_FLUSH_ZERO_ON
          _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#  endif
#  if defined _MM_SET_DENORMALS_ZERO_MODE && defined _MM_DENORMALS_ZERO_ON
          _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#  endif
#endif
        }


    }
}
