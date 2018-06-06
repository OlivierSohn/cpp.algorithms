/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    
#define TYPE_TO_STR(x,y) do { typedef void(*FT)(x); auto tts = type_to_string<FT>(); y = tts(FT()); } while(0)
    
#define COUT_TYPE(x) std::string str_##x; TYPE_TO_STR(x, str_##x); debugging::simplifySymbol(str_##x); std::cout << str_##x;
    
    std::string demangle(const char * type_name, bool remove_namespace = true);
  
    template<typename T>
    struct type_to_string
    {
        template<typename U>
        std::string operator()(void(*)(U))
        {
// TODO support all compilers, see GTEST_HAS_RTTI in gtest-port.h
#ifdef __GXX_RTTI
          return demangle(typeid(U).name(), false);
#else
          return "no-rtti";
#endif
        }
    };
}
