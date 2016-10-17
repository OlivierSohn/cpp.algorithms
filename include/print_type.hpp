#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <cxxabi.h>

namespace imj {

#define TYPE_TO_STR(x,y) do { typedef void(*T)(x); auto tts = type_to_string<T>(); y = tts(T(), #x); } while(0)

#define COUT_TYPE(x) std::string str_##x; TYPE_TO_STR(x, str_##x); std::cout << str_##x;

    
    template<typename T>
    struct type_to_string
    {
        template<typename U>
        std::string operator()(void(*)(U), const std::string& str)
        {
            std::string ret;
            int status;
            auto name = abi::__cxa_demangle(typeid(U).name(), 0, 0, &status);
            if(status == 0) {
                ret = name ? name : typeid(U).name();
            } else {
                ret = std::string("!!! abi::__cxa_demangle error status ") + std::to_string(status);
            }
            free(name);
            return ret;
        }
    };
}
