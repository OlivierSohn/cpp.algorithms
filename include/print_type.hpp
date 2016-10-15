#pragma once

// credit: http://stackoverflow.com/questions/4484982/how-to-convert-typename-t-to-string-in-c

#include <string>
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <cxxabi.h>

namespace imj {

    #define PRINT_TYPE(x) do { typedef void(*T)(x); print_type<T>(T(), #x); } while(0)

    template<typename T>
    struct print_type
    {
        template<typename U>
        print_type(void(*)(U), const std::string& str)
        {
            std::cout << str << " = ";
            int status;
            auto name = abi::__cxa_demangle(typeid(U).name(), 0, 0, &status);
            if(status == 0) {
                std::cout << (name ? name : typeid(U).name()) << std::endl;
            } else {
                std::cout << "!!! abi::__cxa_demangle error status " << status << std::endl;
            }
            free(name);
        }
    };

}
