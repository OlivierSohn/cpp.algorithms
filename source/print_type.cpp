/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    std::string demangle(const char * type_name, bool remove_namespace) {
        std::string ret;
        int status;
        auto name = std::unique_ptr<char, decltype(free)*>(abi::__cxa_demangle(type_name,
                                                                               0, 0, &status),
                                                           free);
        if(status == 0) {
            if(name) {
                ret = name.get();
                if(remove_namespace) {
                    auto pos_namespace = ret.find_last_of("::");
                    if(pos_namespace != std::string::npos) {
                        ret = ret.substr(pos_namespace+1);
                    }
                }
            }
            else {
                ret = type_name;
            }
        } else {
            ret = std::string("abi::__cxa_demangle error status ") + std::to_string(status);
        }
        return ret;
    }
}
