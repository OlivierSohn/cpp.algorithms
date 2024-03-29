/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    namespace debugging {

#ifndef _WIN32 // dbg stack
        std::vector<std::string> getProgramStack(int n_removed)
        {
            std::vector<std::string> all_traces;
#ifdef IMJ_WITH_DEBUG_STACK
            constexpr auto max_n_frames = 100;

            void* addrlist[ max_n_frames + 1 ] = {};

            int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));
            if(!addrlen)
            {
                std::cerr << "logStack: backtrace_symbols returned nullptr" << std::endl;
                return {{"error"}};
            }

            if( addrlen == 0 ) {
                std::cerr << "logStack: trace empty, possibly corrupt" << std::endl;
                return {{"error"}};
            }

            all_traces.reserve(addrlen);

            // resolve addresses into strings containing "filename(function+address)",
            // this array must be free()-ed
            char** symbollist = backtrace_symbols(addrlist, addrlen);

            // iterate over the returned symbol lines. skip the first ones
            for (int i = n_removed + 1; i < addrlen; i++)
            {
                std::string trace(symbollist[i]);

                // parse trace to find symbol, unmangle it
                auto loc_ox = trace.find("0x");
                if(loc_ox != std::string::npos)
                {
                    auto loc_space = trace.find(" ", loc_ox);
                    if(loc_space != std::string::npos)
                    {
                        loc_space++;
                        auto loc_space2 = trace.find(" ", loc_space);
                        if(loc_space2 != std::string::npos)
                        {
                            auto subTrace = trace.substr(loc_space, loc_space2 - loc_space);
                            int demangleStatus;

                            std::unique_ptr<char, decltype(free)*> res {
                                abi::__cxa_demangle(subTrace.c_str(), nullptr, nullptr, &demangleStatus),
                                std::free
                            };

                            if( res && demangleStatus == 0)
                            {
                                subTrace = res.get();

                                simplifySymbol(subTrace);

                                trace.replace(loc_space, loc_space2 - loc_space, subTrace);
                            }
                        }
                    }
                }

                all_traces.push_back(std::move(trace));
            }

            free(symbollist);
#endif
            return std::move(all_traces);
        }
#endif
    } // ns
} // ns
