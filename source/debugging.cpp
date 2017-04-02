
namespace imajuscule {
    namespace debugging {

        static void simplifySymbol(std::string & subTrace)
        {
            std::string s = "imajuscule::";
            
            while(1)
            {
                auto i = subTrace.find(s);
                if (i == std::string::npos) {
                    break;
                }
                subTrace.erase(i, s.length());
            };
        }
        

#ifndef _WIN32
        std::vector<std::string> getProgramStack(int n_removed)
        {
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
            
            std::vector<std::string> all_traces;
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
            
            return std::move(all_traces);
        }
#endif
    } // ns
} // ns
