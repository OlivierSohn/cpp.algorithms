
namespace imajuscule
{
    namespace debugging {

#ifndef _WIN32
        std::vector<std::string> getProgramStack(int n_removed = 0);
#endif
        
        template<typename T>
        struct DumpObjectsOrigins {
            void add(T*o) {
                auto s = getProgramStack(1);
                if(map.count(s)) {
                    map[s].emplace(o);
                }
                else {
                    std::set<T*> set;
                    set.emplace(o);
                    map.emplace(s, std::move(set));
                }
                writeFile();
            }
            
            void remove(T*o) {
                for(auto & s : map) {
                    if(s.second.erase(o)) {
                        writeFile();
                        return;
                    }
                }
                assert(!"an element was not found");
            }
            
        private:
            std::map<std::vector<std::string>, std::set<T*>> map;
            
            void writeFile() const {
                ScopedFileWrite f("stacks.txt", false);
                for(auto const & p : map) {
                    f << "[ " << p.second.size() << " ]" << std::endl << std::endl;
                    for(auto e : p.second) {
                        f << e;
                    }
                    f << "" << std::endl << std::endl;
                    for(auto const & s : p.first) {
                        f << s << std::endl;
                    }
                }
            }
        };        
    }
}
