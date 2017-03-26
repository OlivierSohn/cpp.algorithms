
namespace imajuscule
{
    struct ScopedLog {
        ScopedLog(const char * action, const char * str) : action(action) {
            std::cout << "--> " << action << " " << str << " ... " << std::endl;
        }
        ~ScopedLog() {std::cout << "--> " << action << " Done" << std::endl;}
    private:
        const char * action;
    };
    
    struct ScopedFileWrite {
        ScopedFileWrite(std::string const & str) :
        file(str.c_str()),
        log("Writing", str.c_str())
        {}
        
        ~ScopedFileWrite() { file.close(); }
        
        template<typename T>
        std::ostream & operator << (T v) {
            file << v;
            return file;
        }
        
        template<typename T>
        void write_vec(T const & vec, const char * name) {
            file << name << " = [";
            bool first = true;
            for(auto v : vec) {
                if(first) {
                    first = false;
                }
                file << std::endl;
                file << v;
            }
            file << std::endl << "]'" << std::endl << std::endl;
        }
        
    private:
        ScopedLog log;
        std::ofstream file;
    };
    
    template<typename T>
    class ScopedIndent {
    public:
        ScopedIndent(T & ref) : ref(ref) {
            ref.addIndent();
        }
        ~ScopedIndent() {
            ref.removeIndent();
        }
    private:
        T & ref;
    };
}