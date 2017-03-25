
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
