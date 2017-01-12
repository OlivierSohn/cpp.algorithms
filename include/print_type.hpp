
namespace imajuscule {

#define TYPE_TO_STR(x,y) do { typedef void(*T)(x); auto tts = type_to_string<T>(); y = tts(T()); } while(0)

#define COUT_TYPE(x) std::string str_##x; TYPE_TO_STR(x, str_##x); std::cout << str_##x;

    std::string demangle(const char * type_name, bool remove_namespace = true);
    
    template<typename T>
    struct type_to_string
    {
        template<typename U>
        std::string operator()(void(*)(U))
        {
            return demangle(typeid(U).name(), false);
        }
    };
}
