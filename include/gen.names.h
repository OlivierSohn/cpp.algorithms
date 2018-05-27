
namespace imajuscule {
    class shortIds {
    public:
        static std::string generate();
    private:
        static std::string toString();
        static void increment();
        static int count;
        static char c;
    };
    
    template<char T> struct ascii;
    // with 'enum { value = ... }', the value at runtime was always 0 (don't know why) so I use constexpr instead
    template<> struct ascii<'0'> { static constexpr unsigned char value = 48; };
    template<> struct ascii<'9'> { static constexpr unsigned char value = 57; };
    template<> struct ascii<'A'> { static constexpr unsigned char value = 65; };
    template<> struct ascii<'Z'> { static constexpr unsigned char value = 90; };
    template<> struct ascii<'_'> { static constexpr unsigned char value = 95; };
    template<> struct ascii<'a'> { static constexpr unsigned char value = 97; };
    template<> struct ascii<'z'> { static constexpr unsigned char value = 122; };

    template <char T> constexpr unsigned int ord = ascii<T>::value;
    
    enum Case { UpperCase, LowerCase, AnyCase };

    template<Case C>
    struct charNameIterator {
        
        static unsigned char first() {
            return ord<'0'>;
        }
        static unsigned char last() {
            return (C!=UpperCase) ? ord<'z'> : ord<'_'>;
        }
        
        static char next(unsigned char v) {
            Assert(v >= first());
            Assert(v <= last());
            if(v == ord<'9'>) {
                return (C == LowerCase) ? ord<'_'> : ord<'A'>;
            }
            if((C!=LowerCase) && (v == ord<'Z'>)) {
                return ord<'_'>;
            }
            if(v == ord<'_'>) {
                if(C!=UpperCase) {
                    return ord<'a'>;
                }
                return ord<'0'>;
            }
            if((C!=UpperCase) && (v == ord<'z'>)) {
                return ord<'0'>;
            }
            return v+1;
        }
    };

    /*
     * Generates unique names with a particular case.
     */
    template<Case C = UpperCase>
    struct UniqueNames {
        using char_gen = charNameIterator<C>;

        UniqueNames()
        : name(1, safe_cast<char>(char_gen::first()))
        {}
        
        static UniqueNames & getInstance() {
            static UniqueNames i;
            return i;
        }
        
        std::string next() {
            auto str = name;
            increment();
            return str;
        }
        
    private:
        void increment(){
            for(auto it = name.rbegin(), end = name.rend(); it!=end; ++it) {
                auto & c = *it;
                c = char_gen::next(c);
                if(c != char_gen::first()) {
                    return;
                }
            }
            name.assign(name.size()+1, char_gen::first());
        }
        std::string name;
    };
}
