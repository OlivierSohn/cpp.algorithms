
namespace imajuscule
{
    enum asciiAlphaLower { Min = ord<'a'>, Max = ord<'z'> };
    
    int shortIds::count = 0;
    char shortIds::c = asciiAlphaLower::Min;
    
    std::string shortIds::generate() {
        std::string ret = toString();
        increment();
        return ret;
    }
    std::string shortIds::toString() {
        std::string ret(1,c);
        if(count) {
            ret += std::to_string(count);
        }
        return ret;
    }
    void shortIds::increment() {
        if(c == asciiAlphaLower::Max) {
            c = asciiAlphaLower::Min;
            count ++;
        } else {
            c++;
        }
    }
}
