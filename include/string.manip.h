
namespace imajuscule
{
    enum postProcessing {
        TRIMMED,
        NOT_TRIMMED
    };
    void split_in_lines(const std::string &s, char delim, std::vector<std::string> &elems, postProcessing pp = NOT_TRIMMED);
    std::vector<std::string> split_in_lines(const std::string &s, char delim = '\n');

    std::vector<std::string> Tokenize(const std::string& str, const std::string& delimiters = " ", postProcessing pp = NOT_TRIMMED);
    std::vector<std::string> TokenizeMulti(const std::string& str, const std::string& delimiter, postProcessing pp = NOT_TRIMMED);

    static auto getTime() {
        time_t result;
        result = time(nullptr);
        
        struct tm * pTime = nullptr;
#ifdef _WIN32
        static struct tm time;
        pTime = &time;
        localtime_s(pTime, &result);
#else
        pTime = localtime(&result);
#endif
        return pTime;
    }

    void FormatDate(tm*time, std::string&oDate);
    void FormatDateForComparison(std::string & date);
    void AppendTime(tm*time, std::string&str);

    static inline void WriteCurrentDate(std::string & str) {
        FormatDate(getTime(), str);
    }
    
    static inline void AppendTime(std::string & str) {
        AppendTime(getTime(), str);
    }
    
    bool iequals(const std::string& a, const std::string& b, int nChars = -1);
    bool equals(const std::string& a, const std::string& b, int nChars = -1);
    
    std::string alphaNum(std::string s);
    
    int begins_with(std::string const& s, std::string begin);
    int ibegins_with(std::string const& s, std::string begin);
    
    // trim from start (in place)
    int ltrim(std::string &s);
    
    // trim from end (in place)
    void rtrim(std::string &s);

    // trim from start (in place)
    int ltrim(std::string &s, char c, int maxCount = -1);
    
    // trim from end (in place)
    bool rtrim(std::string &s, char c, int maxCount = -1);
    
    // trim from both ends (in place)
    // returns the number of characters removed on the left side
    inline int trim(std::string &s) {
        rtrim(s);
        return ltrim(s);
    }
    
    // trim from start (copying)
    inline std::string ltrimmed(std::string s) {
        ltrim(s);
        return s;
    }
    
    // trim from end (copying)
    inline std::string rtrimmed(std::string s) {
        rtrim(s);
        return s;
    }
    
    // trim from both ends (copying)
    inline std::string trimmed(std::string s) {
        trim(s);
        return s;
    }
    
    enum Correspondance : unsigned char { NOT_CORRESPONDING, CORRESPONDS_BACKWARD, CORRESPONDS_FORWARD, CORRESPONDS_ANY };

    bool findCorrespondantLocation(std::string const & text, const char c1, const int index1, const char c2, Correspondance &, int & index2);
    bool canCorrespond(const char c, char &cCorrespondant, Correspondance & correspondance);
    
    void removeOutterDoubleQuotes(std::string & s);
    
    bool before_after(std::string & input_then_before, std::string delimiter, std::string & after);

    inline bool isACharName(char c) {
        return std::isalnum(c) || (c=='_');
    }
    
    bool isAName(std::string const & name);
    
    template<typename ... Args>
    std::string string_format(const std::string& format, Args ... args){
        auto size = 1 + snprintf(nullptr, 0, format.c_str(), args ...);
        StackVector<char> buf(size);
        snprintf(buf.data(), size, format.c_str(), args ...);
        return {buf.data()}; // snprintf null-terminates
    }
    
    // the number of chars per line includes the prefix
    std::string auto_format(std::string str,
                            int const chars_per_line,
                            std::string const line1prefix,
                            std::string const lineNprefix);
    
    std::string plural_of_(std::string const & text);

    std::string opposite_of_( std::string const & s );
    
    template<typename T, typename ...Args>
    T lines(Args... args) {
        return {{ (std::string(args) + '\n') }...};
    }
}
