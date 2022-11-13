/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    enum class Trim {
        Yes,
        No
    };
    void split_in_lines(const std::string &s, char delim, std::vector<std::string> &elems, Trim postProcess = Trim::No);
    std::vector<std::string> split_in_lines(const std::string &s, char delim = '\n');

    std::vector<std::string> Tokenize(const std::string& str, const std::string& delimiters = " ", Trim postProcess = Trim::No);
    std::vector<std::string> TokenizeMulti(const std::string& str, const std::string& delimiter, Trim postProcess = Trim::No);

    tm * getTime();

    void FormatDate(tm*time, std::string&oDate);
    void FormatDateForComparison(std::string & date);
    void AppendTime(tm*time, std::string&str);

    static inline void WriteCurrentDate(std::string & str) {
        FormatDate(getTime(), str);
    }

    static inline void AppendTime(std::string & str) {
        AppendTime(getTime(), str);
    }

template<typename T>
std::string formatSecondsDurationWithPrecision(T length_in_seconds, int n_decimals_total) {
    auto str_sec = NumTraits<T>::to_string_with_precision(length_in_seconds, n_decimals_total);
    {
        auto i = str_sec.find('.');
        if(i == std::string::npos) {
            str_sec.push_back('.');
        }
    }
    {
        auto i = str_sec.find('.');
        int n_decimals = str_sec.size() - i - 1;
        if(n_decimals < n_decimals_total) {
            str_sec.append(n_decimals_total - n_decimals, '0');
        }
    }
    return str_sec;
}


    bool iequals(const std::string& a, const std::string& b, int nChars = -1);
    bool equals(const std::string& a, const std::string& b, int nChars = -1);

    std::string alphaNum(std::string s);

    int begins_with(std::string const& s, std::string begin);
    int ibegins_with(std::string const& s, std::string begin);

    inline bool ends_with(std::string const & value, std::string const & ending)
    {
        if (ending.size() > value.size()) {
            return false;
        }
        return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
    }

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

inline std::string justifyRight(int n, std::string s) {
    std::string res;
    if(n > static_cast<int>(s.size())) {
        res = std::string(n-s.size(), ' ');
    }
    res += s;
    return res;
}
inline std::string justifyLeft(int n, std::string s) {
    if(static_cast<int>(s.size()) < n) {
        s += std::string(n-s.size(), ' ');
    }
    return s;
}

std::string stripNewLine(std::string const & s);

    enum Correspondance : unsigned char { NOT_CORRESPONDING, CORRESPONDS_BACKWARD, CORRESPONDS_FORWARD, CORRESPONDS_ANY };

    bool findCorrespondantLocation(std::string const & text, const char c1, const int index1, const char c2, Correspondance &, int & index2);
    bool canCorrespond(const char c, char &cCorrespondant, Correspondance & correspondance);

    bool before_after(std::string & input_then_before, std::string delimiter, std::string & after);

    inline bool isACharName(char c) {
        return std::isalnum(c) || (c=='_');
    }

    bool isAName(std::string const & name);

    template<typename ... Args>
    std::string string_format(const std::string& format, Args ... args){
        char buf[1024];
        snprintf(buf, sizeof(buf), format.c_str(), args ...);
        return {buf}; // snprintf null-terminates
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

  void ReplaceStringInPlace(std::string& subject, const std::string& search,
                            const std::string& replace);

}
