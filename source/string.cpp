/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    void AppendTime(tm*time, std::string&str) {
        str.reserve(str.size() + 8);

        if (!time) {
            str += "..:..:..";
            return;
        }

        auto sHour = std::to_string(time->tm_hour);
        auto sMinute = std::to_string(time->tm_min);
        auto sSecond = std::to_string(time->tm_sec);

        if (time->tm_hour < 10) {
            str += '0';
        }
        str += sHour;
        str += ':';
        if (time->tm_min < 10) {
            str += '0';
        }
        str += sMinute;
        str += ':';
        if (time->tm_sec < 10) {
            str += '0';
        }
        str += sSecond;
    }

    void FormatDate(tm*time, std::string&oDate)
    {
        if (!time) {
            oDate.assign("..../../.. ..:..:..");
        }
        else {
            int zero = 0;
            int day = time->tm_mday;
            int month = time->tm_mon + 1;
            int year = time->tm_year + 1900;
            auto sZero = std::to_string(zero);
            auto sDay = std::to_string(day);
            auto sMonth = std::to_string(month);
            auto sYear = std::to_string(year);

            oDate.append(sYear.c_str());

            oDate.append("/");

            if (month < 10)
                oDate.append(sZero.c_str());
            oDate.append(sMonth.c_str());

            oDate.append("/");

            if (day < 10)
                oDate.append(sZero.c_str());
            oDate.append(sDay.c_str());

            oDate.append(" ");

            AppendTime(time, oDate);
        }
    }

    void FormatDateForComparison(std::string & date)
    {
        const char * numbers = "0123456789";

        if (11 >= date.size()) {
            return;
        }
        if (2 != date.find_first_not_of(numbers, 0)) {
            return;
        }
        if (5 != date.find_first_not_of(numbers, 3)) {
            return;
        }
        if (10 != date.find_first_not_of(numbers, 6)) {
            return;
        }
        //date is with format "dd?mm?yyyy?....." : reverse it
        std::string newDate;
        newDate.append(date.substr(6, 4));
        newDate.append("/");
        newDate.append(date.substr(3, 2));
        newDate.append("/");
        newDate.append(date.substr(0, 2));

        newDate.append(date.substr(10));

        assert(newDate.size() == date.size());

        date.swap(newDate);
    }

    void split_in_lines(const std::string &s, char delim, std::vector<std::string> &elems, postProcessing pp) {
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            if(pp==TRIMMED) {
                trim(item);
            }
            elems.push_back(item);
        }
        if ( s.empty() ) {
            return;
        }
        // add an empty line if new line is the last character
        if(s.back()==delim) {
            elems.emplace_back("");
        }
    }
    std::vector<std::string> split_in_lines(const std::string &s, char delim) {
        std::vector<std::string> elems;
        split_in_lines(s, delim, elems);
        return elems;
    }

    inline void add(std::vector<std::string> & vec, std::string && s, postProcessing pp) {
        if(pp == postProcessing::TRIMMED) {
            trim(s);
            if(s.empty()) {
                return;
            }
        }
        vec.emplace_back(s);
    }
    std::vector<std::string> Tokenize(const std::string& str, const std::string& delimiters, postProcessing pp)
    {
        using sz = std::string::size_type;
        auto const end = std::string::npos;

        std::vector<std::string> tokens;

        sz lastPos = str.find_first_not_of(delimiters, 0);
        sz pos = str.find_first_of(delimiters, lastPos);

        while ( pos != end || lastPos != end)
        {
            // Found a token, add it to the vector.
            add(tokens, str.substr(lastPos, pos - lastPos), pp);
            lastPos = str.find_first_not_of(delimiters, pos);
            pos = str.find_first_of(delimiters, lastPos);
        }
        return tokens;
    }

    std::vector<std::string> TokenizeMulti(const std::string& str, const std::string& delimiter, postProcessing pp)
    {
        std::vector<std::string> tokens;
        if(delimiter.empty()) {
            for(auto c : str) {
                tokens.emplace_back(1, c);
            }
            return tokens;
        }

        using sz = std::string::size_type;
        auto const end = std::string::npos;

        sz lastPos = 0;
        while (str.compare(lastPos, delimiter.size(), delimiter) == 0)
        {
            // str begins with delimiter
            lastPos += delimiter.size();
        }
        if(lastPos == str.size()) {
            return {};
        }

        // lastPos is the first position where we have NOT a delimiter

        sz pos = str.find(delimiter, lastPos);
        if(lastPos == 0 && pos == end) {
            return {};
        }

        // pos is the position AFTER lastPos where we have a delimiter

        assert(pos >lastPos);
        while ( pos != end || lastPos != end)
        {
            assert(pos >lastPos);
            // Found a token, add it to the vector.
            add(tokens, str.substr(lastPos, pos - lastPos), pp);
            if(pos == end) {
                return tokens;
            }
            lastPos = pos + delimiter.size();
            while (str.compare(lastPos, delimiter.size(), delimiter) == 0)
            {
                lastPos += delimiter.size();
            }
            if(lastPos == str.size()) {
                return tokens;
            }

            pos = str.find(delimiter, lastPos);
        }
        return tokens;
    }

    bool iequals(const std::string& a, const std::string& b, int nCharacters)
    {
        int sz = (int)a.size();
        if(nCharacters >= 0) {
            if((int)b.size() < nCharacters || sz < nCharacters) {
                return false;
            }
        } else {
            if ((int)b.size() != sz) {
                return false;
            }
        }
        for (auto i = 0; (nCharacters < 0 && i < sz) || i<nCharacters; ++i) {
            if (tolower(a[i]) != tolower(b[i])) {
                return false;
            }
        }
        return true;
    }

    bool equals(const std::string& a, const std::string& b, int nCharacters)
    {
        int sz = (int)a.size();
        if(nCharacters >= 0) {
            if((int)b.size() < nCharacters || sz < nCharacters) {
                return false;
            }
        } else {
            if ((int)b.size() != sz) {
                return false;
            }
        }
        for (auto i = 0; (nCharacters < 0 && i < sz) || i<nCharacters; ++i) {
            if (a[i] != b[i]) {
                return false;
            }
        }
        return true;
    }

    bool findCorrespondantLocation(std::string const & text, const char c1, const int index1, const char c2, Correspondance & correspondance, int & index2)
    {
        if(correspondance == CORRESPONDS_ANY) { // let's resolve the correspondance
            assert(c1 == c2); // ' or "
            auto count = std::count(text.begin(), text.begin() + index1, c1);
            if(count %2) {
                // odd : c1 is the closing one
                correspondance = CORRESPONDS_BACKWARD;
                auto other = text.find_last_of(c2, index1-1);
                assert(other != std::string::npos);
                index2 = safe_cast<int>(other);
            }
            else {
                // even : c1 is the opening one
                correspondance = CORRESPONDS_FORWARD;
                auto other = text.find_first_of(c2, index1+1);
                if(other == std::string::npos) {
                    return false;
                }
                index2 = safe_cast<int>(other);
            }
            return true;
        }

        int length = (int) text.length();
        int i = index1;

        int countInBetween = 1;
        while (true)
        {
            // 1. in(de)crement iterator and break loop if it is out of bounds

            if (correspondance == CORRESPONDS_FORWARD)
            {
                ++i;
                if (i >= length)
                    break;
            }
            else
            {
                if(0==i)
                    break;
                --i;
            }

            // 2. the new iterator is valid

            char c = text[i];

            // 3. update "countInBetween" according to found char

            if(c1 == c2 && c == c1) {
                --countInBetween;
            }
            else if (c == c1) {
                ++countInBetween;
            }
            else if (c == c2) {
                --countInBetween;
            }

            if (0 == countInBetween) {
                // found it
                index2 = i;
                return true;
            }
        }

        return false;
    }

    bool canCorrespond(const char c, char &cCorrespondant, Correspondance & correspondance)
    {
        correspondance = CORRESPONDS_FORWARD;

        switch (c)
        {
            case '(':
                cCorrespondant = ')';
                return true;
            case '[':
                cCorrespondant = ']';
                return true;
            case '{':
                cCorrespondant = '}';
                return true;
            case '\'':
                cCorrespondant = '\'';
                correspondance = CORRESPONDS_ANY;
                return true;
            case '"':
                cCorrespondant = '"';
                correspondance = CORRESPONDS_ANY;
                return true;
        }

        correspondance = CORRESPONDS_BACKWARD;

        switch (c)
        {
            case ')':
                cCorrespondant = '(';
                return true;
            case ']':
                cCorrespondant = '[';
                return true;
            case '}':
                cCorrespondant = '{';
                return true;
        }

        correspondance = NOT_CORRESPONDING;
        return false;
    }

   int begins_with(std::string const& s, std::string begin) {
        auto size_comparison = (int)begin.size();
        return equals(std::move(begin), s, size_comparison) ? size_comparison : 0;
    }

    int ibegins_with(std::string const& s, std::string begin) {
        auto size_comparison = (int)begin.size();
        return iequals(std::move(begin), s, size_comparison) ? size_comparison : 0;
    }

    void rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(),
                             s.rend(),
                             [](auto c) {return !std::isspace(c);}
                             ).base(),
                s.end());
    }

    int ltrim(std::string &s) {
        auto beg = s.begin();
        auto first_non_space = std::find_if(beg, s.end(), [](char c) {return !std::isspace(c);});
        int range = safe_cast<int>(std::distance(beg, first_non_space));
        s.erase(beg, first_non_space);
        return range;
    }

    bool rtrim(std::string &s, char c, int maxCount) {
        int size = (int)s.size();
        int i = size - 1;
        int n=0;
        while(n!= maxCount && i >= 0) {
            if(!std::isspace(s[i]) && s[i] != c) {
                break;
            }
            --i;
            ++n;
        }
        if(!n) {
            return false;
        }
        s.erase(i+1, n);
        return true;
    }

    int ltrim(std::string &s, char c, int maxCount) {
        int size = (int)s.size();
        int i=0;
        while(i!= maxCount && i < size) {
            if(!std::isspace(s[i]) && s[i] != c) {
                break;
            }
            ++i;
        }
        if(!i) {
            return 0;
        }
        s.erase(0, i);
        return i;
    }

    bool before_after(std::string & input_then_before, std::string const delimiter, std::string & after)
    {
        auto v = TokenizeMulti(input_then_before, delimiter, TRIMMED);
        if(v.size() == 1) {
            if(0 == input_then_before.find(delimiter)) {
                input_then_before.clear();
                after = v[0];
            }
            else {
                input_then_before = v[0];
                after.clear();
            }
            return true;
        }
        if(v.size() != 2) {
            return false;
        }
        input_then_before = v[0];
        after = v[1];
        return true;
    }

    static bool isCharContiguous(char c)
    {
        assert(c != '\n'); // precondition
        return c != ' ';
    }

    bool isAName(std::string const & name) {
        if(name.empty()) {
            return false;
        }

        bool first = true;
        for(auto c : name) {
            if(first) {
                first = false;
                if(std::isdigit(c)) {
                    return false;
                }
            }
            if(isACharName(c)) {
                continue;
            }
            return false;
        }

        return true;
    }


    int format_one_line(int const chars_per_line, int i, std::string & str) {
        int line_beginning = i;
        i += chars_per_line-1;
        int left = -1;
        int right = -1;
        // find left candidate
        for(; i >= line_beginning; --i) {
            if(isCharContiguous(str[i])) { continue; }
            left = i;
            break;
        }
        i = line_beginning + chars_per_line;
        // find right candidate
        for(auto s = safe_cast<int>(str.size()); i < s; ++i) {
            if(isCharContiguous(str[i])) { continue; }
            right = i;
            break;
        }
        // pick best
        int ideal = line_beginning + chars_per_line;
        int chosen;
        if(left == -1) {
            if(right == -1) {
                return safe_cast<int>(str.size());
            }
            chosen = right;
        }
        else if(right == -1) {
            chosen = left;
        }
        else {
            if(ideal - left < right - ideal) {
                chosen = left;
            }
            else {
                chosen = right;
            }
        }
        assert(chosen >= 0);
        assert(chosen < safe_cast<int>(str.size()));

        if(' ' == str[chosen]) {
            str[chosen] = '\n';
        }
        else {
            str.insert(chosen, 1, '\n');
        }
        return chosen + 1;
    }

    std::string auto_format(std::string str, int const chars_per_line_total, std::string const line1prefix, std::string const lineNprefix) {
        std::replace(str.begin(), str.end(), '\n', ' ');
        int i{0};
        while(true)
        {
            assert(i == 0 || str[i-1] == '\n'); // assert that i is at the beginning of a line

            if(i == str.size()) {
                return str;
            }
            std::string const * prefix = (i==0) ? &line1prefix : &lineNprefix;

            int chars_per_line = chars_per_line_total - prefix->size();

            str.insert(str.begin() + i, prefix->begin(), prefix->end());
            i += prefix->size();
            // i is after the line prefix
            if(str.size() - i <= chars_per_line) {
                // the last part is small enough to fit on one line
                return str;
            }
            if(chars_per_line < 1 ) {
                // we cannot meet this constraint. should we use 1 for chars_per_line instead of assert + return ?
                assert(0);
                return str;
            }
            i = format_one_line(chars_per_line, i, str);
        }
    }

    std::string plural_of_(std::string const & text) {
        return text.empty() ? text : (text + "s");
    }

    std::string opposite_of_( std::string const & s ) {
        return "-(" + s + ")";
    }

    std::string alphaNum(std::string s) {
        s.erase(std::remove_if(s.begin(),
                               s.end(),
                               std::not_fn(std::function<int(int)>((int(*)(int))std::isalnum))),
                s.end());
        return s;
    }


}
