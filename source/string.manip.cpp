namespace imajuscule
{
    struct tm * getTime() {
        time_t result;
        result = time(nullptr);

        struct tm * pTime = nullptr;
#ifdef _WIN32 // localtime_s
        thread_local struct tm time;
        pTime = &time;
        localtime_s(pTime, &result);
#else
        pTime = localtime(&result);
#endif
        return pTime;
    }

  void ReplaceStringInPlace(std::string& subject, const std::string& search,
                            const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
      subject.replace(pos, search.length(), replace);
      pos += replace.length();
    }
  }

    std::string stripNewLine(std::string const & s) {
        
        int i1 = 0;
        for(auto const & ch:s) {
            if(ch != '\n') {
                break;
            }
            ++i1;
        }

        int i2 = s.size();
        for(auto it = s.rbegin(), end=s.rend(); it != end; ++it) {
            if(*it != '\n') {
                break;
            }
            --i2;
        }
        
        std::string res;
        res.reserve(std::abs(i2-i1));
        for(auto i=i1; i<i2; ++i) {
            res.push_back(s[i]);
        }
        
        return res;
    }


}
