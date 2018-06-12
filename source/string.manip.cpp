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

}
