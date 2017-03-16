
namespace imajuscule
{
    struct tm * getTime() {
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
}
