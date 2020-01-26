/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace profiling {

        // could be optimized using vDSP_meanv
        template<typename Array, typename T = typename Array::value_type>
        float avg(Array const & a) {
            assert(!a.empty());
//            StringPlot p(20, a.size());
//            p.draw(a, default_curve_char, true);
//            p.log();
            return std::accumulate(a.begin(), a.end(), T{}).count() / static_cast<float>(a.size());
        }

        // could be optimized with vDSP_minv
        template<typename Array, typename T = typename Array::value_type>
        float min_(Array const & a) {
            assert(!a.empty());
//            StringPlot p(20, a.size());
//            p.draw(a, default_curve_char, true);
//            p.log();
            return std::min_element(a.begin(), a.end())->count();
        }


    /*
     Pollutes the caches by reading a large amount of data, in a non-sequential manner.
          
     After all calls to pollute(), it is important to use dummySum and performa a side effect
     with it (printto the console for example) else the compiler could see that pollute()
     has no side effect and optimize it away.
     */
    template<typename T>
    struct CachePolluter {
        
        CachePolluter(T & dummySum)
        : dummySum(dummySum)
        {
            // cache sizes on my machine:

            /*
             The fft (and other) costs measured after flushing caches (including L3) yield to the best approximation of real costs.
             */
            constexpr int64_t l3_cache_size = 4194304;
            constexpr int64_t l2_cache_size = 262144;
            constexpr int64_t l1_cache_size = 32768;
            
            // https://software.intel.com/en-us/forums/software-tuning-performance-optimization-platform-monitoring/topic/744272:
            //
            // " At the gross level, on Xeon E5 v3 systems, reading an array that is 4x larger than the L3 cache size
            //   will clear nearly 100% of the prior data from the L1, L2, and L3 caches."

            constexpr int64_t magic_factor = 4;

            constexpr int64_t sz = l3_cache_size * magic_factor;

            pollution.resize(sz);

            std::iota(pollution.begin(), pollution.end(),
                      0);

            // to reduce the likelyhood that the compiler optimizes the vector away:
            std::swap(pollution[0], pollution[sz-1]);
        }

        void operator()()
        {
            for(auto val : pollution) {
                if(val & 8) {
                    dummySum *= 2;
                }
                else {
                    ++dummySum;
                }
            }
        }
        
    private:
        T & dummySum;
        std::vector<uint8_t> pollution;
    };
    
    template<typename Iterator, typename T = typename Iterator::value_type>
    void loadInCache(Iterator it, Iterator end, T & dummy)
    {
        dummy = std::accumulate(it,
                                end,
                                dummy,
                                [](T acc, T val) {
            if(val & 1) {
                acc *= val;
            }
            else {
                acc += val;
            }
            return acc;
        });
    }

    // inspired from https://stackoverflow.com/questions/5919996/how-to-detect-reliably-mac-os-x-ios-linux-windows-in-c-preprocessor
    
    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
       //define something for Windows (32-bit and 64-bit, this part is common)
       #ifdef _WIN64
          //define something for Windows (64-bit only)
       #else
          //define something for Windows (32-bit only)
       #endif
    #elif __APPLE__
    
#ifndef RUSAGE_THREAD
#define RUSAGE_THREAD 6543634 // a random value
#endif
    // https://stackoverflow.com/questions/5652463/equivalent-to-rusage-thread-darwin
    static inline int my_getrusage(int v, struct rusage *rusage)
    {
        if(v != RUSAGE_THREAD) {
            return getrusage(v, rusage);
        }
        int ret = -1;
        thread_basic_info_data_t info{};
        mach_msg_type_number_t info_count = THREAD_BASIC_INFO_COUNT;
        kern_return_t kern_err;

        kern_err = thread_info(mach_thread_self(),
                               THREAD_BASIC_INFO,
                               (thread_info_t)&info,
                               &info_count);
        if (kern_err == KERN_SUCCESS) {
            memset(rusage, 0, sizeof(struct rusage));
            rusage->ru_utime.tv_sec = info.user_time.seconds;
            rusage->ru_utime.tv_usec = info.user_time.microseconds;
            rusage->ru_stime.tv_sec = info.system_time.seconds;
            rusage->ru_stime.tv_usec = info.system_time.microseconds;
            ret = 0;
        } else {
            errno = EINVAL;
        }
        return ret;
    }
#else
    int my_getrusage(int v, struct rusage *rusage)
    {
        return getrusage(v, rusage);
    }
#endif
    
    struct CpuDuration {
        CpuDuration()
        : user{}
        , kernel{}
        {}

        CpuDuration(rusage const & ru)
        : user(ru.ru_utime)
        , kernel(ru.ru_stime)
        {}

        template<typename D>
        static timeval fromDuration (D d) {
            int64_t const m = std::chrono::duration_cast<std::chrono::microseconds>(d).count();
            timeval t{};
            t.tv_sec = m / 1000000;
            t.tv_usec = m - (t.tv_sec * 1000000);
            return t;
        }
        
        static int64_t countMicroseconds(timeval t) {
            return t.tv_usec + t.tv_sec * 1000000;
        }
        
        CpuDuration operator - (CpuDuration const & o) const {
            auto d = CpuDuration();
            timersub(&user, &o.user, &d.user);
            timersub(&kernel, &o.kernel, &d.kernel);
            return d;
        }
        
        CpuDuration operator + (CpuDuration const & o) const {
            auto copy = *this;
            copy += o;
            return copy;
        }
        
        CpuDuration & operator += (CpuDuration const & o) {
            timeradd(&user, &o.user, &user);
            timeradd(&kernel, &o.kernel, &kernel);
            return *this;
        }
        
        bool operator < (CpuDuration const & o) const {
            return count() < o.count();
        }

        int64_t count() const {
            return countMicroseconds(getGlobal());
        }
        
        timeval getGlobal() const {
            timeval res{};
            timeradd(&user, &kernel, &res);
            return res;
        }
        
        timeval const & getUser() const {
            return user;
        }
        timeval const & getKernel() const {
            return kernel;
        }
    private:
        
        timeval user;
        timeval kernel;
    };
    
    struct ThreadCPUTimer {
        ThreadCPUTimer(std::optional<CpuDuration> & d)
        : d(d)
        {
            start = cpu_time_used();
        }
        
        ~ThreadCPUTimer() {
            end = cpu_time_used();
            if(start && end)
            {
                d = *end - *start;
            }
            else {
                d = {};
            }
        }
        
    private:
        std::optional<CpuDuration> start, end;
        std::optional<CpuDuration> & d;
        
        std::optional<CpuDuration> cpu_time_used() {
           rusage ru;
           int res = my_getrusage(RUSAGE_THREAD,&ru);
            if(res != 0) {
                return {};
            }
            return {ru};
        }
    };
    
    using Timer = ThreadCPUTimer;

    template<typename F>
    CpuDuration measure_thread_cpu_one(F f) {
        std::optional<CpuDuration> duration1;
        {
            // this is just to load the corresponding code in the instruction cache.
            Timer t(duration1);
        }
        
        std::optional<CpuDuration> duration2;
        {
            Timer t(duration2);
            f();
        }
        if(!duration1 || !duration2) // we use both so that the compiler doesn't optimize away the first one
        {
            throw std::runtime_error("failed to measure time interval");
        }
        return *duration2;
    }

    template<typename PREP, typename F>
    std::vector<CpuDuration> measure_thread_cpu_n(int n_warmup, int n, PREP preparation, F f) {
        assert(n > 0);
                    
        std::vector<CpuDuration> durations(n_warmup + n);

        for(auto & d : durations)
        {
            preparation();

            std::optional<CpuDuration> duration;
            {
                Timer t(duration);
                f();
            }
            if(!duration) {
                throw std::runtime_error("failed to measure time interval");
            }
            d = *duration;
        }
        return {durations.begin() + n_warmup, durations.end()};
    }
        /*
        * Measures elapsed system time between construction and destruction of the object.
        */

        template<typename Clock>
        struct WallTimer
        {
            using rep = typename Clock::rep;
            using time_point = typename Clock::time_point;
            using resolution = typename Clock::duration;

            WallTimer(resolution& duration) :
            duration(&duration) {
              startTime = Clock::now();
            }
            ~WallTimer() {
                using namespace std::chrono;
                *duration = duration_cast<resolution>(Clock::now() - startTime);
            }
        private:

            time_point startTime;
            resolution* duration;
        };

    template<typename Clock>
    struct IntervalsTimer
    {
        using rep = typename Clock::rep;
        using time_point = typename Clock::time_point;
        using resolution = typename Clock::duration;

        IntervalsTimer() :
        lastTime(Clock::now())
        {}
        
        resolution elapsedSinceLast() {
            auto prev = lastTime;
            lastTime = Clock::now();
            return lastTime - prev;
        }
    private:

        time_point lastTime;
    };
    
    }

}

static inline std::ostream& operator <<(std::ostream& out, timeval const& tv)
{
    char usfill = ' ';
    std::stringstream ss;
    constexpr int szSec = 4;
    if(true || tv.tv_sec) {
        ss << std::setw(szSec-1) << std::setfill(' ');
        ss << tv.tv_sec;
        ss << ".";
        usfill = '0';
    }
    else {
        ss << std::string(szSec, ' ');
    }
    ss << std::setw(6) << std::setfill(usfill) << tv.tv_usec;
    return out << ss.str();
}
