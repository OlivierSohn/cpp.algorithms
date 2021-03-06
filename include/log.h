
#if defined ( WIN32 )
#  define __func__ __FUNCTION__
#endif

namespace imajuscule {

#ifdef IMJ_LOG_MEMORY
  // initial value is true.
  extern bool & shouldLogMemory();
  // initial value is false.
  enum class ThreadNature {
    Normal,
    RealTime_OS, // the audio callback thread, outside the callback
    RealTime_Program // the audio callback thread, inside the callback.
  };
  extern ThreadNature & threadNature();
#endif

    struct ScopedNoMemoryLogs {
#ifdef IMJ_LOG_MEMORY
      ScopedNoMemoryLogs() : backup(shouldLogMemory())
      {
        shouldLogMemory() = false;
      }
      ~ScopedNoMemoryLogs() {
        shouldLogMemory() = backup;
      }
    private:
      bool backup;
#endif
    };

#ifdef IMJ_LOG_MEMORY
  /*
  * Executes the passed callable only if memory logging is active.
  */
  template <typename Log>
  void logMemory(Log l) {
    using namespace imajuscule;
    if(!shouldLogMemory()) {
      return;
    }
    ScopedNoMemoryLogs s; // do not log the memory allocated / deallocated for the log itself.

    if(threadNature() == ThreadNature::RealTime_OS) {
      printf("** dynamic memory allocated / deallocated in a realtime thread, outside program reach:\n");
    }
    else if(threadNature() == ThreadNature::RealTime_Program) {
      printf("******* Warning: dynamic memory allocated / deallocated in a realtime thread: *********\n");
      logStack();
    }
    else {
      /*
      using namespace std::chrono;
      static auto startTime = high_resolution_clock::now();
      auto timeNow = high_resolution_clock::now();
      if(duration_cast<milliseconds>(timeNow-startTime).count() > 10 * 1000) {
        logStack();
      }
       */
    }
    l();
    if(threadNature() == ThreadNature::RealTime_Program) {
      printf("******* ******* ********\n");
    }
  }
#endif

    typedef enum logLevel
    {
        SCRIPT = 1,
        INFO,
        WARN,
        ERR = 0
    }logLevel;

#ifdef NO_LOGS
# define LG(...)
#else
    void LG(logLevel, /*const char* sModule,*/ const char * format, ...);
#endif

  template <class T>
    void logCoords(const char * message, const T & coords) {
        LG(INFO, "%s %.3f %.3f %.3f", message, coords[0], coords[1], coords[2]);
    }

    template<typename T>
    void print_time(std::chrono::time_point<T> time, std::ostream & os) {
        using namespace std;
        using namespace std::chrono;

        time_t curr_time = T::to_time_t(time);
        char sRep[100];
        // if needed use %Y-%m-%d for year / month / date
        strftime(sRep,sizeof(sRep),"%H:%M:%S",localtime(&curr_time));

        typename T::duration since_epoch = time.time_since_epoch();
        seconds s = duration_cast<seconds>(since_epoch);
        since_epoch -= s;
        milliseconds milli = duration_cast<milliseconds>(since_epoch);

        os << sRep << ":";
        auto c = milli.count();
        if(c < 100) {
            os << "0";
        }
        if(c < 10) {
            os << "0";
        }
        os << c;
    }

    static inline void print_system_time(std::ostream & os) {
        print_time(std::chrono::system_clock::now(), os);
    }
}
