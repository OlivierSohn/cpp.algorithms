

#define LOG_TAG "Imj"

// TODO LOG : redirect std::out and std::err to android logging for android platform
namespace imajuscule
{
#ifdef __ANDROID__
	int toAndroid(logLevel i)
	{
		switch (i)
		{
		case INFO:
			return ANDROID_LOG_INFO;
		case WARN:
			return ANDROID_LOG_WARN;
		default:
		case ERR:
			return ANDROID_LOG_ERROR;
		}
	}
#endif


	const char * levelToChar(logLevel level);

    struct ThreadData {
        ThreadData(size_t idx) : idx(idx) {}
        std::string str;
        size_t idx;
    };
  
    std::atomic_bool & logLock() {
      static std::atomic_bool b = false;
      return b;
    }

    static void getThreadData(ThreadData *& data) {
        auto tid = pthread_self();
        static std::map<decltype(pthread_self()), ThreadData> threads;
      {
        // std::map is not thread safe
        LockGuard l(logLock());

        auto it = threads.find(tid);
        if (unlikely(it == threads.end()))
        {
          auto it = threads.emplace(tid, threads.size());
          data = &it.first->second;
        }
        else {
          data = &it->second;
        }
      }
    }

	void LG(logLevel level, /*const char* sModule,*/ const char * format, ...)
	{
        va_list args;
#ifdef __ANDROID__
		va_start(args, format);
		__android_log_vprint(toAndroid(level), LOG_TAG, format, args);
		va_end(args);
#else
		va_start(args, format);
        auto size = vsnprintf(nullptr, 0, format, args);
		va_end(args);

        ThreadData * thread_data;
        getThreadData(thread_data);

        // to use a StackVector here we would need to have
        // one "stack" pool per thread
        std::string & v = thread_data->str;
        v.resize(size+1);

        va_start(args, format);
        vsnprintf(&v[0], size+1, format, args);
        va_end(args);

        print_system_time();

        fprintf(((level == ERR) ? stderr : stdout),
                "%d|%s|%s\n",
                thread_data->idx,
                levelToChar(level),
                v.data());
#endif
	}

	const char * levelToChar(logLevel level)
	{
		switch (level)
		{
		default:
            case SCRIPT:
            return "SCRIPT";
            case INFO:
            return "INFO";
		case ERR:
			return "ERR";
		case WARN:
			return "WARN";
		}
	}
	// FOR LATER:
	/*
	#ifdef ANDROID
	#  define LOGI(...) LG(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
	#  define LOGW(...) LG(WARN, LOG_TAG, __VA_ARGS__)
	#  define LOGE(...) LG(ERR, LOG_TAG, __VA_ARGS__)
	#else
	#  define QUOTEME_(x) #x
	#  define QUOTEME(x) QUOTEME_(x)
	#  define LOGI(...) printf("I/" LOG_TAG " (" __FILE__ ":" QUOTEME(__LINE__) "): " __VA_ARGS__)
	#  define LOGE(...) printf("E/" LOG_TAG "(" ")" __VA_ARGS__)
	#endif
	*/

}
