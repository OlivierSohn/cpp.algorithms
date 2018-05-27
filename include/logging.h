/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    template<typename T>
    void print_time(std::chrono::time_point<T> time) {
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

        cout << sRep << ":";
        auto c = milli.count();
        if(c < 100) {
            std::cout << "0";
        }
        if(c < 10) {
            std::cout << "0";
        }
        std::cout << c << "|";
    }

    static inline void print_system_time() {
        print_time(std::chrono::system_clock::now());
    }

    struct ScopedLog {
        /*
         * Pass nullptr for action do deactivate logging
         */
        ScopedLog(const char * action, const char * str) : action(action) {
            if(!action) {
                return;
            }
            std::cout << "--> " ;print_system_time(); std::cout << action << " " << str << " ... " << std::endl;
        }
        ~ScopedLog() {
            if(!action) {
                return;
            }
            std::cout << "--> " ;print_system_time(); std::cout << action << " Done" << std::endl;
        }
    private:
        const char * action;
    };

    struct ScopedFileWrite {
        ScopedFileWrite(std::string const & str, bool scopedlog = true) :
        file(str.c_str()),
        log(scopedlog?"Writing":nullptr, str.c_str())
        {}

        ~ScopedFileWrite() { file.close(); }

        template<typename T>
        std::ostream & operator << (T v) {
            file << v;
            return file;
        }

        template<typename T>
        void write_vec(T const & vec, const char * name) {
            file << name << " = [";
            bool first = true;
            for(auto v : vec) {
                if(first) {
                    first = false;
                }
                file << std::endl;
                file << v;
            }
            file << std::endl << "]'" << std::endl << std::endl;
        }

    private:
        ScopedLog log;
        std::ofstream file;
    };

    template<typename T>
    class ScopedIndent {
    public:
        ScopedIndent(T & ref) : ref(ref) {
            ref.addIndent();
        }
        ~ScopedIndent() {
            ref.removeIndent();
        }
    private:
        T & ref;
    };

}
