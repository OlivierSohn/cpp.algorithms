/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{

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
      std::ofstream file;
      ScopedLog log;
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
