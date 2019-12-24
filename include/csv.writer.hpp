

namespace imajuscule {

struct CsvFile {
        CsvFile(std::string path) {
            file.open (path);
        }
        ~CsvFile() {
            file.close();
        }
        
        template<typename T>
        void push(T && v) {
            if(n) {
                file << ", ";
            }
            file << v;
            ++n;
        }
        void newline() {
            file << std::endl;
            file.flush();
            n = 0;
        }
    private:
        std::ofstream file;
        int n = 0;
    };
    

} // NS imajuscule

