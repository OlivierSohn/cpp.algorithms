/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */


#include <cassert>
#include <vector>
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <memory>
#include <map>
#include <array>
#include <exception>
#include <locale>
#include <codecvt>

#ifdef BSONPARSER_DEBUGLOGS
# error redefinition of BSONPARSER_DEBUGLOGS
#endif
// 1 to activate debug logging while reading the file, 0 to deactivate
#define BSONPARSER_DEBUGLOGS 0


namespace platform {
    
    struct incomplete_file : public std::runtime_error {
        incomplete_file() : std::runtime_error("unexpected end of file") {}
    };
    
    struct corrupt_file : public std::runtime_error {
        corrupt_file(const char * str) : std::runtime_error(str) {}
    };
}

namespace platform {
    
    /*
     * reads a file one byte at a time
     */
    struct FileReader {
        using ITER = std::istreambuf_iterator<char>;
        
        FileReader(std::string const & file_name) :
        file(file_name.c_str(),
             std::ios::binary)
        {
            it = ITER{file};
            if(!file) {
                throw std::runtime_error("file not found : " + file_name);
            }
        }
        
#ifndef NDEBUG
        bool isDone() const { return it == ITER{}; }
#endif
        unsigned char getNext() {
#ifndef NDEBUG
            throwIfEnd();
#endif
            
            auto val = *it;
            ++it;
#if BSONPARSER_DEBUGLOGS
            std::cout << byte_index << " : " << safe_cast<int>(val) << std::endl;
            ++byte_index;
#endif
            return val;
        }
        
    private:
#if BSONPARSER_DEBUGLOGS
        int byte_index = 0;
#endif
        std::ifstream file;
        ITER it;
        
#ifndef NDEBUG
        void throwIfEnd() {
            if(isDone()) {
                throw incomplete_file();
            }
        }
#endif
    };
    
    
    /*
     * reads a file all at one
     */
    struct OnceFileReader {
        
        OnceFileReader(std::string const & file_name)
        {
            std::ifstream file{file_name.c_str(), std::ios::binary};
            if(!file) {
                throw std::runtime_error("file not found : " + file_name);
            }
            // Determine the file length
            file.seekg(0, std::ios_base::end);
            auto size=file.tellg();
            file.seekg(0, std::ios_base::beg);
            // Create a vector to store the data
            buf.resize(size);
            // Load the data
            file.read((char*) &buf[0], size);
            // Close the file
            file.close();
        }
        
        unsigned char getNext() {
            auto val = buf[i];
            ++i;
            return val;
        }
        
    private:
        
        std::vector<uint8_t> buf;
        int i = 0;
    };
    
    static inline int is_big_endian()
    {
        union {
            uint32_t i;
            char c[4];
        } u = {0x01020304};
        
        return u.c[0] == 1;
    }

    /*
     * CustomStream allows to read from a Stream, and from a buffer which has
     * a higher priority than the Stream.
     * The buffer is controlled by the user of the class through the "setNextByte" method.
     */
    template<typename Stream>
    struct CustomStream {
        CustomStream(Stream & stream) : stream(stream) {
        }
        
        unsigned char read() {
            if(!next_byte.empty()) {
                auto v = next_byte.front();
                next_byte.erase(next_byte.begin());
                return v;
            }
            return stream.getNext();
        }
        
        void setNextByte(unsigned char b) {
            next_byte.push_back(b);
        }
        
    private:
        std::vector<unsigned char> next_byte;
        Stream & stream;
    };
} // NS platform
