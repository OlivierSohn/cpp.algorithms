/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */


#include <cassert>
#include <vector>
#include <exception>
#include <locale>
#include <codecvt>


namespace utf8 {
    
    constexpr unsigned char b1 = 0x01;
    constexpr unsigned char b2 = 0x02;
    constexpr unsigned char b3 = 0x04;
    constexpr unsigned char b4 = 0x08;
    constexpr unsigned char b5 = 0x10;
    constexpr unsigned char b6 = 0x20;
    constexpr unsigned char b7 = 0x40;
    constexpr unsigned char b8 = 0x80;
    
    struct Utf8Decoder {
        Utf8Decoder(std::wstring & res) : res(res) {
            buffer.reserve(4);
        }
        
        void feed(unsigned char c) {
            using namespace platform;
            
            if(end_of_string) {
                throw corrupt_file("invalid utf8 encoding (end of string in the middle)");
            }
            
            if(buffer.empty()) {
                if(0 == (c >> 7)) {
                    if(c) {
                        res.push_back(c);
                    }
                    else {
                        end_of_string = true;
                    }
                    return;
                }
                if(0x02 != (c >> 6)) {
                    throw corrupt_file("invalid utf8 encoding");
                }
                buffer.push_back(c);
                if(0x06 == (c >> 5)) {
                    ntotal = 2;
                }
                else if(0x0E == (c >> 4)) {
                    ntotal = 3;
                }
                else if(0x1E == (c >> 3)) {
                    ntotal = 4;
                }
                else {
                    throw corrupt_file("invalid utf8 encoding");
                }
            }
            else {
                assert(ntotal >= 2);
                assert(ntotal <= 4);

                buffer.push_back(c);
                
                assert(buffer.size() <= ntotal);
                
                if(buffer.size() == ntotal) {
                    // buffer now contains all the information we need to resolve the utf8 character
                    res.push_back(resolve_utf8());
                    buffer.clear();
                }
            }
        }
        
        bool hasFinished() const {
            return buffer.empty() && end_of_string;
        }
        
    private:
        std::wstring & res;
        int ntotal = 0;
        bool end_of_string = false;
        std::vector<unsigned char> buffer;
        
        /* utf8 decoding algo from http://stackoverflow.com/a/34240796 :
         Read one byte
         If the topmost bit is zero, there is nothing else to do: the code point is 0x00-0x7f
         If the topmost three bits are 110, then you need one extra byte. Take five lowest bits of the first byte, shift them left six bits and OR the lowest six bits from the second byte to get the final value
         If the topmost four bits are 1110, then you need two extra bytes. Take four lowest bits of the first one, shift by 12 bits, or in the six lowest bits from the second byte shifted by six, then finally the six lowest bits of the third byte
         If the topmost five bits are 11110, then you need three extra bytes and will read them, shift etc as previously
         If none of those conditions fit, the data is invalid
         Note that when reading extra bytes, those bytes must have 10 as the most significant bits; anything else is invalid.
         */
        
        int resolve_utf8() const {
            assert(ntotal == buffer.size());
            switch(ntotal) {
                case 2:
                    return
                    ((buffer[0] & (b1|b2|b3|b4|b5)) << 6) |
                    ( buffer[1] & (b1|b2|b3|b4|b5|b6));
                case 3:
                    return
                    ((buffer[0] & (b1|b2|b3|b4)) << 12) |
                    ((buffer[1] & (b1|b2|b3|b4|b5)) << 6 ) |
                    ( buffer[2] & (b1|b2|b3|b4|b5|b6));
                case 4:
                    return
                    ((buffer[0] & (b1|b2|b3)) << 18) |
                    ((buffer[1] & (b1|b2|b3|b4)) << 12 ) |
                    ((buffer[2] & (b1|b2|b3|b4|b5)) << 6 ) |
                    ( buffer[3] & (b1|b2|b3|b4|b5|b6));
                default:
                    throw std::runtime_error("logic error");
            }
            return 0;
        }
    };
    
    static inline auto decode(std::vector<unsigned char> buffer)
    {
        using namespace platform;
        
        std::wstring res;
        Utf8Decoder dec(res);
        for(auto c: buffer){
            dec.feed(c);
        }
        if(!dec.hasFinished()) {
            throw corrupt_file("invalid utf8 encoding");
        }
        return res;
    }
    
    static inline std::wstring s2ws(const std::string& str)
    {
        using convert_typeX = std::codecvt_utf8<wchar_t>;
        std::wstring_convert<convert_typeX, wchar_t> converterX;
        
        return converterX.from_bytes(str);
    }
    
    static inline std::string ws2s(const std::wstring& wstr)
    {
        using convert_typeX = std::codecvt_utf8<wchar_t>;
        std::wstring_convert<convert_typeX, wchar_t> converterX;
        
        return converterX.to_bytes(wstr);
    }
}
