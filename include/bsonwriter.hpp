/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
namespace bsonparser {
    
    struct BSonWriter : public WritableStorage {
        BSonWriter(DirectoryPath d, FileName f) : WritableStorage(d, f) {}

        ~BSonWriter() {
            writeByte(bsonparser::Item::x00);
            
            UpdateFileHeader();
            Finalize();
        }

        eResult Initialize() {
            return doSaveBegin();
        }

        void writeInt(int32_t i) {
            
            writeByte(bsonparser::Item::x10);
            
            writeName("");
            
            n_bytes += WriteData(&i,sizeof(int32_t),1);
        }
        
        void writeBool(bool b) {
            
            writeByte(bsonparser::Item::x08);
            
            writeName("");
            
            uint8_t i = b?1:0;
            n_bytes += WriteData(&i,sizeof(uint8_t),1);
        }
        
        void writeString(std::string const & s) {
            writeByte(bsonparser::Item::x02);

            writeName("");
            
            writeAString(s);
        }

        template<typename F>
        void writeBinary(int32_t nBytes, F generateOneByte) {
            writeByte(bsonparser::Item::x05); // binary
            
            writeName("");

            n_bytes += WriteData(&nBytes,sizeof(int32_t),1);
            
            writeByte(bsonparser::Item::x80); // user-defined

            for(int i=0; i<nBytes; ++i) {
                uint8_t const b = generateOneByte();
                n_bytes += WriteData(&b, sizeof(uint8_t), 1);
            }
        }

    protected:
        
        void DoUpdateFileHeader() override {
            WriteData(&n_bytes, sizeof(int32_t), 1);
        }
        
    private:
        static constexpr auto header_size = 4;
        int32_t n_bytes = header_size;

        void writeAString(std::string const & name) {
            int32_t sz = 1 + name.size();
            n_bytes += WriteData(&sz,sizeof(int32_t),1);
            writeCString(name);
        }
        
        void writeCString(std::string const & name) {
            n_bytes += WriteData(reinterpret_cast<void*>(const_cast<char*>(name.c_str())),
                                 name.size() + 1,
                                 1);
        }
        
        void writeName(std::string name) {
            writeCString(std::move(name));
        }
        
        void writeByte(bsonparser::Item i) {
            unsigned char c = to_underlying(i);
            n_bytes += WriteData(&c, 1, 1);
        }
    };
} // NS bsonparser
}
