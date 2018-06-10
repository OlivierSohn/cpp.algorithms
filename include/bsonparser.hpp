/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    namespace bsonparser {
        
        enum class Item {
            x00 = 0x00,
            x01,
            x02,
            x03,
            x04,
            x05,
            x06,
            x07,
            x08,
            x09,
            x0A,
            x0B,
            x0C,
            x0D,
            x0E,
            x0F,
            x10,
            x11,
            x12,
            x13,
            
            x80 = 0x80,
            x7F = 0x7F,
            xFF = 0xFF,
            
            BEGIN_MASTER_ITEMS,
            
            BYTE = BEGIN_MASTER_ITEMS,
            INT32,
            INT64,
            UINT64,
            DOUBLE,
            DECIMAL128,
            
            DOCUMENT,
            E_LIST,
            ELEMENT,
            E_NAME,
            STRING,
            CSTRING,
            BINARY,
            SUBTYPE,
            CODE_W_S,
            
            END_MASTER_ITEMS
        };
        
        constexpr auto master_rank(Item i) {
            return to_underlying(i) - to_underlying(Item::BEGIN_MASTER_ITEMS);
        }
        
        // '*' in the grammar can translate either to 'ANY_UNTIL_0x00' or to 'FROM_CONTEXT'
        enum class Quantification {
            ANY_UNTIL_0x00,
            FROM_GRAMMAR, // the grammar says exactly how many
            FROM_CONTEXT, // the number is written in the file
        };
        
        struct QuantifiedItem {
            QuantifiedItem(Item item) : QuantifiedItem(item, 1) {}
            
            QuantifiedItem(Item item, int n) :
            quantification(Quantification::FROM_GRAMMAR),
            n_items(n),
            item(item)
            {
                assert(n > 0);
            }
            
            QuantifiedItem(Item item, Quantification q) :
            quantification(q),
            n_items(-1),
            item(item)
            {
                assert(q != Quantification::FROM_GRAMMAR); // use other constructors in that case
            }
            
            auto getQuantification() const { return quantification; }
            
            Item getItem() const { return item; }
            
            int getCount() const { return n_items; }
            
        private:
            Quantification quantification;
            int n_items;
            Item item;
        };
        
        namespace {
            static inline void logIndent(int indent) {
                std::cout << std::string(indent, ' ');
            }
        }
        
        struct Object {
            virtual ~Object() = default;
            
            virtual void log(int indent) const = 0;
        };
        
        template<typename>
        struct Parser;
        
        using FileReaderT = platform::OnceFileReader;
        struct ObjectCreationFunctions {
            
            using F_ONE = std::function<std::unique_ptr<Object>(Parser<FileReaderT> &,
                                                                std::vector<std::unique_ptr<Object>>)>;
            
            using F_MANY = std::function<std::unique_ptr<Object>(int n,
                                                                 Parser<FileReaderT> &,
                                                                 std::vector<std::unique_ptr<Object>>)>;

            ObjectCreationFunctions() = default;
            ObjectCreationFunctions(F_ONE make_one) : make_one(make_one) {}
            ObjectCreationFunctions(F_ONE make_one, F_MANY make_many) : make_one(make_one), make_many(make_many) {}
            
            F_ONE make_one;
            F_MANY make_many;
        };
        
        struct ObjectSpec {
            std::vector<QuantifiedItem> sequenced_items;
            
            ObjectCreationFunctions create_object;
            
            const char * description;
        };
        
        struct Document;
        struct List : public Object {
            friend struct Document;
            
            using Elems = std::vector<std::unique_ptr<Object>>;
            
            List() = default;
            
            List(Elems elems) : elements(std::move(elems)) {}
            
            void append(std::unique_ptr<Object> elem) {
                assert(elem);
                elements.push_back(std::move(elem));
            }
            
            auto const & operator [](int i) const { return elements[i]; }
            
            int size() const {
                return safe_cast<int>(elements.size());
            }
            
            auto begin() const { return elements.begin(); }
            auto end() const { return elements.end(); }
            auto begin() { return elements.begin(); }
            auto end() { return elements.end(); }
            
            void log(int indent) const override {
                logIndent(indent);
                logWithName(indent);
            }
            
            void logWithName(int indent, const char * name = "List", const char * elt_name = nullptr) const {
                std::cout << "[" << name << " (" << size() << ")]";
                if(elt_name) {
                    std::cout << " \"" << elt_name << "\"";
                }
                std::cout << " :" << std::endl;
                indent += 4;
                for(auto const & e : *this) {
                    e->log(indent);
                }
                indent -= 4;
            }
            
        private:
            Elems elements;
        };
        
        struct BinaryData : public Object {
            enum class SubType { // the order should match the bson spec
                GENERIC,
                FUNCTION,
                BIN_OLD,
                UUID_OLD,
                UUID,
                MD5,
                USER_DEFINED
            };
            
            BinaryData(SubType subtype, std::vector<uint8_t> && data) : subtype(subtype), data(std::move(data)) {}
            
            SubType subtype;
            std::vector<uint8_t> data;
            
            void log(int indent) const override {
                std::cout << "[Binary data (subtype " << (int)subtype << " size " << data.size() << ")] ";

                auto writeEOL = [](){
                    std::cout << std::endl;
                };
                
                Periodically<decltype(writeEOL)> p{16, writeEOL};
                
                for(auto d : data) {
                    std::cout << std::hex << d;
                    p.step();
                }
            }
        };
        
        template<typename T>
        struct Vector : public Object {
            
            Vector() = default;
            
            Vector(std::vector<T> && elems) : elements(std::move(elems)) {}
            
            auto const & operator [](int i) const { return elements[i]; }
            
            int size() const {
                return safe_cast<int>(elements.size());
            }
            
            auto begin() const { return elements.begin(); }
            auto end() const { return elements.end(); }
            auto begin() { return elements.begin(); }
            auto end() { return elements.end(); }
            
            void log(int indent) const override {
                logIndent(indent);
                logWithName(indent);
            }
            
            void logWithName(int indent, const char * name = "List", const char * elt_name = nullptr) const {
                std::cout << "[" << name << " (" << size() << ")]";
                if(elt_name) {
                    std::cout << " \"" << elt_name << "\"";
                }
                std::cout << " :" << std::endl;
                indent += 4;
                for(auto const & e : *this) {
                    logIndent(indent);
                    std::cout << e << std::endl;
                }
                indent -= 4;
            }
            
            std::vector<T> & editElements() { return elements; }
        private:
            std::vector<T> elements;
        };
        
        struct Document : public List {
            Document(std::unique_ptr<List> list) : List(std::move(list->elements)) {}
        };
        
        struct Named : public Object {
            
            Named(std::string name) : name(std::move(name)) {}
            
            std::string const & getName() const { return name; }
            
            void log(int indent) const override {
                logIndent(indent);
                logName();
            }
            
            void logName() const {
                std::cout << "\"" << name << "\" : ";
            }
            
        private:
            std::string name;
        };
        
        
        struct Str : public Object {
            
            Str(std::vector<unsigned char> buffer) : str(buffer.begin(),
                                                         buffer.end() - 1) // remove \x00 (maybe it should not be done in all cases?)
            {}
            
            std::string str;
            
            void log(int indent) const override {
                logIndent(indent);
                std::cout << "[Str] \"" << str << "\" : ";
            }
        };
        using KeyName = Str;
        
        template<typename T>
        struct Native : public Named {
            Native(std::string name, T value) : Named(std::move(name)), value(value) {}
            T value;
            
            void log(int indent) const override {
                logIndent(indent);
                std::cout << "["; COUT_TYPE(T); std::cout << "] "; Named::logName();
                std::cout << value << std::endl;
            }
        };
        
        template<typename T>
        struct NativeUnnamed : public Object {
            NativeUnnamed(T value) : value(value) {}
            T value;
            
            void log(int indent) const override {
                logIndent(indent);
                COUT_TYPE(T); std::cout << " : " << value << std::endl;
            }
        };
        
        struct UTCDateTime : public Named {
            int64_t value;
            
            void log(int indent) const override {
                logIndent(indent);
                std::cout << "[time] ";
                logName();
                std::cout << value << std::endl;
            }
        };
        
        struct UTF8Str : public Named {
            UTF8Str(std::string name, std::string value) : Named(std::move(name)), value(std::move(value)) {}
            
            std::string value;
            
            void log(int indent) const override {
                logIndent(indent);
                std::cout << "[str] ";
                logName();
                std::cout << value << std::endl;
            }
        };
        
        struct EmbeddedDocument : public Named {
            EmbeddedDocument(std::string name, std::unique_ptr<Document> doc) :
            Named(std::move(name)), doc(std::move(doc)) {}
            
            auto & editDoc() { return doc; }
            
            void log(int indent) const override {
                logIndent(indent);
                std::cout << "[Embedded doc] ";
                logName();
                std::cout << std::endl;
                indent += 4;
                doc->log(indent);
                indent -= 4;
                //logIndent(indent);
                //std::cout << "End of [Embedded doc]" << std::endl;
            }
            
        private:
            std::unique_ptr<Document> doc;
        };
        
        struct Binary : public Named {

            Binary(std::string name, std::unique_ptr<BinaryData> data) :
            Named(std::move(name)),
            data(std::move(data)) {}
            
            std::unique_ptr<BinaryData> data;
            
            void log(int indent) const override {
                Named::log(indent);
                data->log(indent);
            }
        };
        
        struct Undefined : public Named {
            void log(int indent) const override {
                Named::log(indent);
                std::cout << "[Undefined]" << std::endl;
            }
        };
        
        struct ObjectId : public Named {
            using OID = std::array<unsigned char, 12>;
            
            ObjectId(std::string name, OID const & oid) : Named(std::move(name)), oid(oid) {}
            
            OID oid;
            
            void log(int indent) const override {
                logIndent(indent);
                std::cout << "[oid] ";
                logName();
                for(auto d : oid) {
                    std::cout << "0x" << std::hex << static_cast<int>(d)<< " ";
                }
                std::cout << std::endl;
            }
        };
        
        struct Array : public Named {
            Array(std::string name, std::vector<std::unique_ptr<Object>> elements) :
            Named(std::move(name)),
            list(std::move(elements)) {
            }
            
            void log(int indent) const override {
                logIndent(indent);
                list.logWithName(indent, "Array", getName().c_str());
            }
            
            List list;
        };
        
        template<typename T, typename Reader>
        std::unique_ptr<Object> read_object(Reader & reader) {
            auto val = reader.template read<T>();
            return std::unique_ptr<Object>{std::make_unique<NativeUnnamed<T>>(val).release()};
        }
        
        template<typename T, typename Reader>
        T read_value(Reader & reader) {
            return reader.template read<T>();
        }
        
        template<typename T, typename Reader>
        std::unique_ptr<Object> make_named_native(std::vector<std::unique_ptr<Object>> objs, Reader & reader) {
            return std::unique_ptr<Object>(std::make_unique<Native<T>>(std::move(safe_cast<Str*>(&*objs[1])->str), safe_cast<NativeUnnamed<T>*>(&*objs[2])->value));
        }
        
        template<typename T>
        auto named_native() {
            return [](auto & reader, auto objs){
                return make_named_native<T>(std::move(objs), reader);
            };
        }
        
        template<typename T>
        auto unnamed_native() {
            return [](auto & reader, auto objs){
                assert(objs.empty());
                return read_object<T>(reader);
            };
        }

        template<typename T, typename ...Args>
        auto make_object(Args&&... args) {
            auto underlying = std::make_unique<T>(std::forward<Args>(args)...);
            return std::unique_ptr<Object>(underlying.release());
        }
        
        template<typename T>
        auto n_unnamed_native() {
            return [](int n, auto & reader, auto objs){
                assert(objs.empty());
                std::vector<T> v;
                v.reserve(n);
                for(int i=0; i<n; ++i) {
                    v.push_back(read_value<T>(reader));
                }
                return make_object<Vector<T>>(std::move(v));
            };
        }
        
        template<typename T, typename U>
        auto getAs(std::unique_ptr<U> ptr) {
            return std::unique_ptr<T>(safe_cast<T*>(ptr.release()));
        }
        
        static inline int indexFromName(std::string const & name) {
            auto throw_invalid_index = [](auto const & name, auto const & e) {
                using namespace platform;
                throw corrupt_file((std::string("array index [") + name + "] is invalid ( " + e.what() + " )").c_str());
            };
            
            std::string::size_type sz;
            int index;
            try {
                index = std::stoi(name, &sz);
                if(sz != name.size()) {
                    throw std::invalid_argument("not all");
                }
            }
            catch(std::invalid_argument e) {
                throw_invalid_index(name, e);
            }
            catch(std::out_of_range e) {
                throw_invalid_index(name, e);
            }
            return index;
        }
        
        constexpr auto n_master_items = static_cast<int>(Item::END_MASTER_ITEMS) - static_cast<int>(Item::BEGIN_MASTER_ITEMS);
        using Grammar = std::array<std::vector<ObjectSpec>, n_master_items>;
        
        constexpr auto n_bytes_oid = 12;
        
        static Grammar computeGrammar() {
            // Not all grammar elements are implemented,
            // in case you run the parser on a new file that contains one
            // of these unimplemented grammar elements, it will throw an exception
            // mentionning which element is missing.
            
            // To make the implementation more straightforward, I transformed the way the grammar is formulated to an equivalent form:
            // I move the \x00 at the end of 'Document' to the 'List' definition
            
            std::map<Item, std::vector<ObjectSpec>> grammar = {
                {{Item::BYTE}, {{{}, {unnamed_native<unsigned char>(),n_unnamed_native<unsigned char>()}, "1 byte (8-bits)"}}},
                {{Item::INT32}, {{{}, {unnamed_native<int32_t>()}, "4 bytes (32-bit signed integer, two's complement)" }}},
                {{Item::INT64}, {{{}, {unnamed_native<int64_t>()}, "8 bytes (64-bit signed integer, two's complement)" }}},
                {{Item::UINT64}, {{{}, {unnamed_native<uint64_t>()}, "8 bytes (64-bit unsigned integer)"}}},
                {{Item::DOUBLE}, {{{}, {unnamed_native<double>()}, "8 bytes (64-bit IEEE 754-2008 binary floating point)"}}},
                {{Item::DECIMAL128}, {{{}, {unnamed_native<long double>()}, "16 bytes (128-bit IEEE 754-2008 decimal floating point)"}}},
                
                {{Item::DOCUMENT}, {{{Item::INT32, Item::E_LIST}, {[](auto & reader, auto objs){
                    auto list = getAs<List>(std::move(objs[1]));
                    return make_object<Document>(std::move(list));
                }}, "BSON Document. int32 is the total number of bytes comprising the document."}}},
                
                // this definition is recursive i.e for lists of N items, a stack of 2*N is needed
                {{Item::E_LIST}, {
                    {{Item::x00}, {[](auto & reader, auto objs) {
                        assert(0); // not used anymore, optimized to avoid recursion
                        throw std::logic_error("list parsing should occur in a non recursive way");
                        return make_object<List>();
                    }}, "An empty list"},
                    
                    {{Item::ELEMENT, Item::E_LIST}, {[](auto & reader, auto objs) {
                        assert(0); // not used anymore, optimized to avoid recursion
                        throw std::logic_error("list parsing should occur in a non recursive way");
                        auto & elem = objs[0];
                        auto list = getAs<List>(std::move(objs[1]));
                        list->append(std::move(elem));
                        return std::move(list);
                    }}, "A list containing at least one element"},
                }},
                
                {{Item::ELEMENT}, {
                    
                    {{Item::x01, Item::E_NAME, Item::DOUBLE}, {named_native<double>()}
                        , "64-bit binary floating point" },
                    
                    {{Item::x02, Item::E_NAME, Item::STRING}, {[](auto & reader, auto objs){
                        return make_object<UTF8Str>(std::move(safe_cast<Str*>(&*objs[1])->str),
                                                    std::move(safe_cast<Str*>(&*objs[2])->str));
                    }}, "UTF-8 string"},
                    
                    {{Item::x03, Item::E_NAME, Item::DOCUMENT}, {[](auto & reader, auto objs){
                        auto name = getAs<Str>(std::move(objs[1]));
                        auto doc = getAs<Document>(std::move(objs[2]));
                        return make_object<EmbeddedDocument>(std::move(name->str),
                                                             std::move(doc));
                    }}, "Embedded document"},
                    
                    {{Item::x04, Item::E_NAME, Item::DOCUMENT}, {[](auto & reader, auto objs){
                        using namespace platform;
                        auto name = getAs<Str>(std::move(objs[1]));
                        auto doc = getAs<Document>(std::move(objs[2]));
                        
                        // the objs of the array are not stored in any particular order so we use a map to sort them ...
                        std::map<int, std::unique_ptr<Object>> map;
                        for(auto & doc_item : *doc) {
                            auto const & wname = safe_cast<Named*>(&*doc_item)->getName();
                            auto i = indexFromName(wname);
                            map.emplace(i, std::move(doc_item));
                        }
                        
                        // and we can create an array of elements in the right order like so:
                        std::vector<std::unique_ptr<Object>> array;
                        array.reserve(map.size());
                        for(auto & e : map) {
                            array.push_back(std::move(e.second));
                        }
                        return make_object<Array>(std::move(name->str),
                                                  std::move(array));
                    }}, "Array"},
                    
                    {{Item::x05, Item::E_NAME, Item::BINARY}, {[](auto & reader, auto objs){
                        auto name = getAs<Str>(std::move(objs[1]));
                        return make_object<Binary>(std::move(name->str),
                                                   getAs<BinaryData>(std::move(objs[2])));
                    }}, "Binary data"},
                    
                    {{Item::x06, Item::E_NAME}, {},
                        "Undefined (value) — Deprecated"},
                    
                    {{Item::x07, Item::E_NAME, {Item::BYTE, 12}}, {[](auto & reader, auto objs){
                        std::array<unsigned char, n_bytes_oid> array;
                        constexpr auto offset = 2;
                        for(int i=offset; i<offset+n_bytes_oid; ++i) {
                            array[i-offset] = safe_cast<NativeUnnamed<unsigned char>*>(&*objs[i])->value;
                        }
                        auto name = getAs<Str>(std::move(objs[1]));
                        return make_object<ObjectId>(std::move(name->str),
                                                     array);
                    }}, "ObjectId"},
                    
                    {{Item::x08, Item::E_NAME, Item::BYTE}, {named_native<uint8_t>()},
                        "Boolean"}, // byte should be 0x00 or 0x01
                    
                    {{Item::x09, Item::E_NAME, Item::INT64}, {},
                        "UTC datetime"},
                    
                    {{Item::x0A, Item::E_NAME}, {},
                        "Null value"},
                    
                    {{Item::x0B, Item::E_NAME, {Item::CSTRING, 2}}, {},
                        "Regular expression (pattern, options)"},
                    
                    {{Item::x0C, Item::E_NAME, Item::STRING, {Item::BYTE, 12}}, {},
                        "DBPointer — Deprecated"},
                    
                    {{Item::x0D, Item::E_NAME, Item::STRING}, {},
                        "JavaScript code"},
                    
                    {{Item::x0E, Item::E_NAME, Item::STRING}, {},
                        "Symbol. Deprecated"},
                    
                    {{Item::x0F, Item::E_NAME, Item::CODE_W_S}, {},
                        "JavaScript code w/ scope"},
                    
                    {{Item::x10, Item::E_NAME, Item::INT32}, {named_native<int32_t>()},
                        "32-bit integer"},
                    
                    {{Item::x11, Item::E_NAME, Item::UINT64}, {},
                        "Timestamp"},
                    
                    {{Item::x12, Item::E_NAME, Item::INT64}, {named_native<int64_t>()},
                        "64-bit integer"},
                    
                    {{Item::x13, Item::E_NAME, Item::DECIMAL128}, {named_native<long double>()},
                        "128-bit decimal floating point"},
                    
                    {{Item::xFF, Item::E_NAME}, {},
                        "Min key"},
                    
                    {{Item::x7F, Item::E_NAME}, {},
                        "Max key"},
                }},
                
                {{Item::E_NAME}, {{{Item::CSTRING}, {[](auto & reader, auto objs){
                    return std::move(objs[0]);
                }}, "Key name"}}},
                
                {{Item::STRING}, {{{Item::INT32, {Item::BYTE, Quantification::FROM_CONTEXT}}, {[](auto & reader, auto objs){
                    auto v = getAs<Vector<uint8_t>>(std::move(objs[1]));
                    return make_object<Str>(std::move(v->editElements()));
                }}, "UTF8 String"}}},
                
                {{Item::CSTRING}, {{{{Item::BYTE, Quantification::ANY_UNTIL_0x00}}, {[](auto & reader, auto objs){
                    return std::move(objs[0]);
                }}, "Partial_UTF8 String"}}},
                
                {{Item::BINARY}, {{{Item::INT32, Item::SUBTYPE, {Item::BYTE, Quantification::FROM_CONTEXT}}, {[](auto & reader, auto objs){
                            auto subtype = (BinaryData::SubType)safe_cast<NativeUnnamed<unsigned char>*>(&*objs[1])->value;
                    auto v = getAs<Vector<uint8_t>>(std::move(objs[2]));
                            return make_object<BinaryData>(subtype, std::move(v->editElements()));
                        }}, "Binary"}}},
                
                {{Item::SUBTYPE}, {
                    {{Item::x00}, {}, "Generic binary subtype"},
                    {{Item::x01}, {}, "Function"},
                    {{Item::x02}, {}, "Binary (Old)"},
                    {{Item::x03}, {}, "UUID (Old)"},
                    {{Item::x04}, {}, "UUID"},
                    {{Item::x05}, {}, "MD5"},
                    {{Item::x80}, {[](auto & reader, auto objs) {
                        return std::move(objs[0]);
                    }}, "User defined"},
                }},
                
                {{Item::CODE_W_S}, {{{Item::INT32, Item::STRING, Item::DOCUMENT}, {},
                    "Code w/ scope"}}},
            };
            
            Grammar res;
            for(auto & elt : grammar) {
                res[master_rank(elt.first)] = std::move(elt.second);
            }
            return res;
        }
        
        static const Grammar & getGrammar() {
            static Grammar g = computeGrammar();
            return g;
        }
        
        /*
         * to prevent corrupt files leading to stack overflow
         */
        struct LimitRecursions {
            LimitRecursions(int32_t & i) : i(i) {
                ++i;
                if(i > 10000) {
                    throw platform::corrupt_file("n recursion exceeds 10000");
                }
            }
            ~LimitRecursions() { --i; }
        private:
            int32_t & i;
        };
        
        template<typename Stream>
        struct Parser {
            
            Parser(Stream & stream) : byte_source(stream) {
                assert(!platform::is_big_endian()); // todo support big endian (bson data is little-endian);
            }
            
            std::unique_ptr<Object> readDocument() {
                n_bytes_remaining = seek<int32_t>();
                return readOne(Item::DOCUMENT);
            }
            
            template<typename T>
            T read() {
                constexpr auto n_reads = sizeof(T);
                n_bytes_remaining -= n_reads;
                if(unlikely(n_bytes_remaining < 0)) {
                    throw platform::corrupt_file("reading past file end");
                }
                std::array<unsigned char, n_reads> array;
                for(auto & a : array) {
                    a = byte_source.read();
                }
                return *reinterpret_cast<T*>(&array[0]);
            }
            
            template<typename T>
            T seek() {
                std::array<unsigned char, sizeof(T)> array;
                for(auto & a : array) {
                    a = byte_source.read();
                }
                for(auto val:array) {
                    byte_source.setNextByte(val);
                }
                return *reinterpret_cast<T*>(&array[0]);
            }
            
        private:
            platform::CustomStream<Stream> byte_source;
            
            int deduceSpecIndex(Item item, int &rank)
            {
                using namespace platform;
                rank = master_rank(item);
                if(rank < 0) {
                    // it is not a master element, so it is a byte
                    rank = master_rank(Item::BYTE);
                }
                auto const & gr = getGrammar()[rank];
                if(gr.empty()) {
                    auto b = read<uint8_t>();
                    if(b != to_underlying(item)) {
                        throw corrupt_file("No match found");
                    }
                    return -1;
                }
                
                int chosen_index = 0;
                if(gr.size() > 1){
                    bool found = false;
                    
                    // from this byte we deduce which spec to choose.
                    // we seek instead of read because it belongs to a sub element
                    auto b = seek<uint8_t>();

                    for(auto const & spec : gr) {
                        assert(!spec.sequenced_items.empty());
                        auto const & quantifiedItem = spec.sequenced_items[0];
                        if(quantifiedItem.getCount() != 1) {
                            throw corrupt_file("A single item was expected here");
                        }
                        if(b == to_underlying(quantifiedItem.getItem())) {
                            found = true;
                            break;
                        }
                        ++chosen_index;
                    }
                    
                    if(!found) {
                        // we don't throw here, to support the List case :
                        --chosen_index;
                    }
                }
                return chosen_index;
            }
            
            enum class ReadType {
                VECTOR,
                ONE
            };
            
            std::unique_ptr<Object> read(Item item, int count, ReadType readType) {                
                LimitRecursions lr(recursion);
                
                using namespace platform;
                
                int rank;
                int chosen_index = deduceSpecIndex(item, rank);
                if(chosen_index<0) {
                    return {};
                }

                // 'List' is defined recursively in the grammar
                // so we handle it specially to avoid stack overflow due to recursion
                if(item == Item::E_LIST) {
                    assert(count == 1);
                    std::vector<std::unique_ptr<Object>> list_elements;
                    while(chosen_index != 0) { // 0 is end of list
                        list_elements.push_back(readOne(Item::ELEMENT));
                        chosen_index = deduceSpecIndex(Item::E_LIST, rank);
                    }
                    return std::make_unique<List>(std::move(list_elements));
                }
                auto const & gr = getGrammar()[rank];
                auto const & spec = gr[chosen_index];
                
                auto objs = readComponents(spec.sequenced_items);
                if(readType == ReadType::ONE) {
                    if(!spec.create_object.make_one) {
                        throw std::runtime_error("making one bson element with description [" + std::string{spec.description} + "] is not yet supported");
                    }
                    return spec.create_object.make_one(*this, std::move(objs));
                }
                else {
                    if(!spec.create_object.make_many) {
                        throw std::runtime_error("making many bson elements with description [" + std::string{spec.description} + "] is not yet supported");
                    }
                    return spec.create_object.make_many(count, *this, std::move(objs));
                }
            }
            
            std::unique_ptr<Object> readOne(Item item) {
                return read(item, 1, ReadType::ONE);
            }
            std::unique_ptr<Object> readVector(Item item, int count) {
                return read(item, count, ReadType::VECTOR);
            }

            auto readComponents(std::vector<QuantifiedItem> const & objs) {
                using namespace platform;
                
                bool has_size = false;
                int size = 0;
                
                std::vector<std::unique_ptr<Object>> v;
                v.reserve(objs.size());
                for(auto const &i : objs) {
                    if(i.getQuantification() == Quantification::ANY_UNTIL_0x00) {
                        assert(!has_size);
                        std::vector<uint8_t> buffer;
                        {
                            uint8_t b(0);
                            do {
                                b = read<uint8_t>();
                                buffer.push_back(b);
                            } while(b);
                        }
                        v.emplace_back(std::make_unique<Str>(std::move(buffer)));
                    }
                    else if(i.getQuantification() == Quantification::FROM_GRAMMAR) {
                        int sz = i.getCount();
                        auto needed = v.size() + sz;
                        if(v.capacity() < needed) {
                            v.reserve(needed);
                        }
                        for(auto j=0; j<sz; ++j) {
                            v.push_back(readOne(i.getItem()));
                        }
                        if(i.getItem()==Item::INT32) {
                            if(has_size) {
                                throw corrupt_file("Two sizes found");
                            }
                            has_size = true;
                            size = safe_cast<NativeUnnamed<int32_t>*>(&*v.back())->value;
                        }
                    }
                    else {
                        assert(i.getQuantification() == Quantification::FROM_CONTEXT);
                        if(!has_size) {
                            throw corrupt_file("Size not found in context");
                        }
                        v.push_back(readVector(i.getItem(), size));
                    }
                }
                return std::move(v);
            }
            
            int32_t n_bytes_remaining;
            int32_t recursion = 0;
        };
        
        static inline auto parse(const char * filepath) {
            using namespace bsonparser;
            using namespace platform;
            
            FileReaderT file_reader(filepath);
            Parser<decltype(file_reader)> parser(file_reader);
            return parser.readDocument();
        }
    } // NS bsonparser
} // NS imajuscule
