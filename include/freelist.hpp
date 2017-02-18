
namespace imajuscule {
    
    namespace freelist {
        
        template<typename Integral>
        struct Links {
            template<typename T>
            static Integral fromIndex(size_t index, T const & elements ) { return index; }
            
            template<typename T, typename VALUE = typename T::value_type>
            static Integral fromPtr( VALUE * ptr, T const & elements ) { return ptr - &elements[0].value; }
            
            template<typename T>
            static auto & follow(Integral index, T & elements) { return elements[index].value; }
            
            static constexpr Integral getNull() { return std::numeric_limits<Integral>::max(); }
        };
        
        template<>
        struct Links<void*> {
            template<typename T>
            static auto fromIndex(size_t index, T & elements ) { return &elements[index]; }
            
            template<typename T>
            static auto fromPtr(decltype(std::declval<T>()[0].value) * ptr, T const & elements ) { return ptr; }
            
            template<typename T>
            static auto & follow(void * ptr, T & elements) { return *reinterpret_cast<decltype(std::declval<T>()[0].value)*>(ptr); }
            
            static constexpr void* getNull() { return nullptr; }
        };
        
        
        template<typename Integral>
        struct MaxSize {
            static constexpr size_t value = std::numeric_limits<Integral>::max() - 1;
        };

        template<>
        struct MaxSize<void*> {
            static constexpr size_t value = std::numeric_limits<uint16_t>::max(); // else an array-based freelist is probably not optimal
        };
        
        template< typename T, size_t SIZE, typename LINK_T >
        struct FreeListImpl {
            using value_type = T;
            using link_type = LINK_T;
            static constexpr auto size = SIZE;
            static constexpr auto max_size = MaxSize<LINK_T>::value - 1; // -1 because one element is used to store head
            static constexpr auto null_() { return Links<LINK_T>::getNull(); };
            
            static_assert(sizeof(value_type) >= sizeof(link_type), "");
            
            FreeListImpl() {
                assert(size < max_size);
                initialize();
            }
            
            value_type* Take() {
                constexpr auto null = null_();
                auto & head_ = head();
                if(head_.link == null) {
                    return nullptr;
                }
                auto & ret = Links<LINK_T>::follow(head_.link, elements);
                head_.value = ret;
                
                assert(isOurs(&ret));
                return &ret;
            }
            
            void Return(T*ptr) {
                assert(isOurs(ptr));
                
                auto link = Links<LINK_T>::fromPtr(ptr, elements);
                auto & head_ = head();
                auto & val = Links<LINK_T>::follow(link, elements);
                *reinterpret_cast<link_type*>(&val) = head_.link;
                head_.link = link;
            }
            
        private:
            struct Union {
                union {
                    T value;
                    LINK_T link;
                };
            };
            std::array<Union, size+1> elements;
            
            // I chose the last element as head 
            auto & head() { return elements.back(); }
            
            void initialize() {
                constexpr auto null = null_();
                head().link = Links<LINK_T>::fromIndex(0, elements);
                size_t end = elements.size()-1;
                for(size_t i=1; i<end; ++i) {
                    auto next = Links<LINK_T>::fromIndex(i, elements);
                    elements[i-1].link = next ;
                }
                elements[end-1].link = null;
            }
            
            bool isOurs(T*ptr) const {
                return (ptr - &elements[0].value) < size;
            }
        };
    }
    
    template< typename T, size_t SIZE, typename LINK_T = void* >
    using FreeList = freelist::FreeListImpl<T, SIZE, LINK_T>;
    
} // NS imajuscule

