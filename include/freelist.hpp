
namespace imajuscule {

    enum class Link {
        Pointer, // for bigger types
        Index16    // for types whose sise is < sizeof(void*)
    };
    
    template<Link>
    struct LinkT;
    
    template<>
    struct LinkT<Link::Pointer> {
        using type = void*;
    };
    
    template<>
    struct LinkT<Link::Index16> {
        using type = uint16_t;
    };
    
    
    template<Link>
    struct Links;
    
    template<>
    struct Links<Link::Pointer> {
        template<typename T>
        static auto fromIndex(size_t index, T & elements ) { return &elements[index]; }

        template<typename T>
        static auto fromPtr(typename T::value_type * ptr, T const & elements ) { return ptr; }
        
        template<typename T>
        static auto follow(void * ptr, T & elements) { return ptr; }

        static constexpr void* getNull() { return nullptr; }
    };
    
    template<>
    struct Links<Link::Index16> {
        template<typename T>
        static uint16_t fromIndex(size_t index, T const & elements ) { return index; }
        
        template<typename T, typename VALUE = typename T::value_type>
        static uint16_t fromPtr( VALUE * ptr, T const & elements ) { return ptr - &elements[0]; }
        
        template<typename T>
        static auto follow(uint16_t index, T & elements) { return &elements[index]; }
        
        static constexpr uint16_t getNull() { return std::numeric_limits<uint16_t>::max(); }
    };
    
    
    template<Link>
    struct MaxSize;
    
    template<>
    struct MaxSize<Link::Pointer> {
        static constexpr size_t value = -1; // no limit
    };
    
    template<>
    struct MaxSize<Link::Index16> {
        static constexpr size_t value = std::numeric_limits<uint16_t>::max() - 1;
    };
    
    
    /*
     * T must be default-constructible
     */
    template< typename T, size_t SIZE, Link LINK_T = Link::Pointer >
    struct FreeList {
        using value_type = T;
        using link_type = typename LinkT<LINK_T>::type;
        static constexpr auto size = SIZE;
        static constexpr auto max_size = MaxSize<LINK_T>::value - 1; // -1 because one element is used to store head
        static constexpr auto null_() { return Links<LINK_T>::getNull(); };

        static_assert(sizeof(value_type) >= sizeof(link_type), "");
        
        FreeList() {
            assert(size < max_size);
            initialize();
        }
        
        value_type* Take() {
            static constexpr auto null = null_();
            auto & head_ = head();
            if(head_ == null) {
                return nullptr;
            }
            auto ret = Links<LINK_T>::follow(head_, elements);
            head_ = *reinterpret_cast<link_type*>(ret);
            auto ptr = reinterpret_cast<value_type*>(ret);
            
            //*ptr = {};

            assert(isOurs(ptr));
            return ptr;
        }
        
        void Return(T*ptr) {
            assert(isOurs(ptr));
            
            auto link = Links<LINK_T>::fromPtr(ptr, elements);
            auto & head_ = head();
            *reinterpret_cast<link_type*>(Links<LINK_T>::follow(link, elements)) = head_;
            head_ = link;
        }
        
    private:
        std::array<T, size+1> elements;
        
        link_type & head() { return reinterpret_cast<link_type&>(elements.back()); }
        
        void initialize() {
            static constexpr auto null = null_();
            head() = Links<LINK_T>::fromIndex(0, elements);
            size_t end = elements.size()-1;
            for(size_t i=1; i<end; ++i) {
                auto next = Links<LINK_T>::fromIndex(i, elements);
                elements[i-1] = reinterpret_cast<value_type&>( next );
            }
            elements[end-1] = reinterpret_cast<value_type const&>( null );
        }
        
        bool isOurs(T*ptr) const {
            return (ptr - &elements[0]) < size;
        }
    };

} // NS imajuscule

