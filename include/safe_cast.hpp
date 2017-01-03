
namespace imajuscule {
    
    template<
    typename T,
    typename U
    >
    struct SafeCast {
        template<
        typename T1 = T
        >
        auto
        operator()(U v) ->
        // may need to be changed to std::is_pod<typename std::remove_pointer<T>::type>::value to support cast to pod*
        typename std::enable_if<!std::is_pointer<T1>::value && std::is_pod<T>::value, T>::type {
            static_assert(!std::is_class<T>::value, "");
            return static_cast<T>(v);
        }
        
        template<
        typename T1 = T
        >
        auto
        operator()(U v) ->
        typename std::enable_if<std::is_class<T1>::value, T>::type
        {
            return static_cast<T>(v);
        }
        
        template<
        typename T1 = T
        >
        auto
        operator()(U v) ->
        typename std::enable_if<std::is_reference<T1>::value && std::is_class<typename std::remove_reference<T1>::type>::value, T>::type
        {
            return dynamic_cast<T>(v); // can throw
        }
        
        template<
        typename T1 = T
        >
        auto
        operator()(U ptr) ->
        typename std::enable_if<std::is_pointer<T1>::value && std::is_class<typename std::remove_pointer<T1>::type>::value, T>::type {
            static_assert(std::is_pointer<T>::value, "");
            assert(ptr);
            auto ret = dynamic_cast<T>(ptr);
            assert(ret);
            return ret;
        }
    };
    
    template<typename T, typename U>
    T safe_cast(U&& v) {
#ifndef NDEBUG
        SafeCast<T, U> s;
        return s(v);
#else
        return static_cast<T>(v);
#endif
    }

    template<typename T>
    T safe_cast(void* ptr) {
#ifndef NDEBUG
        assert(ptr);
#endif
        return static_cast<T>(ptr);
    }

} // NS imajuscule

