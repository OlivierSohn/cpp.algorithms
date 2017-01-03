
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
            static_assert(!std::is_pointer<U>::value, "");
            using UNoRef = typename std::remove_reference<U>::type;
            static_assert(std::is_pod<UNoRef>::value, "");
            return static_cast<T>(v);
        }
        
        template<
        typename T1 = T
        >
        auto
        operator()(U v) ->
        typename std::enable_if<std::is_class<T1>::value, T>::type
        {
            using UNoRef = typename std::remove_reference<U>::type;
            static_assert(std::is_base_of<UNoRef, T>::value || std::is_base_of<T, UNoRef>::value, "");
            return static_cast<T>(v);
        }
        
        template<
        typename T1 = T
        >
        auto
        operator()(U v) ->
        typename std::enable_if<std::is_reference<T1>::value && std::is_class<typename std::remove_reference<T1>::type>::value, T>::type
        {
            using UNoRef = typename std::remove_reference<U>::type;
            using TNoRef = typename std::remove_reference<T>::type;
            static_assert(std::is_class<UNoRef>::value, "");
            static_assert(std::is_base_of<UNoRef, TNoRef>::value || std::is_base_of<TNoRef, UNoRef>::value, "");
            return dynamic_cast<T>(v); // can throw
        }
        
        template<
        typename T1 = T
        >
        auto
        operator()(U ptr) ->
        typename std::enable_if<std::is_pointer<T1>::value && std::is_class<typename std::remove_pointer<T1>::type>::value, T>::type {
            using UPtr = typename std::remove_reference<U>::type;
            static_assert(std::is_pointer<UPtr>::value, "");
            static_assert(std::is_class<typename std::remove_pointer<UPtr>::type>::value, "");
            assert(ptr);
            auto ret = dynamic_cast<T>(ptr);
            assert(ret);
            return ret;
        }
    };
    
    template<typename T, typename U>
    T safe_cast_impl(U&& v) {
#ifdef NDEBUG
        throw;
#endif
        SafeCast<T, U> s;
        return s(v);
    }

    template<typename T>
    T safe_cast_impl(void* ptr) {
#ifdef NDEBUG
        throw;
#endif
        assert(ptr);
        return static_cast<T>(ptr);
    }
    
#ifdef safe_cast
# error "redefinition of safe_cast"
#endif
    
#ifndef NDEBUG
# define safe_cast safe_cast_impl
#else
// I use a macro to avoid overhead in release due to parameter passing
# define safe_cast static_cast
#endif
    
} // NS imajuscule

