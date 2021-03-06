/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */


// inspired by https://gcc.gnu.org/bugzilla/show_bug.cgi?id=79518
#if defined(__clang__)
# define AssumeAligned64Core(p) (((uintptr_t(p) % (64)) == 0) ? (p) : (LLVM_BUILTIN_UNREACHABLE, (p)))
#elif defined(_MSC_VER)
# define AssumeAligned64Core(p) __assume((((char*) p) - ((char*) 0)) % (64) == 0)
#elif defined(__INTEL_COMPILER)
# define AssumeAligned64Core(p) __assume_aligned(p, 64)
#elif defined(__GLIBCXX__)
# define AssumeAligned64Core(p) p = __builtin_assume_aligned(p, 64)
#endif

#define AssumeAligned64(p) Assert((uintptr_t(p) % (64)) == 0); AssumeAligned64Core(p);

namespace imajuscule {

template<size_t alignment, typename T>
struct alignas(alignment) Aligned {
    T value;
};

    constexpr auto cache_line_n_bytes = 64;

    enum class Alignment {
        Normal = sizeof(void*),
        SSE    = 16,
        AVX    = 32,
        CACHE_LINE = cache_line_n_bytes,
        PAGE   = 4096
    };

    /*
     * AlignedAlocator allocates aligned memory
     */
    template <typename T, Alignment Align>
    class AlignedAllocator;


    template <Alignment Align>
    class AlignedAllocator<void, Align>
    {
    public:
        typedef void*             pointer;
        typedef const void*       const_pointer;
        typedef void              value_type;

        template <class U> struct rebind { typedef AlignedAllocator<U, Align> other; };
    };

    namespace detail {

        // c++ way to allocate over-aligned memory is wastefull (currently).
        // so we use platformspecific functions.

        static inline void* allocate_aligned_memory(size_t align, size_t size) {
            assert(align >= sizeof(void*));
            assert(is_power_of_two(align));

            if (size == 0) {
                return nullptr;
            }

          void* ptr;
#ifdef _WIN32  // posix_memalign replacement
            _aligned_malloc(size,align);
#else
            if (int rc = posix_memalign(&ptr, align, size)) {
                ptr = nullptr;
            }
#endif

#ifdef IMJ_LOG_MEMORY
            std::printf("+++ %p, size = %zu, alignment = %lu\n",ptr, size, align);
#endif
            return ptr;
        }

        static inline void deallocate_aligned_memory(void* ptr) noexcept {
#ifdef IMJ_LOG_MEMORY
          std::printf("--- %p, alignment = ?\n", ptr);
#endif

#ifdef _WIN32 // posix_memalign replacement
          return _aligned_free(ptr);
#else
          return free(ptr);
#endif
        }
    }

    template <typename T, Alignment Align>
    class AlignedAllocator
    {
        static_assert(std::alignment_of_v<T>
                      <=
                      static_cast<std::underlying_type_t<Alignment>>(Align));
    public:
        typedef T         value_type;
        typedef T*        pointer;
        typedef const T*  const_pointer;
        typedef T&        reference;
        typedef const T&  const_reference;
        typedef size_t    size_type;
        typedef ptrdiff_t difference_type;

        typedef std::true_type propagate_on_container_move_assignment;

        template <class U>
        struct rebind { typedef AlignedAllocator<U, Align> other; };

    public:
        AlignedAllocator() noexcept
        {}

        template <class U>
        AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept
        {}

        size_type
        max_size() const noexcept
        { return (size_type(~0) - size_type(Align)) / sizeof(T); }

        pointer
        address(reference x) const noexcept
        { return std::addressof(x); }

        const_pointer
        address(const_reference x) const noexcept
        { return std::addressof(x); }

        pointer
        allocate(size_type n, typename AlignedAllocator<void, Align>::const_pointer = 0)
        {
            const size_type alignment = static_cast<size_type>( Align );
            void* ptr = detail::allocate_aligned_memory(alignment , n * sizeof(T));
            if (ptr == nullptr) {
                throw std::bad_alloc();
            }

            return reinterpret_cast<pointer>(ptr);
        }

        void
        deallocate(pointer p, size_type) noexcept
        { return detail::deallocate_aligned_memory(p); }

        template <class U, class ...Args>
        void
        construct(U* p, Args&&... args)
        { ::new(reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...); }

        void
        destroy(pointer p)
        { p->~T(); }
    };

    template <typename T, Alignment TAlign, typename U, Alignment UAlign>
    inline
    bool
    operator== (const AlignedAllocator<T,TAlign>&, const AlignedAllocator<U, UAlign>&) noexcept
    { return TAlign == UAlign; }

    template <typename T, Alignment TAlign, typename U, Alignment UAlign>
    inline
    bool
    operator!= (const AlignedAllocator<T,TAlign>&, const AlignedAllocator<U, UAlign>&) noexcept
    { return TAlign != UAlign; }

} // NS
