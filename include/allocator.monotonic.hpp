
namespace imajuscule {

template<typename BaseAllocator>
struct MonotonicBuffer {
    using T = typename BaseAllocator::value_type;
    
    MonotonicBuffer(MonotonicBuffer const&) = delete;
    MonotonicBuffer& operator=(MonotonicBuffer const&) = delete;
    
    MonotonicBuffer(int capacity)
    {
        reserve(capacity);
    }
    
    MonotonicBuffer() = default;
    
    void reserve(int capacity) {
        buffer.reserve(capacity);
    }
    void clear() {
        buffer.clear();
    }
    
    T * take(int sz) {
        if(unlikely(buffer.capacity() < buffer.size() + sz)) {
            throw std::runtime_error("MonotonicBuffer too small");
        }
        auto res = buffer.end();
        buffer.resize(buffer.size() + sz);
        return &*res;
    }
    
    int used() const {
        return buffer.size();
    }
    int remaining() const {
        return buffer.capacity() - buffer.size();
    }
private:
    std::vector<T, BaseAllocator> buffer;
};

template<typename BaseAllocator>
struct StackedMonotonicBuffers {
    static StackedMonotonicBuffers & getThreadLocalInstance() {
        thread_local StackedMonotonicBuffers t;
        return t;
    }
    
    StackedMonotonicBuffers(StackedMonotonicBuffers const&) = delete;
    StackedMonotonicBuffers& operator=(StackedMonotonicBuffers const&) = delete;
    
    StackedMonotonicBuffers() {
        stack.reserve(8);
    }
    
    void push(MonotonicBuffer<BaseAllocator> & b) {
        stack.push_back(&b);
    }
    
    void pop(MonotonicBuffer<BaseAllocator> const & b) {
        if(unlikely(stack.empty() || stack.back() != &b)) {
            throw std::runtime_error("MonotonicBuffer not found on top of stack");
        }
        stack.pop_back();
    }
    
    MonotonicBuffer<BaseAllocator> & top() {
        if(unlikely(stack.empty())) {
            throw std::runtime_error("stack is empty");
        }
        return *stack.back();
    }


private:
    std::vector<MonotonicBuffer<BaseAllocator>*> stack;
};

template<typename T>
struct UseMonotonicBuffer {
    UseMonotonicBuffer(UseMonotonicBuffer const&) = delete;
    UseMonotonicBuffer& operator=(UseMonotonicBuffer const&) = delete;
    
    UseMonotonicBuffer(MonotonicBuffer<T> & buffer)
    : buffer(buffer)
    {
        StackedMonotonicBuffers<T>::getThreadLocalInstance().push(buffer);
    }
    
    ~UseMonotonicBuffer()
    {
        StackedMonotonicBuffers<T>::getThreadLocalInstance().pop(buffer);
    }

private:
    MonotonicBuffer<T> & buffer;
};

template <class BaseAllocator>
class MonotonicAllocator
{
public:
    using value_type = typename BaseAllocator::value_type;
    using size_type = typename BaseAllocator::size_type;
    using difference_type = typename BaseAllocator::difference_type;
    using pointer = typename BaseAllocator::pointer;
    using const_pointer = typename BaseAllocator::const_pointer;
    using reference = typename BaseAllocator::reference;
    using const_reference = typename BaseAllocator::const_reference;
    
    template <class U> struct rebind { typedef MonotonicAllocator<U> other; };
    
    MonotonicAllocator() noexcept
    {};
    
    MonotonicAllocator( const MonotonicAllocator& other ) noexcept
    : MonotonicAllocator()
    {};
    
    MonotonicAllocator(MonotonicAllocator&&o) noexcept
    {}
    MonotonicAllocator& operator =(MonotonicAllocator&&o) noexcept
    {}
    
    template <class U>
    MonotonicAllocator(const MonotonicAllocator<U>&) noexcept
    : MonotonicAllocator()
    {};
    
    size_type max_size() const noexcept {
        return BaseAllocator().max_size();
    };
    
    pointer
    address(reference r) const noexcept
    { return BaseAllocator().address(r); };
    
    const_pointer
    address(const_reference r) const noexcept
    { return BaseAllocator().address(r); };
    
    pointer
    allocate( size_type n )
    {
        return StackedMonotonicBuffers<BaseAllocator>::getThreadLocalInstance().top().take(n);
    };
    
    void
    deallocate(pointer p, size_type n) noexcept
    {};
    
    template <class U, class... Args>
    void
    construct(U* p, Args&&... args)
    { BaseAllocator().construct(p, std::forward<Args>(args)...); }
    
    void
    destroy(pointer p)
    { BaseAllocator().destroy(p); };
};

// Return true if allocators b and a can be safely interchanged. "Safely interchanged" means that b could be
// used to deallocate storage obtained through a and vice versa.
template <class T, class T2>
bool
operator == ( const MonotonicAllocator<T>& a, const MonotonicAllocator<T2>& b)
{ return true; };

template <class T, class T2>
bool
operator != ( const MonotonicAllocator<T>& a, const MonotonicAllocator<T2>& b)
{ return !(operator==(a,b)); };


struct NothingResource {
    void clear() const {}
    void reserve(int) const {}
    int used() const { return 666; }
    int remaining() const { return 666; }
};

template<typename Allocator>
struct MemResource {
    static constexpr bool limited = false;
    using type = NothingResource;
    using use_type = NothingResource;
};

template<typename BaseAllocator>
struct MemResource<MonotonicAllocator<BaseAllocator>> {
    static constexpr bool limited = true;
    using type = MonotonicBuffer<BaseAllocator>;
    using use_type = UseMonotonicBuffer<BaseAllocator>;
};

} // NS imajuscule
