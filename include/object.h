
namespace imajuscule
{
    struct Object {
        virtual ~Object() = default;
        Object() = default;
    };

    struct NonCopyable : public Object {
        NonCopyable() =default;
        NonCopyable(const NonCopyable &) = delete;
        NonCopyable & operator=(const NonCopyable&) = delete;
        
        NonCopyable(NonCopyable &&) = default;
        NonCopyable& operator = (NonCopyable &&) = default;
    };
}
