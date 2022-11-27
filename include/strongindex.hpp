
namespace imajuscule {

template<typename Index, typename T>
class StrongIndexVector : private std::vector<T>
{
public:
    StrongIndexVector() = default;
    explicit StrongIndexVector(typename std::vector<T>::size_type count, const T& value = T())
    : std::vector<T>(count, value)
    {}

    template< class InputIt >
    StrongIndexVector(InputIt first, InputIt last)
    : std::vector<T>(first, last)
    {}
    StrongIndexVector(std::initializer_list<T> init)
    : std::vector<T>(std::move(init))
    {}
    StrongIndexVector(StrongIndexVector const& other) = default;
    StrongIndexVector(StrongIndexVector&& other) = default;
    
    typename std::vector<T>::reference operator[]( Index pos )
    {
        return std::vector<T>::operator[](pos.unsafe_get());
    }
    
    typename std::vector<T>::const_reference operator[]( Index pos ) const
    {
        return std::vector<T>::operator[](pos.unsafe_get());
    }
    
    using typename std::vector<T>::value_type;
    using typename std::vector<T>::allocator_type;
    using typename std::vector<T>::size_type;
    using typename std::vector<T>::difference_type;
    using typename std::vector<T>::reference;
    using typename std::vector<T>::const_reference;
    using typename std::vector<T>::pointer;
    using typename std::vector<T>::const_pointer;
    using typename std::vector<T>::iterator;
    using typename std::vector<T>::const_iterator;
    using typename std::vector<T>::reverse_iterator;
    using typename std::vector<T>::const_reverse_iterator;
    
    StrongIndexVector& operator=(StrongIndexVector const& other) = default;
    StrongIndexVector& operator=(StrongIndexVector&& other) = default;
    using std::vector<T>::operator=;
    
    using std::vector<T>::assign;
    using std::vector<T>::get_allocator;
    using std::vector<T>::at;
    using std::vector<T>::front;
    using std::vector<T>::back;
    using std::vector<T>::data;
    using std::vector<T>::begin;
    using std::vector<T>::cbegin;
    using std::vector<T>::end;
    using std::vector<T>::cend;
    using std::vector<T>::rbegin;
    using std::vector<T>::crbegin;
    using std::vector<T>::rend;
    using std::vector<T>::crend;
    using std::vector<T>::empty;
    using std::vector<T>::size;
    using std::vector<T>::max_size;
    using std::vector<T>::reserve;
    using std::vector<T>::capacity;
    using std::vector<T>::shrink_to_fit;
    using std::vector<T>::clear;
    using std::vector<T>::insert;
    using std::vector<T>::emplace;
    using std::vector<T>::erase;
    using std::vector<T>::push_back;
    using std::vector<T>::emplace_back;
    using std::vector<T>::pop_back;
    using std::vector<T>::resize;
    using std::vector<T>::swap;
};

} // NS imajuscule
