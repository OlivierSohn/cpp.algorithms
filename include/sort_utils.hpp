#include <list>
#include <random>

namespace imj {
    
    // StdSort //
    
    namespace {
        template< typename Container >
        struct StdSort_ { 
            void operator()(Container & ref) {
                std::sort(ref.begin(), ref.end());
            }
        };

        template<typename T>
        struct StdSort_<std::list<T>> { 
            void operator()(std::list<T> & ref) {
                ref.sort();
            }
        };
    }
    
    template< typename Container >
    void StdSort(Container & c) {
        StdSort_<Container> sort_;
        sort_(c);
    }


    // Shuffle //

    template< typename Container >
    void Shuffle( Container & c );

    namespace {
        template< typename Container >
        struct Shuffle_ { 
            void operator()(Container & ref) {
                static auto engine = std::default_random_engine{};
                std::shuffle(ref.begin(), ref.end(), engine);
            }
        };
        
        template<typename T>
        struct Shuffle_<std::list<T>> { 
            void operator()(std::list<T> & ref) {
                std::vector<T> v{ ref.begin(), ref.end() };
                Shuffle(v);
                ref = std::list<T>(v.begin(), v.end());
            }
        };
    }

    template< typename Container >
    void Shuffle( Container & c ) {
        Shuffle_<Container> shuffle_;
        shuffle_(c);
    }
    

    // IsSorted //

    template< typename Container >
    bool IsSorted(Container const & container) {
        typename Container::value_type val;
        for(auto it = container.begin(); it != container.end(); ++it) {
            if(it == container.begin()) {
                val = *it;
                continue;
            }
            if(val > *it) {
                return false;
            }
            val = *it;
        } 
        return true;
    }

} // NS imj

