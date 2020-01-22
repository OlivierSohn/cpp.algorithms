/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    
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

    std::minstd_rand & shuffle_rng_engine();
        
    
    namespace {
        template< typename Container >
        struct Shuffle_ {
            
            // we have a specific engine because we want to reseed it to a given value
            // in order to make tests repeatable
            void operator()(Container & ref) {
                
                std::shuffle(ref.begin(), ref.end(), shuffle_rng_engine());
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
        std::optional<typename Container::value_type> val;
        for(auto v : container) {
            if(val && (*val > v)) {
                return false;
            }
            val = v;
        }
        return true;
    }

} // NS imajuscule

