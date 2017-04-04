/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {

    template< typename iterator>
    struct InsertionSort {
        using iter_range = imajuscule::iter_range<iterator>;
        using value_type = typename std::iterator_traits<iterator>::value_type;
        
        void operator ()(iterator begin, iterator end) {
            if(begin == end) {
                return;
            }
            sort(std::move(begin), std::move(end));
        }
        
    private:
        
        void sort_key(iterator begin, iterator key) {
            assert(key != begin);
            
            auto before = key;
            std::advance(before, -1);
            
            while( *key < *before ) {
                std::swap(*key, *before);
                std::advance(key, -1);
                if(key == begin) {
                    return;
                }
                std::advance(before, -1);
            }
        }
        void sort(iterator begin, iterator end) {
            assert(begin != end);
            
            auto key = begin;
            std::advance(key, 1);

            while(key != end) {
                sort_key(begin, key);
                
                std::advance(key, 1);
            }
        }
    };

    template< typename iterator > 
    void insertion_sort(iterator begin, iterator end) {
        InsertionSort<iterator> s;
        s(std::move(begin), std::move(end));
    }

} // NS imajuscule

