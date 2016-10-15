#pragma once

#include <iterator>
#include <array>

#include "range.hpp"
#include "bounded_lifo.hpp"


namespace imj {

    template< typename iterator>
    std::string logstr(iterator begin, iterator end) {
        std::string str = "\n{ ";    
        std::for_each(begin, end, [&str](typename std::iterator_traits<iterator>::value_type const & e) {
            str += std::to_string(e) + " ";
        });
        str += (" }\n");
        return str;
    }

    template< typename T>
    std::string logstr(T const & container) {
        logstr(container.begin(), container.end());
    }
    
    template< typename iterator>
    void log(iterator begin, iterator end) {
        std::cout << logstr(begin, end);    
    }
 
    template< typename T>
    void log(T const & container) {
        log(container.begin(), container.end());
    }

    template< typename iterator, typename WorkContainer = std::vector<typename std::iterator_traits<iterator>::value_type> >
    struct mergesort {
        using range = imj::range<iterator>;
        using value_type = typename std::iterator_traits<iterator>::value_type;
        
        mergesort(WorkContainer & work) :
            work(work)
        {
        }
        
        void operator ()(range r) {
            static_assert( 
                std::is_same<
                    typename WorkContainer::value_type,
                    value_type 
                >::value,
                "workcontainer and iterator must have same value_type"
                );

            auto d = r.distance();
            if(d > 1) {
                work.resize(d);
                this->sort(r);
            }
        }
        
    private:
        WorkContainer & work;
        
        std::pair<range, range> split(range r, int distance) {
            assert(r.distance() >= 2);
            assert(distance == r.distance());
            auto pivot = distance / 2;
            assert(pivot > 0);
            auto intermediate = r.begin();
            std::advance(intermediate, pivot);
            return {
                {r.begin(), intermediate},
                {intermediate, r.end()}
            };    
        }
        
        void sort(range r) {     
            auto size = r.distance();       
            auto ranges = split(std::move(r), size);

            if(ranges.first.distance() != 1) {
                this->sort(ranges.first);
            }

            if(ranges.second.distance() != 1) {
                this->sort(ranges.second);
            }
            
            if(size == 2) {
                auto & v1 = *(ranges.first.begin());
                auto & v2 = *(ranges.second.begin());
                if( v1 > v2 ) {
                    std::swap(v1, v2);
                }
            } else if(size < 7 ) { // 7 is optimal on my Ubuntu
                std::array<value_type, 7> small_work;
                join(ranges, small_work.begin());
            }
            else {
                join(ranges, work.begin());
            }
        }

        template <typename work_iterator>
        void join(std::pair<range, range> & ranges, work_iterator work_begin) {
            auto & r1 = ranges.first;
            auto & r2 = ranges.second;
            
            assert(r1.distance() >= 1);
            assert(r2.distance() >= 1);
            
            assert(r1 < r2);
            assert(r2.follows( r1 ));    
            
            auto not_sorted_begin = r1.begin();
            // throughout this function, items from 'r1.begin()'
            // to 'not_sorted_begin' (excluded) are sorted

            auto it2 = r2.begin(); // the next unsorted item in r2
            
            // work is supposed to have a size >= the number 
            // of push_back calls that can be made on unsorted_r1
            assert(work.size() >= r1.distance() + r2.distance());
            
            bounded_lifo<value_type, work_iterator> unsorted_r1(work_begin);
            
            auto r1_end = r1.end();
            
            auto r2_end = r2.end();
            
            do {
                if(*not_sorted_begin >= *it2) {
                    unsorted_r1.push_back(*not_sorted_begin);
                    *not_sorted_begin = *it2;
                    ++not_sorted_begin;
                    ++it2;               
                    if(it2 == r2_end) {
                        std::copy(unsorted_r1.begin(), unsorted_r1.end(), not_sorted_begin);
                        return;
                    }
                }
                else {
                    ++not_sorted_begin;
                }
                
                if(not_sorted_begin == r1_end) {
                    break;
                }
                                         
                if(unsorted_r1.empty()) {
                    continue;
                }
            
                bool done = false;
                while(unsorted_r1.front() <= *it2) {
                    unsorted_r1.push_back(*not_sorted_begin);
                    *not_sorted_begin = unsorted_r1.pop_front();
                    ++not_sorted_begin;
                    if(not_sorted_begin == r1_end) {
                        done = true;
                        break;
                    }
                    assert(!unsorted_r1.empty());
                }
                if(done) {
                    break;
                }
            } while(1);
            
            if(not_sorted_begin == it2) {
                return;
            }

            while(it2 != r2_end) {
                while(unsorted_r1.front() <= *it2) {
                    *not_sorted_begin = unsorted_r1.pop_front();
                    ++not_sorted_begin;
                    if(not_sorted_begin == it2) {
                        return;
                    }
                }
            
                // not_sorted_begin within r2 : the elements of r1 
                // that are not yet sorted are in unsorted_r1. However 
                // we know that if unsorted_r1 is non empty, all elements
                // are smaller that *it2
                *not_sorted_begin = *it2;   
                ++not_sorted_begin;
                ++it2;             
            }

            std::copy(unsorted_r1.begin(), unsorted_r1.end(), not_sorted_begin);
        }
    };

    template< typename iterator, typename WorkContainer > 
    void merge_sort_work(iterator begin, iterator end, WorkContainer & work_vector) {
        mergesort<iterator, WorkContainer> sort_functor(work_vector);
        sort_functor({begin, end});
    }

    template< typename iterator > 
    void merge_sort(iterator begin, iterator end) {
        // to optimize things, we reuse the same work vector 
        // todo make thread safe
        static std::vector<typename std::iterator_traits<iterator>::value_type> work_vector; 
        
        merge_sort_work(begin, end, work_vector);
    }

} // NS imj

