#pragma once

#include <iterator>
#include <array>

#include "range.hpp"
#include "bounded_lifo.hpp"
#include "insertion_sort.hpp"

namespace imj {

    constexpr bool zero_or_powOf2(int x) {
        return (x & (x-1)) == 0;
    }
    
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
        return logstr(container.begin(), container.end());
    }
    
    template< typename iterator>
    void log(iterator begin, iterator end) {
        std::cout << logstr(begin, end);    
    }
 
    template< typename T>
    void log(T const & container) {
        log(container.begin(), container.end());
    }

    enum AlgoType {
        RECURSIVE,
        SEQUENTIAL,
        SEQUENTIAL_CACHE_OPTIMIZED
    };
     

    template< typename iterator
              , typename WorkContainer
              , int insertion_sort_below_size >
    struct MergeSort {
        using range = imj::range<iterator>;
        using value_type = typename std::iterator_traits<iterator>::value_type;
        
        MergeSort(WorkContainer & work) :
            work(work)
        {
        }
        
        void operator ()(range r, AlgoType t = SEQUENTIAL) {
            
            static_assert(
                std::is_same<
                    typename WorkContainer::value_type,
                    value_type 
                >::value,
                "workcontainer and iterator must have same value_type"
                );

            auto d = r.distance();
            if(d <= 1) {
                return;
            }
            
            work.resize(d);
            
            // It doesn't seem to optimize, I need to profile to see why
            if(t==SEQUENTIAL_CACHE_OPTIMIZED) {
                this->sort_seq_cache(r);
            }
            else if(t==SEQUENTIAL) {
                this->sort_seq(r);
            }
            else {
                this->sort(r);
            }
        }
        
    private:
        WorkContainer & work;
        
        // not very usefull for performances...
        void sort_seq_cache(range r) {
            enum { max_distance = 4 };
            
            struct RangeSplit {
                RangeSplit(range const & range_) :
                    end(range_.begin()),
                    remaining_distance(range_.distance()) 
                {}
                
                bool next() {
                    if(remaining_distance == 0) {
                        return false;
                    }
                    
                    it = end;
                    advance_end();                    
                    return true;
                }
                
                auto result() const {
                    return range{it, end};
                }                
            
            private : 
                iterator it, end;             
                int remaining_distance;

                void advance_end() {
                    if(remaining_distance >= max_distance) {
                        std::advance(end, max_distance);
                        remaining_distance -= max_distance;
                    }
                    else
                    {
                        std::advance(end, remaining_distance);
                        remaining_distance = 0;
                    }
                    assert(remaining_distance >= 0);
                }
            };
            
            RangeSplit splitter(r);
            int count_splits=0;
            while(splitter.next()) {
                sort_seq(splitter.result());
                count_splits++;
            }

            if(count_splits >= 2) {
                sort_seq(r, max_distance);                
            }
        }
        
        void sort_seq(range r, int cell_size = 1) {
            auto distance = r.distance();
            
            if( cell_size * 2 <= insertion_sort_below_size) {
                sort_by_insertion(r, distance);
                
                static_assert( zero_or_powOf2(insertion_sort_below_size)
                              , "template parameter 'insertion_sort_below_size' must be a power of 2, or zero");
                cell_size = insertion_sort_below_size;
                if(cell_size >= distance) {
                    return;
                }
            }
            
            do {
                sort_seq_level(r, distance, cell_size);                
                cell_size *= 2;
            } while(cell_size < distance);
        }
        
        
        void sort_by_insertion(const range & r, int remaining_distance) {
            
            auto it1 = r.begin();
            
            while(remaining_distance >= insertion_sort_below_size) {
                auto it2 = it1;
                std::advance(it2, insertion_sort_below_size);
                
                insertion_sort(it1, std::move(it2));
                
                remaining_distance -= insertion_sort_below_size;
                std::advance(it1, insertion_sort_below_size);
            }
            
            if(remaining_distance > 1) {
                auto it2 = it1;
                std::advance(it2, remaining_distance);
                insertion_sort(it1, std::move(it2));
            }
        }

        void sort_seq_level(const range & r, int remaining_distance, int cell_size) {
            
            auto it1 = r.begin();                

            auto two_cells = 2 * cell_size;
            while(remaining_distance >= two_cells) {
                auto it2 = it1;
                std::advance(it2, cell_size);
                auto it3 = it2;                
                std::advance(it3, cell_size);

                join(std::pair<range, range>{{it1, it2}, {it2, it3}}, work.begin());
                
                remaining_distance -= two_cells;
                std::advance(it1, two_cells);
            }
            
            if(remaining_distance > 1) {
                while(cell_size > remaining_distance)
                {
                    cell_size /= 2;
                }
                auto it2 = it1;
                std::advance(it2,cell_size);
                auto it3 = r.end();
                if(it3 == it2) {
                    return;
                }
                join(std::pair<range, range>{{it1, it2}, {it2, it3}}, work.begin());
            }
        }
        
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
            if(size <= insertion_sort_below_size) {
                insertion_sort(r.begin(), r.end());
                return;
            }
            
            auto ranges = split(std::move(r), size);

            if(ranges.first.distance() != 1) {
                this->sort(ranges.first);
            }

            if(ranges.second.distance() != 1) {
                this->sort(ranges.second);
            }
            
            // redundant now that we have insertion sort?
            if(size < 7 ) { // 7 is optimal on my Ubuntu
                std::array<value_type, 7> small_work;
                join(ranges, small_work.begin());
            }
            else {
                join(ranges, work.begin());
            }
        }

        template <typename work_iterator>
        void join(std::pair<range, range> const & ranges, work_iterator work_begin) {
            auto & r1 = ranges.first;
            auto & r2 = ranges.second;

            assert(r1.distance() >= 1);
            assert(r2.distance() >= 1);
            
            assert(r1 < r2);
            assert(r2.follows( r1 ));    
            
            
            bounded_lifo<value_type, work_iterator> unsorted_r1(work_begin);
            
            auto not_sorted_begin = r1.begin();
            // throughout this function, items from 'r1.begin()'
            // to 'not_sorted_begin' (excluded) are sorted
            
            auto it2 = r2.begin(); // the next unsorted item in r2
            
            
            auto r1_end = r1.end();
            auto r2_end = r2.end();
            
            do {
                auto v2 = *it2;
                if(*not_sorted_begin >= v2) {
                    unsorted_r1.push_back(*not_sorted_begin);
                    *not_sorted_begin = v2;
                    ++not_sorted_begin;
                    ++it2;               
                    if(it2 == r2_end) {
                        auto r1it = not_sorted_begin;
                        while( r1it != r1_end ) {
                            unsorted_r1.push_back(*r1it);
                            ++r1it;
                        }
                        std::copy(unsorted_r1.begin(), unsorted_r1.end(), not_sorted_begin);
                        return;
                    }
                    v2 = *it2;
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
                while(unsorted_r1.front() <= v2) {
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

            do {
                auto v2 = *it2;
                while(unsorted_r1.front() <= v2) {
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
                *not_sorted_begin = v2;
                ++not_sorted_begin;
                ++it2;
                if(it2 == r2_end) {
                    break;
                }
            } while(it2 != r2_end);

            std::copy(unsorted_r1.begin(), unsorted_r1.end(), not_sorted_begin);
        }
    };

    template< typename iterator, typename WorkContainer, int insertion_sort_below_size >
    void merge_sort_work(iterator begin, iterator end, WorkContainer & work_vector, AlgoType type) {
        MergeSort<iterator, WorkContainer, insertion_sort_below_size> sort_functor(work_vector);
        sort_functor({begin, end}, type);
    }

    template< typename iterator, int insertion_sort_below_size >
    void merge_sort(iterator begin, iterator end, AlgoType type) {
        // to optimize things, we reuse the same work vector 
        // todo make thread safe
        using Container = std::vector<typename std::iterator_traits<iterator>::value_type>;
        static Container work_vector;
        
        merge_sort_work<iterator, Container, insertion_sort_below_size>(begin, end, work_vector, type);
    }

} // NS imj

