#pragma once

#include <iterator>

#include "range.hpp"

namespace imj {

    /*
     * Allows to view a container as a max heap
     */
    template< typename iterator>
    struct HeapView {
        HeapView(iterator begin, iterator end) :
        size_(std::distance(begin, end)),
            begin_(begin)
        {}
        
        void make_max_heap() const {
            // iterate starting from the base of the heap to the top
            // (excluding leaves which are already max heaps)
            
            auto index = ( size_ - 1 )/2;

            while(index != -1) {
                max_heapify(index);

                --index;
            }
        }
        
        void setEnd(iterator end) {
            assert(size_ > 0);
            size_ = std::distance(begin_, end);
        }
        
        // assumes, in the subtree rooted at root, that the only location where the max heap property
        // could be violated is at the root itself. This method fixes this violation
        void max_heapify(size_t index_root) const {
            // if the root value is bigger than its two children there is nothing to do
            // because the children are assumed to be max_heaps already

            index_root++; // heap indexes start at 1
            
            auto left_index = 2*index_root;
            if(left_index > size_) {
                // no left child
                return;
            }
            auto it_left = begin_;
            std::advance(it_left, left_index-1);  // heap indexes start at 1
            
            auto & left_value = *it_left;
            
            auto it_root = begin_;
            std::advance(it_root, index_root-1); // heap indexes start at 1
            
            auto & root_value = *it_root;
            
            auto right_index = left_index + 1;
            if(right_index > size_) {
                // no right child
                
                if(left_value <= root_value) {
                    return;
                }
                
                std::swap(left_value, root_value);
                max_heapify(left_index - 1); // heap indexes start at 1
            }
            else {
                // right child
                auto it_right = it_left;
                std::advance(it_right, 1);
                
                auto & right_value = *it_right;
                
                if(left_value > right_value) {
                    if(left_value <= root_value) {
                        return;
                    }
                    std::swap(left_value, root_value);
                    max_heapify(left_index - 1); // heap indexes start at 1
                }
                else {
                    if(right_value <= root_value) {
                        return;
                    }
                    std::swap(right_value, root_value);
                    max_heapify(right_index - 1); // heap indexes start at 1
                }
            }            
        }
        
        void heapify_root() const {
            max_heapify(0);
        }
        
        bool empty() const {
            return size_ == 0;
        }
    private:
        iterator begin_;
        size_t size_;
    };

    template< typename iterator>
    struct HeapSort {
        using range = imj::range<iterator>;
        using value_type = typename std::iterator_traits<iterator>::value_type;
        
        void operator ()(iterator begin, iterator end) {
            if(begin == end) {
                return;
            }
            sort(std::move(begin), std::move(end));
        }
        
    private:
        
        void sort(iterator begin, iterator end) {
            assert(begin != end);
            auto distance = std::distance(begin, end);
            if(distance <= 1) {
                return;
            }
            
            auto last = begin;
            std::advance(last, distance - 1);
            
            HeapView<iterator> h(begin, end);
            h.make_max_heap();
            
            while(true) {
                std::swap(*begin, *last);
                h.setEnd(last);
                if(h.empty()) {
                    return;
                }
                h.heapify_root();
                std::advance(last, - 1);
            }
        }
    };

    template< typename iterator > 
    void heap_sort(iterator begin, iterator end) {
        HeapSort<iterator> s;
        s(std::move(begin), std::move(end));
    }

} // NS imj

