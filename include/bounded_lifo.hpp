#pragma once

#include <vector>

namespace imajuscule {

    /*
    * Fast "Last-in First-out" queue, with following restrictions:
    *
    * - During the entire lifecycle of this object, the number of calls to 'push_back'
    *       should not be greater than the size of the underlying container associated 
    *       to the iterator passed in the constructor.
    * - During the entire lifecycle of this object, the underlying container's iterators
    *       must remain valid.
    */
    template< typename T, typename iterator = typename std::vector<T>::iterator >
    struct bounded_lifo {

        // 'begin' is the begin iterator of a container whose iterators should not be 
        //   invalidated during the lifetime of this object.
        //   The said container must have size equal or greater to the number of future
        //   calls that will be made to push_back. 
        bounded_lifo(iterator begin) :
            begin_(begin), 
            end_(begin)
        {}
        
        void push_back(T bycopy /*todo use traits*/) {
            *end_ = bycopy;
            ++end_;
        }
        
        T pop_front() {
            assert(!empty());        
            auto & ret = *begin_;
            ++begin_;
            assert(begin_ <= end_);
            return ret;
        }
        
        T front() const {
            assert(!empty());        
            return *begin_;
        }
        
        bool empty() const {
            assert(begin_ <= end_);
            return begin_ == end_;
        }
        
        iterator begin() const {
            return begin_;
        }
        iterator end() const {
            return end_;
        }
        
    private:
        iterator begin_, end_;
    };

} // NS imajuscule

