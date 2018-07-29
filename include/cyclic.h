/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    template<typename Iterator, typename Container>
    void next(Iterator & it_, Container & buf) {
        ++it_;
        if(it_ == buf.end()) {
            it_ = buf.begin();
        }
    }
    template<typename Iterator, typename Container>
    void prev(Iterator & it_, Container & buf) {
        if(it_ == buf.begin()) {
            it_ = buf.end();
        }
        --it_;
    }
  
    template<class T>
    struct cyclic
    {
        using container = typename a64::vector<T>;
        using iterator = typename container::iterator;
        using const_iterator = typename container::const_iterator;
        using value_type = T;
        using ParameterType = T;

        operator container & () { return buf; }
        operator container const& () const { return buf; }

        auto & back() { return buf.back(); }

        auto begin() const { return buf.begin();}
        auto end() const { return buf.end();}
        auto begin() { return buf.begin();}
        auto end() { return buf.end();}
        auto rbegin() const { return buf.rbegin();}
        auto rend() const { return buf.rend();}
        auto rbegin() { return buf.rbegin();}
        auto rend() { return buf.rend();}
        const_iterator cycleEnd() const { return it;}
        auto cycleEnd() { return it;}

        auto size() const { return buf.size(); }
        auto empty() const { return buf.empty(); }

      cyclic() {
        it = buf.begin();
      }

      cyclic(size_t size) : buf(size) {
        it = buf.begin();
      }

        // not copyable
        cyclic(cyclic const &) = delete;
        cyclic& operator =(cyclic const &) = delete;

        // movable
        cyclic(cyclic&&o)
        {
            auto d =  std::distance(o.buf.begin(), o.it);
            buf = std::move(o.buf);
            it = buf.begin() + d;
        }

        cyclic(container i)
        : buf(std::move(i)) {
            it = buf.begin();
        }

        cyclic & operator =(cyclic && o) {
            if(this != &o) {
                buf = std::move(o.buf);
                it = buf.begin() + std::distance(o.it, o.buf.begin());
            }
            return *this;
        }

        void resize(int sz) {
          assert(sz >= 0);
          buf.resize(sz);
          reset();
        }

      void reset() {
        it = buf.begin();
        std::fill(it, buf.end(), T{});
      }
      
        void grow(ParameterType && val) {
            auto dist = std::distance(buf.begin(), it);
            buf.push_back(std::move(val));
            it = buf.begin() + dist;
        }
      
      
      void feedAndFill(ParameterType val) {
        std::fill(buf.begin(), buf.end(), val);
        advance();
      }

        void feed(ParameterType val) {
            *it = std::move(val);
            advance();
        }

        void setByIndex(int i) {
            it = buf.begin() + i;
            assert(it < buf.end());
        }

        int getIndex() const {
            return std::distance(buf.begin(), cycleEnd());
        }

        void advance() {
            next(it, buf);
        }

        void go_back() {
            prev(it, buf);
        }

        void erase(iterator dit) {
            if(it >= dit) {
                if(it != buf.begin()) {
                    --it;
                }
            }
            buf.erase(dit);
        }

        template<typename F>
        void for_each(F f) const {
            auto start = cycleEnd();
            std::for_each(start, end(), f);
            std::for_each(begin(), start, f);
        }

        template<typename F>
        void for_each_bkwd(F f) const {
            auto start = std::reverse_iterator<const_iterator>(cycleEnd());
            std::for_each(start, rend(), f);
            std::for_each(rbegin(), start, f);
        }

        // width == 0 : only on current
        template<typename F>
        void for_each_left_and_right(int width, F f) const {
            assert(width >= 0);
            assert(it != buf.end());
            auto fwdIt = const_iterator(it);
            auto bwdIt = const_iterator(it);

            f(*it);

            for(int i=0; i<width; ++i) {
                next(fwdIt, buf);
                if(fwdIt == bwdIt) {
                    break;
                }
                f(*fwdIt);
                prev(bwdIt, buf);
                if(fwdIt == bwdIt) {
                    break;
                }
                f(*bwdIt);
            }
        }

        auto const & get_backward(int index_backward) const {
            int end_index = getIndex();
            auto real_index = end_index - 1 - index_backward;
            auto sz = buf.size();
            while(real_index < 0) {
                real_index += sz;
            }
            while(real_index >= sz) {
                real_index -= sz;
            }
            return buf[real_index];
        }

        void toBack() {
            it = buf.end()-1;
        }

    private:
        container buf;
        iterator it;
    };

    template<typename T>
    static std::ostream& operator<<(std::ostream& os, const cyclic<T>& c)
    {
        for(auto const & v : c) {
            os << v << " ";
        }
        return os;
    }

    template<typename C, typename T = typename C::value_type>
    a64::vector<T> to_vector(C const & cyclic_) {
        a64::vector<T> vec;
        vec.reserve(cyclic_.size());
        auto cycleStart = cyclic_.cycleEnd();
        auto end = cyclic_.end();
        auto begin = cyclic_.begin();

        std::move(cycleStart, end, std::back_inserter(vec));
        std::move(begin, cycleStart, std::back_inserter(vec));
        return std::move(vec);
    }
}
