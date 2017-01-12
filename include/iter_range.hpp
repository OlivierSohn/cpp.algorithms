
namespace imajuscule {

template <typename iterator>
struct iter_range {
    using value_type = typename std::iterator_traits<iterator>::value_type;
    
    iter_range(iterator begin, iterator end) : begin_(begin), end_(end) {}

    iterator begin() const {
        return begin_;
    }    

    void setBegin(iterator begin) {
        begin_ = begin;
    }

    iterator end() const {
        return end_;
    }
    
    void setEnd(iterator end) {
        end_ = end;
    }
    
    int distance() const {
        return std::distance(begin_, end_);
    }
    
    template <typename F>
    void foreach(F && f) const {
        std::for_each(begin_, end_, f);
    }
        
    bool operator < (iter_range const & other) const {
        return std::distance(end_, other.begin_) >= 0;
    }
    
    bool follows(iter_range const & other) const {
        return begin_ == other.end_;
    }
private:
    iterator begin_, end_;
};

} // NS imajuscule
