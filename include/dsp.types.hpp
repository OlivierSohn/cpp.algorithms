
namespace imajuscule {


struct Latency {
    constexpr explicit Latency(int latency)
    : latency(latency)
    {}
    
    constexpr int toInteger() const {
        return latency;
    }
    
    Latency operator + (Latency const & o) const {
        return Latency(latency + o.latency);
    }
    Latency operator - (Latency const & o) const {
        return Latency(latency - o.latency);
    }
    
    constexpr bool operator == (Latency const & o) const {
        return latency == o.latency;
    }
    constexpr bool operator != (Latency const & o) const {
        return latency != o.latency;
    }
    constexpr bool operator < (Latency const & o) const {
        return latency < o.latency;
    }
    constexpr bool operator > (Latency const & o) const {
        return latency > o.latency;
    }
    constexpr bool operator <= (Latency const & o) const {
        return latency <= o.latency;
    }
    constexpr bool operator >= (Latency const & o) const {
        return latency >= o.latency;
    }
private:
    int latency;
};

}
