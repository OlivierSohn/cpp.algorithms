
namespace imajuscule {

double computeCostWriteNConsecutiveCachelines(int n);
double computeCostReadNConsecutiveCachelines(int n);

template<typename FPT>
double costWriteNConsecutive(int n) {
    if(n <= 0) {
        return 0.;
    }
    
    int const nTouchedCacheLines = 1 + (((n*sizeof(FPT))-1) / cache_line_n_bytes);

    int const key = floor_power_of_two(nTouchedCacheLines);
    double value;
    {
        static std::mutex mut;
        static std::map<int, double> stats;
        std::lock_guard<std::mutex> l(mut);
        
        auto it = stats.find(key);
        if(it != stats.end()) {
            value = it->second;
        }
        else {
            value = computeCostWriteNConsecutiveCachelines(key);
            stats[key] = value;
        }
    }
    return value * (nTouchedCacheLines / static_cast<double>(key));
}

template<typename FPT>
double costReadNConsecutive(int n) {
    if(n <= 0) {
        return 0.;
    }
    
    int const nTouchedCacheLines = 1 + (((n*sizeof(FPT))-1) / cache_line_n_bytes);

    int const key = floor_power_of_two(nTouchedCacheLines);
    double value;
    {
        static std::mutex mut;
        static std::map<int, double> stats;
        std::lock_guard<std::mutex> l(mut);
        
        auto it = stats.find(key);
        if(it != stats.end()) {
            value = it->second;
        }
        else {
            value = computeCostReadNConsecutiveCachelines(key);
            stats[key] = value;
        }
    }
    return value * (nTouchedCacheLines / static_cast<double>(key));
}

}
