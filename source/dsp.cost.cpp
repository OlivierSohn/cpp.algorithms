
namespace imajuscule {

double computeCostWriteNConsecutiveCachelines(int n) {
    
    int const nbytes = n * cache_line_n_bytes;
    
    a64::vector<uint8_t> bytes;
    
    int const Max = std::max(10000000, nbytes);

    bytes.resize(Max);
    
    
    std::vector<int> indexes;
    indexes.reserve(Max);

    // different tests sit at least 4 cachelines appart
    int offsetBetweenTests = nbytes + 4*cache_line_n_bytes;
    // make it a multiple of 512
    if(offsetBetweenTests < 512) {
        offsetBetweenTests = 512;
    }
    offsetBetweenTests = (offsetBetweenTests/512) * 512;
    for(int index=0;
        index+nbytes <= Max;
        index += offsetBetweenTests)
    {
        indexes.push_back(index);
    }
    std::shuffle(indexes.begin(), indexes.end(), lagged_fibonacci<SEEDED::No>());
    auto duration = imajuscule::profiling::measure_thread_cpu_one([&bytes, &indexes, nbytes](){
        for(auto index : indexes) {
            for(int j=0; j<nbytes; ++j) {
                bytes[index + j] = j;
            }
        }
    });
    auto sum = std::accumulate(bytes.begin(), bytes.end(), 0.);
    
    auto costMicroseconds = duration.count() / static_cast<double>(indexes.size());
    if(sum > 3.) {
        // Just to use the value, making sure it's mot optimized away
        costMicroseconds += 0.0000001;
    }
    return costMicroseconds / 1000000.;
}


double computeCostReadNConsecutiveCachelines(int n) {
    
    int const nbytes = n * cache_line_n_bytes;
    
    a64::vector<uint8_t> bytes;
    
    int const Max = std::max(10000000, nbytes);

    bytes.resize(Max);
    
    std::iota(bytes.begin(), bytes.end(), 0);
    
    std::vector<int> indexes;
    indexes.reserve(Max);

    // different tests sit at least 4 cachelines appart
    int offsetBetweenTests = nbytes + 4*cache_line_n_bytes;
    // make it a multiple of 512
    if(offsetBetweenTests < 512) {
        offsetBetweenTests = 512;
    }
    offsetBetweenTests = (offsetBetweenTests/512) * 512;
    for(int index=0;
        index+nbytes <= Max;
        index += offsetBetweenTests)
    {
        indexes.push_back(index);
    }
    std::shuffle(indexes.begin(), indexes.end(), lagged_fibonacci<SEEDED::No>());
    int sum=0;
    auto duration = imajuscule::profiling::measure_thread_cpu_one([&bytes, &indexes, nbytes, &sum](){
        for(auto index : indexes) {
            for(int j=0; j<nbytes; ++j) {
                sum += bytes[index + j];
            }
        }
    });
    
    auto costMicroseconds = duration.count() / static_cast<double>(indexes.size());
    if(sum > 3.) {
        // Just to use the value, making sure it's mot optimized away
        costMicroseconds += 0.0000001;
    }
    return costMicroseconds / 1000000.;
}

} // NS imajuscule

