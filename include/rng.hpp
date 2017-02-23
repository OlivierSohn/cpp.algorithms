
namespace imajuscule {
    struct rng {
        static std::mt19937 & mersenne() {
            thread_local std::mt19937 mersenne_engine_(std::random_device{}());
            return mersenne_engine_;
        }
    };
    
} // NS imajuscule

