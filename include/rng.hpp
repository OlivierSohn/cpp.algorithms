
namespace imajuscule {
    struct rng {
        static std::mt19937 & mersenne() {
            static std::mt19937 mersenne_engine_(std::random_device{}());
            return mersenne_engine_;
        }
    };
    
} // NS imajuscule

