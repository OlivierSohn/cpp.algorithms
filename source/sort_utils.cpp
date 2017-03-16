
namespace imajuscule {
    std::minstd_rand & shuffle_rng_engine() {
        static std::minstd_rand engine;
        return engine;
    }
    
} // NS imajuscule

