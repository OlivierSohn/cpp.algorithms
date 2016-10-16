#include <random>

namespace imj {
    
    struct RNG {

        static RNG & instance() {
            static RNG instance_;
            return instance_;
        }
        
        void seed(int s) {
            engine_.seed(s);
        }
        
        auto & engine() {
            return engine_;
        }
    private:
        std::default_random_engine engine_;

        RNG() = default;
    };

} // NS imj

