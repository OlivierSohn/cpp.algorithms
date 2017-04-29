/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace profiling {
        
        // could be optimized using vDSP_meanv
        template<typename Array, typename T = typename Array::value_type>
        float avg(Array const & a) {
            assert(!a.empty());
//            StringPlot p(20, a.size());
//            p.draw(a, default_curve_char, true);
//            p.log();
            return std::accumulate(a.begin(), a.end(), T{}) / static_cast<float>(a.size());
        }
        
        // could be optimized with vDSP_minv
        template<typename Array, typename T = typename Array::value_type>
        float min_(Array const & a) {
            assert(!a.empty());
//            StringPlot p(20, a.size());
//            p.draw(a, default_curve_char, true);
//            p.log();
            return *std::min_element(a.begin(), a.end());
        }
        
        struct MakeRealTime {
            // deactivated until we prove that this is better...
            //thread::ScopedPriorityChange s{SCHED_OTHER, thread::Priority::Max};
        };
        
        void pollute_cache();
        
        template<typename Clock>
        struct Timer
        {
            using rep = typename Clock::rep;
            using time_point = typename Clock::time_point;
            using resolution = typename Clock::duration;
            
            Timer(rep* duration) :
            duration(duration) {
                startTime = Clock::now();
            }
            ~Timer() {
                using namespace std::chrono;
                *duration = duration_cast<resolution>(Clock::now() - startTime).count();
            }
        private:

            time_point startTime;
            rep* duration;
        };

        template<typename Clock, typename F, typename rep = typename Clock::rep>
        rep measure_one(F f) {
            rep duration;
            
            pollute_cache();

            {
                MakeRealTime rt;
                Timer<Clock> t(&duration);
                
                f();
            }
            
            return duration;
        }
        
        template<typename Clock, typename PREP, typename F, typename rep = typename Clock::rep>
        std::vector<rep> measure_n(int n, PREP preparation, F f) {
            assert(n > 0);
            std::vector<rep> durations(n);
            
            for(auto & duration : durations)
            {
                preparation();
                
                pollute_cache();
                
                {
                    MakeRealTime rt;
                    Timer<Clock> t(&duration);
                    f();
                }
            }
            return std::move(durations);
        }
    }
}
