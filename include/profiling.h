/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace profiling {
        
        template<typename Array, typename T = typename Array::value_type>
        T avg(Array const & a) {
            assert(!a.empty());
            return std::accumulate(a.begin(), a.end(), T{}) / a.size();
        }
        
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
                Timer<Clock> t(duration);
                f();
            }
            
            return duration;
        }
        
        template<typename Clock, typename PREP, typename F, typename rep = typename Clock::rep>
        std::vector<rep> measure_n(int n, PREP preparation, F f) {
            assert(n > 0);
            std::vector<rep> durations(n);
            durations.reserve(n);
            
            for(auto & duration : durations)
            {
                preparation();
                
                pollute_cache();
                
                {
                    Timer<Clock> t(&duration);
                    f();
                }
            }
            return std::move(durations);
        }
    }
}
