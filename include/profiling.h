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
            return std::accumulate(a.begin(), a.end(), T{}).count() / static_cast<float>(a.size());
        }

        // could be optimized with vDSP_minv
        template<typename Array, typename T = typename Array::value_type>
        float min_(Array const & a) {
            assert(!a.empty());
//            StringPlot p(20, a.size());
//            p.draw(a, default_curve_char, true);
//            p.log();
            return std::min_element(a.begin(), a.end())->count();
        }

        void pollute_cache(std::ostream &);

    enum class TimerOption {
        YieldBeforeStart,
        None
    };
        /*
        * Measures elapsed system time between construction and destruction of the object.
        * The constructor first calls 'std::this_thread::yield()' so that the measured action
        * starts at the beginning of a time slice.
        */
        template<typename Clock, TimerOption tt = TimerOption::None>
        struct Timer
        {
            using rep = typename Clock::rep;
            using time_point = typename Clock::time_point;
            using resolution = typename Clock::duration;

            Timer(resolution& duration) :
            duration(&duration) {
                if constexpr (tt == TimerOption::YieldBeforeStart) {
                    // to make sure that we have a full time slice for the subsequent
                    // action we want to measure:
                    std::this_thread::yield();
                }
              startTime = Clock::now();
            }
            ~Timer() {
                using namespace std::chrono;
                *duration = duration_cast<resolution>(Clock::now() - startTime);
            }
        private:

            time_point startTime;
            resolution* duration;
        };
    

        template<typename Clock, typename F, typename rep = typename Clock::rep>
        rep measure_one(F f) {
            rep duration;

            {
                Timer<Clock, TimerOption::YieldBeforeStart> t(duration);

                f();
            }

            return duration;
        }

        template<typename Clock, typename PREP, typename F, typename duration = typename Clock::duration>
        std::vector<duration> measure_n(int n_warmup, int n, PREP preparation, F f) {
            assert(n > 0);
                        
            std::vector<duration> durations(n_warmup + n);

            for(auto & d : durations)
            {
                preparation();

                {
                    Timer<Clock, TimerOption::YieldBeforeStart> t(d);
                    f();
                }
            }
            return {durations.begin() + n_warmup, durations.end()};
        }
    }
}
