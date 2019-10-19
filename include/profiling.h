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

        void pollute_cache(std::ostream &);

        /*
        * Measures elapsed system time between construction and destruction of the object.
        * The constructor first calls 'std::this_thread::yield()' so that the measured action
        * starts at the beginning of a time slice.
        */
        template<typename Clock>
        struct Timer
        {
            using rep = typename Clock::rep;
            using time_point = typename Clock::time_point;
            using resolution = typename Clock::duration;

            Timer(rep& duration) :
            duration(&duration) {
              // TODO it would be nice to also raise the thread priority here,
              // so that the time slices are bigger.

              // to make sure that we have a full time slice for the subsequent
              // action we want to measure:
              std::this_thread::yield();
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

            {
                Timer<Clock> t(duration);

                f();
            }

            return duration;
        }

        template<typename Clock, typename PREP, typename F, typename rep = typename Clock::rep>
        std::vector<rep> measure_n(int n_warmup, int n, PREP preparation, F f) {
            assert(n > 0);
                        
            std::vector<rep> durations(n_warmup + n);

            for(auto & duration : durations)
            {
                preparation();

                {
                    Timer<Clock> t(duration);
                    f();
                }
            }
            return {durations.begin() + n_warmup, durations.end()};
        }
    }
}
