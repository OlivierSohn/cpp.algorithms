/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace thread {
        void logSchedParams();

        enum class Priority {
            Max
        };

        static inline int toNumber(Priority p, int policy) {
            if(p != Priority::Max) {
                throw std::logic_error("unhandled priority");
            }
            return sched_get_priority_max(policy);
        }

        struct PosixSchedParams {
            bool read();
            bool write() const;

            void set(int policy, Priority priority) {
                this->policy = policy;
                param.sched_priority = toNumber(priority, policy);
            }

            void log() const;

        private:
            int policy;
            struct sched_param param;
        };

#if __APPLE__
        struct MachSchedParams {
            bool read();
            bool write() const;

            bool setRealTime();
            bool setNonRealTime();

            void log() const;

        private:
            bool ok = false;
            //thread_extended_policy_data_t ext_policy;
            //thread_precedence_policy_data_t precedence_policy;
            thread_time_constraint_policy_data_t time_constraint_policy;
        };
#endif

        struct SchedParams {
            bool read() {
                bool ok = posix.read();
#if __APPLE__
                ok = mach.read() && ok;
#endif
                return ok;
            }

            bool write() const {
                bool ok = posix.write();
#if __APPLE__
                ok = mach.write() && ok;
#endif
                return ok;
            }

            void log() const {
                posix.log();
#if __APPLE__
                mach.log();
#endif
            }

            PosixSchedParams posix;
#if __APPLE__
            MachSchedParams mach;
#endif
        private:
        };

        struct ScopedPriorityChange {
            ScopedPriorityChange(int policy, Priority priority) {
                using namespace std;
                ok = prev.read();
                if(!ok) {
                    cerr << "could not get priority" << endl;
                    return;
                }
                //prev.log();
                cur = prev;

#if __APPLE__
                ok = cur.mach.setRealTime();
#else
                cur.posix.set(policy, priority);
                //cur.log(); // log Before switching to real time
                ok = cur.write();
#endif

                if(!ok) {
                    cerr << "could not set real time" << endl;
                    return;
                }
            }

            ~ScopedPriorityChange() {
                if(!ok) {
                    return;
                }
#if __APPLE__
                if(!prev.mach.setNonRealTime()) {
                    std::cerr << "could not set non realtime" << std::endl;
                }
#else
                if(!prev.write()) {
                    std::cerr << "could not set previous priority" << std::endl;
                    return;
                }
#endif
                //prev.log();
            }

        private:
            bool ok;
            SchedParams prev, cur;
        };
    }
}
