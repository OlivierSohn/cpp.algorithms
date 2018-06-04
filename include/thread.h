/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace thread {
        void logSchedParams();

        struct PosixSchedParams {
            bool read();
            bool write() const;

            void setMaxPriority(int policy) {
                this->policy = policy;
                param.sched_priority = sched_get_priority_max(policy);
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
#if __APPLE__
                return mach.read();
#else
                return posix.read();
#endif
            }

            bool write() const {
#if __APPLE__
                return mach.write();
#else
                return posix.write();
#endif
            }

            void log() const {
#if __APPLE__
                mach.log();
#else
                posix.log();
#endif
            }

#if __APPLE__
            MachSchedParams mach;
#else
            PosixSchedParams posix;
#endif
        };
        
        static inline bool priorityIsReadOnly() {
#if __APPLE__
            return false;
#else
            SchedParams p;
            if(!p.read()) {
                return false;
            }
            p.posix.setMaxPriority(SCHED_RR);
            return !p.write();
#endif
        }

        struct ScopedRTPriority {
            ScopedRTPriority() {
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
                cur.posix.setMaxPriority(SCHED_RR);
                //cur.log(); // log Before switching to real time
                ok = cur.write();
#endif
            }

            ~ScopedRTPriority() {
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
