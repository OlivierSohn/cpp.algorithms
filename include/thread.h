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
            //return 63; // osx
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
        
        struct MachSchedParams {
            bool read();
            bool write() const;
            
            bool setRealTime();
            bool setNonRealTime();
            
            void log() const;
            
        private:
            bool ok = false;
            thread_extended_policy_data_t ext_policy;
            thread_precedence_policy_data_t precedence_policy;
            thread_time_constraint_policy_data_t time_constraint_policy;
        };
        
        struct SchedParams {
            bool read() {
                bool ok = mach.read();
                ok = posix.read() && ok;
                return ok;
            }
            
            bool write() const {
                bool ok = mach.write();
                ok = posix.write() && ok;
                return ok;
            }
            
            void log() const {
                mach.log();
                posix.log();
            }
            
            MachSchedParams mach;
            PosixSchedParams posix;
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
                ok = cur.mach.setRealTime();
                if(!ok) {
                    cerr << "could not set real time" << endl;
                    return;
                }
//                cur.posix.set(policy, priority);
                //cur.log(); // log Before switching to real time
                /*if(!cur.write()) {
                    cerr << "could not set priority" << endl;
                    return;
                }*/
            }
            
            ~ScopedPriorityChange() {
                if(!ok) {
                    return;
                }
  /*              if(!prev.write()) {
                    std::cerr << "could not set previous priority" << std::endl;
                    return;
                }*/
                if(!prev.mach.setNonRealTime()) {
                    std::cerr << "could not set non realtime" << std::endl;
                }
                //prev.log();
            }
            
        private:
            bool ok;
            SchedParams prev, cur;
        };
    }
}
