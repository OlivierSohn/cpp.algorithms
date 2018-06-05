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

      struct CtrlRTPriority {
        CtrlRTPriority() {
          using namespace std;
          ok = prev.read();
          if(!ok) {
            cerr << "could not get priority" << endl;
            return;
          }
          cur = prev;
#ifndef __APPLE__
          cur.posix.setMaxPriority(SCHED_RR);
#endif
        }
        
        void lock() {
          if(!ok) {
            return;
          }
#if __APPLE__
          ok = cur.mach.setRealTime();
#else
          ok = cur.write();
#endif
        }
        
        void unlock() {
          if(!ok) {
            return;
          }
#if __APPLE__
          ok = prev.mach.setNonRealTime();
          if(!ok) {
            std::cerr << "could not set non realtime" << std::endl;
          }
#else
          ok = prev.write();
          if(!ok) {
            std::cerr << "could not set previous priority" << std::endl;
          }
#endif
        }
        
      private:
        bool ok;
        SchedParams prev, cur;
      };


      struct ScopedRTPriority {
        ScopedRTPriority() {
          ctrl.lock();
        }
        
        ~ScopedRTPriority() {
          ctrl.unlock();
        }
        
      private:
        CtrlRTPriority ctrl;
      };
      
    }
}
