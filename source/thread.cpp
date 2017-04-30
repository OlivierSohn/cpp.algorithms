/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    namespace thread {
        
        void handle_error(int en, const char * msg) {
            errno = en;
            perror(msg);
        }
        
#if __APPLE__
        bool MachSchedParams::read() {
             ok = true;
            /*
             thread_port_t mach_thread_id = pthread_mach_thread_np(pthread_self());
             {
             mach_msg_type_number_t count = THREAD_EXTENDED_POLICY_COUNT;
             kern_return_t result = thread_policy_get(
             mach_thread_id, THREAD_EXTENDED_POLICY,
             reinterpret_cast<thread_policy_t>(&ext_policy), &count, &ext_default);
             if(result != KERN_SUCCESS) {
             ok = false;
             }
             }
             {
             mach_msg_type_number_t count = THREAD_PRECEDENCE_POLICY_COUNT;
             kern_return_t result = thread_policy_get(
             mach_thread_id, THREAD_PRECEDENCE_POLICY,
             reinterpret_cast<thread_policy_t>(&precedence_policy), &count, &precedence_default);
             if(result != KERN_SUCCESS) {
             ok = false;
             }
             }
             {
             mach_msg_type_number_t count = THREAD_TIME_CONSTRAINT_POLICY_COUNT;
             kern_return_t result = thread_policy_get(
             mach_thread_id, THREAD_TIME_CONSTRAINT_POLICY,
             reinterpret_cast<thread_policy_t>(&time_constraint_policy), &count, &time_constraint_default);
             if(result != KERN_SUCCESS) {
             ok = false;
             }
             }*/
            return ok;
        }
        
        bool MachSchedParams::setRealTime() {
            using namespace std;

#if TARGET_OS_IPHONE
            return false;
#else
            static int bus_speed = 0;
            if(!bus_speed) {
                int ret;
                int mib[2] = { CTL_HW, HW_BUS_FREQ };
                size_t len;
                
                len = sizeof(bus_speed);
                ret = sysctl (mib, 2, &bus_speed, &len, NULL, 0);
                if (ret < 0) {
                    bus_speed = 0;
                    throw std::logic_error("sysctl failed");
                    return false;
                }
            }
            
            auto HZ = bus_speed;
            time_constraint_policy.period=0;
            time_constraint_policy.computation=HZ/500;
            time_constraint_policy.constraint=HZ/300;
            time_constraint_policy.preemptible=0;
            int ret;
            if ((ret=thread_policy_set(pthread_mach_thread_np(pthread_self()),
                                       THREAD_TIME_CONSTRAINT_POLICY, (thread_policy_t)&time_constraint_policy,
                                       THREAD_TIME_CONSTRAINT_POLICY_COUNT)) != KERN_SUCCESS) {
                cerr << "write(time_constraint_policy) failed" << endl;
                return false;
            }
            return true;
#endif
        }
        
        bool MachSchedParams::setNonRealTime() {
#if TARGET_OS_IPHONE
            return false;
#else
            kern_return_t ret;
            thread_standard_policy_data_t pt;
            mach_msg_type_number_t cnt = THREAD_STANDARD_POLICY_COUNT;
            boolean_t get_default = TRUE;
            
            ret = thread_policy_get(pthread_mach_thread_np(pthread_self()),
                                    THREAD_STANDARD_POLICY,
                                    (thread_policy_t)&pt,
                                    &cnt, &get_default);
            if (KERN_SUCCESS != ret)
                return false;
            
            ret = thread_policy_set(pthread_mach_thread_np(pthread_self()),
                                    THREAD_STANDARD_POLICY,
                                    (thread_policy_t)&pt,
                                    THREAD_STANDARD_POLICY_COUNT);
            if (KERN_SUCCESS != ret)
                return false;
            
            return true;
#endif
        }
        
        bool MachSchedParams::write() const {
            /*
            using namespace std;
            thread_port_t threadport = pthread_mach_thread_np(pthread_self());
            
            int ret;
            if(ext_default) {
                if ((ret=thread_policy_set(threadport,
                                           THREAD_EXTENDED_POLICY, (thread_policy_t)&ext_policy,
                                           THREAD_EXTENDED_POLICY_COUNT)) != KERN_SUCCESS) {
                    cerr << "write(ext_policy) failed" << endl;
                    return false;
                }
            }
            if(precedence_default) {
                if ((ret=thread_policy_set(threadport,
                                           THREAD_PRECEDENCE_POLICY, (thread_policy_t)&precedence_policy,
                                           THREAD_PRECEDENCE_POLICY_COUNT)) != KERN_SUCCESS) {
                    cerr << "write(precedence_policy) failed" << endl;
                    return false;
                }
            }
            if(time_constraint_default) {
                if ((ret=thread_policy_set(threadport,
                                           THREAD_TIME_CONSTRAINT_POLICY, (thread_policy_t)&time_constraint_policy,
                                           THREAD_TIME_CONSTRAINT_POLICY_COUNT)) != KERN_SUCCESS) {
                    cerr << "write(time_constraint_policy) failed" << endl;
                    return false;
                }
            }
             */
            return true;
        }
        
        void MachSchedParams::log() const {
            /*
            using namespace std;
            if(ext_default) {
                cout << "extended timeshare: " << (int)ext_policy.timeshare;
            }
            if(precedence_default) {
                cout << "precendence priority: " << precedence_policy.importance;
            }
            if(time_constraint_default) {
                cout << "time constraint";
            }
            cout << endl;
             */
        }
#endif
        
        bool PosixSchedParams::read() {
            
            auto s = pthread_getschedparam(pthread_self(), &policy, &param);
            if (s) {
                handle_error(s, "pthread_getschedparam");
                return false;
            }
            return true;
        }
        
        bool PosixSchedParams::write() const {
            auto s = pthread_setschedparam(pthread_self(), policy, &param);
            if (s) {
                handle_error(s, "pthread_setschedparam");
                return false;
            }
            return true;
        }
        
        void PosixSchedParams::log() const {
            using namespace std;
            cout << "policy ";
            switch(policy) {
                case SCHED_FIFO:
                    cout << "FIFO";
                    break;
                case SCHED_OTHER:
                    cout << "OTHER";
                    break;
                case SCHED_RR:
                    cout << "RR";
                    break;
                default:
                    cout << policy;
                    break;
            }
            cout << endl;
            cout
            << "priority " << param.sched_priority
            << " (" << sched_get_priority_min(policy) << " / " << sched_get_priority_max(policy) << ")" << endl;
        }
        
        void logSchedParams() {
            SchedParams p;
            bool ok = p.read();
            if(!ok) {
                std::cerr << "could not get priority" << std::endl;
                return;
            }
            p.log();
        }
    }
}
