

// we should assert that once it decreases, it does so until 0

namespace imajuscule {
    
    struct ControlledPoolGrowth;
    struct ControlledPoolGrowthObject;

    struct Pool {
#ifndef NDEBUG
        friend struct ControlledPoolGrowth;
        friend struct ControlledPoolGrowthObject;
#endif
        static constexpr auto N = 4000;
        
        static Pool & getInstance();
        
        void * GetNext(size_t const alignment, size_t const n_bytes, size_t const n_elems ) {
#ifndef NDEBUG
            assert(state == Growing);
            // if the previous assert breaks it means either
            // - the container using the Pool is reallocating to grow. If we let it happen, there
            // will be memory fragmentation. To fix it, use vector::reserve for example
            // - or the pool is used by more than one container and they don't respect the '2 phase'
            // principle : pool is not used, then all allocation happen, then all deallocation happen
            // to return to the state where pool is not used
#endif
            auto index = buffer.size() - space_left;
            if(space_left) {
                void* ptr = &buffer[index];
                auto prev = space_left;
                if(std::align(alignment, n_bytes, ptr, space_left)) {
                    assert(space_left <= prev);
                    align_wasted += prev - space_left;
                    space_left -= n_bytes;
                    count_elems += n_elems;
                    max_size_used = std::max(max_size_used, buffer.size() - space_left);
                    return ptr;
                }
            }
            if(index != 0) {
                // our buffer is fully utilized, we cannot reallocate it because there are allocated
                // elements in the program. We will delay reallocation at the time when we are the primary pool
                // and the number of elements is back to 0.
                {
                    // we do that in order to be sure that future allocations will not use this pool
                    align_wasted += space_left;
                    space_left = 0;
                }
                if(!overflow) {
                    allocate_overflow();
                }
                return overflow->GetNext(alignment, n_bytes, n_elems);
            }
            return nullptr;
        }
        
        // this method is called only for the primary pool
        void Free(void * p, size_t n_bytes, size_t n_elems) {
            auto res = doFree(p, n_bytes, n_elems);
            assert(res);
            if(count_elems != 0) {
                return;
            }
            // resize the vector to fit all the data in our pool next time
            auto npools = countPools();
            auto new_size = npools * buffer.size();
            if(npools > 1) {
                std::cout << "the pool grows to " << new_size << ", " << npools << " pools were used" << std::endl;
            }
            init(new_size);
        }

        size_t get_max_size_used() const { return max_size_used; }
        void reset_max_size_used() { max_size_used = 0; }
        
        size_t size() const { return buffer.size() + (overflow ? overflow->size() : 0); }
        size_t count() const { return count_elems + (overflow ? overflow->count() : 0); }
        bool used() const { return buffer.size() != space_left; }
        
        int countPools() const {
            if(overflow) {
                return 1 + overflow->countPools();
            }
            return 1;
        }
        
        // not taking into account alignment constraints here...
        size_t maxElemSize() const {
            return buffer.size();
        }
        
        void resize(size_t sz = N) {
            init(sz);
        }
    private:
        
        Pool(size_t const sz=N) {
            resize(sz);
        }
        
        void allocate_overflow();
        
        void init(size_t const n)
        {
            overflow.release();
            if(n != buffer.size()) {
                buffer.clear();
                buffer.resize(n);
            }
            space_left = n;
            align_wasted = 0;
#ifndef NDEBUG
            state = Growing;
#endif
        }
        
        bool doFree(void * p, size_t & n_bytes, size_t & n_elems) {
#ifndef NDEBUG
            state = Shrinking;
#endif
            if(overflow && overflow->doFree(p, n_bytes, n_elems)) {
#ifndef NDEBUG
                if(overflow->count_elems == 0) {
                    if(controlled_pool_count > 0 && count_elems == controlled_pool_count) {
                        state = Growing;
                    }                    
                }
#endif
                return true;
            }
            if(count_elems == 0) {
                return false;
            }
            count_elems -= n_elems;
            space_left += n_bytes;
            if(count_elems >= 0) {
#ifndef NDEBUG
                if(controlled_pool_count > 0 && count_elems == controlled_pool_count) {
                    state = Growing;
                }
#endif
                assert(space_left <= buffer.size());
                if(count_elems == 0) {
                    assert(space_left + align_wasted == buffer.size());
                    space_left = buffer.size();
                }
                return true;
            }
            n_elems = -count_elems;
            space_left = buffer.size();
            assert(space_left > buffer.size());
            n_bytes = space_left - buffer.size();
            count_elems = 0;
            align_wasted = 0;
            assert(0); // are we supposed to have this scenario??
            // it means the order of allocations and the order of deallocations is not consistent
            return false;
        }
        
        static Pool * instance;

        // we could make that an array instead and leave it as first member so that
        // big elements find a good alignment at the beginning
        // If we do that the pool cannot grow by itself anymore when the last element
        // is freed, instead we need an extra level of indirection ('PoolManager')
        // that handles that : the allocator would have a pointer to the pool given
        // by the pool manager and when Free returns '0' (last element is freed), it
        // would need to tell the pool manager that he can resize it now
        std::vector<unsigned char> buffer;
#ifndef NDEBUG
        enum State { Growing, Shrinking };
        static State state;
        int32_t controlled_pool_count = -1; // -1 means this pool has no controlled level
#endif
        size_t space_left;
        int32_t align_wasted = 0;
        int32_t count_elems = 0;
        std::unique_ptr<Pool> overflow;
        size_t max_size_used = 0;

#ifndef NDEBUG
        struct ControlPoint {
            int32_t pool_index = -1, count_elems;
            bool operator ==(const ControlPoint & o) {
                return pool_index == o.pool_index && count_elems == o.count_elems;
            }
        };
        
        struct Control {
            ControlPoint cur, prev;
        };
        
        void saveCurrentControl(ControlPoint & ctrl_point) {
            ++ ctrl_point.pool_index;
            if(controlled_pool_count >= 0) {
                ctrl_point.count_elems = controlled_pool_count;
                controlled_pool_count = -1;
            }
            else if(overflow) {
                overflow->saveCurrentControl(ctrl_point);
            }
            else {
                ctrl_point.count_elems = controlled_pool_count;
            }
        }
        
        void setCurrentControl(ControlPoint & ctrl_point) {
            ++ ctrl_point.pool_index;
            if(overflow && overflow->count_elems) {
                overflow->setCurrentControl(ctrl_point);
            }
            else {
                ctrl_point.count_elems = count_elems;
                controlled_pool_count = count_elems;
            }
        }
        
        void applyControl(ControlPoint & ctrl_point) {
            if( 0 == ctrl_point.pool_index) {
                controlled_pool_count = ctrl_point.count_elems;
            }
            else {
                --ctrl_point.pool_index;
                assert(overflow);
                overflow->applyControl(ctrl_point);
            }
        }
        
        Control pushControl() {
            Control ctrl;
            saveCurrentControl(ctrl.prev);
            setCurrentControl(ctrl.cur);
            return ctrl;
        }
        
        void popControl(Control sz) {
            ControlPoint check;
            saveCurrentControl(check);
            // verify that current control is the one we set
            assert( check == sz.cur );
            applyControl(sz.prev);
        }
#endif
    };

    /*
     * Use when you are using the pool for a container,
     * you reserved this container fully,
     * and you want to be able to use the pool for other
     * containers that have a smaller lifecycle than your
     * container
     */
    struct ControlledPoolGrowthObject {
#ifndef NDEBUG
        ControlledPoolGrowthObject()
        :used(false)
        {}
        
        ~ControlledPoolGrowthObject() {
            if(used) {
                release();
            }
        }
        
        // movability added to be able to move StaticVector
        ControlledPoolGrowthObject(ControlledPoolGrowthObject && o):
        used(o.used),
        ctrl(o.ctrl) {
            o.used = false;
        }
        ControlledPoolGrowthObject& operator =(ControlledPoolGrowthObject &&o)  {
            if(used) {
                release();
            }
            used = o.used;
            o.used = false;
            ctrl = o.ctrl;
            return *this;
        }
#endif
        
        void acquire() {
#ifndef NDEBUG
            assert(!used);
            ctrl = Pool::getInstance().pushControl();
            used = true;
#endif
        }
        
        void release() {
#ifndef NDEBUG
            assert(used);
            Pool::getInstance().popControl(ctrl);
            used = false;
#endif
        }
        
    private:
#ifndef NDEBUG
        bool used : 1;
        Pool::Control ctrl;
#endif
    };

    
    struct ControlledPoolGrowth {
#ifndef NDEBUG
        ControlledPoolGrowth() { cpgo.acquire(); }
        ~ControlledPoolGrowth() { cpgo.release(); }
    private:
        ControlledPoolGrowthObject cpgo;
#endif
    };

} // ns imajuscule
