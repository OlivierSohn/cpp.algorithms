

// we should assert that once it decreases, it does so until 0

namespace imajuscule {
    struct Pool {

        static constexpr auto N = 4000;
        
        static Pool & getInstance();
        
        void * GetNext(size_t alignment, size_t n_bytes, size_t n_elems ) {
            assert(state == Growing);
            auto index = elems.size() - space_left;
            void* ptr = &elems[index];
            if(std::align(alignment, n_bytes, ptr, space_left)) {
                space_left -= n_bytes;
                i += n_elems;
                return ptr;
            }
            if(index != 0) {
                // our buffer is fully utilized, we cannot reallocate it because there are allocated
                // elements in the program. We will delay reallocation at the time when we are the primary pool
                // and the number of elements is back to 0.
                if(!overflow) {
                    allocate_overflow();
                }
                return overflow->GetNext(alignment, n_bytes, n_elems);
            }
            return nullptr;
        }
        
        // this method is called only for the primary pool
        void Free(void * p, size_t n_elems) {
            auto res = doFree(p, n_elems);
            assert(res);
            if(i != 0) {
                return;
            }
            // resize the vector to fit all the data in our pool next time
            auto npools = countPools();
            auto new_size = npools * elems.size();
            if(npools > 1) {
                std::cout << "the pool grows to " << new_size << ", " << npools << " pools were used" << std::endl;
            }
            init(new_size);
        }

        size_t size() const { return elems.size() + (overflow ? overflow->size() : 0); }
        size_t count() const { return i + (overflow ? overflow->count() : 0); }
        bool empty() const { return elems.size() == space_left; }
        
        int countPools() const {
            if(overflow) {
                return 1 + overflow->countPools();
            }
            return 1;
        }
        
        // not taking into account alignment constraints here...
        size_t maxElemSize() const {
            return elems.size();
        }
        
        void resize(size_t sz = N) {
            init(sz);
        }
    private:
        
        Pool(size_t sz=N) {
            resize(sz);
        }
        
        void allocate_overflow();
        
        void init(size_t n)
        {
            overflow.release();
            if(n != elems.size()) {
                elems.clear();
                elems.resize(n);
            }
            space_left = n;
            state = Growing;
        }
        
        bool doFree(void * p, size_t & n_elems) {
            state = Shrinking;
            if(overflow && overflow->doFree(p, n_elems)) {
                return true;
            }
            if(i == 0) {
                return false;
            }
            i -= n_elems;
            if(i >= 0) {
                return true;
            }
            n_elems = -i;
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
        std::vector<unsigned char> elems;
        enum State { Growing, Shrinking } state : 1;
        size_t space_left;
        int i = 0;
        std::unique_ptr<Pool> overflow;
    };

} // ns imajuscule
