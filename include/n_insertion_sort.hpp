namespace imajuscule {
    
    template<int n, typename Iterator>
    void nInsertSort(Iterator b, Iterator e)
    {
        static_assert(n >= 1);
        
        int sz = std::distance(b, e);
        if(sz == 0) {
            return;
        }
        
        // insertion sort on 'n' first elements of v
        for(auto itJ=b, itFin=b+(std::min(n, sz)-1); itJ != itFin; ++itJ) {
            // v is sorted on [0, j], and we make it sorted on [0, j+1]
            for(auto itK=itJ;;--itK) {
                if(*(itK+1) >= *itK) {
                    break;
                }
                std::swap(*(itK+1), *itK);
                if(itK == b) {
                    break;
                }
            }
        }
        
        using T = std::remove_reference_t<decltype(*b)>;
        T work[n];
        
        for(int i=n-1;;)
        {
            int end = std::min(i+n, sz-1);
            int m = end-i;
            if(m<=0) {
                return;
            }
            assert(m >= 1);
            
            // v is sorted on [0,i]
            // we make it sorted on [0,end]
            
            // copy next 'm' elements of v to work
            for(int k=0; k<m; ++k) {
                work[k] = *(b+(i+k+1));
            }
            
            // insertion sort on 'm' first elements of work
            for(int j=1; j != m; ++j) {
                // work is sorted on [0, j-1], and we make it sorted on [0, j]
                for(int k=j; k != 0; --k) {
                    if(work[k] >= work[k-1]) {
                        break;
                    }
                    std::swap(work[k], work[k-1]);
                }
            }
            
            // insert 'm' first elements of work in v, from biggest to smallest
            
            int j=end;
            for(int k=m; k != 0; --k) {
                *(b+j) = work[k-1];
                
                // compare elements that are k indexes apart until we find the insertion point
                
                for(;j >= k; --j) {
                    if(*(b+j) >= *(b+(j-k))) {
                        break;
                    }
                    if(k==1) {
                        std::swap(*(b+j), *(b+(j-k)));
                    }
                    else {
                        *(b+(j-1)) = *(b+j);
                        *(b+j) = *(b+(j-k));
                    }
                }
                --j;
            }
            i=end;
        }
    }
    
    template<int n, typename VEC>
    void nInsertSort(VEC & v)
    {
        static_assert(n >= 1);
        
        int sz = v.size();
        
        // insertion sort on 'n' first elements of v
        for(int j=1, fin=std::min(n, sz); j < fin; ++j) {
            // v is sorted on [0, j-1], and we make it sorted on [0, j]
            for(int k=j; k != 0; --k) {
                if(v[k] >= v[k-1]) {
                    break;
                }
                std::swap(v[k], v[k-1]);
            }
        }
        
        using T = typename VEC::value_type;
        T work[n];
        
        for(int i=n-1;;)
        {
            int end = std::min(i+n, sz-1);
            int m = end-i;
            if(m<=0) {
                return;
            }
            assert(m >= 1);
            
            // v is sorted on [0,i]
            // we make it sorted on [0,end]
            
            // copy next 'm' elements of v to work
            for(int k=0; k<m; ++k) {
                work[k] = v[i+k+1];
            }
            
            // insertion sort on 'm' first elements of work
            for(int j=1; j != m; ++j) {
                // work is sorted on [0, j-1], and we make it sorted on [0, j]
                for(int k=j; k != 0; --k) {
                    if(work[k] >= work[k-1]) {
                        break;
                    }
                    std::swap(work[k], work[k-1]);
                }
            }
            
            // insert 'm' first elements of work in v, from biggest to smallest
            
            int j=end;
            for(int k=m; k != 0; --k) {
                v[j] = work[k-1];
                
                // compare elements that are k indexes apart until we find the insertion point
                
                for(;j >= k; --j) {
                    if(v[j] >= v[j-k]) {
                        break;
                    }
                    if(k==1) {
                        std::swap(v[j], v[j-k]);
                    }
                    else {
                        v[j-1] = v[j];
                        v[j] = v[j-k];
                    }
                }
                --j;
            }
            i=end;
        }
    }
}
