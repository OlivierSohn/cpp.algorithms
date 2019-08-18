#include "cplx_FPT.c"

#define N_LOCAL_BUTTERFLIES       replace_N_LOCAL_BUTTERFLIES // must be a power of 2
#define N_GLOBAL_BUTTERFLIES      replace_N_GLOBAL_BUTTERFLIES // must be a power of 2, and >= N_LOCAL_BUTTERFLIES
#define LOG2_N_GLOBAL_BUTTERFLIES replace_LOG2_N_GLOBAL_BUTTERFLIES
#define MINUS_PI_over_N_GLOBAL_BUTTERFLIES replace_MINUS_PI_over_N_GLOBAL_BUTTERFLIES
#define FFT_SIZE                  replace_FFT_SIZE

inline int expand(int idxL, int log2N1, int mm) {
  return ((idxL-mm) << 1) + mm;
}

inline int h_index(int i, // index of input
                int j, // input offset
                int n  // number of partitions
                ) {
  int res = i-j;
  return res < 0? (res+n) : res;
}

/*
 1. performs the fft of 'inputoutput'
 2. copy the result to the appropriate location in 'fft_of_xs'
 3. accumulate the multiplications of 'fft_of_xs' with 'fft_of_hs'
 4. performs the ifft of the result
 5. returns the result in 'inputoutput'.
 ('fft_of_xs' are the fft of a real signal,
 hence the result is real)
 */
__kernel void kernel_func(__global FPT * inputoutput, // input and return values, size fft_size
                          __global const struct cplx_FPT * fft_of_hs, // size n*fft_size
                          __global struct cplx_FPT * fft_of_xs, // size n*fft_size
                          __local struct cplx_FPT* pingpong,
                          int i_offset, // input offset
                          int const n // number of partitions
                          ) {
  
  int const k = get_global_id(0);
  int const base_idx = k * N_LOCAL_BUTTERFLIES;

  __local struct cplx_FPT *prev = pingpong;
  __local struct cplx_FPT *next = pingpong + 2*N_GLOBAL_BUTTERFLIES;

  for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
    int const m = get_global_size(0) * j + k;
    // coalesced global memory read, local memory write with no bank conflict.
    prev[m] = complexFromReal(inputoutput[m]);
  }

  //  1. performs the fft of 'inputoutput'
  for(int i=1, LOG2_N_GLOBAL_BUTTERFLIES_over_i = LOG2_N_GLOBAL_BUTTERFLIES, log2i = 0;
      i <= N_GLOBAL_BUTTERFLIES;
      i <<= 1, --LOG2_N_GLOBAL_BUTTERFLIES_over_i, ++log2i)
  {
    barrier(CLK_LOCAL_MEM_FENCE);
    
    for(int j=0; j<N_LOCAL_BUTTERFLIES; ++j)
    {
      int const m = base_idx + j;
      int const mm = m & (i-1);
      
      int idxD = expand(m, log2i, mm); // 'm << 1' for first butterfly
      int const tIdx = mm << LOG2_N_GLOBAL_BUTTERFLIES_over_i; // 0 for first butterfly
      
      // butterfly,
      //  from : m, m+N_GLOBAL_BUTTERFLIES
      //  to   : idxD, idxD+i
      butterfly_outofplace(m,idxD,prev,next, N_GLOBAL_BUTTERFLIES, i,
                           polar(tIdx * MINUS_PI_over_N_GLOBAL_BUTTERFLIES));
    }

    // swap(prev,next)
    {
      __local struct cplx_FPT * tmp = prev;
      prev = next;
      next = tmp;
    }
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);
  
  // by now, prev contains the fft of the input.
  
  // 2. copy the result to the appropriate location in 'fft_of_xs'
  if(n>1) // else we don't need it
  {
    __global struct cplx_FPT * target_fft_of_x = fft_of_xs + i_offset*FFT_SIZE;
   
    for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
      int const m = get_global_size(0) * j + k;
      // coalesced global memory write
      target_fft_of_x[m] = prev[m];
    }
  }
  
  // we don't use a global memory barrier here, because we won't read from the location
  // we just wrote to: instead, we use 'prev' which contains the same information.
  
  // 3. accumulate the multiplications of 'fft_of_xs' with 'fft_of_hs'
  
  // We make a special case for the first iteration, for input at index i_offset, and h0
  //   because prev already contains fft_of_xj, so we don't need to read it from
  //   global memory.
  for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
    int const m = get_global_size(0) * j + k;
    // coalesced global memory read, local memory write with no bank conflict.
    prev[m] = cplxMult(prev[m], fft_of_hs[m]);
  }
  
  // For other iterations, we read fft_of_xl from global memory:
  {
    int h_index = n-1;
    for(int l=i_offset+1; l != i_offset;) {
      if(l==n) {
        l = 0;
        continue;
      }
      __global const struct cplx_FPT * fft_of_h = fft_of_hs + h_index * FFT_SIZE;
      __global struct cplx_FPT * fft_of_x = fft_of_xs + l * FFT_SIZE;
      for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
        int const m = get_global_size(0) * j + k;
        // coalesced global memory reads, local memory write with no bank conflict.
        prev[m] = cplxAdd(prev[m], cplxMult(fft_of_x[m], fft_of_h[m]));
      }
      ++l;
      --h_index;
    }
  }
  
  // 4. performs the ifft of the result

  // we conjugate the result, to do the first step of the ifft.

  for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
    int const m = get_global_size(0) * j + k;
    // as explained in https://www.dsprelated.com/showarticle/800.php
    // we could instead swap real and imaginary values
    prev[m].imag = -prev[m].imag;
  }
  
  // 2nd step of the ifft (it is simply the same code as the previous fft)
  // TODO factorize this code
  for(int i=1, LOG2_N_GLOBAL_BUTTERFLIES_over_i = LOG2_N_GLOBAL_BUTTERFLIES, log2i = 0;
      i <= N_GLOBAL_BUTTERFLIES;
      i <<= 1, --LOG2_N_GLOBAL_BUTTERFLIES_over_i, ++log2i)
  {
    barrier(CLK_LOCAL_MEM_FENCE);
    
    for(int j=0; j<N_LOCAL_BUTTERFLIES; ++j)
    {
      int const m = base_idx + j;
      int const mm = m & (i-1);
      
      int idxD = expand(m, log2i, mm); // 'm << 1' for first butterfly
      int const tIdx = mm << LOG2_N_GLOBAL_BUTTERFLIES_over_i; // 0 for first butterfly
      
      // butterfly,
      //  from : m, m+N_GLOBAL_BUTTERFLIES
      //  to   : idxD, idxD+i
      butterfly_outofplace(m,idxD,prev,next, N_GLOBAL_BUTTERFLIES, i,
                           polar(tIdx * MINUS_PI_over_N_GLOBAL_BUTTERFLIES));
    }
    
    // swap(prev,next)
    {
      __local struct cplx_FPT * tmp = prev;
      prev = next;
      next = tmp;
    }
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);
  
  // 3rd step of the ifft : we should conjugate the result, but since
  // we expect the result to be real, we skip this step.
  
  // 5. returns the result in 'inputoutput'.

  for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
    int const m = get_global_size(0) * j + k;
    // coalesced global memory write
    inputoutput[m] = prev[m].real;
    // prev[m].imag is expected to be 0
  }
}
