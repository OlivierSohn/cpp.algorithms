#include "cplx.c"

#define N_LOCAL_BUTTERFLIES       replace_N_LOCAL_BUTTERFLIES // must be a power of 2
#define N_GLOBAL_BUTTERFLIES      replace_N_GLOBAL_BUTTERFLIES // must be a power of 2, and >= N_LOCAL_BUTTERFLIES
#define LOG2_N_GLOBAL_BUTTERFLIES replace_LOG2_N_GLOBAL_BUTTERFLIES
#define MINUS_PI_over_N_GLOBAL_BUTTERFLIES replace_MINUS_PI_over_N_GLOBAL_BUTTERFLIES


inline int expand(int idxL, int log2N1, int mm) {
  return ((idxL-mm) << 1) + mm;
}

/*
 1. performs the fft of 'inputoutput'
 2. multiplies the result with 'multFreq'
 3. performs the ifft of the result
 4. returns the result in 'inputoutput'.
 ('multLeft' is the fft of a real signal,
 hence the result is real)
 */
__kernel void kernel_func(__global float *inputoutput,
                          __global struct cplx * multFreq,
                          __local struct cplx* pingpong
                          //,__global float *justoutput
                          ) {
  int const k = get_global_id(0);
  int const base_idx = k * N_LOCAL_BUTTERFLIES;

  __local struct cplx *prev = pingpong;
  __local struct cplx *next = pingpong + 2*N_GLOBAL_BUTTERFLIES;

  for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
    int const m = get_global_size(0) * j + k;
    // coalesced global memory read, local memory write with no bank conflict.
    prev[m] = complexFromReal(inputoutput[m]);
  }

  // fft
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
      __local struct cplx * tmp = prev;
      prev = next;
      next = tmp;
    }
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);
  
  // by now, prev contains the fft of the input.
  
  // we multiply the fft of the input by multFreq,
  // and we conjugate the result, to do the first step of the ifft.

  for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
    int const m = get_global_size(0) * j + k;
    // as explained in https://www.dsprelated.com/showarticle/800.php
    // we could instead swap real and imaginary values
    prev[m] = cplxConjugate(cplxMult(prev[m], multFreq[m]));
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
      __local struct cplx * tmp = prev;
      prev = next;
      next = tmp;
    }
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);
  
  // 3rd step of the ifft : we should conjugate the result, but since
  // we expect the result to be real, we skip this step.

  for(int j=0; j<2*N_LOCAL_BUTTERFLIES; ++j) {
    int const m = get_global_size(0) * j + k;
    // coalesced global memory write
    inputoutput[m] = prev[m].real;
    // prev[m].imag is expected to be 0
  }
}
