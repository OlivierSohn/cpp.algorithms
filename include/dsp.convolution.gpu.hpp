// Using 'FIRFilterGPUAsync', it seems reasonable to think that
// the size of an fft should be bigger or equal to the size of a callback buffer,
// so that when the CPU waits for a GPU kernel to finish, it's a kernel that was launched
// in a previous callback call. This way, the CPU doesn't have to wait too much (TODO measure this).

// implement the dual of "Finegrained partition" on the GPU:
//    fft of recent input
//    copy to delayed inputs ffts
//    accumulativeMult of "fft_delayed_inputs" * "fft_h"
//    ifft

// Implement ffts of larger size on the CPU (where the inputs won't fit
// in the shared memory of a work group)

// TODO (FIRFilterGPUAsyncN)
// use one io_mem_obj and one command queue per asynchronicity level
// so that the GPU can parallelize tasks.
// Measure if this is faster (maybe the overhead or task parallelization is
// detrimental to the latency. Also there could be a bug in the device implementation
// where a queue that has more recent work to do is prioritary vs. another one
// that has old work in it).
//
// Also, an intermediate solution would be to use a pool of pair<io_mem_obj, command_queue>

// TODO (FIRFilterGPUAsyncN, FIRFilterGPUAsync)
// use clGetEventInfo(CL_EVENT_COMMAND_ EXECUTION_STATUS) to measure the
// number of times that the CPU waits for the GPU.
// We could optimize the number of levels of asynchronicity based on that metric
// to encure that the CPU never waits.

namespace imajuscule
{
  inline void kill() {
    throw std::runtime_error("program error");
  }
  inline void CHECK_CL_ERROR(int res) {
    if(res==CL_SUCCESS) {
      return;
    }
    fprintf(stderr, "OpenCL Error %d\n", res);
    kill();
  }
  
  inline int shared_memory_size(cl_device_id device_id) {
    cl_ulong local_mem_sz;
    cl_int ret = clGetDeviceInfo(device_id,
                                 CL_DEVICE_LOCAL_MEM_SIZE,
                                 sizeof(local_mem_sz), &local_mem_sz, NULL);
    CHECK_CL_ERROR(ret);
    return local_mem_sz;
  }
  
  struct OpenCLContext {
    
    OpenCLContext() {
      cl_uint ret_num_devices;
      cl_uint ret_num_platforms;
      cl_platform_id platform_id = NULL;
      cl_int ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
      CHECK_CL_ERROR(ret);
      ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_DEFAULT, 1,
                           &device_id, &ret_num_devices);
      CHECK_CL_ERROR(ret);
      
      context = clCreateContext( NULL, 1, &device_id, NULL, NULL, &ret);
      CHECK_CL_ERROR(ret);
      
      command_queue = clCreateCommandQueue(context, device_id,
                                           CL_QUEUE_PROFILING_ENABLE, &ret);
      CHECK_CL_ERROR(ret);
    }
    
    ~OpenCLContext() {
      cl_int ret = clFlush(command_queue);
      CHECK_CL_ERROR(ret);
      ret = clFinish(command_queue);
      CHECK_CL_ERROR(ret);
      ret = clReleaseCommandQueue(command_queue);
      CHECK_CL_ERROR(ret);
      ret = clReleaseContext(context);
      CHECK_CL_ERROR(ret);
    }
    
    cl_device_id device_id = NULL;
    cl_context context;
    cl_command_queue command_queue;
  };
  
  static auto const & getOpenCLContext() {
    static OpenCLContext c;
    return c;
  }

  /*
   Returns the maximum number of convolution coefficients
   that can be used in a single partition when using the GPU,
   and the pingponging stockham implementation.
   */
  template<typename T>
  int maxCountGPUStockhamFFTCoeffs() {
    int const max_local_memory = shared_memory_size(getOpenCLContext().device_id);
    // We use floor_power_of_two because our fft algorithms work on sizes that are powers of 2.
    // The factor 2 is because we are pingponging between 2 bufffers in shared memory.
    int const max_fft_length = floor_power_of_two(max_local_memory / (2*2*sizeof(T)));
    int const max_n_coeffs = max_fft_length/2;
    
    return max_n_coeffs;
  }
  

#define imj_xstr(s) imj_str(s)
#define imj_str(s) #s
  
  constexpr const char * src_root() {
    return imj_xstr(SRC_ROOT); // 'SRC_ROOT' is defined in CMakeFiles.txt
  }

  inline std::string ReplaceString(std::string subject, const std::string& search,
                            const std::string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != std::string::npos) {
      subject.replace(pos, search.length(), replace);
      pos += replace.length();
    }
    return subject;
  }
  
  inline bool device_supports_double(cl_device_id device_id) {
    cl_uint preferred_vec_width_double;
    cl_int ret;
    ret = clGetDeviceInfo(device_id, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,
                          sizeof(preferred_vec_width_double),
                          &preferred_vec_width_double,
                          NULL);
    CHECK_CL_ERROR(ret);
    return preferred_vec_width_double > 0;
  }

  template<typename T>
  struct ScopedKernel {
    
    cl_program program = 0;
    cl_kernel kernel = 0;
    int nButterfliesPerThread;
    
    void setup(cl_context context, cl_device_id device_id, std::string const & kernel_src, size_t const input_size) {
      release();
      
      using namespace imajuscule;
      if(input_size < 2) {
        throw std::invalid_argument("input_size");
      }
      
      int const nButterflies = input_size/2;
      
      if(std::is_same<double,T>::value && !device_supports_double(device_id)) {
        throw std::runtime_error("The GPU device has no native support for double (cl_khr_fp64 extension is missing). Please use 'float' instead of 'double'.");
      }
         
      // TODO if local memory can hold the output, use the kernel with local memory,
      // else use the kernel with global memory
      // We could think of a mixed approach where we compute the fft by parts:
      //   do the first levels by blocks, using local memory + writeback,
      //   omit the last writeback, use the local memory + global memory for
      //     the other levels
      //   do the writeback of the omitted portion
      
      for(nButterfliesPerThread = 1;;) {
        char buf[256];
        memset(buf, 0, sizeof(buf));
        snprintf(buf, sizeof(buf), "%a", (T)(-M_PI/nButterflies));

        std::vector<std::pair<std::string,std::string>> replacements{
          {"replace_MINUS_PI_over_N_GLOBAL_BUTTERFLIES", buf},
          {"replace_FFT_SIZE",             std::to_string(2*nButterflies)},
          {"replace_N_GLOBAL_BUTTERFLIES", std::to_string(nButterflies)},
          {"replace_LOG2_N_GLOBAL_BUTTERFLIES", std::to_string(power_of_two_exponent(nButterflies))},
          {"replace_N_LOCAL_BUTTERFLIES", std::to_string(nButterfliesPerThread)},
          {"FPT", std::is_same<float,T>::value ? "float" : "double"},
        };
        
        std::string replaced_str = kernel_src;
        for(auto [from,to] : replacements) {
          replaced_str = ReplaceString(replaced_str, from, to);
        }
        if(!std::is_same<float,T>::value) {
          replaced_str = std::string("#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n") + replaced_str;
        }


        size_t const replaced_source_size = replaced_str.size();
        const char * rep_src = replaced_str.data();
        
        // Create a program from the kernel source
        cl_int ret;
        program = clCreateProgramWithSource(context, 1,
                                            (const char **)&rep_src, (const size_t *)&replaced_source_size, &ret);
        CHECK_CL_ERROR(ret);
        
        auto options = std::string("-I ") + src_root() + " -cl-denorms-are-zero -cl-strict-aliasing -cl-fast-relaxed-math";
        
        // Build the program
        ret = clBuildProgram(program, 1, &device_id,
                             // -cl-fast-relaxed-math makes the twiddle fators computation a little faster
                             // but a little less accurate too.
                             options.c_str(),
                             NULL, NULL);
        CHECK_CL_ERROR(ret);
        
        // Create the OpenCL kernel
        kernel = clCreateKernel(program, "kernel_func", &ret);
        CHECK_CL_ERROR(ret);
        
        size_t workgroup_max_sz;
        ret = clGetKernelWorkGroupInfo(kernel,
                                       device_id,
                                       CL_KERNEL_WORK_GROUP_SIZE,
                                       sizeof(workgroup_max_sz), &workgroup_max_sz, NULL);
        CHECK_CL_ERROR(ret);
        std::cout << "workgroup max size: " << workgroup_max_sz << " for " << nButterfliesPerThread << " butterfly per thread." << std::endl;
        
        if(nButterflies > nButterfliesPerThread * workgroup_max_sz) {
          release();
          // To estimate the next value of 'nButterfliesPerThread',
          // we make the reasonnable assumption that "work group max size"
          // won't be bigger if we increase 'nButterfliesPerThread':
          nButterfliesPerThread = nButterflies / workgroup_max_sz;
          continue;
        }
        break;
      }
    }
    ScopedKernel() = default;
    
    ~ScopedKernel() {
      release();
    }
    
  private:
    void release() {
      if(kernel) {
        cl_int ret = clReleaseKernel(kernel);
        CHECK_CL_ERROR(ret);
        kernel = 0;
      }
      if(program) {
        cl_int ret = clReleaseProgram(program);
        CHECK_CL_ERROR(ret);
        program = 0;
      }
    }
    
    ScopedKernel(const ScopedKernel&) = delete;
    ScopedKernel& operator=(const ScopedKernel&) = delete;
    ScopedKernel(ScopedKernel&&) = delete;
    ScopedKernel& operator=(ScopedKernel&&) = delete;
  };
  
  
  inline std::string fullpath(std::string const & file) {
    return std::string(src_root()) + "/" + file;
  }
  
  inline void get_file_contents2(const std::string &filename, std::string & str )
  {
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in)
    {
      throw std::runtime_error("file not found");
    }
    in.seekg(0, std::ios::end);
    int res = in.tellg();
    if(-1 == res) {
      throw std::runtime_error("tellg error");
    }
    Assert(res >= 0);
    str.resize(static_cast<size_t>(res));
    in.seekg(0, std::ios::beg);
    in.read(&str[0], str.size());
    in.close();
  }
  
  inline std::string read_kernel(const std::string & kernel) {
    std::string ret;
    get_file_contents2(fullpath(kernel), ret);
    return ret;
  }
  
  namespace detail {
    template<typename T>
    struct GPUWork {
      void waitForCompletion() const {
        if(0 == finished) {
          return;
        }
        cl_int ret = clWaitForEvents(1, &finished);
        CHECK_CL_ERROR(ret);
      }
      cl_event finished = 0;
      a64::vector<T> result;
    };
  }
  
  
  /*
   A FIR Filter using the GPU.
   
   Memory transfers to/from the GPU and GPU kernel execution
   overlap with CPU work. When using several levels of asynchronicity,
   it can be possible to ensure that the CPU never waits for the GPU.
   
   The impulse response is split in smaller parts such that
   the sizes of ffts are small enough to be computed on the gpu
   using the "in shared memory only" kernel.
   
   TODO use a separate command queue for reading back the data, and n output memory objects
   so that reading back the data can be done in parallel with writing the input for the next kernel.
   Or is the device already doing that optimization?
   */
  template <typename T>
  struct PartitionnedFIRFilterGPUAsyncN {
    static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value);
    static constexpr auto kernel_file = "fft_mult_ifft_partitionned.cl";
    using Tag = imj::Tag;
    using FPT = T;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;
    using Signal_value_type = typename fft::RealSignal_<Tag, FPT>::value_type;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    static constexpr auto mult_assign = fft::RealFBins_<Tag, FPT>::mult_assign;
    static constexpr auto scale = fft::RealFBins_<Tag, FPT>::scale;
    
    using Algo = typename fft::Algo_<Tag, FPT>;
    using Contexts = fft::Contexts_<Tag, FPT>;
    
    struct SetupParam {
      int n_levels; // max number of simultaneously queued kernels
    };
    
    void applySetup(SetupParam const & s) {
      waitForGPUJobsCompletion();
      
      gpu_works.clear();
      if(s.n_levels <= 0) {
        throw std::invalid_argument("number of levels must be > 0");
      }
      gpu_works.resize(s.n_levels);
    }
    
    void setCoefficients(a64::vector<T> coeffs_) {
      waitForGPUJobsCompletion();
      reset();
      
      if(coeffs_.size() <= 2) {
        coeffs_.resize(2);
      }
      
      int const max_n_coeffs_per_partition = maxCountGPUStockhamFFTCoeffs<T>();
      
      n_partitions = 1 + (coeffs_.size()-1) / max_n_coeffs_per_partition;
      int const n_coeffs_per_partition =
      (n_partitions == 1) ?
      ceil_power_of_two(coeffs_.size()) :
      max_n_coeffs_per_partition;
      int const fft_length = n_coeffs_per_partition * 2;
      
      // pad with zeros
      coeffs_.resize(n_coeffs_per_partition * n_partitions);
      
      N = fft_length/2;
      
      // the signal will be 0-padded.
      // Note that to optimize GPU memory transfer, we could copy just the first half,
      // and change the kernel to initialize shared memory with zeroes for the 2nd half.
      signal.resize(fft_length);
      previous_result.resize(fft_length);
      
      if(gpu_works.empty()) {
        throw std::logic_error("0 levels. did you forget to call applySetup?");
      }
      for(auto & w : gpu_works) {
        w.result.clear();
        w.result.resize(fft_length);
      }
      gpu_work = &gpu_works[0];
      
      fft.setContext(Contexts::getInstance().getBySize(fft_length));
      
      
      std::vector<CplxFreqs> fft_of_hs;
      {
        fft_of_hs.resize(n_partitions);
        auto it_coeffs = coeffs_.begin();
        
        RealSignal coeffs_slice(fft_length, Signal_value_type(0)); // initialize with zeros (second half is padding)
        for(auto & fft_of_h : fft_of_hs) {
          fft_of_h.resize(fft_length);
          auto end_coeffs = it_coeffs + n_coeffs_per_partition;
          assert(end_coeffs <= coeffs_.end());
          auto slice_it = coeffs_slice.begin();
          for(;it_coeffs != end_coeffs; ++it_coeffs, ++slice_it) {
            using RealT = typename RealSignal::value_type;
            *slice_it = RealT(*it_coeffs);
          }
          
          // coeffs_slice is padded with 0, because it is bigger than partition_size
          // and initialized with zeros.
          fft.forward(coeffs_slice.begin(), fft_of_h, fft_length);
          
          // to avoid a division in the last step of the ifft, we scale the frequencies:
          scale(fft_of_h, 1./fft_length);
        }
        assert(it_coeffs == coeffs_.end());
      }
      
      
      // create kernel
      auto kernel_src = read_kernel(kernel_file);
      kern.setup(getOpenCLContext().context,
                 getOpenCLContext().device_id,
                 kernel_src,
                 fft_length);
      cl_int ret;
      
      auto const bufSz = fft_of_hs.size() * fft_length * sizeof(fft_of_hs[0][0]);
      // Create the memory buffer on the device to hold fft_of_h
      fft_of_h_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_ONLY,
                                        bufSz, NULL, &ret);
      CHECK_CL_ERROR(ret);
      if(n_partitions > 1) // if we have a single partition, we don't store fft of x in global memory.
      {
        fft_of_x_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_WRITE,
                                          bufSz, NULL, &ret);
        CHECK_CL_ERROR(ret);
        T zero(0);
        ret = clEnqueueFillBuffer(getOpenCLContext().command_queue, fft_of_x_mem_obj,
                                  &zero, sizeof(zero),
                                  0,
                                  bufSz,
                                  0, NULL, NULL);
        CHECK_CL_ERROR(ret);
      }
      // Create the memory buffer on the device to hold the input and the output
      io_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_WRITE,
                                  signal.size() * sizeof(signal[0]), NULL, &ret);
      CHECK_CL_ERROR(ret);
      
      // assign kernel arguments (io, fft_of_h)
      ret = clSetKernelArg(kern.kernel, 0, sizeof(cl_mem), (void *)&io_mem_obj);
      CHECK_CL_ERROR(ret);
      ret = clSetKernelArg(kern.kernel, 1, sizeof(cl_mem), (void *)&fft_of_h_mem_obj);
      CHECK_CL_ERROR(ret);
      ret = clSetKernelArg(kern.kernel, 2, sizeof(cl_mem), (void *)&fft_of_x_mem_obj);
      CHECK_CL_ERROR(ret);
      const int n_bytes_local_memory = 2*sizeof(T) * 2*signal.size();
      if(n_bytes_local_memory > shared_memory_size(getOpenCLContext().device_id)) {
        // this should never happen, the very notion of partitions was introduced to
        // cope with that issue.
        throw std::logic_error("not enough memory on the gpu device");
      }
      ret = clSetKernelArg(kern.kernel, 3, n_bytes_local_memory, NULL); // local memory
      CHECK_CL_ERROR(ret);
      ret = clSetKernelArg(kern.kernel, 5, sizeof(n_partitions), &n_partitions);
      CHECK_CL_ERROR(ret);
      
      
      int i=0;
      for(auto const & fft_of_h : fft_of_hs) {
        // Copy fft_of_h to its memory buffers.
        // This can crash if the GPU has not enough memory.
        int const sz = fft_of_h.size() * sizeof(fft_of_h[0]);
        ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, fft_of_h_mem_obj,
                                   (i==n_partitions-1) ? CL_TRUE:CL_FALSE, // block only at the last iteration
                                   i * sz, // offset
                                   sz,
                                   &fft_of_h[0],
                                   0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        ++i;
      }
    }
    
    ~PartitionnedFIRFilterGPUAsyncN() {
      waitForGPUJobsCompletion();
      if(fft_of_h_mem_obj) {
        cl_int ret = clReleaseMemObject(fft_of_h_mem_obj);
        CHECK_CL_ERROR(ret);
      }
      if(fft_of_x_mem_obj) {
        cl_int ret = clReleaseMemObject(fft_of_x_mem_obj);
        CHECK_CL_ERROR(ret);
      }
      if(io_mem_obj) {
        cl_int ret = clReleaseMemObject(io_mem_obj);
        CHECK_CL_ERROR(ret);
      }
    }
    
    bool isValid() const {
      return std::is_same<float,T>::value || device_supports_double(getOpenCLContext().device_id);
    }
    
    int getLatency() const {
      return
      N-1 + // we need N inputs before we can start the GPU kernel
      (gpu_works.size()-1)*N; // we have some levels of asynchronicity
    }
    
    T step(T val) {
      signal[curIndex++] = val;
      if(curIndex == N) {
        // signal now contains N elements, so we can process it.
        
        // copy signal to its memory buffer.
        // This can crash if the GPU has not enough memory.
        cl_int ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                          CL_FALSE,
                                          0,
                                          signal.size() * sizeof(signal[0]), // TODO copy less than that because at least half of the vector is 0
                                          &signal[0],
                                          0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        
        // Execute the OpenCL kernel
        
        size_t global_item_size = get_fft_length(N)/(2*kern.nButterfliesPerThread);
        size_t local_item_size = global_item_size;
        ret = clSetKernelArg(kern.kernel, 4, sizeof(input_index), &input_index);
        ++input_index;
        if(input_index == n_partitions) {
          input_index = 0;
        }
        CHECK_CL_ERROR(ret);
        ret = clEnqueueNDRangeKernel(getOpenCLContext().command_queue, kern.kernel, 1, NULL,
                                     &global_item_size,
                                     &local_item_size,
                                     0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        // Read the memory buffer 'io_mem_obj' on the device to the local variable output
        ret = clEnqueueReadBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                  CL_FALSE,
                                  0,
                                  (N*2) * sizeof(decltype(gpu_work->result[0])),
                                  &gpu_work->result[0],
                                  0, NULL,
                                  &gpu_work->finished);
        CHECK_CL_ERROR(ret);
        
        advanceGPUWork();
        
        // wait for the job results to be available
        gpu_work->waitForCompletion();
        
        // previous_result : old new
        // gpu_work        : new future
        //
        // will become:
        //
        // previous_result : (new+new) future
        // gpu_work        : unchanged
        
        for(int i=0; i<N; ++i) {
          previous_result[i] = previous_result[i+N] + gpu_work->result[i];
          previous_result[i+N] = gpu_work->result[i+N];
        }
        
        curIndex = 0;
      }
      return previous_result[curIndex];
    }
    
    template<typename FPT2>
    void stepAddVectorized(FPT2 const * const input_buffer,
                           FPT2 * output_buffer,
                           int nSamples)
    {
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] += step(input_buffer[i]);
        }
    }
      
    double getEpsilon() const {
      return n_partitions * fft::getFFTEpsilon<FPT>(get_fft_length(N)) + 2 * std::numeric_limits<FPT>::epsilon();
    }
    
  private:
    int N = 0;
    int curIndex = 0;
    a64::vector<T> signal;
    a64::vector<T> previous_result;
    
    detail::GPUWork<T> * gpu_work = nullptr;
    std::vector<detail::GPUWork<T>> gpu_works;
    int input_index = 0, n_partitions = 0;
    
    ScopedKernel<T> kern;
    cl_mem io_mem_obj = 0;
    Algo fft;
    cl_mem fft_of_h_mem_obj = 0;
    cl_mem fft_of_x_mem_obj = 0;
    
    static auto get_fft_length(int n) {
      auto N_nonzero_y = 2 * n;
      return ceil_power_of_two(N_nonzero_y);
    }
    
    void advanceGPUWork() {
      if(gpu_work == &gpu_works[gpu_works.size()-1]) {
        gpu_work = &gpu_works[0];
      }
      else {
        ++gpu_work;
      }
    }
    
    void waitForGPUJobsCompletion() const {
      for (auto const & w : gpu_works) {
        w.waitForCompletion();
      }
    }
    
    void reset() {
      input_index = 0;
    }
  };
  
  
  
  // PartitionnedFIRFilterGPUAsyncN is the final step of an "evolution",
  // the classes below are, by chronological order, the different steps of that
  // evolution, introducing one new element at a time.
  
  
  
  // A naïve implementation of a FIR Filter using the GPU.
  // It is (very) unoptimal because we block during memory transfers and gpu kernel execution.
  // Due to the GPU kernel limitation, only small number of coefficients are allowed.
  /*template <typename T>
  struct FIRFilterGPU {
    static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value);
    static constexpr auto kernel_file = "fft_mult_ifft.cl";
    using Tag = imj::Tag;
    using FPT = T;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    static constexpr auto mult_assign = fft::RealFBins_<Tag, FPT>::mult_assign;
    static constexpr auto scale = fft::RealFBins_<Tag, FPT>::scale;
    
    using Algo = typename fft::Algo_<Tag, FPT>;
    using Contexts = fft::Contexts_<Tag, FPT>;
    
    struct SetupParam {};
    
    void applySetup(SetupParam const &) const {}
    
    void setCoefficients(a64::vector<T> coeffs_) {
      if(coeffs_.size() <= 2) {
        coeffs_.resize(2);
      }
      
      int fft_length = get_fft_length(coeffs_.size());
      N = fft_length/2;
      
      // the signal will be 0-padded.
      // Note that to optimize GPU memory transfer, we could copy just the first half,
      // and change the kernel to initialize shared memory with zeroes for the 2nd half.
      signal.resize(fft_length);
      previous_result.resize(fft_length);
      gpu_work.resize(fft_length);
      
      fft.setContext(Contexts::getInstance().getBySize(fft_length));
      
      CplxFreqs fft_of_h;
      fft_of_h.resize(fft_length);
      
      // pad impulse response with 0
      
      coeffs_.resize(fft_length, {});
      
      // compute fft of padded impulse response
      auto coeffs = makeRealSignal(std::move(coeffs_));
      fft.forward(coeffs.begin(), fft_of_h, fft_length);
      
      // to avoid a division in the last step of the ifft, we scale the frequencies:
      scale(fft_of_h, 1./fft_length);
      
      // create kernel
      auto kernel_src = read_kernel(kernel_file);
      kern.setup(getOpenCLContext().context,
                 getOpenCLContext().device_id,
                 kernel_src,
                 fft_length);
      cl_int ret;
      // Create the memory buffer on the device to hold fft_of_h
      fft_of_h_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_ONLY,
                                        fft_of_h.size() * sizeof(typename decltype(fft_of_h)::value_type), NULL, &ret);
      CHECK_CL_ERROR(ret);
      // Create the memory buffer on the device to hold the input and the output
      io_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_WRITE,
                                  signal.size() * sizeof(signal[0]), NULL, &ret);
      CHECK_CL_ERROR(ret);
      
      // assign kernel arguments (io, fft_of_h)
      ret = clSetKernelArg(kern.kernel, 0, sizeof(cl_mem), (void *)&io_mem_obj);
      CHECK_CL_ERROR(ret);
      ret = clSetKernelArg(kern.kernel, 1, sizeof(cl_mem), (void *)&fft_of_h_mem_obj);
      CHECK_CL_ERROR(ret);
      const int n_bytes_local_memory = 2*sizeof(T) * 2*fft_length;
      ret = clSetKernelArg(kern.kernel, 2, n_bytes_local_memory, NULL); // local memory
      CHECK_CL_ERROR(ret);
      
      // Copy fft_of_h to its memory buffers.
      // This can crash if the GPU has not enough memory.
      ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, fft_of_h_mem_obj,
                                 CL_TRUE,
                                 0,
                                 fft_of_h.size() * sizeof(typename decltype(fft_of_h)::value_type),
                                 &fft_of_h[0],
                                 0, NULL, NULL);
    }
    
    ~FIRFilterGPU() {
      if(fft_of_h_mem_obj) {
        cl_int ret = clReleaseMemObject(fft_of_h_mem_obj);
        CHECK_CL_ERROR(ret);
      }
      if(io_mem_obj) {
        cl_int ret = clReleaseMemObject(io_mem_obj);
        CHECK_CL_ERROR(ret);
      }
    }
    
    bool isValid() const {
      return std::is_same<float,T>::value || device_supports_double(getOpenCLContext().device_id);
    }
    
    int getLatency() const { return N-1; }
    
    T step(T val) {
      signal[curIndex++] = val;
      if(curIndex == N) {
        // signal now contains N elements, so we can process it.
        
        // copy signal to its memory buffer.
        // This can crash if the GPU has not enough memory.
        cl_int ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                          CL_FALSE,
                                          0,
                                          signal.size() * sizeof(signal[0]), // TODO copy less than that because at least half of the vector is 0
                                          &signal[0],
                                          0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        
        // Execute the OpenCL kernel
        
        size_t global_item_size = get_fft_length(N)/(2*kern.nButterfliesPerThread);
        size_t local_item_size = global_item_size;
        ret = clEnqueueNDRangeKernel(getOpenCLContext().command_queue, kern.kernel, 1, NULL,
                                     &global_item_size,
                                     &local_item_size,
                                     0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        // Read the memory buffer 'io_mem_obj' on the device to the local variable output
        // This call is blocking. See also 'FIRFilterGPUAsync' that does not block here,
        // at the cost of a higher latency.
        ret = clEnqueueReadBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                  CL_TRUE,
                                  0,
                                  gpu_work.size() * sizeof(decltype(gpu_work[0])),
                                  &gpu_work[0],
                                  0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        
        // previous_result : old new
        // gpu_work        : new future
        //
        // will become:
        //
        // previous_result : (new+new) future
        // gpu_work        : unchanged
        
        for(int i=0; i<N; ++i) {
          previous_result[i] = previous_result[i+N] + gpu_work[i];
          previous_result[i+N] = gpu_work[i+N];
        }
        
        curIndex = 0;
      }
      return previous_result[curIndex];
    }
    
    double getEpsilon() const {
      return fft::getFFTEpsilon<FPT>(get_fft_length(N)) + 2 * std::numeric_limits<FPT>::epsilon();
    }
    
  private:
    int N = 0;
    int curIndex = 0;
    a64::vector<T> signal;
    a64::vector<T> previous_result;
    a64::vector<T> gpu_work;
    ScopedKernel<T> kern;
    cl_mem io_mem_obj = 0;
    Algo fft;
    cl_mem fft_of_h_mem_obj = 0;
    
    static auto get_fft_length(int n) {
      auto N_nonzero_y = 2 * n;
      return ceil_power_of_two(N_nonzero_y);
    }
    
  };*/
  

   // A naïve implementation of a FIR Filter using the GPU.
   // Memory transfers to/from the GPU and GPU kernel execution
   // can overlap with CPU work.
   // Due to the GPU kernel limitation, only small number of coefficients are allowed.
   // See also FIRFilterGPUAsyncN which is a generalization.
  /*
  template <typename T>
  struct FIRFilterGPUAsync {
    static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value);
    static constexpr auto kernel_file = "fft_mult_ifft.cl";
    using Tag = imj::Tag;
    using FPT = T;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    static constexpr auto mult_assign = fft::RealFBins_<Tag, FPT>::mult_assign;
    static constexpr auto scale = fft::RealFBins_<Tag, FPT>::scale;
    
    using Algo = typename fft::Algo_<Tag, FPT>;
    using Contexts = fft::Contexts_<Tag, FPT>;
    
    struct SetupParam {};
    
    void applySetup(SetupParam const &) const {}
    
    void setCoefficients(a64::vector<T> coeffs_) {
      waitForGPUJobsCompletion();
      
      if(coeffs_.size() <= 2) {
        coeffs_.resize(2);
      }
      
      int fft_length = get_fft_length(coeffs_.size());
      N = fft_length/2;
      
      // the signal will be 0-padded.
      // Note that to optimize GPU memory transfer, we could copy just the first half,
      // and change the kernel to initialize shared memory with zeroes for the 2nd half.
      signal.resize(fft_length);
      previous_result.resize(fft_length);
      gpu_work_a.result.resize(fft_length);
      gpu_work_b.result.resize(fft_length);
      gpu_work = &gpu_work_a;

      fft.setContext(Contexts::getInstance().getBySize(fft_length));
      
      CplxFreqs fft_of_h;
      fft_of_h.resize(fft_length);
      
      // pad impulse response with 0
      
      coeffs_.resize(fft_length, {});
      
      // compute fft of padded impulse response
      auto coeffs = makeRealSignal(std::move(coeffs_));
      fft.forward(coeffs.begin(), fft_of_h, fft_length);
      
      // to avoid a division in the last step of the ifft, we scale the frequencies:
      scale(fft_of_h, 1./fft_length);
      
      // create kernel
      auto kernel_src = read_kernel(kernel_file);
      kern.setup(getOpenCLContext().context,
                 getOpenCLContext().device_id,
                 kernel_src,
                 fft_length);
      cl_int ret;
      // Create the memory buffer on the device to hold fft_of_h
      fft_of_h_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_ONLY,
                                        fft_of_h.size() * sizeof(typename decltype(fft_of_h)::value_type), NULL, &ret);
      CHECK_CL_ERROR(ret);
      // Create the memory buffer on the device to hold the input and the output
      io_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_WRITE,
                                  signal.size() * sizeof(signal[0]), NULL, &ret);
      CHECK_CL_ERROR(ret);
      
      // assign kernel arguments (io, fft_of_h)
      ret = clSetKernelArg(kern.kernel, 0, sizeof(cl_mem), (void *)&io_mem_obj);
      CHECK_CL_ERROR(ret);
      ret = clSetKernelArg(kern.kernel, 1, sizeof(cl_mem), (void *)&fft_of_h_mem_obj);
      CHECK_CL_ERROR(ret);
      ret = clSetKernelArg(kern.kernel, 2, 2*sizeof(T) * 2*signal.size(), NULL); // local memory
      CHECK_CL_ERROR(ret);
      
      // Copy fft_of_h to its memory buffers.
      // This can crash if the GPU has not enough memory.
      ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, fft_of_h_mem_obj,
                                 CL_TRUE,
                                 0,
                                 fft_of_h.size() * sizeof(typename decltype(fft_of_h)::value_type),
                                 &fft_of_h[0],
                                 0, NULL, NULL);
    }
    
    ~FIRFilterGPUAsync() {
      waitForGPUJobsCompletion();
      if(fft_of_h_mem_obj) {
        cl_int ret = clReleaseMemObject(fft_of_h_mem_obj);
        CHECK_CL_ERROR(ret);
      }
      if(io_mem_obj) {
        cl_int ret = clReleaseMemObject(io_mem_obj);
        CHECK_CL_ERROR(ret);
      }
    }
    
    bool isValid() const {
      return std::is_same<float,T>::value || device_supports_double(getOpenCLContext().device_id);
    }
    
    int getLatency() const {
      return
      N-1 + // we need N inputs before we can start the GPU kernel
      N; // we have one level of asynchronicity
    }
    
    T step(T val) {
      signal[curIndex++] = val;
      if(curIndex == N) {
        // signal now contains N elements, so we can process it.
        
        // copy signal to its memory buffer.
        // This can crash if the GPU has not enough memory.
        cl_int ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                          CL_FALSE,
                                          0,
                                          signal.size() * sizeof(signal[0]), // TODO copy less than that because at least half of the vector is 0
                                          &signal[0],
                                          0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        
        // Execute the OpenCL kernel
        
        size_t global_item_size = get_fft_length(N)/(2*kern.nButterfliesPerThread);
        size_t local_item_size = global_item_size;
        ret = clEnqueueNDRangeKernel(getOpenCLContext().command_queue, kern.kernel, 1, NULL,
                                     &global_item_size,
                                     &local_item_size,
                                     0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        // Read the memory buffer 'io_mem_obj' on the device to the local variable output
        ret = clEnqueueReadBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                  CL_FALSE,
                                  0,
                                  (N*2) * sizeof(decltype(gpu_work->result[0])),
                                  &gpu_work->result[0],
                                  0, NULL,
                                  &gpu_work->finished);
        CHECK_CL_ERROR(ret);
        
        swapGPUWork();
        
        // wait for the forelast job results to be available
        gpu_work->waitForCompletion();

        // previous_result : old new
        // gpu_work        : new future
        //
        // will become:
        //
        // previous_result : (new+new) future
        // gpu_work        : unchanged
        
        for(int i=0; i<N; ++i) {
          previous_result[i] = previous_result[i+N] + gpu_work->result[i];
          previous_result[i+N] = gpu_work->result[i+N];
        }
        
        curIndex = 0;
      }
      return previous_result[curIndex];
    }
    
    double getEpsilon() const {
      return fft::getFFTEpsilon<FPT>(get_fft_length(N)) + 2 * std::numeric_limits<FPT>::epsilon();
    }
    
  private:
    int N = 0;
    int curIndex = 0;
    a64::vector<T> signal;
    a64::vector<T> previous_result;
    
    detail::GPUWork<T> * gpu_work = nullptr;
    detail::GPUWork<T> gpu_work_a, gpu_work_b;

    ScopedKernel<T> kern;
    cl_mem io_mem_obj = 0;
    Algo fft;
    cl_mem fft_of_h_mem_obj = 0;
    
    static auto get_fft_length(int n) {
      auto N_nonzero_y = 2 * n;
      return ceil_power_of_two(N_nonzero_y);
    }
    
    void swapGPUWork() {
      if(gpu_work == &gpu_work_a) {
        gpu_work = &gpu_work_b;
      }
      else {
        gpu_work = &gpu_work_a;
      }
    }
    
    void waitForGPUJobsCompletion() const {
      gpu_work_a.waitForCompletion();
      gpu_work_b.waitForCompletion();
    }
  };*/

  
   // A naïve implementation of a FIR Filter using the GPU.
   // Memory transfers to/from the GPU and GPU kernel execution
   // overlap with CPU work. We use several levels of asynchronicity,
   // so that it can be possible to ensure that the CPU never waits for the GPU.
   // Due to the GPU kernel limitation, only small number of coefficients are allowed.
  /*
  template <typename T>
  struct FIRFilterGPUAsyncN {
    static_assert(std::is_same<T, float>::value || std::is_same<T, double>::value);
    static constexpr auto kernel_file = "fft_mult_ifft.cl";
    using Tag = imj::Tag;
    using FPT = T;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    static constexpr auto makeRealSignal = fft::RealSignal_<Tag, FPT>::make;
    
    using CplxFreqs = typename fft::RealFBins_<Tag, FPT>::type;
    static constexpr auto mult_assign = fft::RealFBins_<Tag, FPT>::mult_assign;
    static constexpr auto scale = fft::RealFBins_<Tag, FPT>::scale;
    
    using Algo = typename fft::Algo_<Tag, FPT>;
    using Contexts = fft::Contexts_<Tag, FPT>;
    
    struct SetupParam {
      int n_levels; // max number of simultaneously queued kernels
    };
    
    void applySetup(SetupParam const & s) {
      waitForGPUJobsCompletion();

      gpu_works.clear();
      if(s.n_levels <= 0) {
        throw std::invalid_argument("number of levels must be > 0");
      }
      gpu_works.resize(s.n_levels);
    }
    
    void setCoefficients(a64::vector<T> coeffs_) {
      waitForGPUJobsCompletion();

      if(coeffs_.size() <= 2) {
        coeffs_.resize(2);
      }
      
      int fft_length = get_fft_length(coeffs_.size());
      N = fft_length/2;
      
      // the signal will be 0-padded.
      // Note that to optimize GPU memory transfer, we could copy just the first half,
      // and change the kernel to initialize shared memory with zeroes for the 2nd half.
      signal.resize(fft_length);
      previous_result.resize(fft_length);

      if(gpu_works.empty()) {
        throw std::logic_error("0 levels. did you forget to call applySetup?");
      }
      for(auto & w : gpu_works) {
        w.result.clear();
        w.result.resize(fft_length);
      }
      gpu_work = &gpu_works[0];
      
      fft.setContext(Contexts::getInstance().getBySize(fft_length));
      
      CplxFreqs fft_of_h;
      fft_of_h.resize(fft_length);
      
      // pad impulse response with 0
      
      coeffs_.resize(fft_length, {});
      
      // compute fft of padded impulse response
      auto coeffs = makeRealSignal(std::move(coeffs_));
      fft.forward(coeffs.begin(), fft_of_h, fft_length);
      
      // to avoid a division in the last step of the ifft, we scale the frequencies:
      scale(fft_of_h, 1./fft_length);
      
      // create kernel
      auto kernel_src = read_kernel(kernel_file);
      kern.setup(getOpenCLContext().context,
                 getOpenCLContext().device_id,
                 kernel_src,
                 fft_length);
      cl_int ret;
      // Create the memory buffer on the device to hold fft_of_h
      fft_of_h_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_ONLY,
                                        fft_of_h.size() * sizeof(typename decltype(fft_of_h)::value_type), NULL, &ret);
      CHECK_CL_ERROR(ret);
      // Create the memory buffer on the device to hold the input and the output
      io_mem_obj = clCreateBuffer(getOpenCLContext().context, CL_MEM_READ_WRITE,
                                  signal.size() * sizeof(signal[0]), NULL, &ret);
      CHECK_CL_ERROR(ret);
      
      // assign kernel arguments (io, fft_of_h)
      ret = clSetKernelArg(kern.kernel, 0, sizeof(cl_mem), (void *)&io_mem_obj);
      CHECK_CL_ERROR(ret);
      ret = clSetKernelArg(kern.kernel, 1, sizeof(cl_mem), (void *)&fft_of_h_mem_obj);
      CHECK_CL_ERROR(ret);
      ret = clSetKernelArg(kern.kernel, 2, 2*sizeof(T) * 2*signal.size(), NULL); // local memory
      CHECK_CL_ERROR(ret);
      
      // Copy fft_of_h to its memory buffers.
      // This can crash if the GPU has not enough memory.
      ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, fft_of_h_mem_obj,
                                 CL_TRUE,
                                 0,
                                 fft_of_h.size() * sizeof(typename decltype(fft_of_h)::value_type),
                                 &fft_of_h[0],
                                 0, NULL, NULL);
    }
    
    ~FIRFilterGPUAsyncN() {
      waitForGPUJobsCompletion();
      if(fft_of_h_mem_obj) {
        cl_int ret = clReleaseMemObject(fft_of_h_mem_obj);
        CHECK_CL_ERROR(ret);
      }
      if(io_mem_obj) {
        cl_int ret = clReleaseMemObject(io_mem_obj);
        CHECK_CL_ERROR(ret);
      }
    }
    
    bool isValid() const {
      return std::is_same<float,T>::value || device_supports_double(getOpenCLContext().device_id);
    }
    
    int getLatency() const {
      return
      N-1 + // we need N inputs before we can start the GPU kernel
      (gpu_works.size()-1)*N; // we have some levels of asynchronicity
    }
    
    T step(T val) {
      signal[curIndex++] = val;
      if(curIndex == N) {
        // signal now contains N elements, so we can process it.
        
        // copy signal to its memory buffer.
        // This can crash if the GPU has not enough memory.
        cl_int ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                          CL_FALSE,
                                          0,
                                          signal.size() * sizeof(signal[0]), // TODO copy less than that because at least half of the vector is 0
                                          &signal[0],
                                          0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        
        // Execute the OpenCL kernel
        
        size_t global_item_size = get_fft_length(N)/(2*kern.nButterfliesPerThread);
        size_t local_item_size = global_item_size;
        ret = clEnqueueNDRangeKernel(getOpenCLContext().command_queue, kern.kernel, 1, NULL,
                                     &global_item_size,
                                     &local_item_size,
                                     0, NULL, NULL);
        CHECK_CL_ERROR(ret);
        // Read the memory buffer 'io_mem_obj' on the device to the local variable output
        ret = clEnqueueReadBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                  CL_FALSE,
                                  0,
                                  (N*2) * sizeof(decltype(gpu_work->result[0])),
                                  &gpu_work->result[0],
                                  0, NULL,
                                  &gpu_work->finished);
        CHECK_CL_ERROR(ret);
        
        advanceGPUWork();
        
        // wait for the job results to be available
        gpu_work->waitForCompletion();
        
        // previous_result : old new
        // gpu_work        : new future
        //
        // will become:
        //
        // previous_result : (new+new) future
        // gpu_work        : unchanged
        
        for(int i=0; i<N; ++i) {
          previous_result[i] = previous_result[i+N] + gpu_work->result[i];
          previous_result[i+N] = gpu_work->result[i+N];
        }
        
        curIndex = 0;
      }
      return previous_result[curIndex];
    }
    
    double getEpsilon() const {
      return fft::getFFTEpsilon<FPT>(get_fft_length(N)) + 2 * std::numeric_limits<FPT>::epsilon();
    }
    
  private:
    int N = 0;
    int curIndex = 0;
    a64::vector<T> signal;
    a64::vector<T> previous_result;
    
    detail::GPUWork<T> * gpu_work = nullptr;
    std::vector<detail::GPUWork<T>> gpu_works;
    
    ScopedKernel<T> kern;
    cl_mem io_mem_obj = 0;
    Algo fft;
    cl_mem fft_of_h_mem_obj = 0;
    
    static auto get_fft_length(int n) {
      auto N_nonzero_y = 2 * n;
      return ceil_power_of_two(N_nonzero_y);
    }
    
    void advanceGPUWork() {
      if(gpu_work == &gpu_works[gpu_works.size()-1]) {
        gpu_work = &gpu_works[0];
      }
      else {
        ++gpu_work;
      }
    }
    
    void waitForGPUJobsCompletion() const {
      for (auto const & w : gpu_works) {
        w.waitForCompletion();
      }
    }
  };
   */
  
  
}
