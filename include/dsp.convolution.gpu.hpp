
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
      
      cl_int ret;
      
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
        snprintf(buf, sizeof(buf), "%a", (float)(-M_PI/nButterflies));
        
        std::string const replaced_str = ReplaceString(ReplaceString(ReplaceString(ReplaceString(kernel_src,
                                                                                                 "replace_MINUS_PI_over_N_GLOBAL_BUTTERFLIES",
                                                                                                 buf),
                                                                                   "replace_N_GLOBAL_BUTTERFLIES",
                                                                                   std::to_string(nButterflies)),
                                                                     "replace_LOG2_N_GLOBAL_BUTTERFLIES",
                                                                     std::to_string(power_of_two_exponent(nButterflies))),
                                                       "replace_N_LOCAL_BUTTERFLIES",
                                                       std::to_string(nButterfliesPerThread));
        size_t const replaced_source_size = replaced_str.size();
        const char * rep_src = replaced_str.data();
        
        // Create a program from the kernel source
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
  
  /*
   A naïve implementation of a FIR Filter using the GPU.
   
   It is (very) unoptimal because we block during memory transfers and gpu kernel execution.
   
   Due to the GPU kernel limitation, only small number of coefficients are allowed.
   */
  template <typename T>
  struct FIRFilterGPU {
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
      ret = clSetKernelArg(kern.kernel, 2, 2*sizeof(float) * 2*signal.size(), NULL); // local memory
      CHECK_CL_ERROR(ret);
      
      // Copy fft_of_h to its memory buffers.
      // This can crash if the GPU has not enough memory.
      ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, fft_of_h_mem_obj,
                                 CL_TRUE /* blocking call */,
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
      return true;
    }
    
    int getLatency() const { return N-1; }
    
    T step(T val) {
      signal[curIndex++] = val;
      if(curIndex == N) {
        // signal now contains N elements, so we can process it.
        
        // copy signal to its memory buffer.
        // This can crash if the GPU has not enough memory.
        cl_int ret = clEnqueueWriteBuffer(getOpenCLContext().command_queue, io_mem_obj,
                                          CL_TRUE /* blocking call */,
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
                                  CL_TRUE /* blocking call */,
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
    ScopedKernel kern;
    cl_mem io_mem_obj = 0;
    Algo fft;
    cl_mem fft_of_h_mem_obj = 0;

    static auto get_fft_length(int n) {
      auto N_nonzero_y = 2 * n;
      return ceil_power_of_two(N_nonzero_y);
    }
    
  };

}
