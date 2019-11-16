
namespace imajuscule
{
  /*
   Asynchronous convolution on a separate thread.
   */
  template <typename Async>
  struct AsyncCPUConvolution {
      using FPT = typename Async::FPT;

      static constexpr int nCoefficientsFadeIn = Async::nCoefficientsFadeIn;
      static constexpr bool has_subsampling = Async::has_subsampling;

      // We leave room for one more element in worker_2_rt queue so that the worker
      //   can push immediately if its computation is faster than the real-time thread
      //   (which is very unlikely, though).
      static constexpr int room_sz = 1;

      AsyncCPUConvolution() = default;
      
      // not movable, because the worker lambda captures this
      AsyncCPUConvolution(AsyncCPUConvolution && o) = delete;
      AsyncCPUConvolution& operator=(AsyncCPUConvolution &&) = delete;

      ~AsyncCPUConvolution() {
          terminateAsyncJobs();
      }

    struct SetupParam : public Cost {
        using InnerParams = typename Async::SetupParam;
        SetupParam(int inputSubmissionPeriod,
                   int queueSize,
                   InnerParams const & innerParams)
        : inputSubmissionPeriod(inputSubmissionPeriod)
        , queueSize(queueSize)
        , innerParams(innerParams)
        {}
        
        int inputSubmissionPeriod;
        int queueSize;
        InnerParams innerParams;
        
        void logSubReport(std::ostream & os) override {
            os << "submission_period : " << inputSubmissionPeriod << std::endl;
            os << "queue_size : " << queueSize << std::endl;
            innerParams.logSubReport(os);
        }
    };
    
    void applySetup(SetupParam const & s) {
      terminateAsyncJobs();
        
        N = s.inputSubmissionPeriod;
        queueSize = s.queueSize;
        
        algo.applySetup(s.innerParams);
    }
    
    void setCoefficients(a64::vector<FPT> coeffs) {
      terminateAsyncJobs();

      algo.setCoefficients(std::move(coeffs));

        previous_result = std::make_unique<vec>();
        previous_result->resize(N, {});
        
        signal = std::make_unique<vec>();
        signal->resize(N); // no need to 0-initialize

        worker_vec = std::make_unique<vec>();
        worker_vec->resize(N); // no need to 0-initialize
        
        Assert(jobs==0);

        worker_2_rt = std::make_unique<Queue>(queueSize);
        rt_2_worker = std::make_unique<Queue>(queueSize);

        for(int i=0; i<queueSize - room_sz; ++i) {
            auto zeros = std::make_unique<vec>();
            zeros->resize(N, {});
            worker_2_rt->push(zeros);
            ++jobs;
        }

        curIndex = 0;
        
        worker = std::make_unique<std::thread>([this]() {
            while(true) {
                vec_ptr work;
                while(!rt_2_worker->try_pop(work)) {
                    std::this_thread::yield();
                }
                if(unlikely(!work)) {
                    // sentinel
                    vec_ptr sentinel;
                    worker_2_rt->push(std::move(sentinel));
                    return;
                }
                Assert(N == worker_vec->size());
                Assert(N == work->size());
                algo.stepAssignVectorized(work->data(),
                                          worker_vec->data(),
                                          worker_vec->size());
                worker_2_rt->push(std::move(worker_vec));
                worker_vec = std::move(work);
            }
        });
    }
      
      void reset() {
          terminateAsyncJobs();
          algo.reset();
      }
      
      void flushToSilence() {
          terminateAsyncJobs();
          algo.flushToSilence();
      }
      
    bool isZero() const {
        return algo.isZero();
    }
      
    bool isValid() const {
        return N > 0 && queueSize > 0 && algo.isValid();
    }
    
    int getLatency() const {
      return
      N-1 + // we need N inputs before we can submit
      (queueSize-room_sz)*N + // we have some levels of asynchronicity
      algo.getLatency();
    }
      
      int countCoefficients() const {
          return algo.countCoefficients();
      }
    
    FPT step(FPT val) {
      (*signal)[curIndex++] = val;
      if(unlikely(curIndex == N)) {
          submit_signal(std::move(signal));
          signal.swap(previous_result);
          if(unlikely(!try_receive_result(previous_result))) {
              ++error_worker_too_slow;
              previous_result = receive_result();
          }
          
          curIndex = 0;
      }
      return (*previous_result)[curIndex];
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
    template<typename FPT2>
    void stepAssignVectorized(FPT2 const * const input_buffer,
                              FPT2 * output_buffer,
                              int nSamples)
    {
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] = step(input_buffer[i]);
        }
    }
    template<typename FPT2>
    void stepAddInputZeroVectorized(FPT2 * output_buffer,
                                    int nSamples)
    {
        for(int i=0; i<nSamples; ++i) {
            output_buffer[i] += step({});
        }
    }
      
    double getEpsilon() const {
        return algo.getEpsilon();
    }
    
      auto const & getAsyncAlgo() const {
          return algo;
      }
      
      int getResultQueueSize() const {
          if(!worker_2_rt) {
              return 0;
          }
          return worker_2_rt->unsafe_num_elements();
      }

  private:
    int N = 0;
    int curIndex = 0;
      Async algo;
      using vec = std::vector<FPT>;
      using vec_ptr = std::unique_ptr<vec>;
      
      // Let's assume single producer single consumer for the moment.
      // But in the end, when we have multiple convolutions for multiple channels, we may want to
      // do multiple producer, single consumer.
      
      using Queue = atomic_queue::AtomicQueueB2<
      /* T = */ vec_ptr,
      /* A = */ std::allocator<vec_ptr>,
      /* MAXIMIZE_THROUGHPUT */ true,
      /* TOTAL_ORDER = */ true,
      /* SPSC = */ true
      >;
      
      std::unique_ptr<Queue> rt_2_worker, worker_2_rt;
      int jobs = 0;
      int queueSize = 0;
      int error_worker_too_slow = 0;

        vec_ptr signal, previous_result;
      vec_ptr worker_vec;
      std::unique_ptr<std::thread> worker;

    void terminateAsyncJobs() {
        if(!worker) {
            Assert(jobs==0);
            return;
        }
        submit_signal({}); // sentinel
        while(true) {
            flush_results();
            if(jobs == 0) {
                break;
            }
            std::this_thread::yield();
        }
        worker->join();
        worker.reset();
        rt_2_worker.reset();
        worker_2_rt.reset();
    }
      
      void submit_signal(vec_ptr s) {
          rt_2_worker->push(std::move(s));
          ++jobs;
      }
      
      void flush_results() {
          vec_ptr res;
          while(try_receive_result(res)) {
          }
      }
      
      bool try_receive_result(vec_ptr & res) {
          bool success = worker_2_rt->try_pop(res);
          if(success) {
              --jobs;
          }
          return success;
      }
      
      vec_ptr receive_result() {
          vec_ptr res = std::move(worker_2_rt->pop());
          --jobs;
          return res;
      }
          
  };
  
}
