
namespace imajuscule
{
  /*
   Asynchronous convolution on a separate thread.
   */
  template <typename Async>
  struct AsyncCPUConvolution {
      using FPT = typename Async::FPT;

      static constexpr int nComputePhaseable = Async::nComputePhaseable;
      static constexpr int nCoefficientsFadeIn = Async::nCoefficientsFadeIn;
      static constexpr bool has_subsampling = Async::has_subsampling;

      // We leave room for one more element in worker_2_rt queue so that the worker
      //   can push immediately if its computation is faster than the real-time thread
      //   (which is very unlikely, though).
      static constexpr int queue_room_sz = 1;

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
    
    void setup(SetupParam const & s) {
      terminateAsyncJobs();
        
        N = s.inputSubmissionPeriod;
        queueSize = s.queueSize;
        
        algo.setup(s.innerParams);
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

        for(int i=0; i<queueSize - queue_room_sz; ++i) {
            auto zeros = std::make_unique<vec>();
            zeros->resize(N, {});
            worker_2_rt->push(std::move(zeros));
            ++jobs;
        }
        
        // by now:
        // #rt_2_worker 0
        // #worker_2_rt queueSize - queue_room_sz
        //
        // later, we will submit and receive (in that order) in the same iteration.
        //
        // so in the worst case, if "the worker is doing nothing" (not reading from rt_2_worker),
        // we will empty worker_2_rt:
        //
        // #rt_2_worker queueSize - queue_room_sz
        // #worker_2_rt 0
        //
        // and then, we will submit one more to rt_2_worker, then block on reading from worker_2_rt:
        // #rt_2_worker queueSize
        // #worker_2_rt 0
        //
        // Note that the blocking could also occur ** before **, when writing to rt_2_worker
        // if a sentinel has been written to rt_2_worker just before.
        //
        // the situation will unblock when the worker reads from rt_2_worker:
        // #rt_2_worker queueSize - queue_room_sz
        // #worker_2_rt 0
        // and then writes to worker_2_rt:
        // #rt_2_worker queueSize - queue_room_sz
        // #worker_2_rt 1

        curIndex = 0;
        
        worker = std::make_unique<std::thread>([this]()
        {
            std::mutex dummyMut; // used for rt_2_worker_cond.
            auto popWork = [&dummyMut, this]() -> vec_ptr
            {
                vec_ptr work;
                
                auto try_pop = [this, &work] {
                    return rt_2_worker->try_pop(work);
                };
                
                // We look in the queue before yielding, for performance in case
                // the queue is not empty.
                
                if(try_pop()) {
                    return std::move(work);
                }
                
                // !!! We yield here to reduce the likelyhood of a preemption
                // at the location of the race condition hereunder.

                std::this_thread::yield();
                if(try_pop()) {
                    return std::move(work);
                }
                
                // !!! Race condition if the producer notifies HERE:
                // the consumer is not yet waiting so it won't be woken up.
                //
                // To reduce the likelihood of that race condition, we yield just before try_pop()
                // so the scheduler will probably not preempt our thread now because it has just begun
                // its time slice.
                //
                // Still, the race exists even without preemption, in case of multi-core concurrence.
                //
                // If the race condition occurs, the worker will miss the wake up,
                // and will wake up at the next notification (one audio callback later)
                //
                // To ensure that this will not generate an audio dropout,
                // we augment the size of the queue by the number of frames in an audio callback.
                
                while(true)
                {
                    std::unique_lock l(dummyMut);
                    rt_2_worker_cond.wait(l);

                    // This could be :
                    // - a real wake-up :
                    //     - the producer has pushed to the queue
                    // - a race-condition wake-up :
                    //     - the producer has pushed to the queue, but since the producer didn't own the mutex 'dummyMut'
                    //       when pushing to the queue, there is a race and from the point of view of the consumer thread,
                    //       the queue may not have changed yet.
                    // - a spurious wake-up :
                    //     - the producer has not pushed to the queue
                    //
                    // To account for race-condition wake-up, we retry to read from the queue for some time, and then
                    // to account for spurious wake-ups, we return to condition-variable wait.

                    if(likely(try_pop())) {
                        // This was a real wake-up
                        return std::move(work);
                    }
                    for(int busy_wait=10; busy_wait; --busy_wait) {
                        if(try_pop()) {
                            // this was a race-condition wake-up
                            return std::move(work);
                        };
                    }
                    for(auto dt = std::chrono::nanoseconds(100);
                        dt < std::chrono::microseconds(200);
                        dt *= 10) {
                        std::this_thread::sleep_for(dt);
                        if(try_pop()) {
                            // this was a race-condition wake-up (the information took way longer to arrive)
                            return std::move(work);
                        }
                    }
                    // this was a spurious wake-up,
                    // or a race-condition wake-up and it takes more than 200 microseconds for cpu stuff to synchronize?
                }
            };
            
            while(true) {
                
                auto work = popWork();
                
                if(unlikely(!work)) {
                    // sentinel
                    vec_ptr sentinel;
                    auto res = worker_2_rt->try_push(std::move(sentinel));
                    Assert(res);
                    return;
                }
                Assert(N == worker_vec->size());
                Assert(N == work->size());
                algo.stepAssignVectorized(work->data(),
                                          worker_vec->data(),
                                          worker_vec->size());
                auto res = worker_2_rt->try_push(std::move(worker_vec));
                Assert(res);
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
          N*((queueSize-queue_room_sz) // we have some levels of asynchronicity
             +1) -1 // we need N inputs before we can submit
          + algo.getLatency();
      }
      
      std::array<int, nComputePhaseable> getComputePeriodicities() const {
          return algo.getComputePeriodicities();
      }
      // in [0, getComputePeriodicity())
      std::array<int, nComputePhaseable> getComputeProgresses() const {
          return algo.getComputeProgresses();
      }
      void setComputeProgresses(std::array<int, nComputePhaseable> const & progresses) {
          algo.setComputeProgresses(progresses);
      }

      int countCoefficients() const {
          return algo.countCoefficients();
      }
    
    FPT step(FPT val) {
      (*signal)[curIndex++] = val;
      if(unlikely(curIndex == N))
      {
          while(unlikely(!try_submit_signal(std::move(signal)))) {
              ++error_worker_too_slow;
              rt_2_worker_cond.notify_one(); // in case this is because of the race condition
          }
          signal.swap(previous_result);
          while(unlikely(!try_receive_result(previous_result))) {
              ++error_worker_too_slow;
              rt_2_worker_cond.notify_one(); // in case this is because of the race condition
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
      int getSignalQueueSize() const {
          if(!rt_2_worker) {
              return 0;
          }
          return rt_2_worker->unsafe_num_elements();
      }
      
      int countErrorsWorkerTooSlow() const {
          return error_worker_too_slow;
      }
      
      int getQueueSize() const {
          return queueSize;
      }
      int getSubmissionPeriod() const {
          return N;
      }

  private:
    int N = 0;
    int curIndex = 0;
      Async algo;
      using vec = a64::vector<FPT>;
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
      std::condition_variable rt_2_worker_cond;
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
            rt_2_worker_cond.notify_one(); // in case this is because of the race condition
            std::this_thread::yield();
        }
        worker->join();
        worker.reset();
        rt_2_worker.reset();
        worker_2_rt.reset();
    }

      void submit_signal(vec_ptr s) {
          rt_2_worker->push(std::move(s));

          rt_2_worker_cond.notify_one();
          ++jobs;
      }

      bool try_submit_signal(vec_ptr && s) {
          bool res = rt_2_worker->try_push(std::forward<vec_ptr>(s));
          if(likely(res)) {
              rt_2_worker_cond.notify_one();
              ++jobs;
          }
          return res;
      }
      
      void flush_results() {
          vec_ptr res;
          while(try_receive_result(res)) {
          }
      }
      
      bool try_receive_result(vec_ptr & res) {
          bool success = worker_2_rt->try_pop(res);
          if(likely(success)) {
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
