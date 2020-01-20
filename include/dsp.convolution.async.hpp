
namespace imajuscule
{
enum class PolicyOnWorkerTooSlow {
    PermanentlySwitchToDry,
    Wait // for testing purposes
};

struct AsyncCPUConvolutionConstants {
    // We leave room for one more element in worker_2_rt queue so that the worker
    //   can push immediately if its computation is faster than the real-time thread
    //   (which is very unlikely, though).
    static constexpr int queue_room_sz = 1;
};

template<typename InnerParams>
struct AsyncSetupParam : public Cost {
    static constexpr int queue_room_sz = AsyncCPUConvolutionConstants::queue_room_sz;

    AsyncSetupParam(int inputSubmissionPeriod,
               int queueSize,
               InnerParams const & asyncParams)
    : inputSubmissionPeriod(inputSubmissionPeriod)
    , queueSize(queueSize)
    , asyncParams(asyncParams)
    {}
    
    int inputSubmissionPeriod;
    int queueSize;
    InnerParams asyncParams;
    
    bool handlesCoefficients() const {
        return queueSize > 0 && inputSubmissionPeriod > 0 && asyncParams.handlesCoefficients();
    }

    Latency getImpliedLatency() const {
        Assert(handlesCoefficients());
        return asyncParams.getImpliedLatency() +
        Latency( inputSubmissionPeriod*((queueSize-queue_room_sz) + 1) - 1);
    }
    
    void adjustWork(int targetNCoeffs) {
        asyncParams.adjustWork(targetNCoeffs);
        if(!asyncParams.handlesCoefficients()) {
            queueSize = 0; // not sure if that is needed
        }
        if(!handlesCoefficients()) {
            setCost(0.);
        }
    }
    
    void logSubReport(std::ostream & os) const override {
        os << "Async, period : " << inputSubmissionPeriod <<  " size : " << queueSize << std::endl;
        {
            IndentingOStreambuf i(os);
            asyncParams.logSubReport(os);
        }
    }
    
    template<typename F>
    void forEachUsingSameContext(F f) const {
        f(*this);
        // do not recurse into asyncParams because the async part uses another context.
    }
};


  /*
   Asynchronous convolution on a separate thread.
   */
  template <typename Async, PolicyOnWorkerTooSlow OnWorkerTooSlow>
  struct AsyncCPUConvolution {
      using FPT = typename Async::FPT;
      
      using AsyncPart = Async;

      static constexpr int nComputePhaseable = Async::nComputePhaseable;
      static constexpr int nCoefficientsFadeIn = Async::nCoefficientsFadeIn;
      static constexpr bool has_subsampling = Async::has_subsampling;
      static constexpr bool step_can_error = true;
      static constexpr int queue_room_sz = AsyncCPUConvolutionConstants::queue_room_sz;

      AsyncCPUConvolution() = default;
      
      // not movable, because the worker lambda captures this
      AsyncCPUConvolution(AsyncCPUConvolution && o) = delete;
      AsyncCPUConvolution& operator=(AsyncCPUConvolution &&) = delete;

      ~AsyncCPUConvolution() {
          terminateAsyncJobs();
      }

    using SetupParam = AsyncSetupParam<typename Async::SetupParam>;

    void logComputeState(std::ostream & os) const {
        os << "Async [" << curIndex << "/" << N << "], queueSize : " << queueSize << std::endl;
        {
            IndentingOStreambuf i(os);
            algo.logComputeState(os);
        }
    }
      
    void setup(SetupParam const & s) {
      terminateAsyncJobs();
        
        N = s.inputSubmissionPeriod;
        queueSize = s.queueSize;
        
        algo.setup(s.asyncParams);
    }
    
    void setCoefficients(a64::vector<FPT> coeffs) {
      terminateAsyncJobs();

      algo.setCoefficients(std::move(coeffs));
        
        buffer.clear();
        buffer.resize(
                      N * // size of a block
                      (1+ // previous_result
                       1+ // signal
                       1+ // worker_vec
                       queueSize-queue_room_sz)
                      , {});
        
        int i=0;
        previous_result = i;
        i += N;
        signal = i;
        i += N;
        worker_vec = i;
        i += N;
        
        Assert(jobs==0);

        worker_2_rt = std::make_unique<Queue>(queueSize);
        rt_2_worker = std::make_unique<Queue>(queueSize);

        for(int j=0; j<queueSize - queue_room_sz; ++j) {
            worker_2_rt->push(i);
            i += N;
            ++jobs;
        }
        Assert(i==buffer.size());
        
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
        
        worker = std::make_unique<std::thread>([this]()
        {
            std::mutex dummyMut; // used for rt_2_worker_cond.
            auto popWork = [&dummyMut, this]() -> block_index
            {
                block_index work;
                
                auto try_pop = [this, &work] {
                    return rt_2_worker->try_pop(work);
                };
                
                // We look in the queue before yielding, for performance in case
                // the queue is not empty.
                
                if(try_pop()) {
                    return work;
                }
                
                // !!! We yield here to reduce the likelyhood of a preemption
                // at the location of the race condition hereunder.

                std::this_thread::yield();
                if(try_pop()) {
                    return work;
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
                        return work;
                    }
                    for(int busy_wait=10; busy_wait; --busy_wait) {
                        if(try_pop()) {
                            // this was a race-condition wake-up
                            return work;
                        };
                    }
                    for(auto dt = std::chrono::nanoseconds(100);
                        dt < std::chrono::microseconds(200);
                        dt *= 10) {
                        std::this_thread::sleep_for(dt);
                        if(try_pop()) {
                            // this was a race-condition wake-up (the information took way longer to arrive)
                            return work;
                        }
                    }
                    // this was a spurious wake-up,
                    // or a race-condition wake-up and it takes more than 200 microseconds for cpu stuff to synchronize?
                }
            };
            
            int const N = this->N;
            auto const bufferData = this->buffer.data();
            
            while(true) {
                
                auto work = popWork();
#ifdef IMJ_WITH_ASYNCCONV_STATS
                std::optional<profiling::CpuDuration> threadCPUDuration;
                {
                    profiling::ThreadCPUTimer tt(threadCPUDuration);
#endif // IMJ_WITH_ASYNCCONV_STATS
                    if(unlikely(work==sentinel)) {
                        auto res = worker_2_rt->try_push(sentinel);
                        Assert(res);
                        return;
                    }
                    algo.stepAssignVectorized(bufferData + work,
                                              bufferData + worker_vec,
                                              N);
                    auto res = worker_2_rt->try_push(worker_vec);
                    Assert(res);
                    worker_vec = work;
#ifdef IMJ_WITH_ASYNCCONV_STATS
                }
                if(!threadCPUDuration) {
                    throw std::runtime_error("cannot measure time");
                }
                async_durations.emplace_back(*threadCPUDuration);
#endif // IMJ_WITH_ASYNCCONV_STATS
            }
        });
    }
      
      void reset() {
          terminateAsyncJobs();
          Assert(jobs == 0);
          algo.reset();
          queueSize = 0;
          error_worker_too_slow = false;
          N = 0;
          buffer = {};
      }
      
      void flushToSilence() {
          terminateAsyncJobs();
          algo.flushToSilence();
      }
      
    bool isZero() const {
        return algo.isZero();
    }
      
      bool isValid() const {
          if(!algo.isValid()) {
              return false;
          }
          if(algo.handlesCoefficients()) {
              return N > 0 && queueSize > 0;
          }
          return queueSize == 0;
      }

      bool handlesCoefficients() const {
          return N > 0 && queueSize > 0 && algo.handlesCoefficients();
      }
      
      Latency getLatency() const {
          Assert(handlesCoefficients());
          return
          Latency(
                  N*((queueSize-queue_room_sz) // we have some levels of asynchronicity
             +1) -1 // we need N inputs before we can submit
          ) + algo.getLatency();
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
    
    FPT step(FPT val) {
      if constexpr (OnWorkerTooSlow == PolicyOnWorkerTooSlow::PermanentlySwitchToDry) {
          if(unlikely(error_worker_too_slow)) {
              return val;
          }
      }
      buffer[signal+(curIndex++)] = val;
      if(unlikely(curIndex == N))
      {
          curIndex = 0;
          while(unlikely(!try_submit_signal(signal))) {
              error_worker_too_slow = true;
              if constexpr (OnWorkerTooSlow == PolicyOnWorkerTooSlow::PermanentlySwitchToDry) {
                  return val;
              }
              else {
                  rt_2_worker_cond.notify_one(); // in case this is because of the race condition
              }
          }
          std::swap(signal, previous_result);
          while(unlikely(!try_receive_result(previous_result))) {
              error_worker_too_slow = true;
              if constexpr (OnWorkerTooSlow == PolicyOnWorkerTooSlow::PermanentlySwitchToDry) {
                  return val;
              }
              else {
                  rt_2_worker_cond.notify_one(); // in case this is because of the race condition
              }
          }
      }
      return buffer[previous_result + curIndex];
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
      
      bool hasStepErrors() const {
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

      using block_index = int;
      static constexpr block_index sentinel = -1;
      
      // Let's assume single producer single consumer for the moment.
      // But in the end, when we have multiple convolutions for multiple channels, we may want to
      // do multiple producer, single consumer.
      
      using Queue = atomic_queue::AtomicQueueB2<
      /* T = */ block_index,
      /* A = */ std::allocator<block_index>,
      /* MAXIMIZE_THROUGHPUT */ true,
      /* TOTAL_ORDER = */ true,
      /* SPSC = */ true
      >;
      
      std::unique_ptr<Queue> rt_2_worker, worker_2_rt;
      std::condition_variable rt_2_worker_cond;
      int jobs = 0;
      int queueSize = 0;
      bool error_worker_too_slow = false;

      block_index signal, previous_result;
      block_index worker_vec;
      a64::vector<FPT> buffer; // block_index are indices into this buffer.
      std::unique_ptr<std::thread> worker;

#ifdef IMJ_WITH_ASYNCCONV_STATS
  public:
      std::vector<profiling::CpuDuration> async_durations;
  private:
#endif

    void terminateAsyncJobs() {
        curIndex = 0;
        if(!worker) {
            Assert(jobs==0);
            return;
        }
        submit_signal(sentinel);
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

      void submit_signal(block_index s) {
          rt_2_worker->push(s);

          rt_2_worker_cond.notify_one();
          ++jobs;
      }

      bool try_submit_signal(block_index s) {
          bool res = rt_2_worker->try_push(s);
          if(likely(res)) {
              rt_2_worker_cond.notify_one();
              ++jobs;
          }
          return res;
      }
      
      void flush_results() {
          block_index res;
          while(try_receive_result(res)) {
          }
      }
      
      bool try_receive_result(block_index & res) {
          bool success = worker_2_rt->try_pop(res);
          if(likely(success)) {
              --jobs;
          }
          return success;
      }
      
      block_index receive_result() {
          block_index res = std::move(worker_2_rt->pop());
          --jobs;
          return res;
      }
          
  };
  
}
