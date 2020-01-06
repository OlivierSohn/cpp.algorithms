template <typename Async>
struct DescAsyncCPUConvolution
{
    static constexpr int nCoefficientsFadeIn = Async::nCoefficientsFadeIn;
    static constexpr bool has_subsampling = Async::has_subsampling;
    static constexpr bool step_can_error = true;
};

namespace imajuscule {
template <typename Async, PolicyOnWorkerTooSlow OnWorkerTooSlow>
struct AlgoAsyncCPUConvolution;

template <typename Async, PolicyOnWorkerTooSlow OnWorkerTooSlow>
struct StateAsyncCPUConvolution {
    using Algo = AlgoAsyncCPUConvolution<typename Async::Algo, OnWorkerTooSlow>;
    using Desc = DescAsyncCPUConvolution<typename Async::Desc>;
    
    using FPT = typename Async::FPT;
    using Tag = typename Async::Tag;

    StateAsyncCPUConvolution() = default;
    
    // not movable, because the worker lambda captures this
    StateAsyncCPUConvolution(StateAsyncCPUConvolution && o) = delete;
    StateAsyncCPUConvolution& operator=(StateAsyncCPUConvolution &&) = delete;
    
    ~StateAsyncCPUConvolution() {
        terminateAsyncJobs();
    }

    MinSizeRequirement setCoefficients(Algo const & algo,
                                       a64::vector<FPT> coeffs) {
        terminateAsyncJobs();
        
        auto async_conv = std::make_unique<Convolution<typename Async::Algo>>();
        async_conv->setup(algo.asyncParams);
        Assert(async_conv->isValid());
        async_conv->setCoefficients(std::move(coeffs));
        epsilon = async_conv->getEpsilon();
        
        int const N = algo.getSubmissionPeriod();
        int const queueSize = algo.getQueueSize();
        
        buffer.clear();
        buffer.resize(N * // size of a block
                      (1+ // previous_result
                       1+ // signal
                       1+ // worker_vec
                       queueSize-Algo::queue_room_sz)
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
        
        for(int j=0; j<queueSize - Algo::queue_room_sz; ++j) {
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
        
        worker = std::make_unique<std::thread>([this,
                                                async_conv=std::move(async_conv),
                                                N]()
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
                    for(int i=0; i<N; ++i) {
                        bufferData[worker_vec + i] = async_conv->step(bufferData[work + i]);
                    }
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
        
        return {
            1, // x size : pour l'instant on copie x par x donc un seul suffit.
            1, // y size : on écrit un par un pour l'instant
            0, // y anticipé : on n'écrit pas dans le futur
            {} // ffts : besoin de rien pour l'instant puisque les fft sont recalculées dans la partie async
        };
    }
    
    void reset() {
        terminateAsyncJobs();
        Assert(jobs == 0);
        error_worker_too_slow = false;
        buffer = {};
    }
    
    void flushToSilence() {
        terminateAsyncJobs();
    }
        
    double getEpsilon(Algo const & algo) const {
        return epsilon;
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
    
    int curIndex = 0; // todo redondant avec x.progress
private:
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
public:
    std::condition_variable rt_2_worker_cond;
    bool error_worker_too_slow = false;
    
    block_index signal, previous_result;
private:
    int jobs = 0;
    block_index worker_vec;
public:
    a64::vector<FPT> buffer; // block_index are indices into this buffer.
private:
    std::unique_ptr<std::thread> worker;
    
#ifdef IMJ_WITH_ASYNCCONV_STATS
public:
    std::vector<profiling::CpuDuration> async_durations;
private:
#endif
    double epsilon = 0.;

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
    
public:
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


template <typename Async, PolicyOnWorkerTooSlow OnWorkerTooSlow>
struct AlgoAsyncCPUConvolution {
    using State = StateAsyncCPUConvolution<typename Async::State, OnWorkerTooSlow>;
    using Desc = DescAsyncCPUConvolution<typename Async::Desc>;

    using FPT = typename Async::FPT;
    using Tag = typename Async::Tag;
    
    using RealSignal = typename fft::RealSignal_<Tag, FPT>::type;
    
    static constexpr auto get_signal = fft::RealSignal_<Tag, FPT>::get_signal;

    // We leave room for one more element in worker_2_rt queue so that the worker
    //   can push immediately if its computation is faster than the real-time thread
    //   (which is very unlikely, though).
    static constexpr int queue_room_sz = 1;
    
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
        
        int getImpliedLatency() const {
            return
            inputSubmissionPeriod*((queueSize-queue_room_sz) + 1)
            - 1
            + innerParams.getImpliedLatency();
        }
        
        void logSubReport(std::ostream & os) const override {
            os << "Async, period : " << inputSubmissionPeriod <<  " size : " << queueSize << std::endl;
            {
                IndentingOStreambuf i(os);
                innerParams.logSubReport(os);
            }
        }
    };
    
    void setup(SetupParam const & s) {
        N = s.inputSubmissionPeriod;
        queueSize = s.queueSize;
        asyncParams = s.innerParams;
    }
    
    void reset() {
        asyncParams = {};
        queueSize = 0;
        N = 0;
    }
    
    bool isValid() const {
        return N > 0 && queueSize > 0 && asyncParams.isValid();
    }
    
    int getLatency() const {
        return
        N*((queueSize-queue_room_sz) // we have some levels of asynchronicity
           +1) -1 // we need N inputs before we can submit
        + asyncParams.getLatency();
    }
    
    void step(State & s,
              XAndFFTS<FPT, Tag> const & x_and_ffts,
              Y<FPT, Tag> & y) const
    {
        /*
         si une nouvelle fft est dispo, mettre une demande dans la file vers le worker avec les infos sur la fft.
         Il faudra que cette fft reste dispo au moins pendant le délai de traitement theorique maximum.
         */
        if constexpr (OnWorkerTooSlow == PolicyOnWorkerTooSlow::PermanentlySwitchToDry) {
            if(unlikely(s.error_worker_too_slow)) {
                y.y[y.progress] += x_and_ffts.x[x_and_ffts.progress-1];
                return;
            }
        }
        s.buffer[s.signal+(s.curIndex++)] = get_signal(x_and_ffts.x[x_and_ffts.progress-1]);
        if(unlikely(s.curIndex == N))
        {
            s.curIndex = 0;
            while(unlikely(!s.try_submit_signal(s.signal))) {
                s.error_worker_too_slow = true;
                if constexpr (OnWorkerTooSlow == PolicyOnWorkerTooSlow::PermanentlySwitchToDry) {
                    y.y[y.progress] += x_and_ffts.x[x_and_ffts.progress-1];
                    return;
                }
                else {
                    s.rt_2_worker_cond.notify_one(); // in case this is because of the race condition
                }
            }
            std::swap(s.signal, s.previous_result);
            while(unlikely(!s.try_receive_result(s.previous_result))) {
                s.error_worker_too_slow = true;
                if constexpr (OnWorkerTooSlow == PolicyOnWorkerTooSlow::PermanentlySwitchToDry) {
                    y.y[y.progress] += x_and_ffts.x[x_and_ffts.progress-1];
                    return;
                }
                else {
                    s.rt_2_worker_cond.notify_one(); // in case this is because of the race condition
                }
            }
        }
        y.y[y.uProgress] += typename RealSignal::value_type(s.buffer[s.previous_result + s.curIndex]);
    }
    
    int getQueueSize() const {
        return queueSize;
    }
    int getSubmissionPeriod() const {
        return N;
    }
    
    typename Async::SetupParam asyncParams;
private:
    int N = 0;
    int queueSize = 0;
};

}
