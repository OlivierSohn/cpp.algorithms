namespace imajuscule {

template <typename Async>
struct DescAsyncCPUConvolution
{
    static constexpr int nCoefficientsFadeIn = Async::nCoefficientsFadeIn;
    static constexpr bool has_subsampling = Async::has_subsampling;
    static constexpr bool step_can_error = true;
    
    /*
     Will be usefull if, in the future, ffts are shared between synchronous and asynchronous parts.
     */
    static int const getMinFftRingbufferSize(int const queue_size,
                                             int const submission_period,
                                             int const fft_length) {
        // Ffts are stored in a ring buffer:
        // - The  synchronous thread writes and reads the buffer
        // - The asynchronous thread reads the buffer
        //
        // Since there is no locking mechanism (to avoid priority inversions) between synchronous and asynchronous threads,
        // the ring buffer needs to be big enough to ensure that "in the worst case" (see below),
        // the ffts read by the asynchronous worker have not yet been overwritten.
        //
        // The worst case is when:
        // - async worker has popped buffer 'W' from rt_2_worker queue,
        // - rt_2_worker queue is full (it contains buffers 'W+1' to 'W+queue_size', included), and
        // - the synchronous thread has already computed ffts corresponding to buffer 'W+queue_size+1'
        //
        // Hence, the minimum ring buffer size is the maximum number of ffts of length fft_length generated for
        //  'S = (queue_size+2)*submission_period' consecutive samples:
        //  minRingSize = 1 + (S-1) / fft_compute_period
        //  (where fft_compute_period = fft_length/2)
        
        Assert(is_power_of_two(fft_length));
        
        int const fft_compute_period = fft_length/2;
        
        int const S = (queue_size+2) * submission_period;
        
        Assert(fft_compute_period > 0);
        
        int const minRingSize = 1 + (S-1) / fft_compute_period;
        
        return minRingSize;
    }
};

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

private:
    using block_index = int;
    static constexpr block_index worker_vec_initial = 0;
public:
    
    MinSizeRequirement setCoefficients(Algo const & algo,
                                       a64::vector<FPT> coeffs) {
        reset();
        
        async_conv = std::make_unique<Convolution<typename Async::Algo>>();
        Assert(algo.asyncParams);
        async_conv->setup(*algo.asyncParams);
        Assert(async_conv->isValid());
        async_conv->setCoefficients(std::move(coeffs));
        
        int const N = algo.getSubmissionPeriod();
        int const queueSize = algo.getQueueSize();
        
        buffer.resize(N * // size of a block
                      (1+ // previous_result
                       1+ // signal
                       1+ // worker_vec
                       queueSize-Algo::queue_room_sz)
                      , {});
        
        int i=0;
        // worker_vec is now a variable in the worker thread to avoid false sharing
        //worker_vec = i;
        Assert(i == worker_vec_initial);
        i += N;
        signal = i;
        i += N;
        previous_result = i;
        i += N;
        
        Assert(jobs==0);
        
        worker_2_rt = Queue(queueSize);
        rt_2_worker = Queue(queueSize);
        
        for(int j=0; j<queueSize - Algo::queue_room_sz; ++j) {
            worker_2_rt.push(i);
            i += N;
            ++jobs;
        }
        Assert(i==buffer.size());
        
        //    ... by now, buffer is partitionned like this:
        //
        // |worker_vec     |signal         |previous_result|work_res_1     |work_res_2     |...|work_res_q-1   |
        //                  >               >
        //                  sync write      sync read
        //
        //    ... after one rt period:
        //
        // |worker_vec     |work_input_x   |signal         |previous_result|work_res_2     |...|work_res_q-1   |
        //  !               !               >               >
        //  worker write    worker read     sync write      sync read
        //
        //    ... after another rt period (assuming worker has not been scheduled):
        //
        // |worker_vec     |work_input_x   |work_input_x+1 |signal         |previous_result|...|work_res_q-1   |
        //  !               !                               >               >
        //  worker write    worker read                     sync write      sync read
        //
        //    ... and after worker has finished his job:
        //
        // |work_res_q     |worker_vec     |work_input_x+1 |signal         |previous_result|...|work_res_q-1   |
        //                  !               !               >               >
        //                  worker write    worker read     sync write      sync read
        //
        //
        // TODO to avoid false sharing (worker and sync writing to the same cache line),
        //   it is important that a block has a size of at least one cache line.


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
        // (depending on PolicyOnWorkerTooSlow)
        // the situation will unblock when the worker reads from rt_2_worker:
        // #rt_2_worker queueSize - queue_room_sz
        // #worker_2_rt 0
        // and then writes to worker_2_rt:
        // #rt_2_worker queueSize - queue_room_sz
        // #worker_2_rt 1
        
        worker = std::make_unique<std::thread>([this,
                                                N]()
                                               {
            std::mutex dummyMut; // used for rt_2_worker_cond.
            auto popWork = [&dummyMut, this]() -> block_index
            {
                block_index work;
                
                auto try_pop = [this, &work] {
                    return rt_2_worker.try_pop(work);
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
                // we have augmented the size of the queue by the number of frames in an audio callback.
                
                while(true)
                {
                    std::unique_lock l(dummyMut);
                    rt_2_worker_cond.wait(l);
                    
                    // It could be :
                    // - a real wake-up :
                    //     - the producer has pushed to the queue
                    // - a race-condition wake-up :
                    //     - the producer has pushed to the queue, but since the producer didn't own the mutex 'dummyMut'
                    //       when pushing to the queue, there is a race and from the point of view of the consumer thread,
                    //       the queue may not have changed yet.
                    // - a spurious wake-up :
                    //     - the producer has not pushed to the queue
                    
                    if(likely(try_pop())) {
                        // It was a real wake-up
                        return work;
                    }
                    // To account for race-condition wake-up, we retry to read from the queue for some time, and then
                    for(int busy_wait=10; busy_wait; --busy_wait) {
                        if(try_pop()) {
                            // It was a race-condition wake-up
                            return work;
                        };
                    }
                    for(auto dt = std::chrono::nanoseconds(100);
                        dt < std::chrono::microseconds(200);
                        dt *= 10) {
                        std::this_thread::sleep_for(dt);
                        if(try_pop()) {
                            // It was a race-condition wake-up (the information took a long time to arrive)
                            return work;
                        }
                    }
                    // It was either a spurious wake-up,
                    // or a race-condition wake-up and it took more than 200 microseconds for the cpu to synchronize

                    // to account for spurious wake-ups, we return to condition-variable wait.
                }
            };
            
            auto const bufferData = this->buffer.data();
            block_index worker_vec = worker_vec_initial;
            
            while(true) {
                
                auto work = popWork();
#ifdef IMJ_WITH_ASYNCCONV_STATS
                std::optional<profiling::CpuDuration> threadCPUDuration;
                {
                    profiling::ThreadCPUTimer tt(threadCPUDuration);
#endif // IMJ_WITH_ASYNCCONV_STATS
                    if(unlikely(work==sentinel)) {
                        auto res = worker_2_rt.try_push(sentinel);
                        Assert(res);
                        return;
                    }
                    // Verify that worker_vec and work represent adjacent blocks in the ring:
                    Assert(worker_vec == work-N || (work==0 && worker_vec==(buffer.size()-N)));
                    for(int i=0; i<N; ++i) {
                        bufferData[worker_vec + i] = async_conv->step(bufferData[work + i]);
                    }
                    auto res = worker_2_rt.try_push(worker_vec);
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
    
    bool isZero() const {
        if(!async_conv) {
            return true;
        }
        return async_conv->isZero();
    }
    
    void reset() {
        terminateAsyncJobs();
        Assert(jobs == 0);
        async_conv.reset();
        error_worker_too_slow = false;
        buffer.clear();
        previous_result = signal = 0; // when both values are equal the convolution is not active
    }
        
    double getEpsilon(Algo const & algo) const {
        if(!async_conv) {
            return 0;
        }
        return async_conv->getEpsilon();
    }
    
    void logComputeState(Algo const & algo, std::ostream & os) const {
        int const N = algo.getSubmissionPeriod();
        int const queueSize = algo.getQueueSize();

        os << "Async [_/" << N << "], queueSize : " << queueSize << std::endl;

        IndentingOStreambuf i(os);
        if(async_conv) {
            async_conv->logComputeState(os);
        }
        else {
            os << "nothing" << std::endl;
        }
    }

    int getResultQueueSize() const {
        return worker_2_rt.unsafe_num_elements();
    }
    int getSignalQueueSize() const {
        return rt_2_worker.unsafe_num_elements();
    }
    
    bool hasStepErrors() const {
        return error_worker_too_slow;
    }
    
    int curIndex = 0; // todo redondant avec x.progress
    block_index signal, previous_result;
private:
    using vec = a64::vector<FPT>;
    
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
    
    int jobs = 0;
public:
    std::condition_variable rt_2_worker_cond;
private:
    Queue rt_2_worker = Queue(0);
    Queue worker_2_rt = Queue(0);
    std::unique_ptr<Convolution<typename Async::Algo>> async_conv;
public:
    bool error_worker_too_slow = false;
    
    a64::vector<FPT> buffer; // block_index are indices into this buffer.
private:
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
        rt_2_worker = Queue(0);
        worker_2_rt = Queue(0);
    }
    
public:
    void submit_signal(block_index s) {
        rt_2_worker.push(s);
        
        rt_2_worker_cond.notify_one();
        ++jobs;
    }
    
    bool try_submit_signal(block_index s) {
        bool res = rt_2_worker.try_push(s);
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
        bool success = worker_2_rt.try_pop(res);
        if(likely(success)) {
            --jobs;
        }
        return success;
    }
    
    void flushToSilence(Algo const & algo)
    {
        int queueSize = algo.getQueueSize();
        if((previous_result != signal) && worker) {
            // wait for worker to stop processing
            while(worker_2_rt.unsafe_num_elements() != (queueSize-Algo::queue_room_sz))
            {
                rt_2_worker_cond.notify_one(); // in case this is because of the race condition
                std::this_thread::yield();
            }
        }
        
        std::fill(buffer.begin(),
                  buffer.end(),
                  FPT());
        if(async_conv)
        {
            async_conv->flushToSilence();
        }
        curIndex = 0;
        // it is important to _not_ change signal and previous_result
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

    static constexpr int queue_room_sz = AsyncCPUConvolutionConstants::queue_room_sz;

    using AsyncParams = typename Async::SetupParam;
    using SetupParam = AsyncSetupParam<AsyncParams>;

    void setup(SetupParam const & s) {
        N = s.inputSubmissionPeriod;
        queueSize = s.queueSize;
        asyncParams = s.asyncParams;
    }
    
    bool handlesCoefficients() const {
        return N > 0 && queueSize > 0 && asyncParams && asyncParams->handlesCoefficients();
    }
    bool isValid() const {
        if(asyncParams) {
            if(!asyncParams->isValid()) {
                return false;
            }
            if(asyncParams->handlesCoefficients()) {
                return N > 0 && queueSize > 0;
            }
        }
        return queueSize == 0;
    }

    Latency getLatency() const {
        Assert(handlesCoefficients());
        Assert(asyncParams);
        return
        asyncParams->getImpliedLatency() +
        Latency(
                N*((queueSize-queue_room_sz) // we have some levels of asynchronicity
                   +1) -1 // we need N inputs before we can submit
        );
    }
    
    void step(State & s,
              XAndFFTS<FPT, Tag> const & x_and_ffts,
              Y<FPT, Tag> & y) const
    {
        if(unlikely(s.signal == s.previous_result)) {
            return;
        }
        if constexpr (OnWorkerTooSlow == PolicyOnWorkerTooSlow::PermanentlySwitchToDry) {
            if(unlikely(s.error_worker_too_slow)) {
                y.y[y.uProgress] += x_and_ffts.x[x_and_ffts.progress-1];
                return;
            }
        }
        s.buffer[s.signal+(s.curIndex)] = get_signal(x_and_ffts.x[x_and_ffts.progress-1]);
        ++s.curIndex;
        if(unlikely(s.curIndex == N))
        {
            s.curIndex = 0;
            while(unlikely(!s.try_submit_signal(s.signal))) {
                s.error_worker_too_slow = true;
                if constexpr (OnWorkerTooSlow == PolicyOnWorkerTooSlow::PermanentlySwitchToDry) {
                    y.y[y.uProgress] += x_and_ffts.x[x_and_ffts.progress-1];
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
                    y.y[y.uProgress] += x_and_ffts.x[x_and_ffts.progress-1];
                    return;
                }
                else {
                    s.rt_2_worker_cond.notify_one(); // in case this is because of the race condition
                }
            }
            // Verify that s.signal and s.previous_result represent adjacent blocks in the ring:
            Assert(s.signal == s.previous_result-N || (s.previous_result==0 && s.signal==(s.buffer.size()-N)));

        }
        y.y[y.uProgress] += typename RealSignal::value_type(s.buffer[s.previous_result + s.curIndex]);
    }
    
    void flushToSilence(State & s) const
    {
        s.flushToSilence(*this);
    }
    
    int getQueueSize() const {
        return queueSize;
    }
    int getSubmissionPeriod() const {
        return N;
    }
    
    std::optional<AsyncParams> asyncParams;
private:
    int N = 0;
    int queueSize = 0;
};

template <typename Async, PolicyOnWorkerTooSlow OnWorkerTooSlow>
struct corresponding_legacy_dsp<AlgoAsyncCPUConvolution<Async, OnWorkerTooSlow>> {
    using type = AsyncCPUConvolution<corresponding_legacy_dsp_t<Async>, OnWorkerTooSlow>;
};

}
