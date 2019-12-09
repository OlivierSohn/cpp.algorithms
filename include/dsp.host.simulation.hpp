
namespace imajuscule
{

struct Periodically {
    Periodically(std::chrono::steady_clock::duration time_step,
                 int n_max_periods = 1000000)
    : time_step(time_step)
    {
        time_calls_start.reserve(n_max_periods);
        time_calls_end.reserve(n_max_periods);
        dur_sleep.reserve(n_max_periods);
    }
    
    template<typename F>
    bool exec(F f) {
        t_zero = std::chrono::steady_clock::now();
        while(true) {
            if(!nextPeriod()) {
                // the previous call to f took too long to complete
                return false;
            }
            if(!f()) {
                // end of the simulation
                return true;
            }
        }
    }
    
    int i=1; // start at one, so that the first iteration is not seen as a missed deadline
    std::chrono::steady_clock::time_point t_zero;
    std::chrono::steady_clock::duration time_step;
    std::vector<std::chrono::steady_clock::time_point> time_calls_start, time_calls_end;
    std::vector<std::pair<std::chrono::steady_clock::duration, std::chrono::steady_clock::duration>> dur_sleep;
private:
    bool nextPeriod() {
        auto now = std::chrono::steady_clock::now();
        time_calls_end.push_back(now);
        if(likely(!time_calls_start.empty()))
        {
            auto dtLastCall = time_calls_end.back()-time_calls_start.back();
            if(dtLastCall > time_step) {
                // either f had too much work to do, or the thread has been preempted
                
                // we will uncomment this code when the thread will be configured for realtime audio
                /*
                std::cout << "all costs:" << std::endl;
                for(int i=0; i<time_calls_start.size(); ++i) {
                    auto dtCall = time_calls_end[i+1] - time_calls_start[i];
                    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(dtCall).count() << std::endl;
                }
                std::cout << "all sleeps:" << std::endl;
                for(auto const & [need, slept] : dur_sleep) {
                    std::cout <<
                    std::chrono::duration_cast<std::chrono::microseconds>(need).count() << " " <<
                    std::chrono::duration_cast<std::chrono::microseconds>(slept).count() << " " <<
                    std::chrono::duration_cast<std::chrono::microseconds>(slept-need).count() << std::endl;
                }

                std::cout << "The last f call took too long at i=" << i << std::endl;
                std::cout << "time step:" << std::chrono::duration_cast<std::chrono::microseconds>(time_step).count() << std::endl;
                std::cout << "cost f   :" << std::chrono::duration_cast<std::chrono::microseconds>(dtLastCall).count() << std::endl;
                return false;//*/
            }
        }
        std::chrono::steady_clock::duration needSleep{};
        auto next_time = time_step * i + t_zero;
        ++i;
        if(now < next_time) {
            needSleep = next_time - now;
            std::this_thread::sleep_for(needSleep);
        }
        else {
            // the last sleep was too long, we can igonre it.
        }
        auto wakeup = std::chrono::steady_clock::now();
        auto slept = wakeup-now;
        dur_sleep.emplace_back(needSleep, slept);
        time_calls_start.push_back(wakeup);
        return true;
    };
};

struct AudioHostSimulator {
    double const frame_rate;
    int const n_audio_frames_per_cb;
private:
    std::vector<a64::vector<float>> inputs, outputs;
    int curFrame = 0;
    int maxSz;
public:
    
    AudioHostSimulator(double const frame_rate,
                       int const n_audio_frames_per_cb,
                       int const n_input_channels,
                       int const n_output_channels)
    : frame_rate(frame_rate)
    , n_audio_frames_per_cb(n_audio_frames_per_cb)
    , maxSz(n_audio_frames_per_cb)
    {
        inputs.resize(n_input_channels);
        outputs.resize(n_output_channels);
        
        for(auto & i : inputs) {
            i.resize(n_audio_frames_per_cb);
        }
        for(auto & o : outputs) {
            o.resize(n_audio_frames_per_cb);
        }
    }
    
    AudioHostSimulator(double const frame_rate,
                       int const n_audio_frames_per_cb,
                       std::vector<a64::vector<float>> inputs,
                       std::vector<a64::vector<float>> outputs)
    : frame_rate(frame_rate)
    , n_audio_frames_per_cb(n_audio_frames_per_cb)
    , inputs(std::move(inputs))
    , outputs(std::move(outputs))
    {
        maxSz = 0;
        for(auto const & i : this->inputs) {
            if(i.size() < n_audio_frames_per_cb) {
                throw std::logic_error("input too short");
            }
            maxSz = std::max(maxSz, static_cast<int>(i.size()));
        }
        for(auto const & o : this->outputs) {
            if(o.size() < n_audio_frames_per_cb) {
                throw std::logic_error("output too short");
            }
            maxSz = std::max(maxSz, static_cast<int>(o.size()));
        }
    }

    template<typename F>
    std::optional<Periodically> simulate(F f) {
        
        auto f2 = [this, f](){
            if(curFrame + n_audio_frames_per_cb > maxSz) {
                curFrame = 0;
            }
            auto res = f(inputs, outputs, curFrame, n_audio_frames_per_cb);
            curFrame += n_audio_frames_per_cb;
            return res;
        };

        Periodically p(std::chrono::nanoseconds(static_cast<int>(1e9 * n_audio_frames_per_cb / frame_rate)));
        if(!p.exec(f2)) {
            // one call took too long to execute
            return {};
        }
        return p;
    }
};

} // namespace imajuscule
