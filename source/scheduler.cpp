

using namespace imajuscule;

void Scheduler::update(float dt) {
    for(auto it = functions.begin(); it != functions.end(); ) {
        auto & f = *it;
        f.remaining -= dt;
        if( f.remaining <= 0.f ) {
            auto res = f.f();
            if(f.when == Periodicity::ONCE) {
                it = functions.erase(it);
                continue;
            } else {
                f.remaining = (res > 0.f)? res : f.interval;
            }
        }
        ++it;
    }
}

Scheduler::Function::Function( Periodicity p, float interval, func && f) :
interval(interval)
, remaining((p == Periodicity::PERIODIC) ? 0.f : interval)
, when(p)
, f(std::move(f)) {}
