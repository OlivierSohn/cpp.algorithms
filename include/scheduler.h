
namespace imajuscule {

    enum class Periodicity {
        PERIODIC,
        ONCE
    };

    struct Scheduler {
        using func = std::function<float(void)>;
        
        struct Function {
            // when function is periodic:
            //   if the lambda returns 0, then the initial interval is used for the next period
            //   else the returned value is used
            Function( Periodicity, float interval, func &&);
            float interval;
            float remaining;
            Periodicity when;
            func f;
        };

        void update(float);
        void schedule( Function && f ) {
            functions.emplace_back( std::move(f) );
        }
        
    private:
        // we use a list to allow adding a scheduled function from within a scheduled function
        std::list< Function > functions;
    };
}
