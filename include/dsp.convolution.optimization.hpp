

namespace imajuscule
{

template<typename T>
struct PartitionAlgo;

template<typename SetupParam>
struct PartitionningSpec
{    
    Optional<SetupParam> optimal_setup;
    float getCost() const { return optimal_setup ? optimal_setup->getCost() : std::numeric_limits<float>::max(); }
    GradientDescent<SetupParam> gd;
    
    void logReport(int n_channels,
                   double theoretical_max_avg_time_per_frame,
                   std::ostream & os)
    {
        using namespace std;
        if(optimal_setup) {
            optimal_setup->logReport(n_channels,
                                     theoretical_max_avg_time_per_frame,
                                     os);
        }
        
        constexpr auto debug_gradient_descent = false;
        if constexpr (debug_gradient_descent) {
            os << "Gradient descent report :" << endl;
            gd.debug(true, os);
        }
    }
};

template<typename SetupParam>
struct PartitionningSpec2
{
    Optional<SetupParam> optimal_setup;
    float getCost() const { return optimal_setup ? optimal_setup->getCost() : std::numeric_limits<float>::max(); }
  
    void logReport(int n_channels,
                   double theoretical_max_avg_time_per_frame,
                   std::ostream & os)
    {
        using namespace std;
        if(optimal_setup) {
            optimal_setup->logReport(n_channels,
                                     theoretical_max_avg_time_per_frame,
                                     os);
        }
    }
};

struct Cost {
    virtual ~Cost() = default;

    Cost() = default;
    
    Cost(int phase)
    : phase(phase)
    {}
    
  void setCost(float c) { cost = c; }

  float getCost() const { return cost; }
    
    void setPhase(int ph) { phase = ph; }
    std::optional<int> getPhase() const { return phase; }

    virtual void logSubReport(std::ostream & os) const = 0;

    void logReport(int n_channels, double theoretical_max_avg_time_per_frame, std::ostream & os)
    {
        using namespace std;
        auto actual = getCost();
        auto theoretical = theoretical_max_avg_time_per_frame / static_cast<float>(n_channels);
        auto ratio = actual / theoretical;
        
        /*
         os << "Dropouts probability    : ";
         static_assert(ratio_soft_limit < ratio_hard_limit);
         if(ratio >= ratio_soft_limit) {
         if(ratio > ratio_hard_limit) {
         os << "100 %";
         }
         else {
         os << " 50 %";
         }
         }
         else {
         os << " 0 %";
         }
         os << endl;
         */
        
        auto nActual = countDecimalNumbersBeforeTheDot(actual);
        auto nTheoretical = countDecimalNumbersBeforeTheDot(theoretical);
        std::string prefixActual =
        std::string(std::max(0, nTheoretical-nActual), ' ');
        std::string prefixTheoretical =
        std::string(std::max(0, nActual-nTheoretical), ' ');
        auto const precision = os.precision();
        os << "Foreseen CPU load (1 core)          : " << std::fixed << std::setprecision(2) << 100.*ratio << "%" << endl;
        os << "Average computation time per sample : " << std::fixed << std::setprecision(0) << actual << " ns" << endl;
        //os << "- Actual                : " << prefixActual      << actual      << " ns" << endl;
        //os << "- Allowed (theoretical) : " << prefixTheoretical << theoretical << " ns" << endl;
        os.precision(precision);
        if(phase) {
            os << "Optimized computation phase : " << *phase << endl;
        }

        {
            os << "Composition:" << std::endl;
            IndentingOStreambuf i(os);
            logSubReport(os);
        }
        
        //os << endl;
        //os << *partitionning.cost << endl;
    }
private:
    // worst computation time over one callback, averaged per sample.
  float cost = std::numeric_limits<float>::max();
    // phase (of the late handler), in frames (in case 2 or more convolutions run at the same time, we can dephase them)
    std::optional<int> phase;
};

}
