

namespace imajuscule
{

template<typename Param, typename FPT, typename FFTTag>
struct PartitionAlgo;

struct Cost {
    virtual ~Cost() = default;

    Cost() = default;
    
    Cost(int phase)
    : phase(phase)
    {}
    
  void setCost(float c) { cost = c; }

  float getCost() const { return cost; }
    
    void setPhase(float ph) { phase = ph; }
    std::optional<float> getPhase() const { return phase; }

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

        os << "Algorithm:" << std::endl;
        {
            IndentingOStreambuf i(os);
            logSubReport(os);
        }
        
        //os << endl;
        //os << *partitionning.cost << endl;
    }
private:
    // worst computation time over one callback, averaged per sample.
  float cost = std::numeric_limits<float>::max();
    // phase in frames (in case 2 or more convolutions run at the same time, we can dephase them)
    std::optional<float> phase;
};

}
