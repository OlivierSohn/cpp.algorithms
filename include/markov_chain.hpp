
namespace imajuscule {
    
    enum Move {
        ENTER,
        LEAVE
    };
    
    struct MarkovNode {
        using f = std::function<void(Move, MarkovNode&me, MarkovNode&from_to)>;

        MarkovNode(f func) : on_state_change(std::move(func)) {}

        f on_state_change;
        std::vector<std::pair<float,MarkovNode*>> destinations;
    };
    
    /*
     * a markov chain where each node has a probability to transition towards a neighbour,
     * and a lambda to let the user program stuff during transitions
     */
    
    struct MarkovChain {
        void initialize(int n) {
            current = nodes[n].get();
        }
        
        MarkovNode * step() {
            assert(current);
            if(current->destinations.empty()) {
                return current;
            }
            auto proba = std::uniform_real_distribution<float>{0.f, 1.f}(rng::mersenne());
            auto accum = 0.f;
            for(auto const & dest : current->destinations) {
                accum += dest.first;
                if(accum < proba) {
                    continue;
                }
                current->on_state_change(LEAVE, *current, *dest.second);
                dest.second->on_state_change(ENTER, *dest.second, *current);
                current = dest.second;
                break;
            }
            return current;
        }
        
        MarkovNode * getCurrent() const {
            return current;
        }

        // why doesn't this compile??
        /*template <class... Args>
        MarkovNode * emplace_(Args&&... args) {
            nodes.push_back(std::make_unique<MarkovNode>(std::forward<Args>(args)...));
            return nodes.back().get();
        }*/
        
        template <class... Args>
        MarkovNode * emplace(Args&&... args) {
            nodes.emplace_back(std::forward<Args>(args)...);
            return nodes.back().get();
        }
        
        std::vector<std::unique_ptr<MarkovNode>> nodes;
    private:
        MarkovNode * current = nullptr;
    };
    
    
    static inline void def_markov_transition(MarkovNode&n1, MarkovNode&n2, float p) {
        n1.destinations.emplace_back(p, &n2);
    }
    static inline void def_markov_transition(MarkovNode*n1, MarkovNode*n2, float p) {
        assert(n1 && n2);
        def_markov_transition(*n1, *n2, p);
    }
    
    // why doesn't this compile??
    /*
        template <class... Args>
        static inline MarkovNode * emplace(MarkovChain & mc, Args&&... args) {
            using namespace std;
            mc.nodes.push_back(make_unique<MarkovNode>(std::forward<Args>(args)...));
            return mc.nodes.back().get();
        }
    */

} // NS imajuscule

