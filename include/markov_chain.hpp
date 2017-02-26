
namespace imajuscule {
    
    enum class Move {
        CREATE,
        DELETE,
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
    
    enum Probabilities {
        NORMALIZE,
        UNCHANGED
    };
    
    struct MarkovChain {
        void initialize(int n) {
            assert(n < nodes.size());
            current = nodes[n].get();
        }
        
        template<bool EXEC_TRANSITIONS>
        MarkovNode * step_normalized() {
            return do_step<Probabilities::NORMALIZE, EXEC_TRANSITIONS>();
        }
        
        template<bool EXEC_TRANSITIONS>
        MarkovNode * step() {
            return do_step<Probabilities::UNCHANGED, EXEC_TRANSITIONS>();
        }
        
        template<Probabilities PROBA, bool EXEC_TRANSITIONS>
        MarkovNode * do_step() {
            assert(current);
            if(current->destinations.empty()) {
                return current;
            }
            float sum;
            if(PROBA == NORMALIZE) {
                sum = 0.f;
                for(auto const & dest : current->destinations) {
                    sum += dest.first;
                }
            }
            else {
                sum = 1.f;
            }
            
            auto proba = std::uniform_real_distribution<float>{0.f, sum}(rng::mersenne());
            auto accum = 0.f;
            MarkovNode * to(nullptr);
            for(auto const & dest : current->destinations) {
                accum += dest.first;
                if(accum > proba) {
                    to = dest.second;
                    break;
                }
            }
            if((PROBA == NORMALIZE) && !to) {
                // handle numerical errors
                to = current->destinations.back().second;
            }
            if((PROBA == UNCHANGED) && !to) {
                return current;
            }
            assert(to);
            if(EXEC_TRANSITIONS) {
                current->on_state_change(Move::LEAVE, *current, *to);
                to->on_state_change(Move::ENTER, *to, *current);
            }
            current = to;
            return current;
        }
        
        MarkovNode * getCurrent() const {
            return current;
        }

        template <class... Args>
        MarkovNode * emplace(Args&&... args) {
            nodes.push_back(std::make_unique<MarkovNode>(std::forward<Args>(args)...));
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

} // NS imajuscule

