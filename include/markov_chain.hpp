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
    
    enum class ExecuteLambdas {
        Yes,
        No
    };
    
    struct MarkovChain {
        void initialize(int n) {
            assert(n < nodes.size());
            current = nodes[n].get();
        }
        
        bool empty() const { return nodes.empty(); }
        
        template<ExecuteLambdas exec>
        MarkovNode * step_normalized(float p) {
            return do_step<Probabilities::NORMALIZE, exec>(p);
        }
        
        template<ExecuteLambdas exec>
        MarkovNode * step(float p) {
            return do_step<Probabilities::UNCHANGED, exec>(p);
        }
        
        template<Probabilities PROBA, ExecuteLambdas exec>
        MarkovNode * do_step(float proba) {
            assert(proba >= 0.f);
            assert(proba <= 1.f);
            assert(current);
            if(!current) {
                return nullptr;
            }
            if(current->destinations.empty()) {
                return current;
            }
            if(PROBA == NORMALIZE) {
                float sum{};
                for(auto const & dest : current->destinations) {
                    sum += dest.first;
                }
                proba *= sum;
            }
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
            if(exec == ExecuteLambdas::Yes) {
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
        
        void clear() {
            nodes.clear();
            current = nullptr;
        }
        
    private:
        std::vector<std::unique_ptr<MarkovNode>> nodes;
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

