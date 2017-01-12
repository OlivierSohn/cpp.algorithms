
// this doesn't compile in markov_chain.hpp header
namespace imajuscule {
    template <class... Args>
    static inline MarkovNode * emplace(MarkovChain & mc, Args&&... args) {
        using namespace std;
        mc.nodes.push_back(make_unique<MarkovNode>(forward<Args>(args)...));
        return mc.nodes.back().get();
    }
} // NS imajuscule

