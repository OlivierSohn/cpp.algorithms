/*
 
 Find the minimal distance between two strings, for example:
 
 - hello
 - yellow
 
 should lead to : one substitution (h/y), one addition (w)
 
 */

namespace imajuscule {
    enum class EditDistanceAction {
        Insert,
        Delete,
        Substitute,
        Match
    };
    std::string toString(EditDistanceAction e) {
        switch(e) {
            case EditDistanceAction::Insert:
                return "insert";
            case EditDistanceAction::Delete:
                return "Delete";
            case EditDistanceAction::Substitute:
                return "Substitute";
            case EditDistanceAction::Match:
                return "Match";
        }
    }
    
    struct EditDistanceStep {
        EditDistanceStep(EditDistanceAction a, std::size_t i1, std::size_t i2) : index_1(i1), index_2(i2), action(a) {}
        EditDistanceStep(EditDistanceStep const & start, EditDistanceAction move_) : action(move_) {
            std::tie(index_1, index_2) = start.getNextIndices();
        }

        std::size_t index_1, index_2;
        EditDistanceAction action;

        bool operator == (EditDistanceStep const & o) const {
            return index_1 == o.index_1 && index_2 == o.index_2 && action == o.action;
        }
        
        std::pair<std::size_t, std::size_t> getNextIndices() const {
            switch(action) {
                case EditDistanceAction::Insert:
                    return {index_1  , index_2+1};
                    
                case EditDistanceAction::Delete:
                    return {index_1+1, index_2  };
                    
                case EditDistanceAction::Match:
                case EditDistanceAction::Substitute:
                    return {index_1+1, index_2+1};
            }
            assert(0);
        }
        
        EditDistanceStep succ(EditDistanceAction m) const {
            return {*this, m};
        }
    };
    
    struct EditDistanceResult {
        EditDistanceResult(std::vector<EditDistanceAction> v) : EditDistanceResult() {
            res.reserve(v.size());
            for(auto e:v) { add(e); }
        }
        
        EditDistanceResult() = default;
        
        void Substitute() { add(EditDistanceAction::Substitute); }
        void Insert() { add(EditDistanceAction::Insert); }
        void Delete() { add(EditDistanceAction::Delete); }
        void Match() { add(EditDistanceAction::Match); }
        
        void add(EditDistanceAction a) {
            std::size_t i1,i2;
            std::tie(i1,i2) = getIndices();
            res.emplace_back(a, i1, i2);
        }

        std::pair<std::size_t, std::size_t> getIndices() const {
            if(res.empty()) {
                return {0,0};
            }
            return res.back().getNextIndices();
        }

        auto const & result() const { return res; }
        auto size() const { return res.size(); }

        double cost(std::function<double(std::optional<EditDistanceAction> const, EditDistanceAction const)> const & weights,
                   std::size_t start_i) const {
            double c = 0.;
            std::optional<EditDistanceAction> prevAction;
            for(auto i=start_i; i<res.size(); ++i) {
                c += weights(prevAction, res[i].action);
                prevAction = res[i].action;
            }
            return c;
        }
        
        std::optional<EditDistanceAction> getMaybeLastMove() const {
            if(res.empty()) {
                return {};
            }
            return res.back().action;
        }
        
        EditDistanceResult dropLeft(std::size_t n) const {
            EditDistanceResult ret;
            ret.res.reserve(res.size()-n);
            for(int i=n; i<res.size(); ++i) {
                ret.res.push_back(res[i]);
            }
            return ret;
        }
        
        void withSuffix(EditDistanceResult const & suffix) {
            res.reserve(res.size() + suffix.res.size());
            res.insert(res.end(), suffix.res.begin(), suffix.res.end());
        }
        
        
        void diagnose(std::function<double(std::optional<EditDistanceAction> const, EditDistanceAction const)> const & weights) {
            std::cout << std::endl;
            if(res.empty()) {
                std::cout << "empty" << std::endl;
            }
            for(auto & e:res) {
                std::cout << toString(e.action) << "\t" << e.index_1 << "\t" << e.index_2 << std::endl;
            }
            std::cout << "cost:" << cost(weights, 0) << std::endl;
            std::cout << std::endl;
        }
        

    private:
        std::vector<EditDistanceStep> res;
    };
    
    namespace ed_naiveRecursive {
        /*
         The naive implementation using recursion
         */

        EditDistanceResult
        editDistance(std::string const & s1,
                     std::string const & s2,
                     std::function<double(std::optional<EditDistanceAction> const, EditDistanceAction const)> const & weights,
                     EditDistanceResult res = {}) {
            std::size_t i1,i2;
            std::tie(i1,i2) = res.getIndices();
            
            if(i1 == s1.size() || i2 == s2.size()) {
                // at least one of the two strings have been fully processed:
                // we have either only deletions or only insertions left.
                
                for(; i2 != s2.size(); ++i2) {
                    res.Insert();
                }
                for(; i1 != s1.size(); ++i1) {
                    res.Delete();
                }
            }
            else {
                if(s1[i1] == s2[i2]) {
                    // it's a match!
                    // TODO this assumes that weight(match) <= weight(any other)
                    res.Match();
                    res = editDistance(s1, s2, weights, res);
                }
                else {
                    // The minimal-cost solution could either be a substitution, an insert or a delete.
                    // So we try those three possibilities:
                    std::optional<double> bestCost;
                    auto prefixSize = res.result().size();
                    EditDistanceResult bestRes;
                    for(int i=0; i<3; ++i) {
                        auto resOption = res;
                        resOption.add(static_cast<EditDistanceAction>(i));
                        resOption = editDistance(s1, s2, weights, resOption);
                        auto costOption = resOption.cost(weights,
                                                         prefixSize  // the prefix is the same for all options so we can ignore it in the cost
                                                         );
                        if(!bestCost || *bestCost > costOption) {
                            bestCost = costOption;
                            bestRes = resOption;
                        }
                    }
                    return bestRes;
                }
            }
            
            return res;
        }
    } // NS ed_naiveRecursive
    
    struct MemoKey {
        MemoKey(EditDistanceStep const & step) {
            std::tie(i1,i2) = step.getNextIndices();
            prevAction = step.action;
        }
        
        MemoKey():
        i1(0), i2(0)
        {}

        std::optional<EditDistanceAction> prevAction;
        std::size_t i1;
        std::size_t i2;
        
        bool operator == (const MemoKey & o) const {
            return prevAction==o.prevAction && i1==o.i1 && i2==o.i2;
        }
    };

    template <typename...> struct hash;
    
    template<typename T>
    struct hash<T>
    : public std::hash<T>
    {
        using std::hash<T>::hash;
    };
    
    
    template <typename T, typename... Rest>
    struct hash<T, Rest...>
    {
        inline std::size_t operator()(const T& v, const Rest&... rest) {
            std::size_t seed = hash<Rest...>{}(rest...);
            seed ^= hash<T>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };
}

namespace std
{
    template<>
    struct hash< ::imajuscule::MemoKey >
    {
        typedef ::imajuscule::MemoKey argument_type;
        std::size_t operator()( const argument_type& arg ) const
        {
            ::imajuscule::hash<
                std::optional<::imajuscule::EditDistanceAction>,
                std::size_t,
                std::size_t
            > hasher;
            return hasher(arg.prevAction, arg.i1, arg.i2);
        }
    };
}

namespace imajuscule {
    namespace ed_memoizedRecursion {
        /*
         Using memoization to prune the recursion tree
         */
        /*
        using MemoMap = std::unordered_map<
        MemoKey,
        EditDistanceResult // suffix only
        >;
        struct Memo {
            MemoMap m;
            mutable std::size_t n_errorLookups = 0;
            mutable std::size_t n_successLookups = 0;
            
            void memoize(MemoKey const & k, EditDistanceResult const & fullRes, std::size_t ignoredPrefix) {
                assert(m.count(k) == 0);
                m[k] = fullRes.dropLeft(ignoredPrefix);
            }
            auto find(MemoKey const & k) const {
                auto it = m.find(k);
                if(it == m.end()) {
                    ++n_errorLookups;
                }
                else {
                    ++n_successLookups;
                }
                return it;
            }
            
            auto end() const { return m.end(); }
            
            void diagnose() const {
                std::cout << "memo size = " << m.size() <<
                " lookups success = " << n_successLookups <<
                " error = " << n_errorLookups << std::endl;
            }
        };
        EditDistanceResult editDistance(std::string const & s1,
                                        std::string const & s2,
                                        std::function<double(std::optional<EditDistanceAction> const, EditDistanceAction const)> const & weights,
                                        Memo & memo,
                                        EditDistanceResult res = {}) {
            std::size_t i1,i2;
            std::tie(i1,i2) = res.getIndices();
            auto const prefixSize = res.result().size();
            auto key = std::make_tuple(res.getMaybeLastMove(), i1, i2);
            auto it = memo.find(key);
            if(it != memo.end()) {
                // we already have a result for these parameters, so instead of computing it again
                // we use the result:
                res.withSuffix(it->second);
                return res;
            }
            
            if(i1 == s1.size() || i2 == s2.size()) {
                // at least one of the two strings have been fully processed:
                // we have either only deletions or only insertions left.
                
                for(; i2 != s2.size(); ++i2) {
                    res.Insert();
                }
                for(; i1 != s1.size(); ++i1) {
                    res.Delete();
                }
            }
            else {
                if(s1[i1] == s2[i2]) {
                    // it's a match!
                    // TODO this assumes that weight(match) <= weight(any other)
                    res.Match();
                    res = editDistance(s1, s2, weights, memo, res);
                }
                else {
                    // The minimal-cost solution could either be a substitution, an insert or a delete.
                    // So we try those three possibilities:
                    std::optional<double> bestCost;
                    EditDistanceResult bestRes;
                    for(int i=0; i<3; ++i) {
                        auto resOption = res;
                        resOption.add(static_cast<EditDistanceAction>(i));
                        resOption = editDistance(s1, s2, weights, memo, resOption);
                        auto costOption = resOption.cost(weights,
                                                         prefixSize  // the prefix is the same for all options so we can ignore it in the cost
                                                         );
                        if(!bestCost || *bestCost > costOption) {
                            bestCost = costOption;
                            bestRes = resOption;
                        }
                    }
                    memo.memoize(key, bestRes, prefixSize);
                    return bestRes;
                }
            }
            
            memo.memoize(key, res, prefixSize);
            return res;
        }
         */
    } // NS ed_memoizedRecursion
    
    namespace ed_optimizedMemoizedRecursion {
        /*
         Optimizing the memoization data structure to use less memory
         */
        
        struct PartialPathSummary {
            std::size_t size = 0;
            double cost = {};
            std::optional<EditDistanceAction> firstMove;
        };
    }
}

namespace std
{
    template<>
    struct hash< ::imajuscule::ed_optimizedMemoizedRecursion::PartialPathSummary >
    {
        typedef ::imajuscule::ed_optimizedMemoizedRecursion::PartialPathSummary argument_type;
        std::size_t operator()( const argument_type& arg ) const
        {
            ::imajuscule::hash<
            std::size_t,
            double,
            std::optional<::imajuscule::EditDistanceAction>
            > hasher;
            return hasher(arg.size, arg.cost, arg.firstMove);
        }
    };
}

namespace imajuscule {
    namespace ed_optimizedMemoizedRecursion {

        PartialPathSummary concat(PartialPathSummary const & p1, PartialPathSummary const & p2) {
            PartialPathSummary res;
            res.size = p1.size + p2.size;
            res.cost = p1.cost + p2.cost;
            res.firstMove = p1.firstMove;
            return res;
        }
        
        using MemoMap = std::unordered_map<
        MemoKey,
        PartialPathSummary
        >;
        struct Memo {
            MemoMap m;
            mutable std::size_t n_errorLookups = 0;
            mutable std::size_t n_successLookups = 0;
            
            void memoize(MemoKey const & k, PartialPathSummary const & v) {
                assert(m.count(k) == 0);
                m[k] = v;
            }
            auto find(MemoKey const & k) const {
                auto it = m.find(k);
                if(it == m.end()) {
                    ++n_errorLookups;
                }
                else {
                    ++n_successLookups;
                }
                return it;
            }
            
            auto end() const { return m.end(); }
            
            void diagnose() const {
                std::cout << "memo size = " << m.size() <<
                " lookups success = " << n_successLookups <<
                " error = " << n_errorLookups << std::endl;
            }
        };
        struct Computer {
            Computer(std::string const & s1,
                     std::string const & s2,
                     Memo & memo,
                     std::function<double(std::optional<EditDistanceAction> const, EditDistanceAction const)> const & weights) :
            s1(s1), s2(s2), memo(memo), weights(weights)
            {}
            
            PartialPathSummary operator ()(MemoKey const k = {}) {
                auto it = memo.find(k);
                if(it != memo.end()) {
                    return it->second;
                }
                auto memoize = [this, &k](PartialPathSummary const & res){
                    memo.memoize(k, res);
                    return res;
                };
                
                std::optional<EditDistanceAction> action;
                bool active1 = k.i1 != s1.size();
                bool active2 = k.i2 != s2.size();
                if(!active1 || !active2) { // at least one of the two strings have been fully processed
                    if (active1) {
                        action = EditDistanceAction::Delete;
                    }
                    else if(active2) {
                        action = EditDistanceAction::Insert;
                    }
                }
                else if(s1[k.i1] == s2[k.i2]) {
                    action = EditDistanceAction::Match;
                }
                else {
                    std::optional<PartialPathSummary> best;
                    for(int i=0; i<3; ++i) {
                        auto option = recurse(k, static_cast<EditDistanceAction>(i));
                        if(!best || best->cost > option.cost) {
                            best = option;
                        }
                    }
                    return memoize(*best);
                }
                
                return memoize(action ? recurse(k, *action) : PartialPathSummary());
            }
            
        private:
            std::string const & s1, s2;
            Memo & memo;
            std::function<double(std::optional<EditDistanceAction> const, EditDistanceAction const)> const & weights;
            
            PartialPathSummary recurse(MemoKey const & k,
                                       EditDistanceAction const move) {
                PartialPathSummary prefix;
                
                prefix.size = 1;
                prefix.firstMove = move;
                prefix.cost = weights(k.prevAction, move);
                
                return concat(prefix,
                              (*this)(MemoKey(EditDistanceStep(move, k.i1, k.i2))));
            }
        };
        
        PartialPathSummary computeShortestDistance(std::string const & s1,
                                                   std::string const & s2,
                                                   std::function<double(std::optional<EditDistanceAction> const, EditDistanceAction const)> const & weights,
                                                   Memo & memo) {
            return Computer(s1, s2, memo, weights)();
        }


        bool almostEq(double d1, double d2, double eps = 1e-8) {
            if(d1 == 0. || d2 == 0.) {
                return std::abs(d1) < eps && std::abs(d2) < eps;
            }
            return (std::abs(d1-d2) / std::max(std::abs(d1), std::abs(d2))) < eps;
        }
        EditDistanceResult extractResultFromSummary(std::optional<EditDistanceAction> firstMove,
                                                    Memo const & memo) {
            EditDistanceResult res;
            
            auto move = firstMove;
            while(move) {
                res.add(*move);
                auto it = memo.find({res.result().back()});
                assert(it != memo.end());
                move = it->second.firstMove;
            }
            return res;
        }
        EditDistanceResult editDistance(std::string const & s1,
                                        std::string const & s2,
                                        std::function<double(std::optional<EditDistanceAction> const, EditDistanceAction const)> const & weights,
                                        Memo & memo) {
            auto s = computeShortestDistance(s1,s2, weights, memo);
            auto res = extractResultFromSummary(s.firstMove, memo);
            assert(res.size() == s.size);
            return res;
        }

    } // NS ed_optimizedMemoizedRecursion
    
}

namespace editDistanceTest {
    void test(std::string s1,
              std::string s2,
              std::vector<imajuscule::EditDistanceAction> as) {
        using namespace imajuscule;
        
        //using namespace imajuscule::ed_naiveRecursive;
        //using namespace imajuscule::ed_memoizedRecursion;
        using namespace imajuscule::ed_optimizedMemoizedRecursion;

        auto weights = [](std::optional<EditDistanceAction> const prev, EditDistanceAction const e) {
            switch(e) {
                case EditDistanceAction::Insert:
                case EditDistanceAction::Delete:
                    if(prev && *prev == EditDistanceAction::Substitute) {
                        // to avoid identical costs for
                        // substitute -> insert and insert -> substitute
                        // or
                        // substitute -> delete and delete -> substitute
                        // we favor the case where substitution occurs first.
                        return 0.9;
                    }
                    return 1.;
                case EditDistanceAction::Substitute:
                    return 1.;
                case EditDistanceAction::Match:
                    return 0.;
            }
        };

        Memo memo;
        auto res = editDistance(s1, s2, weights, memo);
        EditDistanceResult er(std::move(as));
        
        std::cout << std::endl;
        std::cout << s1 << std::endl;
        std::cout << s2 << std::endl;
        res.diagnose(weights);
        er.diagnose(weights);
        
        memo.diagnose();
        
        ASSERT_EQ(res.result(), er.result());
    }
}
TEST(EditDistance, ed) {
    using namespace editDistanceTest;
    using namespace imajuscule;
    
    test("","",{});
    test(" ","",{EditDistanceAction::Delete});
    test(""," ",{EditDistanceAction::Insert});
    test(" "," ",{EditDistanceAction::Match});
    test("a","b",{EditDistanceAction::Substitute});

    test("ax","bx",{
        EditDistanceAction::Substitute,
        EditDistanceAction::Match
    });

    test("ayx","bx",{
        EditDistanceAction::Substitute,
        EditDistanceAction::Delete,
        EditDistanceAction::Match
    });

    test("ax","byx",{
        EditDistanceAction::Substitute,
        EditDistanceAction::Insert,
        EditDistanceAction::Match
    });

    test("hello you",
         "yellow you",
         {
             EditDistanceAction::Substitute, // h <> y
             EditDistanceAction::Match, // e
             EditDistanceAction::Match, // l
             EditDistanceAction::Match, // l
             EditDistanceAction::Match, // o
             EditDistanceAction::Insert, // + w
             EditDistanceAction::Match, // ' '
             EditDistanceAction::Match, // y
             EditDistanceAction::Match, // o
             EditDistanceAction::Match // u
         });
}
