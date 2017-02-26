
using namespace imajuscule;


TEST(MarkovChain_, one_node) {
    
    MarkovChain mc;
    
    int n = 0;
    int balance = 0;
    
    auto unique_node = mc.emplace([&](Move const m, MarkovNode&me, MarkovNode&from_to){
        ++n;
        if(m==Move::ENTER) {
            balance++;
        }
        else {
            balance--;
        }
        EXPECT_EQ(&me, &from_to);
    });

    def_markov_transition(unique_node, unique_node, 0.5f);
    
    mc.initialize(0);
    
    for(int i=0; i<100; i++) {
        EXPECT_EQ(unique_node, mc.step<true>());
        EXPECT_EQ(unique_node, mc.getCurrent());
    }
    
    ASSERT_NE(0, n);
    ASSERT_EQ(0, balance);
}

void f_markov(MarkovChain m) {
    std::cout << &m << std::endl;
}

TEST(MarkovChain_, two_nodes) {
    
    MarkovChain mc;
    
    auto n = 0;
    auto node1 = mc.emplace([&](Move const m, MarkovNode&me, MarkovNode&from_to){
        ++n;
    });
    auto node2 = mc.emplace([&](Move const m, MarkovNode&me, MarkovNode&from_to){
    });

    constexpr auto transition_p = 0.1f;

    def_markov_transition(node1, node2, transition_p);
    def_markov_transition(node2, node1, transition_p);
    
    mc.initialize(0);
    
    constexpr auto nSteps = 100000;
    auto count_1 = 0;
    for(int i=0; i<nSteps; i++) {
        mc.step<true>();
        if(mc.getCurrent() == node1) {
            ++count_1;
        }
    }
    
    auto expected_n_trans = transition_p*nSteps;
    EXPECT_LT(.75f * expected_n_trans, n);
    EXPECT_GT(1.25f * expected_n_trans, n);

    auto expected_count_1 = nSteps/2;
    EXPECT_LT(.75f * expected_count_1, count_1);
    EXPECT_GT(1.25f * expected_count_1, count_1);
    
    f_markov(std::move(mc));
}


