TEST(StaticVectorSingleThreadSCMP, singleThread_0_elt) {
  using namespace imajuscule;
  
  singlethread::static_vector<int> a(0);
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n);
  }
  
  EXPECT_FALSE(a.tryInsert(1));
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n);
  }
}

TEST(StaticVectorSingleThreadSCMP, singleThread_1_elt) {
  using namespace imajuscule;
  
  singlethread::static_vector<int> a(1);
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n);
  }
  
  EXPECT_TRUE(a.tryInsert(1)); // add an element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 1) {
        elt = 8;
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 8) {
        ++n;
        return false; // delete the element
      }
      return true;
    });
    EXPECT_EQ(1,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n); // because the element was deleted
  }
  
  EXPECT_TRUE(a.tryInsert(4)); // add an element
  EXPECT_FALSE(a.tryInsert(5)); // add an element
  EXPECT_FALSE(a.tryInsert(6)); // add an element
  EXPECT_FALSE(a.tryInsert(7)); // add an element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 4) {
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(1,n);
  }
  
}

TEST(StaticVectorSingleThreadSCMP, singleThread_2_elts) {
  using namespace imajuscule;
  
  singlethread::static_vector<int> a(2);
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return true;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(0,n);
  }
  
  EXPECT_TRUE(a.tryInsert(1)); // add an element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 1) {
        elt = 8;
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  EXPECT_TRUE(a.tryInsert(2)); // add another element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 2) {
        elt = 9;
        ++n;
        return true;
      }
      if(elt == 8) {
        elt = 10;
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(2,n);
  }
  
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 10) {
        ++n;
        return false; // delete the first element
      }
      return true;
    });
    EXPECT_EQ(1,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 9) {
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(1,n);
  }
  
  EXPECT_TRUE(a.tryInsert(3)); // add another element
  
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      if(elt == 9) {
        ++n;
        return true;
      }
      if(elt == 3) {
        ++n;
        return true;
      }
      return false;
    });
    EXPECT_EQ(0,nRemoved);
    EXPECT_EQ(2,n);
  }
  {
    int n = 0;
    int nRemoved = a.forEach([&n](int & elt) {
      ++n;
      return false;
    });
    EXPECT_EQ(2,nRemoved);
    EXPECT_EQ(2,n);
  }
}

TEST(StaticVectorSingleThreadSCMP, OnRemovalAssignFromDefault_UP) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    singlethread::static_vector<std::unique_ptr<Destructible>, singlethread::OnRemoval::AssignFromDefault> a(10);
    
    a.tryInsert(std::make_unique<Destructible>());
    EXPECT_EQ(1,nLiveObjects());
    a.tryInsert(std::make_unique<Destructible>());
    EXPECT_EQ(2,nLiveObjects());
    
    a.forEach([](auto &){ return false; });
    EXPECT_EQ(0,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}

TEST(StaticVectorSingleThreadSCMP, OnRemovalAssignFromDefault) {
  using namespace imajuscule;

  EXPECT_EQ(0,nLiveObjects());
  {
    singlethread::static_vector<Destructible, singlethread::OnRemoval::AssignFromDefault> a(10);
    
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    
    // remove all elements
    a.forEach([](auto){ return false; });
    EXPECT_EQ(10,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}


TEST(StaticVectorSingleThreadSCMP, OnRemovalDoNothing_UP) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    singlethread::static_vector<std::unique_ptr<Destructible>, singlethread::OnRemoval::DoNothing> a(10);
    
    a.tryInsert(std::make_unique<Destructible>());
    EXPECT_EQ(1,nLiveObjects());
    a.tryInsert(std::make_unique<Destructible>());
    EXPECT_EQ(2,nLiveObjects());
    
    a.forEach([](auto &){ return false; });
    EXPECT_EQ(2,nLiveObjects()); // the removal hasn't changed the number of live objects
    
    // the next 2 inserts will succeed and not change the number of live objects
    EXPECT_TRUE(a.tryInsert(std::make_unique<Destructible>()));
    EXPECT_EQ(2,nLiveObjects());
    EXPECT_TRUE(a.tryInsert(std::make_unique<Destructible>()));
    EXPECT_EQ(2,nLiveObjects());

    // subsequent inserts will change the number of live objects
    EXPECT_TRUE(a.tryInsert(std::make_unique<Destructible>()));
    EXPECT_EQ(3,nLiveObjects());
    EXPECT_TRUE(a.tryInsert(std::make_unique<Destructible>()));
    EXPECT_EQ(4,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}

TEST(StaticVectorSingleThreadSCMP, OnRemovalDoNothing) {
  using namespace imajuscule;
  
  EXPECT_EQ(0,nLiveObjects());
  {
    singlethread::static_vector<Destructible, singlethread::OnRemoval::DoNothing> a(10);
    
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    
    // remove all elements
    a.forEach([](auto){ return false; });
    EXPECT_EQ(10,nLiveObjects());

    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
    a.tryInsert({});
    EXPECT_EQ(10,nLiveObjects());
  }
  EXPECT_EQ(0,nLiveObjects());
}


