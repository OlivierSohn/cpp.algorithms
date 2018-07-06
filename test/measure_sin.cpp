
template<typename T, typename F>
void testSin(F sin_func) {
  using namespace imajuscule;
  using namespace imajuscule::profiling;
  using namespace std;
  using namespace std::chrono;
  high_resolution_clock::rep duration = 0;
  T sum {};
  {
    Timer<high_resolution_clock> t(duration);
    for(int i=0; i<100000000; ++i) {
      sum += sin_func(static_cast<T>(i));
    }
  }
  cout << duration << endl;
  cout << "                          " << sum << endl; // so that the vector values are used.
}


void runSinTests() {
  testSin<float> ([] (float  v) { return std::cos(v); });
  testSin<float> ([] (float  v) { return cosf(v); });
  testSin<double>([] (double v) { return std::cos(v); });
  testSin<double>([] (double v) { return cos(v); });
}
TEST(MeasureSin,measure) {
  // run the tests 3 times in a row
  for(int i=0; i<3; ++i) {
    runSinTests();
  }
}
