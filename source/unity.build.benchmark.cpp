#include "unity.build.cpp"

namespace imajuscule::bench::pow {
double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}
}

void bench_pow() {
  using namespace imajuscule;
  using namespace imajuscule::profiling;
  
  std::vector<double> exps;
  exps.reserve(10000000);
  
  for(int i = 0; i<10000000; ++i) {
    exps.push_back(static_cast<double>(1 + i/10000000.));
  }
  
  {
    double res{};
    
    auto dur = measure_thread_cpu_one([&res, &exps](){
      for(auto exp : exps) {
        res += std::pow(2., exp);
      }
    });
    std::cout << dur.count() << std::endl;
    std::cout << res << std::endl;
  } // 44000
  {
    double res{};
    
    auto dur = measure_thread_cpu_one([&res, &exps](){
      for(auto exp : exps) {
        res += imajuscule::bench::pow::fastPow(2., exp);
      }
    });
    std::cout << dur.count() << std::endl;
    std::cout << res << std::endl;
  } // 11000
  
  
  std::optional<int> idx;
  double sumErrors = 0;
  double sumRes = 0;
  double maxDiff;
  int i = -1;
  for(auto exp : exps) {
    ++i;
    double r1 = imajuscule::bench::pow::fastPow(2., exp);
    double r2 = std::pow(2., exp);
    sumRes += r2;
    auto diff = std::abs(r1 - r2);
    sumErrors += diff;
    if (!idx || diff > maxDiff) {
      maxDiff = diff;
      idx = i;
    }
  }
  std::cout << *idx << std::endl;
  double r1 = imajuscule::bench::pow::fastPow(2., exps[*idx]);
  std::cout << r1 << std::endl;
  double r2 = std::pow(2., exps[*idx]);
  std::cout << r2 << std::endl;
  std::cout << "avg res " << sumRes / exps.size() << std::endl;
  std::cout << "avg error " << sumErrors / exps.size() << std::endl;

  std::optional<double> prev;
  int same = 0;
  std::optional<int> maxSame;
  for(auto exp : exps) {
    double p = imajuscule::bench::pow::fastPow(2., exp);
    if (prev && *prev == p) {
      ++same;
    } else {
      same = 0;
    }
    if (!maxSame || *maxSame < same) {
      maxSame = same;
    }
    prev = p;
  }
  if (!maxSame || *maxSame < same) {
    maxSame = same;
  }
  std::cout << "max same " << *maxSame << std::endl;
}

void bench_log() {
  using namespace imajuscule;
  using namespace imajuscule::profiling;
  
  std::vector<double> exps;
  exps.reserve(10000000);
  
  for(int i = 0; i<10000000; ++i) {
    exps.push_back(static_cast<double>(1 + i/10000000.));
  }
  
  {
    double res{};
    
    auto dur = measure_thread_cpu_one([&res, &exps](){
      for(auto exp : exps) {
        res += std::log2(exp);
      }
    });
    std::cout << dur.count() << std::endl;
    std::cout << res << std::endl;
  } // 44000
  {
    double res{};
    
    auto dur = measure_thread_cpu_one([&res, &exps](){
      for(auto exp : exps) {
        res += sprout::log2(exp);
      }
    });
    std::cout << dur.count() << std::endl;
    std::cout << res << std::endl;
  } // 11000

}

int main() {
  std::cout << "pow" << std::endl
  bench_pow();
  std::cout << "log" << std::endl
  // shows that constexpr sprout::log is much more cpu intensive that non-constexpr std::log
  bench_log();
  return 0;
}
