namespace imajuscule {

template<typename T>
bool areNear(T a, T b, double eps) {
  if(std::abs(a) <= eps && std::abs(b) <= eps) {
    return true;
  }
  if(std::abs(a-b) <= eps) {
    return true;
  }
  return std::abs(a-b)/std::abs(std::max(a,b)) <= eps;
}

template<typename T>
std::unique_ptr<T> extractFromEnd(std::vector<std::unique_ptr<T>> &c, T*o) {
  for(auto i = c.rbegin(), e = c.rend(); i!=e; ++i) {
    auto & p = *i;
    if(p.get() != o) {
      continue;
    }
    std::unique_ptr<T> res;
    res.swap(p);
    c.erase(std::next(i).base());
    return res;
  }
  return {};
}

template<typename T>
static inline std::vector<T> mkDirac(int sz, T amplitude = 1.) {
  std::vector<T> res;
  res.resize(sz);
  if(res.empty()) {
    throw std::logic_error("0-size");
  }
  res[0] = amplitude;
  return res;
}

template<typename T>
auto mkTestCoeffs(int const sz){
  a64::vector<T> v(sz);
  auto index = 0;
  for(auto & value: v) {
    value = (sz - index) / static_cast<T>(sz);
    ++index;
  }
  return std::move(v);
}

// the coefficients should tend to 0, to mimic real responses
// (else, with scaling, the end of the dirac's response will be wrong)
static inline a64::vector<double> mkCoefficientsTriangle(int sz) {
    a64::vector<double> res;
    res.reserve(sz);
    constexpr int period = 100;
    for(int i=0; i<sz; ++i) {
        auto j = i%(2*period);
        if(j < period) {
            res.push_back(1.-(j/(double)period));
        }
        else {
            res.push_back((j-period)/(double)period);
        }
    }
    
    constexpr int fadeout = 40;
    if(sz > fadeout) {
        for(int i=0; i<fadeout; ++i) {
            res[res.size()-1-i] *= (i+1) / static_cast<double>(fadeout);
        }
    }
    return res;
}

std::vector<Scaling> mkNaiveScaling(int firstSz, int const countCoeffs);
std::vector<Scaling> mkBetterScaling(int firstSz, int const countCoeffs);


}
