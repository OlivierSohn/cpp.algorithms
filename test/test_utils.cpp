namespace imajuscule {

std::vector<Scaling> mkNaiveScaling(int firstSz, int const countCoeffs) {
    std::vector<Scaling> scalings;
    int remainingCoeffs = countCoeffs;
    int sz = firstSz;
    while(remainingCoeffs > 0) {
        scalings.emplace_back(sz /* sz */,
                              1 /* repeat */);
        remainingCoeffs -= sz;
        sz *= 2;
    }
    return scalings;
}

static inline float randf(float high = 1.f, float low = 0.f)
{
  return low + ((high-low) * ((float)std::rand() / (float)(RAND_MAX + 1.f)));
}

}
