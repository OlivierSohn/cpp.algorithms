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

}
