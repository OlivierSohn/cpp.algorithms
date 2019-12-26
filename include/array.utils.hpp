
namespace imajuscule {

template<std::size_t i, std::size_t j>
auto arrayConcat(std::array<int, i> const & a,
                 std::array<int, j> const & b) -> std::array<int, i+j>
{
    std::array<int, i+j> r;
    for(int k =0; k<i; ++k) {
        r[k] = a[k];
    }
    for(int k =0; k<j; ++k) {
        r[i+k] = b[k];
    }
    return r;
}

template<std::size_t i, std::size_t j>
auto arraySplit(std::array<int, i+j> const & r) -> std::pair<std::array<int, i>, std::array<int, j>>
{
    std::array<int, i> a;
    std::array<int, j> b;
    for(int k =0; k<i; ++k) {
        a[k] = r[k];
    }
    for(int k =0; k<j; ++k) {
        b[k] = r[i+k];
    }
    return {a,b};
}

}
