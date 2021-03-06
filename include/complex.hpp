/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */


namespace imajuscule
{
    // stl's complex is slow because it supports NaN / Inf, so I write my own
    // that doesn't support Nan / Inf
    
    template<typename T>
    struct complex {
        using FPT = T;
        using Tr = NumTraits<T>;
        
        explicit complex()
        : re{}, im{} {};
        
        explicit complex(T const re)
        : re(re), im(0) {}

        complex(T const re, T const im)
        : re(re), im(im) {}
        
        T real() const { return re; }
        T imag() const { return im; }
        
        T * realAddr() { return &re; }
        T * imagAddr() { return &im; }
        
        complex & operator *=(complex const & z) {
            T c = re;
            T d = im;
            re = z.real() * c - z.imag() * d;
            im = z.real() * d + z.imag() * c;
            return *this;
        }
        
        complex & operator +=(complex const & z) {
            re += z.real();
            im += z.imag();
            return *this;
        }
        
        complex & operator *=(T const v) {
            re *= v;
            im *= v;
            return *this;
        }
        
        void convert_to_conjugate() {
            im = -im;
        }
        
    private:
        T re, im;
    };
    
    template<typename T>
    complex<T> operator -(complex<T> const & a, complex<T> const & b) {
        return {
            a.real() - b.real(),
            a.imag() - b.imag()
        };
    }
    
    
    template<typename T>
    complex<T> operator +(complex<T> const & a, complex<T> const & b) {
        return {
            a.real() + b.real(),
            a.imag() + b.imag()
        };
    }
    
    template<typename T>
    complex<T> operator *(complex<T> const & a, complex<T> const & b) {
        auto ret = a;
        ret *= b;
        return ret;
    }
    
    /*
     '/' operator is not provided to encourage the user to replace divisions by multiplications
     */
    
    template<typename T>
    auto operator *(T f, complex<std::remove_const_t<T>> const & b) {
        return decltype(b) {
            f*b.real(),
            f*b.imag()
        };
    }
    
    template<typename T, typename C = complex<std::remove_const_t<T>>>
    auto polar(T theta) -> C {
        return {
            std::cos(theta),
            std::sin(theta)
        };
    }
    
    template<typename T>
    T norm(complex<T> const & z) {
        return z.real()*z.real() + z.imag()*z.imag();
    }
    
    template<typename T>
    T abs(complex<T> const & z) {
        return std::sqrt(norm(z));
    }
    
    template<typename T>
    auto conj(complex<T> const & z) -> complex<T> {
        return {z.real(), -z.imag()};
    }
    
    template<typename T>
    T arg(complex<T> const & z) {
        return std::atan2(z.imag(), z.real());
    }
    
    
    template<typename T, typename ITER = typename std::vector<T>::iterator>
    std::vector<complex<T>> complexify(ITER it, ITER end) {
        std::vector<complex<T>> ret;
        ret.reserve(std::distance(it, end));
        for(; it!=end; ++it) {
            ret.push_back({*it, 0});
        }
        return ret;
    }
}
