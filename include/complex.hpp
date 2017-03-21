

namespace imajuscule
{
    // stl's complex is slow because it supports NaN / Inf, so I write my own
    // that doesn't support Nan / Inf
    
    template<typename T>
    struct complex {
        using FPT = T;
        using Tr = NumTraits<T>;
        
        explicit complex() = default;
        complex(T re, T im)
        : re(re), im(im) {}
        
        T real() const { return re; }
        T imag() const { return im; }
        
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
        
        complex & operator *=(T v) {
            re *= v;
            im *= v;
            return *this;
        }
        
        complex & operator /=(T v) {
            assert(v != Tr::zero());
            return operator*=(Tr::one()/v);
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
    complex<T> operator *(float f, complex<T> const & b) {
        return {
            f*b.real(),
            f*b.imag()
        };
    }
    
    template<typename T>
    complex<T> polar(T theta) {
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
    T arg(complex<T> const & z) {
        return std::atan2(z.imag(), z.real());
    }
    
    
}
