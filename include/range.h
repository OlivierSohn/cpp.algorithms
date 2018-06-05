/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    // don't change these values, see body_chunk_positinoned::store method
    enum where : signed char { lower = -1, iN = 0, UPPER = 1 };

    enum splitMethod {
        PixelWidth = 0,
        NoPixelWidth = 1
    };

    template <class T>
    struct range
    {
        using Tr = NumTraits<T>;

        range() : min_ (Tr::one()), max_(-Tr::one()) {}

        range(const range & other ) : min_(other.min_), max_(other.max_) {}
        range(range && o) : min_(std::move(o.min_)), max_(std::move(o.max_)) {}

        range & operator=(const range & other) {
            min_ = other.min_;
            max_ = other.max_;
            return *this;
        }

        range& operator=(range&& o) {
            min_ = std::move(o.min_);
            max_ = std::move(o.max_);
            return *this;
        }

        range(T min, T max) {
            assert(min <= max);
            set(min, max);
        }

        bool operator == ( const range & other ) const {
            return !( (*this)!= other );
        }

        bool operator != ( const range & other ) const;

        bool before( range const & other ) const {
            assert( !empty() );

            if( std::is_integral<T>() ) {
                return ( max_ < other.min_ );
            } else {
                auto tmp = other.avg_ - avg_;
                tmp -= ( other.halfSpan_ + halfSpan_ );
                return tmp > 0.f;
            }
        }

        T minTo( T val ) const {
            assert( !empty() );

            if( std::is_integral<T>() ) {
                return val - min_;
            } else {
                val -= halfSpan_;
                return val + halfSpan_;
            }
        }
        T maxTo( T val ) const {
            assert( !empty() );

            if( std::is_integral<T>() ) {
                return val - max_;
            } else {
                val -= halfSpan_;
                return val - halfSpan_;
            }
        }

        T getAt(T ratio) const {
            assert(ratio >= 0.f);
            assert(ratio <= 1.f);
            assert(!empty());
            assert(!std::is_integral<T>());
            return getMin() + ratio * getSpan();
        }

        bool empty() const {
            if( std::is_integral<T>() ) {
                return max_ < min_;
            } else {
                return halfSpan_ < Tr::zero();
            }
        }

        void make_empty() {
            if( std::is_integral<T>() ) {
                min_ = Tr::one();
                max_ = -Tr::one();
            } else {
                halfSpan_ = -Tr::one();
            }
        }

        T delta() const {
            assert(!empty());
            if( std::is_integral<T>() ) {
                return max_ - min_;
            } else {
                return 2.f * halfSpan_;
            }
        }

        void set(T val) {
            if( std::is_integral<T>() ) {
                min_ = max_ = val;
            } else {
                assert(!std::isnan((float)val ));
                avg_ = val;
                halfSpan_ = 0.f;
            }
        }

        void set(T Min, T Max) {
            assert(Min <= Max);
            if( std::is_integral<T>() ) {
                min_ = Min;
                max_ = Max;
            } else {
                assert(!std::isnan((float)Min ));
                assert(!std::isnan((float)Max ));
                avg_ = 0.5f * (Min + Max);
                halfSpan_ = 0.5f * ( Max - Min );
            }
        }

        template <class U = T, typename = std::enable_if<!std::is_integral<U>::value>>
        void setByAvgHalfSpan( T avg, T halfSpan) {
            static_assert(!std::is_integral<T>(), "wrong type");
            avg_ = avg;
            assert( halfSpan >= Tr::zero() );
            halfSpan_ = halfSpan;
        }

        bool extend(T val);
        void extend(range<T> const & r);

        void includeMargin(T margin) {
            assert(margin <= Tr::zero() || !empty()); // adding a margin to an empty range makes no sense
            if( std::is_integral<T>() ) {
                min_ -= margin;
                max_ += margin;
            } else {
                halfSpan_ += margin;
            }
        }

        bool contains(T val) const {
            if( empty() ) {
                return false;
            }
            if( std::is_integral<T>() ) {
                return (val <= max_ && val >= min_);
            } else {
                val -= avg_;
                return( (val <= halfSpan_) && ( val >= - halfSpan_ ) );
            }
        }

        bool contains(range const & r) const {
            return contains(r.getMin()) && contains(r.getMax());
        }

        where whereIs(T val) const {
            assert(!empty());
            if ( std::is_integral<T>() ) {
                if ( val < min_ ) {
                    return lower;
                } else if ( val > max_ ) {
                    return UPPER;
                }
            } else {
                val -= avg_;
                if ( val < -halfSpan_ ) {
                    return lower;
                } else if ( val > halfSpan_ ) {
                    return UPPER;
                }
            }

            return iN;
        }

        bool overlaps( range const & r ) const;

        range intersection( range const & r ) const;

        template <class U = T, typename = std::enable_if<!std::is_integral<U>::value>>
        void trim( T val ) {
            static_assert( !std::is_integral<T>(), "wrong type" );
            assert( val >= Tr::zero() );

            // don't test for that : it's ok if an empty range is reduced, it is still empty after ( halfSpan_ < 0.f )
            //if ( halfSpan_ >= 0.f ) {
                halfSpan_ *= val;
            //}
        }

        T getSpan() const {
            if( std::is_integral<T>() ) {
                return max_ - min_;
            } else {
                return 2.f * halfSpan_;
            }
        }

        T getMax() const {
            if( std::is_integral<T>() ) {
                return max_;
            } else {
                return avg_ + halfSpan_;
            }
        }
        T getMin() const {
            if( std::is_integral<T>() ) {
                return min_;
            } else {
                return avg_ - halfSpan_;
            }
        }

        /*
         Waring: for integral types, rounds the value.
         */
        T getCenter() const {
            if( std::is_integral<T>() ) {
                return ( min_ + max_ ) / 2;
            } else {
                return avg_;
            }
        }

        T getExpCenter() const {
            return exp_mean(getMin(), getMax());
        }

        void setMin(T val) {
            if( std::is_integral<T>() ) {
                min_ = val;
            } else {
                if( empty() ) {
                    set( val );
                } else {
                    val -= avg_;
                    auto span = halfSpan_ - val;
                    if( span < 0.f ) {
                        // newMin > oldMax
                        set( val );
                    } else {
                        auto diffMin = val + halfSpan_;
                        avg_ += diffMin * 0.5f;
                        halfSpan_ = 0.5f * span;
                    }
                }
            }
        }

        void setMax(T val) {
            if( std::is_integral<T>() ) {
                max_ = val;
            } else {
                if( empty() ) {
                    set( val );
                } else {
                    val -= avg_;
                    auto span = val + halfSpan_;
                    if( span < 0.f ) {
                        // newMax < oldMin
                        set( val );
                    } else {
                        auto diffMax = val - halfSpan_;
                        avg_ += diffMax * 0.5f;
                        halfSpan_ = 0.5f * span;
                    }
                }
            }
        }

        T minDist( T val ) const;

        void translate( T val ) {
            if( std::is_integral<T>() ) {
                min_ += val;
                max_ += val;
            } else {
                avg_ += val;
            }
        }

        range & scaleNormalized(float minScale, float maxScale) {
            assert(!empty());
            auto d = delta();
            auto m = getMin();
            set(m + minScale * d,
                m + maxScale * d);
            return *this;
        }

        void homothety( T origin, T factor ) {
            assert( factor >= Tr::zero() );

            if( std::is_integral<T>() ) {
                min_ = origin + (min_ - origin) * factor;
                max_ = origin + (max_ - origin) * factor;
            } else {
                avg_ = origin + (avg_ - origin) * factor;
                halfSpan_ *= factor;
            }
        }

        template < splitMethod METHOD >
        void split( float ratio, range & r1, range & r2) const {

            assert( !empty() );

            if( std::is_integral<T>() ) {
                if( PixelWidth == METHOD) {
                    // min == max == 0 means a pixel of width 1
                    auto const pixel_size = Tr::one();

                    T length = max_ + pixel_size - min_;

                    T ratio_length = (T)((((float)length) * ratio) + 0.5f);

                    r2.min_ = min_ + ratio_length;
                    r1.max_ = r2.min_ - pixel_size;
                } else {
                    T length = max_ - min_;

                    T ratio_length = (T)((((float)length) * ratio) + 0.5f);

                    r2.min_ = min_ + ratio_length;
                    r1.max_ = r2.min_;
                }
            } else {
                r1.halfSpan_ = halfSpan_ * ratio;
                auto diff = halfSpan_ - r1.halfSpan_;
                assert( diff >= 0.f);
                r2.halfSpan_ = diff;

                r1.avg_ = avg_ - diff;
                r2.avg_ = avg_ + r1.halfSpan_;
            }
        }


    private:
        union {
            T min_; // integral types
            T avg_; // floating types
        };
        union {
            T max_; // integral types
            T halfSpan_; // floating types
        };
    };

    extern template struct range<double>;
    extern template struct range<float>;
    extern template struct range<int32_t>;
    extern template struct range<int16_t>;

    template<typename T>
    void logRange(range<T> r) {
        if(r.empty()) {
            LG(INFO, "[empty]");
        }
        else {
            LG(INFO, "[%f .. %f]", r.getMin(), r.getMax());
        }
    }
}
