/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
template <class T>
bool range<T>::operator != ( const range & other ) const {
    if( other.empty() ) {
        return !empty();
    }
    if( empty() ) {
        return true;
    }
    if( std::is_integral<T>() ) {
        return (max_ != other.max_) || (min_ != other.min_);
    } else {
        return (avg_ != other.avg_) || (halfSpan_ != other.halfSpan_);
    }
}

template <class T>
bool range<T>::extend(T val)
{
    if(empty()) {
        set(val);
        return true;
    }
    
    if( std::is_integral<T>() ) {
        if(val > max_) {
            max_ = val;
            return true;
        }
        
        if(val < min_) {
            min_ = val;
            return true;
        }
    }
    else {
        assert(!std::isnan( (float)val ));
        
        // compute wrt avg_ for precision
        
        val -= avg_;
        auto diffMax = val - halfSpan_;
        if( diffMax > 0.f ) {
            halfSpan_ = 0.5f * ( halfSpan_ + val );
            avg_ += diffMax * 0.5f;
            return true;
        } else {
            auto diffMin = val + halfSpan_;
            if( diffMin < 0.f ) {
                halfSpan_ = 0.5f * ( halfSpan_ - val );
                avg_ += diffMin * 0.5f;
                return true;
            }
        }
    }
    
    return false;
}

template <class T>
void range<T>::extend(range<T> const & r)
{
    if(r.empty()) {
        return;
    }
    
    if( std::is_integral<T>() ) {
        extend( r.min_ );
        extend( r.max_ );
    } else {
        
        // compute wrt avg_ for precision
        
        auto otherRelAvg = r.avg_ - avg_;
        auto relMax = otherRelAvg + r.halfSpan_;
        auto relMin = otherRelAvg - r.halfSpan_;
        
        auto diffMax = relMax - halfSpan_;
        if( diffMax >= 0.f ) {
            if( relMin <= -halfSpan_ ) {
                avg_ = r.avg_;
                halfSpan_ = r.halfSpan_;
            } else {
                halfSpan_ = 0.5f * ( halfSpan_ + relMax);
                avg_ += diffMax * 0.5f;
            }
        } else {
            auto diffMin = relMin + halfSpan_;
            if( diffMin < 0.f ) {
                halfSpan_ = 0.5f * ( halfSpan_ - relMin );
                avg_ += diffMin * 0.5f;
            }
        }
    }
}

template <class T>
T range<T>::minDist( T val ) const
{
    T distA, distB;
    
    if( std::is_integral<T>() ) {
        distA = min_ - val;
    } else {
        val -= avg_;
        distA = - halfSpan_ - val;
    }
    
    if( distA >= Tr::zero() ) {
        return distA;
    }
    
    if( std::is_integral<T>() ) {
        distB = val - max_;
    } else {
        distB = val - halfSpan_;
    }
    
    if( distB >= Tr::zero() ) {
        return distB;
    }
    
    assert( distA < Tr::zero() );
    assert( distB < Tr::zero() );
    return std::max( distA, distB );
}

template <class T>
bool range<T>::overlaps( range<T> const & r ) const {
    if( r.empty() ) {
        return false;
    }
    if( empty() ) {
        return false;
    }
    
    if( std::is_integral<T>() ) {
        if( min_ > r.max_ ) { return false; }
        if( max_ < r.min_ ) { return false; }
    } else {
        // compute wrt avg_ for precision
        auto otherRelAvg = r.avg_ - avg_;
        
        auto relMax = otherRelAvg + r.halfSpan_;
        if( - halfSpan_ > relMax ) { return false; }
        
        auto relMin = otherRelAvg - r.halfSpan_;
        if( halfSpan_ < relMin ) { return false; }
    }
    
    return true;
}

template <class T>
range<T> range<T>::intersection( range<T> const & r ) const {
    range inter;
    
    if( !r.empty() && !empty() ) {
        if( std::is_integral<T>() ) {
            if( ( min_ <= r.max_ ) && ( r.min_ <= max_ ) ) {
                auto Min = std::max( min_, r.min_ );
                auto Max = std::min( max_, r.max_ );
                assert(Max >= Min);
                inter.set( Min, Max );
            }
        } else {
            auto otherRelAvg = r.avg_ - avg_;
            
            auto relMax = otherRelAvg + r.halfSpan_;
            auto span = relMax + halfSpan_;
            
            if( span >= 0.f ) {
                auto relMin = otherRelAvg - r.halfSpan_;
                auto span2 = halfSpan_ - relMin;
                
                if( span2 >= 0.f ) {
                    auto diffMax = relMax - halfSpan_;
                    
                    if( -halfSpan_ > relMin ) {
                        inter = *this;
                        
                        if( diffMax < 0.f ) {
                            inter.avg_ += diffMax * 0.5f;
                            inter.halfSpan_ = 0.5f * span;
                        }
                    } else {
                        inter = r;
                        
                        if( diffMax > 0.f ) {
                            inter.avg_ -= diffMax * 0.5f;
                            inter.halfSpan_ = 0.5f * span2;
                        }
                    }
                }
            }
        }
    }
    
    return inter;
}

template struct range<double>;
template struct range<float>;
template struct range<int32_t>;
template struct range<int16_t>;
template struct range<uint16_t>;
}
