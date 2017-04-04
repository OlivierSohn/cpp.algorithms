/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

using namespace imajuscule;

enumTraversal::enumTraversal(unsigned int low, unsigned int up, std::function<const char*(int)> real_value_desc, std::function<bool(int)> value_real) :
    mLowerBound(low),
    mUpperBound(up),
    mValueIsReal(value_real),
    mRealValueDesc(real_value_desc),
    mRealValuesComputed(false)
{
    assert(low < up);
}

enumTraversal::enumTraversal(unsigned int low, unsigned int up, std::function<const char*(int)> real_value_desc) :
mLowerBound(low),
mUpperBound(up),
mValueIsReal([this](int v)->bool{return (v>=mLowerBound && v<mUpperBound);}),
mRealValueDesc(real_value_desc),
mRealValuesComputed(false)
{
    assert(low < up);
}

std::string enumTraversal::valToString(int val) const
{
    if (!mValueIsReal(val)) {
        return "invalid";
    }
    return mRealValueDesc(val);
}

bool enumTraversal::valToRealValueIndex(int val, int & index) const {
    int i=0;
    for(auto v : realValues()) {
        if(val==v) {
            index = i;
            return true;
        }
        ++i;
    }
    return false;
}

bool enumTraversal::valFromString(const std::string & s, int & i) const
{
    for(auto val : realValues()) {
        if (iequals(s, std::string(mRealValueDesc(val))))
        {
            i = val;
            return true;
        }
    }
    
    return false;
}

const std::vector<int> & enumTraversal::realValues() const
{
    if (!mRealValuesComputed) {
        mRealValuesComputed = true;

        for (int i = mLowerBound; i < mUpperBound; i++)
        {
            if (mValueIsReal(i)) {
                mRealValues.push_back(i);
            }
        }
    }
    return mRealValues;
}

bool enumTraversal::valIsReal(int val) const
{
    return mValueIsReal(val);
}
