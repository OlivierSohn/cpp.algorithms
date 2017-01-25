
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
mValueIsReal([this](int v)->bool{return (v>mLowerBound && v<mUpperBound);}),
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

bool enumTraversal::valFromString(const std::string & s, int & i) const
{
    const std::vector<int> & vec = realValues();
    for (auto val : vec)
    {
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

        for (int i = mLowerBound + 1; i < mUpperBound; i++)
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
int enumTraversal::low() const
{
    return mLowerBound;
}

int enumTraversal::high() const
{
    return mUpperBound;
}
