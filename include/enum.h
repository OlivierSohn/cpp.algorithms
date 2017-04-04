/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule
{
    struct enumTraversal
    {
        enumTraversal(unsigned int low, unsigned int upperBound, std::function<const char*(int)> realValueDesc);
        enumTraversal(unsigned int low, unsigned int upperBound, std::function<const char*(int)> realValueDesc, std::function<bool(int)> valueIsReal);
        virtual ~enumTraversal() {}
        
        std::string valToString(int) const;
        bool valFromString(const std::string &, int &) const;
        bool valToRealValueIndex(int, int &) const;
        const std::vector<int> & realValues() const;

        bool valIsReal(int) const;

    private:
        std::function<bool(int)> mValueIsReal;
        std::function<const char*(int)> mRealValueDesc;
        uint8_t mLowerBound:7, mUpperBound:7; // lower bound is included, upper bound is excluded
        mutable bool mRealValuesComputed:1;
        mutable std::vector<int> mRealValues;
    };
}
