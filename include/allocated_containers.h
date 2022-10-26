/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {

namespace adaptive_stack_allocated {

template<typename T>
using vector = std::vector<T, StackAllocator<T>>;

template<typename K, typename V, class Compare = std::less<K>>
using map = std::map<K, V, Compare, StackAllocator<std::pair<const K,V>>>;

template<typename K, typename V, class Hash = std::hash<K>, class Pred = std::equal_to<K>>
using unordered_map = std::unordered_map<K, V, Hash, Pred, StackAllocator<std::pair<const K,V>>>;

} // NS adaptive_stack_allocated

namespace aP {

template<typename T>
using Alloc = AlignedAllocator<T, Alignment::PAGE>;

} // NS aP

namespace a64 {

template<typename T>
using Alloc = AlignedAllocator<T, Alignment::CACHE_LINE>;

template<typename T>
using vector = std::vector<T, Alloc<T>>;

} // NS a64

namespace monotonic {

namespace aP {
template<typename T>
using Alloc = MonotonicAllocator<AlignedAllocator<T, Alignment::PAGE>>;

} // NS monotonic::aP

template<typename BaseAllocator>
using vector = std::vector<typename BaseAllocator::value_type, MonotonicAllocator<BaseAllocator>>;
} // NS monotonic

}
