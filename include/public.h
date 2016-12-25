#pragma once

#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>
#include <map>

#include "pool.h"
#include "allocator.hpp"

namespace imajuscule {
    namespace pool {
        template<typename T>
        using vector = std::vector<T, Allocator<T>>;
        
        template<typename K, typename V, class Compare = std::less<K>>
        using map = std::map<K, V, Compare, Allocator<std::pair<const K,V>>>;
    }
}
