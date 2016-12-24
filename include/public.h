#pragma once

#include <cassert>
#include <cstddef>
#include <iostream>
#include <memory>
#include <vector>

#include "pool.h"
#include "allocator.hpp"

namespace imajuscule {
    
    template<typename T>
    using pool_vector = std::vector<T, Allocator<T>>;
    
}
