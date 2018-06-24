
namespace imajuscule {
  template<typename T>
  std::unique_ptr<T> extractFromEnd(std::vector<std::unique_ptr<T>> &c, T*o) {
    for(auto i = c.rbegin(), e = c.rend(); i!=e; ++i) {
      auto & p = *i;
      if(p.get() != o) {
        continue;
      }
      std::unique_ptr<T> res;
      res.swap(p);
      c.erase(std::next(i).base());
      return res;
    }
    return {};
  }

  namespace detail {
    template<typename C, typename Lock>
    void makeSomeRoom(int n, int32_t cap1, C & container, Lock & l)
    {
      {
        using V = std::vector<typename C::value_type>;
        static constexpr auto growth_factor = 2;
        V v;

        // allocation happens outside the lock scope:
        {
          auto newCap = std::min(cap1 + 1, growth_factor * cap1);
          Assert(newCap < 10000000);
          v.reserve(newCap);
        }

        l.lock();
        container.trySwapUnderlyingContainer(v);
        l.unlock();

        // deallocation happens outside the lock scope:
        // we make vector destruction explicit for clarity but we could also
        // remove this line and let the vector be destructed when going out of scope.
        v = V{};
      }

      l.lock();

      if(auto amount = container.shouldGrow(n)) {
        auto cap3 = container.underlyingContainerCapacity();
        l.unlock();

        makeSomeRoom(n, cap3+amount, container, l);
      }
    }
  }

  enum class CanRealloc {
    Yes,
    No
  };

  /*
   * Reallocates, if needed, the underlying container so that 'n'
   * additional elements can be pushed. When this function returns, the lock 'l'
   * is being taken.
   *
   * 'n' : the number of elements that we want to be able to add to the container.
   *
   * returns false if a reallocation is needed but canRealloc == CanRealloc::No
   * hence it could not be done.
   */
  template<CanRealloc canRealloc, typename C, typename Lock>
  [[nodiscard]] bool reserveAndLock(int n, C & container, Lock & l) {
    l.lock();

    if(auto amount = container.shouldGrow(n)) {
      if constexpr (canRealloc == CanRealloc::No) {
        return false;
      }
      auto cap1 = container.underlyingContainerCapacity();

      l.unlock();

      detail::makeSomeRoom(n, cap1+amount, container, l);
      // the lock was taken by the previous call.
    }
    return true;
  }

  template<typename T>
  struct vectorWrapper {
    using container = std::vector<T>;
    using value_type = typename container::value_type;
    vectorWrapper(container & v) : v(v) {}

    /*
     * Returns the additional capacity that should be added to the current underlying container
     * in order to support adding 'nAdds' element.
     *
     * 0 is returned when we can add the elements without growing, else
     * a strictly positive number is returned.
     */
    int shouldGrow(int nAdds) const noexcept {
      int nElems = v.size();
      int maxCount = v.capacity();
      int diff = nElems + nAdds - maxCount;
      return std::max(0,diff);
    }

    int32_t underlyingContainerCapacity() const noexcept {
      return v.capacity();
    }

    /*
     * This method will swap the underlying container if the passed container has
     * a strictly bigger capacity than the original one.
     *
     * @Returns true if the swap has been done, false otherwise.
     *
     * If the swap was done, the original storage is returned in c.
     */
    bool trySwapUnderlyingContainer(container & c) {
      Assert(c.size() == 0);

      if(c.capacity() <= v.capacity()) {
        return false;
      }

      // move the elements to the new container
      std::move(v.begin(),v.end(),std::back_insert_iterator<container>(c));
      std::swap(c,v);
      return true;
    }

  private:
    container & v;
  };

  template<typename T>
  auto mkVectorWrapper(std::vector<T>&v) {
    return vectorWrapper<T>{v};
  }

} // NS imajuscule
