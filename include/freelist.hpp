/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {

    namespace freelist {

        template<typename Integral>
        struct Links {
            template<typename T>
            static Integral fromIndex(size_t index, T const & elements ) { return index; }

            template<typename T, typename VALUE = typename T::value_type>
            static Integral fromPtr( VALUE * ptr, T const & elements ) { return ptr - &elements[0].value; }

            template<typename T>
            static auto & follow(Integral index, T & elements) { return elements[index].value; }

            static constexpr Integral getNull() { return std::numeric_limits<Integral>::max(); }
        };

        template<>
        struct Links<void*> {
            template<typename T>
            static auto fromIndex(size_t index, T & elements ) { return &elements[index]; }

            template<typename T>
            static auto fromPtr(decltype(std::declval<T>()[0].value) * ptr, T const & elements ) { return ptr; }

            template<typename T>
            static auto & follow(void * ptr, T & elements) { return *reinterpret_cast<decltype(std::declval<T>()[0].value)*>(ptr); }

            static constexpr void* getNull() { return nullptr; }
        };


        template<typename Integral>
        struct MaxSize {
            static constexpr size_t value = std::numeric_limits<Integral>::max() - 1;
        };

        template<>
        struct MaxSize<void*> {
            static constexpr size_t value = std::numeric_limits<uint16_t>::max(); // else an array-based freelist is probably not optimal
        };

        /* SIZE_BLOCK is the count of T objects stored contiguously.
        The free list will recursively instantiate sub free lists to
        handle overflow. */
        template< typename T, size_t SIZE_BLOCK, typename LINK_T >
        struct FreeListImpl {
            using value_type = T;
            using link_type = LINK_T;

            static_assert(sizeof(value_type) >= sizeof(link_type)); // else the free list is not optimal

            static constexpr auto sizeBlock = SIZE_BLOCK;

          private:
            using MyT = FreeListImpl<T,SIZE_BLOCK,LINK_T>;
            static constexpr auto max_size = MaxSize<LINK_T>::value - 1; // -1 because one element is used to store head
            static constexpr auto null_() { return Links<LINK_T>::getNull(); };

            static_assert(sizeBlock < max_size);

            struct Union {
                union {
                    T value;
                    LINK_T link;
                };
            };

          public:
            FreeListImpl() {
                initialize();
            }

            value_type* Take() {
                auto & head_ = head();
                if(head_.link == null_()) {
                    return useNext().Take();
                }
                auto & ret = Links<LINK_T>::follow(head_.link, elements);
                head_.value = ret;

                assert(isOurs(&ret));
                return &ret;
            }

            void Return(T*ptr) {
                if(isOurs(ptr)) {
                    auto link = Links<LINK_T>::fromPtr(ptr, elements);
                    auto & head_ = head();
                    auto & val = Links<LINK_T>::follow(link, elements);
                    *reinterpret_cast<link_type*>(&val) = head_.link;
                    head_.link = link;
                }
                else {
                    useNext().Return(ptr);
                }
            }

        private:
            std::array<Union, sizeBlock+1> elements;
            std::unique_ptr<MyT> next;

            // I chose the last element as head
            auto & head() { return elements.back(); }

            void initialize() {
                head().link = Links<LINK_T>::fromIndex(0, elements);
                size_t end = elements.size()-1;
                for(size_t i=1; i<end; ++i) {
                    auto next = Links<LINK_T>::fromIndex(i, elements);
                    elements[i-1].link = next ;
                }
                elements[end-1].link = null_();
            }

            bool isOurs(T*ptr) const {
                std::ptrdiff_t i = ptr - &elements[0].value;
                return i >= 0 && i < sizeBlock;
            }

            auto & useNext() {
                if(!next) {
                    next = std::make_unique<MyT>();
                }
                return *next;
            }
        };
    }

    template< typename T, size_t SIZE, typename LINK_T = void* >
    using FreeList = freelist::FreeListImpl<T, SIZE, LINK_T>;

} // NS imajuscule
