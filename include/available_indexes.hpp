
namespace imajuscule
{
    // for vector
    template<class T, class Elem>
    void make_room_for_index(T index, std::vector<Elem> & c)
    {
        if(c.size() < index+1) {
            c.resize(index+1);
        }
    }

    // for maps, lists, etc...
    template<class T, class Container>
    void make_room_for_index(T index, Container & c)
    {
    }
    
    /*
     * To keep track of which indexes are in-use in a container.
     *
     * Works only if :
     * - the container adds its elements at the index resulting of a call to Take
     * - the indexes of removed elements are Returned by calling Return
     * - THere is no overflow
     */
    template<typename T>
    struct AvailableIndexes {
        using index_t = T;
        
        void reserve(size_t size) {
            availables.reserve(size);
        }
        
        template<typename Container>
        T Take(Container & container) {
            if (!availables.empty()) {
                auto key = availables.back();
                availables.pop_back();
                return key;
            }
            assert(container.size() < std::numeric_limits<T>::max());
            T index = safe_cast<T>(container.size());
            make_room_for_index(index, container);
            return index;
        }

        void Return(T index) {
            availables.emplace_back(index);
        }
        
        auto size() {
            return availables.size();
        }
    private:
        std::vector<T> availables;
    };

}
