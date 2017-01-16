
namespace imajuscule {
    namespace adaptive_stack_allocated {
        template<typename T>
        using vector = std::vector<T, StackAllocator<T>>;
        
        template<typename K, typename V, class Compare = std::less<K>>
        using map = std::map<K, V, Compare, StackAllocator<std::pair<const K,V>>>;
        
        template<typename K, typename V, class Hash = std::hash<K>, class Pred = std::equal_to<K>>
        using unordered_map = std::unordered_map<K, V, Hash, Pred, StackAllocator<std::pair<const K,V>>>;
    }
}
