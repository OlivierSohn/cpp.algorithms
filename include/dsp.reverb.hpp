
namespace imajuscule
{

/*
 
 */
template<typename C>
void dephase(int const total_instances,
             int const index_instance,
             C & rev)
{
    if constexpr(!C::has_subsampling) {
        Assert(total_instances);
        float ratio = index_instance / static_cast<float>(total_instances);
        rev.dephaseByGroupRatio(ratio);
    }
}


  
}
