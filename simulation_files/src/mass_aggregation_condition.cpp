#include <mass_aggregation_condition.hh>

namespace simulation{

template<typename type>
mass_aggregation_condition<type>::mass_aggregation_condition(system<type> *ref_sys)
{

    sys_state = ref_sys;

}

template<typename type>
void mass_aggregation_condition<type>::show_out()
{ std::cout<<"object"<<std::endl;}

template<typename type>
bool mass_aggregation_condition<type>::agg_condition(constituent<type> *c_1, constituent<type> *c_2)
{

    id_1 = c_1->get_id();
    id_2 = c_2->get_id();

    agg_1 = sys_state->get_id_map(id_1);
    agg_2 = sys_state->get_id_map(id_2);

    cluster_1 = sys_state->get_aggregate(agg_1);
    cluster_2 = sys_state->get_aggregate(agg_2);

    if ((cluster_1->get_mass() + cluster_2->get_mass()) >= sys_state->get_seedmass()){
        return true;
    }

    else
        return false; 


}

template class mass_aggregation_condition<int>;
template class mass_aggregation_condition<double>;


}