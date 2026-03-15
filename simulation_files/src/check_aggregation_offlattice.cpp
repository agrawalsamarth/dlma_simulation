#include <check_aggregation_offlattice.hh>

namespace simulation{

template<typename type>
check_aggregation_offlattice<type>::check_aggregation_offlattice(system<type> *system_state, normal_bind<type> *bind_system, aggregation_condition<type> *ref_condition, double tolerance){

    sys_state     = system_state;
    bind_sys      = bind_system;
    condition     = ref_condition;
    agg_tolerance = tolerance;
}

template <typename type>
void check_aggregation_offlattice<type>::check_for_aggregation(constituent<type> *c_1){

    bool is_checked = false;
    int  particle_id;

    cluster_id = c_1->get_id();

    for (int i = 0; i < c_1->get_size(); i++){

        neighbours  = c_1->get_neighbour_list(i);
        particle_id = c_1->get_element_id(i); 

        for (int j = 0; j < neighbours.size(); j++) { 

            neighbour_id = neighbours[j];

            neighbour_cluster_id = sys_state->get_id_map(neighbour_id);
            temp = sys_state->get_aggregate(neighbour_cluster_id);

            if (!temp){
                std::cout<<"NULL pointer encountered"<<std::endl;
                exit(EXIT_FAILURE);
            }

            if (neighbour_cluster_id != cluster_id){

                particle_1 = sys_state->get_particle_by_index(particle_id);
                particle_2 = sys_state->get_particle_by_index(neighbour_id);

                distance = sys_state->get_interparticle_distance(particle_1, particle_2);

                if ((distance) <= ((1.+agg_tolerance) * (0.5 * (particle_1->get_diameter()+particle_2->get_diameter()))))
                {
                    if (condition->agg_condition(particle_1, particle_2)){
                        //sys_state->add_attachment(particle_id, neighbour_id);
                        check_for_aggregation(bind_sys->bind_aggregates(c_1, temp));
                        is_checked = true;
                    }
                }
                
            }

            if (is_checked)
                break;
            
        }

        if (is_checked)
            break;

    }


}

template class check_aggregation_offlattice<int>;
template class check_aggregation_offlattice<double>;

}