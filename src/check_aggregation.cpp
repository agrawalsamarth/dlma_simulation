#include <check_aggregation_onlattice.hh>

namespace simulation{

template<typename type>
check_aggregation_onlattice<type>::check_aggregation_onlattice(system<type> *system_state, normal_bind<type> *bind_system, aggregation_condition<type> *ref_condition){

    sys_state = system_state;
    bind_sys  = bind_system;
    condition = ref_condition;
}

template<typename type>
void check_aggregation_onlattice<type>::display_compute_times()
{
    std::cout<<"time 1 = "<<time_1*1e-9<<std::endl;
    std::cout<<"time 2 = "<<time_2*1e-9<<std::endl;
    std::cout<<"time 3 = "<<time_3*1e-9<<std::endl;
}

template <typename type>
void check_aggregation_onlattice<type>::check_for_aggregation(constituent<type> *c_1){

    //cp_1 = std::chrono::steady_clock::now();
    bool is_checked = false;
    int  particle_id;

    cluster_id = c_1->get_id();

    //cp_2 = std::chrono::steady_clock::now();
    //time_1 += std::chrono::duration_cast<std::chrono::nanoseconds>(cp_2 - cp_1).count();

    for (int i = 0; i < c_1->get_size(); i++){

        //cp_3 = std::chrono::steady_clock::now();

        neighbours  = c_1->get_neighbour_list(i);
        particle_id = c_1->get_element_id(i);

        //cp_4 = std::chrono::steady_clock::now();

        for (int j = 0; j < neighbours.size(); j++) { 

            neighbour_id = neighbours[j];

            if (neighbour_id != -1){

                neighbour_cluster_id = sys_state->get_id_map(neighbour_id);
                

                /*if (!temp){
                    std::cout<<"NULL pointer encountered"<<std::endl;
                    exit(EXIT_FAILURE);
                }*/


                if (neighbour_cluster_id != cluster_id){

                    //particle_1 = sys_state->get_particle_by_id(particle_id);
                    //particle_2 = sys_state->get_particle_by_id(neighbour_id);

                    particle_1 = sys_state->get_particle_by_index(particle_id);
                    particle_2 = sys_state->get_particle_by_index(neighbour_id);

                    if (condition->agg_condition(particle_1, particle_2)){
                        temp = sys_state->get_aggregate(neighbour_cluster_id);
                        sys_state->add_attachment(c_1);
                        check_for_aggregation(bind_sys->bind_aggregates(c_1, temp));
                        is_checked = true;
                    }
                }


            }

            if (is_checked)
                break;

            
        }

        //cp_5 = std::chrono::steady_clock::now();
        
        //time_2 += std::chrono::duration_cast<std::chrono::nanoseconds>(cp_4 - cp_3).count();
        //time_3 += std::chrono::duration_cast<std::chrono::nanoseconds>(cp_5 - cp_4).count();

        if (is_checked)
            break;

    }


}

template class check_aggregation_onlattice<int>;
template class check_aggregation_onlattice<double>;

}