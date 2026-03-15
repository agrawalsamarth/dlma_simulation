#include <normal_bind.hh>

namespace simulation{

template<typename type>
normal_bind<type>::normal_bind(system<type> *system_ptr){

    sys_state = system_ptr;

}

template<typename type>
constituent<type>* normal_bind<type>::bind_aggregates(constituent<type> *c_1, constituent<type> *c_2){

    //std::cout<<"1"<<std::endl;

    std::string name_type = "cluster";
    constituent<type> *temp;

    temp = factory.create_constituent(sys_state->get_latest_cluster_id(), sys_state->get_lattice(), sys_state->get_dim(), name_type, sys_state->get_box());

    for (int i = 0; i < c_1->get_size(); i++){
        temp->add_constituent(c_1->get_element(i));
    }

    for (int i = 0; i < c_2->get_size(); i++){
        temp->add_constituent(c_2->get_element(i));
    }

    temp->calculate_aggregate_mass();
    temp->set_current_seed_status(1);

    sys_state->remove_aggregate(c_1->get_id());
    sys_state->remove_aggregate(c_2->get_id());
    delete c_1;
    delete c_2;
    sys_state->add_aggregate(temp);
    sys_state->build_id_map();
    sys_state->build_idx_map_for_agg();
    sys_state->calculate_propensity();

    return temp;

}


template class normal_bind<int>;
template class normal_bind<double>;

}