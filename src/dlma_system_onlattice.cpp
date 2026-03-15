#include <dlma_system_onlattice.hh>

namespace simulation{

template <typename type>
dlma_system_onlattice<type>::dlma_system_onlattice(char *params_name)
{
    this->read_params_parser(params_name);
    build_L_arrs();

    if (this->system_type == "dlma"){
        build_site_vector();
        initialize_system();
    }

    else if (this->system_type == "random_site_percolation"){
        std::srand(this->rng_seed);
        build_site_vector_for_rsp();
        initialize_system_for_percolation();
    }

    else{
        std::cout<<"system not defined"<<std::endl;
        exit(EXIT_FAILURE);
    }

}

template <typename type>
void dlma_system_onlattice<type>::build_L_arrs()
{

    L_eff    = (int*)malloc(sizeof(int) * this->D);
    temp_arr = (type*)malloc(sizeof(type) * this->D);
    axis_arr = (int*)malloc(sizeof(int) * this->D);
    neighbour_pos = (type*)malloc(sizeof(type) * this->D);

    for (int axis = 0; axis < this->D; axis++)
        L_eff[axis] = 1;

    for (int axis = 0; axis < this->D; axis++)
    {
        for (int itr = axis+1; itr < this->D; itr++)
            L_eff[axis] *= this->L[itr];
    }

    //for (int axis = 0; axis < this->D; axis++)
        //std::cout<<"L_eff = "<<L[axis]<<std::endl;

    L_total = 1;

    for (int axis = 0; axis < this->D; axis++)
        L_total *= this->L[axis];

}

template <typename type>
void dlma_system_onlattice<type>::initialize_system()
{

    bool is_placed;
    constituent<type> *temp;

    std::string name_type = "particle";

    type temp_pos[this->D];

    for (int i = 0; i < this->N; i++){

        temp = this->factory.create_constituent(i, this->lattice, this->D, name_type, this->box);

        is_placed = false;

        temp->set_diameter(1.);

        if (i < this->N_s){
            temp->set_mass(this->seed_mass);
            temp->set_original_seed_status(1);
            temp->set_current_seed_status(1);
        }

        else{
            temp->set_mass(1.);
            temp->set_original_seed_status(0);
            temp->set_current_seed_status(0);
        }

        this->all_particles.push_back(temp);

        /*while (is_placed == false){

            for (int axis = 0; axis < this->D; axis++)
                temp_pos[axis] = (int)(this->dis(this->generator) * this-> L[axis]);

            if (this->box->get_particle_id(temp_pos) == -1){

                for (int axis = 0; axis < this->D; axis++)
                    temp->pos(axis) = temp_pos[axis];

                temp->add_constituent_to_cell();
                this->all_particles.push_back(temp);
                is_placed = true;
            }

        }*/

    }

    int  counter_pc = 0;
    int  random_index;
    int  pos_index;
    //int* temp_pos_arr;

    while (counter_pc < this->N) {

        random_index = available_sites.size() * this->dis(this->generator);
        pos_index    = available_sites[random_index];
        temp_pos_arr = get_array_from_index(pos_index);

        for (int axis = 0; axis < this->D; axis++)
            this->all_particles[counter_pc]->pos(axis) = temp_pos_arr[axis];

        this->all_particles[counter_pc]->add_constituent_to_cell();

        remove_site_from_vector(pos_index);

        counter_pc++;
    }

    name_type = "cluster";

    for (int i = 0; i < this->N; i++){

        temp = this->factory.create_constituent(this->get_latest_cluster_id(), this->lattice, this->D, name_type, this->box);
        temp->add_constituent(this->all_particles[i]);
        temp->calculate_aggregate_mass();
        //std::cout<<"seed mass = "<<get_seedmass()<<std::endl;

        this->aggregates.push_back(temp);

    }

    this->build_id_map();
    this->calculate_propensity();
    this->build_idx_map_for_agg();
    
}

template<typename type>
void dlma_system_onlattice<type>::initialize_system_for_percolation()
{

    constituent<type> *temp;
    std::string name_type = "particle";

    int div, rem;

    this->N = L_total * this->phi;

    int random_index;
    int pos_index;

    //std::cout<<"1 init"<<std::endl;

    for (int i = 0; i < this->N; i++){

        //if (this->dis(this->generator) < this->phi){

            temp = this->factory.create_constituent(i, this->lattice, this->D, name_type, this->box);
            //this->N += 1;
            temp->set_diameter(1.);
            temp->set_mass(this->seed_mass);
            temp->set_original_seed_status(1);
            temp->set_current_seed_status(1);

            //random_index = available_sites.size() * this->dis(this->generator);
            //std::cout<<"size = "<<available_sites.size()<<"\t random_index = "<<random_index<<std::endl;
            pos_index    = available_sites[i];
            temp_pos_arr = get_array_from_index(pos_index);
            //remove_site_from_vector(pos_index);

            for (int axis = 0; axis < this->D; axis++)
                temp->pos(axis) = temp_pos_arr[axis];

            temp->add_constituent_to_cell();
            this->all_particles.push_back(temp);

        //}

    }

    //std::cout<<"size of = "<<sizeof(&temp)<<std::endl;

    //std::cout<<"2 init"<<std::endl;

    name_type = "cluster";

    for (int i = 0; i < this->N; i++){

        temp = this->factory.create_constituent(this->get_latest_cluster_id(), this->lattice, this->D, name_type, this->box);
        temp->add_constituent(this->all_particles[i]);
        temp->calculate_aggregate_mass();
        //std::cout<<"seed mass = "<<get_seedmass()<<std::endl;

        this->aggregates.push_back(temp);

    }

    //std::cout<<"3 init"<<std::endl;

    /*for (int i = 0; i < this->aggregates.size(); i++){
        this->aggregates[i]->add_agg_to_cell();
    }*/

    this->attachments.resize(this->N);
    this->build_id_map();

    for (int i = 0; i < this->N; i++)
        this->add_attachment(this->aggregates[i]);

    //std::cout<<"4 init"<<std::endl;
    

}

template <typename type>
void dlma_system_onlattice<type>::build_site_vector()
{

    available_sites.clear();
    std::vector<int> temp_neighbours;
    int neighbour_flag;

    //for (int i = 0; i < L_total; i++)
        //available_sites.push_back(i);

    //constituent<type> *temp;
    //temp = this->factory.create_constituent(0, this->lattice, this->D, name_type, this->box);

    for (int i = 0; i < L_total; i++){

        temp_pos_arr    = get_array_from_index(i);
        temp_neighbours = this->box->get_neighbour_list(temp_pos_arr);

        neighbour_flag = 0;

        //if (temp_neighbours.size() == 0){

        for (int j = 0; j < temp_neighbours.size(); j++){

            if (temp_neighbours[j] != -1){
                neighbour_flag = 1;
                break;
            }
        }
        //}

        //if (neighbour_flag == 0){
            available_sites.push_back(i);
            this->box->add_particle_to_cell(i, temp_pos_arr);
        //}
                
    }

    this->box->clear_cell_field();

}

template <typename type>
void dlma_system_onlattice<type>::build_site_vector_for_rsp()
{
    available_sites.clear();

    for (int i = 0; i < L_total; i++)
        available_sites.push_back(i);

    std::random_shuffle(available_sites.begin(), available_sites.end());

}

template <typename type>
void dlma_system_onlattice<type>::remove_site_from_vector(int index)
{
    available_sites.erase(std::remove(available_sites.begin(), available_sites.end(), index), available_sites.end());
}

template <typename type>
void dlma_system_onlattice<type>::exclude_neighbours(int random_index)
{

    temp_pos_arr = get_array_from_index(random_index);

    for (int axis = 0; axis < this->D; axis++){
        for (int i = 0; i < 2; i++){

            for (int axis_2 = 0; axis_2 < this->D; axis_2++)
                axis_arr[axis_2] = (axis_2 == axis) * (-1 + 2*i);

            for (int axis_2 = 0; axis_2 < this->D; axis_2++){
                neighbour_pos[axis_2]  = (temp_pos_arr[axis_2] + axis_arr[axis_2]);
                neighbour_pos[axis_2] += (this->L[axis_2] * (neighbour_pos[axis_2]==-1)) - (this->L[axis_2] * (neighbour_pos[axis_2]==this->L[axis_2]));
            }

            neighbour_index = get_index_from_array(neighbour_pos);
            remove_site_from_vector(neighbour_index);

        }
    }



}

template <typename type>
type* dlma_system_onlattice<type>::get_array_from_index(int index)
{
    for (int axis = 0; axis < this->D; axis++){
        div = index/L_eff[axis];
        temp_arr[axis] = div;
        index = index % L_eff[axis];
    }

    return temp_arr;
}

template <typename type>
int dlma_system_onlattice<type>::get_index_from_array(type *arr)
{
    counter = 0;

    for (int axis = 0; axis < this->D; axis++)
            counter += arr[axis] * L_eff[axis];

    return counter;
}


template <typename type>
bool dlma_system_onlattice<type>::check_viability(constituent<type> *c_1, type *dr)
{
    type temp_pos[this->D];

    for (int axis = 0; axis < this->D; axis++)
        temp_pos[axis] = 0;

    int cluster_id = c_1->get_id();
    int neighbour_cluster_id;
    int neighbour_id;

    constituent<type> *temp;


    this->is_viable = true;

    for (int i = 0; i < c_1->get_size(); i++){

        temp = c_1->get_element(i);

        for (int axis = 0; axis < this->D; axis++){
            temp_pos[axis] = this->box->get_refill(temp->pos(axis)+dr[axis], axis);
        }

        neighbour_id = this->box->get_particle_id(temp_pos);

        if (neighbour_id != -1){
            neighbour_cluster_id = this->get_id_map(neighbour_id);

            if (neighbour_cluster_id != cluster_id)
                this->is_viable=false;
        }

    }

    return this->is_viable;

}

template <typename type>
void dlma_system_onlattice<type>::move_aggregate(int i, type *dr)
{

    if (check_viability(this->aggregates[i], dr)){
        this->aggregates[i]->remove_constituent_from_cell();
        this->aggregates[i]->move(dr);
        this->aggregates[i]->add_constituent_to_cell();
    }


}

template <typename type>
void dlma_system_onlattice<type>::print_grid()
{
    type temp_pos[this->D];
    int temp_id;

    for (int i = 0; i < this->L[0]; i++){
        for (int j = 0; j < this->L[1]; j++){

            temp_pos[0] = i;
            temp_pos[1] = j;

            temp_id = this->box->get_particle_id(temp_pos);

            std::cout<<temp_id<<"\t";

        }
        std::cout<<"\n"<<std::endl;
    }

}

template<typename type>
dlma_system_onlattice<type>::~dlma_system_onlattice()
{
    free(L_eff);
    free(temp_arr);
    free(axis_arr);
    free(neighbour_pos);
    free(this->L_flag);
    free(this->periodic_flag);
    free(this->L);
    free(this->halfL);
    
    delete this->box;

    for (int i = 0; i < this->N; i++)
        delete this->all_particles[i];

    for (int i = 0; i < this->aggregates.size(); i++)
        delete this->aggregates[i];

}

template class dlma_system_onlattice<int>;
template class dlma_system_onlattice<double>;

}