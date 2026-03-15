#include <dlma_system_offlattice.hh>

namespace simulation{

template <typename type>
dlma_system_offlattice<type>::dlma_system_offlattice(char *params_name)
{ 
    //std::cout<<"offlattice 1"<<std::endl;
    this->read_params_parser(params_name);
    //std::cout<<"offlattice 2"<<std::endl;

    if (this->system_type == "dlma"){
        initialize_system();
        int N_pair = (this->N * (this->N - 1))/2; 
        rij = (double*)malloc(sizeof(double)*N_pair);
    }

    else if (this->system_type == "erdos_renyi"){
        init_erdos_renyi();
    }

    else {
        std::cout<<"please check system type"<<std::endl;
    }


}

template <typename type>
void dlma_system_offlattice<type>::initialize_system()
{
    bool is_placed;
    constituent<type> *temp;

    std::string name_type = "particle";

    type distance;

    pos_arr_for_particle = (int*)malloc(sizeof(int)*this->D);
    delta_L_temp         = (double*)malloc(sizeof(double)*this->D);
    L_eff_temp           = (int*)malloc(sizeof(int)*this->D);
    pos_pop_temp         = (int*)malloc(sizeof(int)*this->D);
    neigh_pop_temp       = (int*)malloc(sizeof(int)*this->D);
    pop_arr              = (int*)malloc(sizeof(int)*this->D);

    build_site_vector();

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

        /*while (is_placed == false){

            for (int axis = 0; axis < this->D; axis++)
                temp->pos(axis) = (this->dis(this->generator) * this-> L[axis]);

                is_placed = true;
                
                for (int j = 0; j < i; j++){
            
                    distance = this->get_interparticle_distance(temp, this->all_particles[j]);

                    if (distance < ((1.+this->tolerance) * 0.5 * (temp->get_diameter() + this->all_particles[j]->get_diameter()))){
                        is_placed = false;
                        break;
                    }

                }            

        }*/

        get_array_from_index_map(available_sites_temp[i], pos_arr_for_particle);

        for (int axis = 0; axis < this->D; axis++)
            temp->pos(axis) = (1.*pos_arr_for_particle[axis]*delta_L_temp[axis]) + (0.5*delta_L_temp[axis]);

        temp->add_constituent_to_cell();
        this->all_particles.push_back(temp);

    }

    name_type = "cluster";

    for (int i = 0; i < this->N; i++){

        temp = this->factory.create_constituent(this->get_latest_cluster_id(), this->lattice, this->D, name_type, this->box);
        temp->add_constituent(this->all_particles[i]);
        temp->calculate_aggregate_mass();
        this->aggregates.push_back(temp);

    }

    this->build_id_map();

    name_type = "particle";
    image = this->factory.create_constituent(this->N, this->lattice, this->D, name_type, this->box);

    this->calculate_propensity();
    this->build_idx_map_for_agg();

    
}

template <typename type>
void dlma_system_offlattice<type>::init_erdos_renyi()
{
    constituent<type> *temp;
    std::string name_type = "particle";

    type distance;
    bool is_placed;

    double temp_r;
    double temp_phi;

    for (int i = 0; i < this->N; i++){

        temp = this->factory.create_constituent(i, this->lattice, this->D, name_type, this->box);
        temp->set_diameter(1.);    
        temp->set_mass(1.);
        temp->set_original_seed_status(1);
        temp->set_current_seed_status(1);

        is_placed = false;
    
        while (!is_placed){

            //for (int axis = 0; axis < this->D; axis++)
            temp_r   = this->dis(this->generator);
            temp_phi = 2. * M_PI *  this->dis(this->generator);
            /*temp->pos(0)  = this->L[0] * temp_r * sin(temp_phi) * 0.5;
            temp->pos(1)  = this->L[1] * temp_r * cos(temp_phi) * 0.5;
            temp->pos(0) += 0.5 * this->L[0];
            temp->pos(1) += 0.5 * this->L[1];*/

            for (int axis = 0; axis < this->D; axis++)
                temp->pos(axis) = this->L[axis] * this->dis(this->generator);


            //std::cout<<temp->pos(0)<<"\t"<<temp->pos(1)<<std::endl;

            is_placed = true;
            
            /*for (int j = 0; j < i; j++){
        
                distance = this->get_interparticle_distance(temp, this->all_particles[j]);

                if (distance < ((1.+this->tolerance) * 0.5 * (temp->get_diameter() + this->all_particles[j]->get_diameter()))){
                    is_placed = false;
                    break;
                }

            }*/           

        }

        //temp->add_constituent_to_cell();
        this->all_particles.push_back(temp);

    }

    //for (int i = 0; i < this->N; i++)
        //std::cout<<this->all_particles[i]->pos(0)<<"\t"<<this->all_particles[i]->pos(1)<<std::endl;

    //std::cout<<"size = "<<this->attachments.size()<<std::endl;
    this->attachments.resize(this->N);

    double rcg_distance = 0.;

    if (this->distance_metric_rgg == 1){

        for (int i = 0; i < this->N; i++){
            for (int j = i+1; j < this->N; j++){

                
                    rcg_distance = this->get_interparticle_distance(this->all_particles[i], this->all_particles[j]);

                    if (rcg_distance < this->phi)
                        this->add_attachment(i, j);
                

            }
        }

    }

    else {

        for (int i = 0; i < this->N; i++){
            for (int j = i+1; j < this->N; j++){

                    rcg_distance = this->get_manhattan_distance(this->all_particles[i], this->all_particles[j]);
                    //std::cout<<"rcg_distance = "<<rcg_distance<<std::endl;
                
                    if (rcg_distance < this->phi)
                        this->add_attachment(i, j);
                
            }
        }

    }

}

/*template <typename type>
void dlma_system_offlattice<type>::initialize_system()
{

    bool is_placed;
    constituent<type> *temp;

    std::string name_type = "particle";

    type distance;

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

        //while (is_placed == false){

            for (int axis = 0; axis < this->D; axis++)
                temp->pos(axis) = 1. * i;

        //}

        temp->add_constituent_to_cell();
        this->all_particles.push_back(temp);

    }

    temp = this->get_particle_by_id(49);
    temp->pos(0) = 5.;

    name_type = "cluster";

    for (int i = 0; i < this->N; i++){

        temp = this->factory.create_constituent(this->get_latest_cluster_id(), this->lattice, this->D, name_type, this->box);
        temp->add_constituent(this->all_particles[i]);
        temp->calculate_aggregate_mass();
        this->aggregates.push_back(temp);

    }

    this->build_id_map();

    name_type = "particle";
    image = this->factory.create_constituent(this->N, this->lattice, this->D, name_type, this->box);

    this->calculate_propensity();
    
}*/

/*template <typename type>
void dlma_system_offlattice<type>::move_aggregate(int i, type *dr)
{

    type scale = fix_overlap(i, dr);

    for (int axis = 0; axis < this->D; axis++)
        dr[axis] = scale * dr[axis];

    this->aggregates[i]->remove_constituent_from_cell();
    this->aggregates[i]->move(dr);
    this->aggregates[i]->add_constituent_to_cell();

    //std::cout<<"i="<<this->aggregates[i]->get_id()<<std::endl;
    //calc_rij();

}*/

template <typename type>
void dlma_system_offlattice<type>::move_aggregate(int i, type *dr)
{

    type scale = fix_overlap(i, dr);

    for (int axis = 0; axis < this->D; axis++)
        dr[axis] = scale * dr[axis];

    this->aggregates[i]->remove_constituent_from_cell();
    this->aggregates[i]->move(dr);
    this->aggregates[i]->add_constituent_to_cell();


}

template <typename type>
void dlma_system_offlattice<type>::build_site_vector()
{
    double  L_total_temp;
    L_disc   = (int)(this->L[0]/((1.+this->tolerance)));

    if (L_disc%2 == 1){
        L_disc -= 1;
    }

    L_total_temp = 1;

    for (int axis = 0; axis < this->D; axis++)
        L_total_temp *= L_disc;

    //int    *nx_temp;
    nx_temp = (int*)malloc(sizeof(int)*L_total_temp);

    for (int i = 0; i < L_total_temp; i++)
        nx_temp[i] = -1;

    for (int axis = 0; axis < this->D; axis++){
        delta_L_temp[axis] = (this->L[axis]/(1.*L_disc));
    }  

    for (int axis = 0; axis < this->D; axis++){
        L_eff_temp[axis] = 1;
    }

    for (int axis = 0; axis < this->D; axis++)
    {
        for (int itr = axis+1; itr < this->D; itr++)
            L_eff_temp[itr] *= L_disc;
    }

    available_sites_temp.clear();

    //int temp_counter;
    int neighbour_counter_temp;
    bool adj_flag;

    for (int i = 0; i < L_total_temp; i++){

        get_array_from_index_map(i, pos_pop_temp);
        adj_flag = true;

        for (int axis = 0; axis < this->D; axis++){
            for (int j = 0; j < 2; j++){

                neighbour_counter_temp = i + (L_eff_temp[axis] * (2*j - 1));
                get_array_from_index_map(neighbour_counter_temp, neigh_pop_temp);

                for (int axis_2 = 0; axis_2 < this->D; axis_2++){
                    neigh_pop_temp[axis_2] += L_disc*(neigh_pop_temp[axis_2]<0) - L_disc*(neigh_pop_temp[axis_2] >= L_disc); 
                }

                neighbour_counter_temp = get_index_from_array_map(neigh_pop_temp);
                
                if (nx_temp[neighbour_counter_temp] != -1)
                    adj_flag = false;

            }
        }

        if (adj_flag)
            available_sites_temp.push_back(i);

    }

    if (this->N > available_sites_temp.size()){
        std::cout<<"phi,N,L is not appropriate for initialization"<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::random_shuffle(available_sites_temp.begin(), available_sites_temp.end());
   
}

template<typename type>
void dlma_system_offlattice<type>::get_array_from_index_map(int pop_i, int* for_pop_arr){

    for (int axis = 0; axis < this->D; axis++){
        pop_div           = pop_i/L_eff_temp[axis];
        for_pop_arr[axis] = pop_div%L_disc; 
    }

}

template<typename type>
int dlma_system_offlattice<type>::get_index_from_array_map(int* temp_arr_for_pop){

    pop_counter = 0;

    for (int axis = 0; axis < this->D; axis++){
        pop_counter += temp_arr_for_pop[axis] * L_eff_temp[axis];
    }

    return pop_counter;

}


/*template<typename type>
type dlma_system_offlattice<type>::fix_overlap(int i, type *dr)
{
    hs = build_collision_list(i, 0., dr);

    double t_c = 1.;
    double fraction = 1.;

    bool hs_coll = false;

    for (int index = 0; index < hs.size(); index++) {

        if (hs[index].min_time < fraction) {
            t_c = hs[index].min_time;
            hs_coll = true;
        }

    }

    bonds = build_collision_list(i, this->tolerance, dr);

    double start_zone = 1.;
    double end_zone   = 1e8;

    int index_i;
    int index_j;

    for (int index = 0; index < bonds.size(); index++) {
        
        if (bonds[index].min_time < start_zone) {
            start_zone = bonds[index].min_time;
            index_i = bonds[index].i;
            index_j = bonds[index].j;
        }


        if (bonds[index].max_time < end_zone)
            end_zone = bonds[index].max_time;

    }

    bool tbc = false;
    double tbc_min = 1.;
    double tbc_max = t_c;

    double tbc_time_start;
    double tbc_time_end;

    for (int index = 0; index < bonds.size(); index++) {

        if ((bonds[index].i == index_i) && (bonds[index].j != index_j))
        {

            tbc_time_start = bonds[index].min_time;
            tbc_time_end   = bonds[index].max_time;

            tbc = true;

            if (tbc_time_start < end_zone) {

                

            }



        }


    }
    

    return fraction;

}

template<typename type>
type dlma_system_offlattice<type>::fix_overlap(int i, type *dr)
{
    hs = build_collision_list(i, 0., dr);

    double t_c = 1.;
    double fraction = 1.;

    bool hs_coll = false;

    int hs_index_i;
    int hs_index_j;

    for (int index = 0; index < hs.size(); index++) {

        if (hs[index].min_time < t_c) {
            t_c = hs[index].min_time;
            hs_index_i = hs[index].i;
            hs_index_j = hs[index].j;

            hs_coll = true;
        }

    }

    bonds = build_collision_list(i, this->tolerance, dr);

    double start_zone = 1.;
    double end_zone   = 1e8;

    int index_i;
    int index_j;

    bool coll = false;

    for (int index = 0; index < bonds.size(); index++) {
        
        if (bonds[index].min_time < start_zone) {
            start_zone = bonds[index].min_time;
            end_zone   = bonds[index].max_time;
            index_i    = bonds[index].i;
            index_j    = bonds[index].j;

            coll = true;
        }

    }

    bool tbc = false;
    double tbc_min = 1.;
    double tbc_max = t_c;

    double tbc_time_start;
    double tbc_time_end;

    for (int index = 0; index < bonds.size(); index++) {

        if ((bonds[index].i == index_i) && (bonds[index].j != index_j) && (bonds[index].min_time < end_zone)) {

            tbc = true;
            
            if (bonds[index].min_time < tbc_min)
                tbc_min = bonds[index].min_time;

            if (bonds[index].max_time < tbc_max)
                tbc_max = bonds[index].max_time;

        }


    }

    if (coll){
    
        if (tbc == true) {

            if (tbc_max < end_zone)
                end_zone = tbc_max;

            fraction = start_zone + (end_zone - start_zone) * this->dis(this->generator);

            if (fraction > t_c)
                fraction = t_c;

            //std::cout<<"test"<<std::endl;

        }

        else {

            if ((hs_coll == true) && (hs_index_i == index_i) && (hs_index_j = index_j))
                fraction = t_c;
                
            else {
                
                fraction = 0.5 * (start_zone + end_zone);
                
                if (fraction > t_c)
                    fraction = t_c;

                //std::cout<<"test"<<std::endl;


            }
    

        }

    }


    return fraction;

}*/

template<typename type>
type dlma_system_offlattice<type>::fix_overlap(int i, type *dr)
{
    hs = build_collision_list(i, 0., dr);

    double t_c = 1.;
    double fraction = 1.;

    bool hs_coll = false;

    int hs_index_i;
    int hs_index_j;

    for (int index = 0; index < hs.size(); index++) {

        if ((fabs(hs[index].min_time) < t_c) && (hs[index].min_time >= 0.)) {
            t_c = hs[index].min_time;
            hs_index_i = hs[index].i;
            hs_index_j = hs[index].j;

            hs_coll = true;
        }

    }

    bonds = build_collision_list(i, this->tolerance, dr);

    /*for (auto x:bonds){
        std::cout<<"i = "<<i<<" j = "<<x.j<<" min_time = "<<x.min_time<<" max_time = "<<x.max_time<<"\n";
    }

    std::cout<<"tc = "<<t_c<<std::endl;*/

    //int index_first_i;
    //int index_second_i;
    //int index_first_j;
    //int index_second_j;

    double min_i;
    double max_i;

    double min_j;
    double max_j;

    double start_zone = t_c;
    double end_zone = t_c;

    bool tbc = false;

    for (int i = 0 ; i < bonds.size(); i++){

        min_i = bonds[i].min_time;
        max_i = bonds[i].max_time;

        if ((min_i > 0.) && (min_i < t_c)){
        
            for (int j = i+1; j < bonds.size(); j++){

                min_j = bonds[j].min_time;
                max_j = bonds[j].max_time;

                if ((min_j > 0.)){

                    if ((min_i > min_j) && (min_i < max_j)){
                        tbc = true;

                        if (min_j < start_zone){
                            start_zone = min_j;

                            if (max_i < max_j)
                                end_zone   = max_i;

                            else
                                end_zone   = max_j;
                        }

                    }


                    if ((min_j > min_i) && (min_j < max_i)){
                        tbc = true;

                        if (min_i < start_zone){
                            start_zone = min_i;

                            if (max_i < max_j)
                                end_zone   = max_i;

                            else
                                end_zone   = max_j;
                        }

                    }                
                }
            }

        }

    }

    if (start_zone > t_c){
        tbc = false;
    }

    if (end_zone > t_c){
        end_zone = t_c;
    }

    //std::cout<<"tbc = "<<tbc<<" hs = "<<hs_coll<<std::endl;

    if (tbc){
        fraction = start_zone + (this->dis(this->generator) * (end_zone-start_zone));
        //std::cout<<fraction<<std::endl;
        //std::cout<<"t_c = "<<t_c<<" start = "<<start_zone<<" end_zone = "<<end_zone<<" fraction = "<<fraction<<std::endl;
        return fraction;
    }

    if (hs_coll){
        return t_c;
    }
        
    return 1.;

}


template<typename type>
std::vector<coll_deets> dlma_system_offlattice<type>::build_collision_list(int i, double alpha, type *dr)
{

    type diff[this->D];
    type beta;
    type q;
    type dia_distance;
    type diff_2;

    std::vector<coll_deets> list_coll;

    constituent<type> *c_1 = this->aggregates[i];
    constituent<type> *ref_particle;
    constituent<type> *nb_particle;
    type distance;
    type fraction = 1.;

    int particle_id;
    int cluster_id = c_1->get_id();
    int neighbour_id;
    int neighbour_cluster_id;

    coll_deets t_c;

    for (int i = 0; i < c_1->get_size(); i++){

        ref_particle = c_1->get_element(i);

        for (int axis = 0; axis < this->D; axis++)
            image->pos(axis) = ref_particle->pos(axis);

        image->move(dr);

        actual_list = ref_particle->get_neighbour_list();
        image_list  = image->get_neighbour_list();
        
        neighbours.clear();
        neighbours.reserve(actual_list.size() + image_list.size());
        neighbours.insert( neighbours.end(), actual_list.begin(), actual_list.end() );
        neighbours.insert( neighbours.end(), image_list.begin(), image_list.end() );
           
        for (int j = 0; j < neighbours.size(); j++) { 

            neighbour_id = neighbours[j];
            neighbour_cluster_id = this->get_id_map(neighbour_id);

            if (neighbour_cluster_id != cluster_id){

                nb_particle = this->get_particle_by_id(neighbour_id);
                t_c = calc_quad_eqn(dr, alpha, ref_particle, nb_particle);
                list_coll.push_back(t_c);

            }

        }

    }

    return list_coll;


}

template<typename type>
coll_deets dlma_system_offlattice<type>::calc_quad_eqn(type *dr, double alpha, constituent<type> *ref_particle, constituent<type> *nb_particle)
{

    //nb_particle = this->get_particle_by_id(neighbour_id);

    coll_deets temp_coll;
    temp_coll.i = ref_particle->get_id();
    temp_coll.j = nb_particle->get_id();

    type diff[this->D];

    for (int axis = 0; axis < this->D; axis++)
        diff[axis] = this->box->get_periodic_distance(nb_particle->pos(axis), ref_particle->pos(axis), axis);

    double q = 0.;

    for (int axis = 0; axis < this->D; axis++)
        q += diff[axis] * dr[axis];

    double d_mag_2 = 0.;

    for (int axis = 0; axis < this->D; axis++)
        d_mag_2 += (dr[axis] * dr[axis]);

    double d_mag = sqrt(d_mag_2);

    //std::cout<<"ref = "<<ref_particle->get_id()<<" nb = "<<nb_particle->get_id()<<" d_mag = "<<d_mag<<" d_mag_2 = "<<d_mag_2<<std::endl;

    double diff_2 = 0.;

    for (int axis = 0; axis < this->D; axis++)
        diff_2 += diff[axis] * diff[axis];

    double dia_distance = 0.;

    dia_distance = (1.+ alpha)* 0.5 * (ref_particle->get_diameter() + nb_particle->get_diameter());
    dia_distance = dia_distance * dia_distance;

    double beta = 0.;
    double t_c[2];

    beta = (q * q) - (d_mag_2*(diff_2 - dia_distance));

    if (beta >= 0){
        t_c[0] = (q - sqrt(beta))/d_mag_2;
        t_c[1] = (q + sqrt(beta))/d_mag_2;
    }

    else{
        t_c[0] = -2.;
        t_c[1] = -2.;
    }

    //if (nb_particle->get_id() == 0){
        //std::cout<<"q = "<<q<<" beta = "<<beta<<std::endl;
    //}

    if (t_c[0] >= 0.){

        if (fabs(t_c[0]) < 1e-4){            
            temp_coll.min_time = 0.;
            temp_coll.max_time = 0.;
        }

        else{
            temp_coll.min_time = t_c[0];
            temp_coll.max_time = t_c[1];
        }

    }

    else {
        temp_coll.min_time = -2.;
        temp_coll.max_time = -2.;
    }

    return temp_coll;

}

template<typename type>
void dlma_system_offlattice<type>::add_attachment(const int i, const int j)
{
    this->attachments[i].push_back(j);
    this->attachments[j].push_back(i);
}

template<typename type>
void dlma_system_offlattice<type>::calc_rij()
{

    int count = 0;

    constituent<type> *p_1;
    constituent<type> *p_2;

    for (int i = 0; i < this->N; i++){
        p_1 = this->all_particles[i];
        for (int j = i+1; j < this->N; j++){
            p_2 = this->all_particles[j];
            rij[count] = this->get_interparticle_distance(p_1, p_2);

            if (rij[count] < (1.-1e-4)){
                std::cout<<"particle i = "<<p_1->get_aggregate_id()<<"\t particle j = "<<p_2->get_aggregate_id()<<"\t distance="<<rij[count]<<std::endl;
                exit(EXIT_FAILURE);
            }

            count++;

        }
    }



}

template<typename type>
void dlma_system_offlattice<type>::build_attachment_list()
{

    constituent<type>* particle_1;
    constituent<type>* particle_2;
    double distance;

    std::vector<int> neighbours;

    int particle_id;
    int neighbour_id;

    for (int i = 0; i < this->N; i++){

        particle_1 = this->get_particle_by_id(i);
        neighbours = particle_1->get_neighbour_list();

        particle_id = i;

        for (int j = 0; j < neighbours.size(); j++){

            neighbour_id = neighbours[j];

            if (neighbour_id != particle_id){

                particle_2 = this->get_particle_by_id(neighbour_id);
                distance   = this->get_interparticle_distance(particle_1, particle_2);

                if ((distance-1e-4) < ((1. + this->tolerance) * 0.5 * (particle_1->get_diameter() + particle_2->get_diameter()))){

                    if(!(std::find(this->attachments[particle_id].begin(), this->attachments[particle_id].end(), neighbour_id) != this->attachments[particle_id].end())){

                        this->attachments[particle_id].push_back(neighbour_id);
                        this->attachments[neighbour_id].push_back(particle_id);

                    } 

                }

            }


        }

    }


}






template<typename type>
dlma_system_offlattice<type>::~dlma_system_offlattice()
{
    free(this->L_flag);
    free(this->periodic_flag);
    free(this->L);
    free(rij);
    free(this->halfL);

    delete this->box;

    for (int i = 0; i < this->N; i++)
        delete this->all_particles[i];

    for (int i = 0; i < this->aggregates.size(); i++)
        delete this->aggregates[i];

    delete image;

    free(pos_arr_for_particle);
    free(delta_L_temp);
    free(L_eff_temp);
    free(pos_pop_temp);
    free(neigh_pop_temp);
    free(pop_arr); 
    free(nx_temp);   

}

template class dlma_system_offlattice<int>;
template class dlma_system_offlattice<double>;

}