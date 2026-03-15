#include <dlma_iterator.hh>

namespace simulation{

template <typename type>
std::vector<std::string> dlma_iterator<type>::split_string_by_delimiter(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

template <typename type>
dlma_iterator<type>::dlma_iterator(char *filename)
{
    std::ifstream parser(filename, std::ifstream::in);

    if (parser.fail()){
        std::cout<<"either file does not exist or does not have permissions"<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::string bind_name;
    bool bind_flag=false;
    std::string agg_condition_name;
    bool agg_condition_flag=false;
    std::string check_agg_name;
    bool check_agg_flag=false;
    std::string movement_name;
    bool movement_flag=false;
    std::string save_config_name;
    bool save_config_flag=false;
    int  rng_seed;
    bool rng_seed_flag=false;
    int  lattice;
    bool lattice_flag=false;
    double tolerance;
    bool   tolerance_flag=false;

    int count = 0;

    std::string str;
    std::vector<std::string> results;

    while (getline(parser, str)){

        results = split_string_by_delimiter(str, '=');

        if (results[0] == "bind"){
            bind_name      = results[1];
            bind_flag      = true;
        }

        if (results[0] == "aggregation_condition"){
            agg_condition_name      = results[1];
            agg_condition_flag      = true;
        }

        if (results[0] == "aggregation_type"){
            check_agg_name      = results[1];
            check_agg_flag      = true;
        }


        if (results[0] == "movement"){
            movement_name      = results[1];
            movement_flag      = true;
        }

        if (results[0] == "system"){
            save_config_name      = results[1];
            save_config_flag      = true;
        }

        if (results[0] == "rng_seed"){
            rng_seed      = stoi(results[1]);
            rng_seed++;
            rng_seed_flag = true;
        }

        if (results[0] == "lattice"){
            lattice      = stoi(results[1]);
            lattice_flag = true;
        }

        if (results[0] == "agg_dist_tolerance"){
            tolerance      = stod(results[1]);
            tolerance_flag = true;
        }

        if (results[0] == "final_aggregate_number"){
            final_aggregate_number = stoi(results[1]);
        }        

    }

    parser.close();

    /*if (bind_flag == false){
        std::cout<<"please provide bind"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if (agg_condition_flag == false){
        std::cout<<"please provide aggregation condition"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if (check_agg_flag == false){
        std::cout<<"please provide aggregation type"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if (movement_flag == false){
        std::cout<<"please provide movement"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if (save_config_flag == false){
        std::cout<<"please provide system type"<<std::endl;
        exit(EXIT_FAILURE);
    }

    if (lattice_flag == false){
        std::cout<<"please provide lattice type"<<std::endl;
        exit(EXIT_FAILURE);
    }*/

    if (bind_flag == false){
        bind_name = "normal";
    }

    if (agg_condition_flag == false){
        agg_condition_name = "mass";
    }

    if (check_agg_flag == false){
        check_agg_name = "normal";
    }

    if (movement_flag == false){
        movement_name = "brownian";
    }

    if (save_config_flag == false){
        save_config_name = "dlma";
    }

    if (lattice_flag == false){
        lattice = 1;
    }

    if (rng_seed_flag == false){
        rng_seed=1;
    }

    if (tolerance_flag == false){
        tolerance = 0.;
    }

    //std::cout<<"1"<<std::endl;
    //std::cout<<"creating system"<<std::endl;
    sys_state   = factory->create_new_system(save_config_name, lattice, filename);
    //std::cout<<"system created"<<std::endl;
    //std::cout<<"2"<<std::endl;
    binding_obj = factory->create_bind_system(bind_name, sys_state);
    //std::cout<<"3"<<std::endl;
    agg_condition = factory->create_aggregation_condition(agg_condition_name, sys_state);
    //std::cout<<"4"<<std::endl;
    aggregation_check_obj = factory->create_check_aggregation(check_agg_name, lattice, sys_state, binding_obj, agg_condition, tolerance);
    //std::cout<<"5"<<std::endl;
    movement_test = factory->create_movement(movement_name, sys_state->get_dim(), rng_seed, lattice);
    //std::cout<<"6"<<std::endl;
    save_obj = factory->create_save_config(save_config_name, sys_state, sys_state->get_box());
    //std::cout<<"7"<<std::endl;


    if (save_config_name == "dlma"){

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

        for (int i = 0; i < sys_state->get_latest_cluster_id_without_increment(); i++){
        
            temp_c = sys_state->get_aggregate(i);

            if (temp_c){
                aggregation_check_obj->check_for_aggregation(temp_c);
            }

            //std::cout<<"i = "<<i<<std::endl;

        }

        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        //std::cout << "Time difference for aggregation calc = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    }

    //std::cout<<"total aggregates = "<<sys_state->total_aggregates()<<std::endl;    
}

template <typename type>
void dlma_iterator<type>::iteration_step()
{

    //std::cout<<"1"<<std::endl;
    //cp_1 = std::chrono::steady_clock::now();
    temp   = sys_state->choose_aggregate();
    //cp_2 = std::chrono::steady_clock::now();
    //std::cout<<"2"<<std::endl;
    temp_c = sys_state->get_constituent(temp);
    //cp_3 = std::chrono::steady_clock::now();
    //std::cout<<"3"<<std::endl;
    sys_state->move_aggregate(temp, movement_test->delta_x());
    //cp_4 = std::chrono::steady_clock::now();
    //std::cout<<"4"<<std::endl;
    aggregation_check_obj->check_for_aggregation(temp_c);
    //cp_5 = std::chrono::steady_clock::now();
    //std::cout<<"5"<<std::endl;

    //time_1 += std::chrono::duration_cast<std::chrono::seconds>(cp_2 - cp_1).count();
    //time_2 += std::chrono::duration_cast<std::chrono::seconds>(cp_3 - cp_2).count();
    //time_3 += std::chrono::duration_cast<std::chrono::seconds>(cp_4 - cp_3).count();
    //time_4 += std::chrono::duration_cast<std::chrono::seconds>(cp_5 - cp_4).count();


}

template <typename type>
void dlma_iterator<type>::run_system()
{
    //std::cout<<"1"<<std::endl;

    /*while (sys_state->total_aggregates() > final_aggregate_number){
        //std::cout<<"aggregates="<<sys_state->total_aggregates()<<std::endl;
        iteration_step();
    }*/
    //std::cout<<"aggregates="<<sys_state->total_aggregates()<<std::endl;
    //std::cout<<"here"<<std::endl;
    sys_state->build_attachment_list();
    //std::cout<<"and then here"<<std::endl;

    /*std::cout<<"aggregate choosing time = "<<time_1<<"(s)";
    std::cout<<"get constituent time = "<<time_2<<"(s)";
    std::cout<<"move aggregate time = "<<time_3<<"(s)";
    std::cout<<"check aggregation time = "<<time_4<<"(s)";*/

    //aggregation_check_obj->display_compute_times();
    

    return;
}

template <typename type>
void dlma_iterator<type>::save_config_file()
{
    save_obj->save_configuration();
}


template <typename type>
void dlma_iterator<type>::save_config_file(char *filename)
{
    save_obj->save_configuration(filename);
}

template <typename type>
void dlma_iterator<type>::create_movie_files(char *filename)
{
    int itr = 0;
    char *filename_m;
    int index;
    std::string index_s;
    std::string temp_s;
    std::string ext = ".csv";

    int freq = 25;

    while (sys_state -> total_aggregates() != 1)
    {
        iteration_step();

        if (itr % freq == 0){

            //std::cout<<"here1"<<std::endl;

            index   = itr/freq;
            index_s = std::to_string(index);
            filename_m = (char*)malloc(1 + strlen(filename)+ index_s.length() + ext.length());
            strcpy(filename_m, filename);
            strcat(filename_m, index_s.c_str());
            strcat(filename_m, ext.c_str());

            //std::cout<<"here"<<std::endl;

            save_obj->save_configuration(filename_m);

            free(filename_m);
        }

        itr++;
    }

    index++;
    index_s = std::to_string(index);
    filename_m = (char*)malloc(1 + strlen(filename)+ index_s.length() + ext.length());
    strcpy(filename_m, filename);
    strcat(filename_m, index_s.c_str());
    strcat(filename_m, ext.c_str());

    save_obj->save_configuration(filename_m);

    free(filename_m);


}

template <typename type>
void dlma_iterator<type>::run_system_for_percolation()
{
    sys_state->build_attachment_list();
    return;
}

template <typename type>
void dlma_iterator<type>::run_system_for_erdos_renyi()
{
    return;
}

template <typename type>
dlma_iterator<type>::~dlma_iterator()
{

    delete sys_state;
    delete binding_obj;
    delete agg_condition;
    delete aggregation_check_obj;
    delete movement_test;
    delete save_obj;

}

template class dlma_iterator<int>;
template class dlma_iterator<double>;

}