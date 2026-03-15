#include "system_iterator.hh"
#include <chrono>

#ifndef DLMA_ITERATOR_H
#define DLMA_ITERATOR_H

namespace simulation{

template <typename type>
class dlma_iterator: public system_iterator<type>{

    public:

        dlma_iterator(char *filename);
        ~dlma_iterator();
        void run_system();
        void iteration_step();
        void save_config_file();
        void save_config_file(char *filename);
        void create_movie_files(char *filename);
        void run_system_for_percolation();
        void run_system_for_erdos_renyi();



    private:

        std::vector<std::string> split_string_by_delimiter(const std::string& s, char delimiter);

        system<type> *sys_state;
        normal_bind<type> *binding_obj;
        aggregation_condition<type> *agg_condition;
        check_aggregation<type> *aggregation_check_obj;
        particle_movement<type> *movement_test;
        save_config<type> *save_obj;
        iterator_factory<type> *factory;


        int temp;
        constituent<type> *temp_c;
        int final_aggregate_number = 1;

        double time_1 = 0.;
        double time_2 = 0.;
        double time_3 = 0.;
        double time_4 = 0.;

        std::chrono::steady_clock::time_point cp_1;
        std::chrono::steady_clock::time_point cp_2;
        std::chrono::steady_clock::time_point cp_3;
        std::chrono::steady_clock::time_point cp_4;
        std::chrono::steady_clock::time_point cp_5;

};



}



#endif