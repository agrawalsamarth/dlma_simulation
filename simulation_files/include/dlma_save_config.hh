#include <save_config.hh>

#ifndef DLMA_SAVE_CONFIG_H
#define DLMA_SAVE_CONFIG_H

namespace simulation{

template<typename type>
class dlma_save_config: public save_config<type>{


    public:

        dlma_save_config(system<type> *ref_sys, simulation_box<type> *ref_box);
        ~dlma_save_config() = default;
        void save_configuration();
        void save_configuration(char *filename);
        

    private:

        system<type>         *sys_state;
        simulation_box<type> *box;
        std::vector<std::string> print_data;
        
        void convert_data_to_string();

};


}


#endif