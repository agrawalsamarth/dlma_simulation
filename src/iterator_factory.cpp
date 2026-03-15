#include <iterator_factory.hh>

namespace simulation{

template <typename type>
particle_movement<type>* iterator_factory<type>::create_movement(std::string name_type, int dim, int rng_seed, int lattice)
{
    if ((strcmp(name_type.c_str(), "brownian")==0) && (lattice == 1))
        return new brownian_movement<type>(dim, rng_seed);
    else if ((strcmp(name_type.c_str(), "brownian")==0) && (lattice == 0))
        return new brownian_movement_offlattice<type>(dim, rng_seed);
    else{
        std::cout<<"unknown movement type"<<std::endl;
        exit(EXIT_FAILURE);
    }
}

template <typename type>
aggregation_condition<type>* iterator_factory<type>::create_aggregation_condition(std::string name_type, system<type> *system_state)
{
    if (strcmp(name_type.c_str(), "mass")==0)
        return new mass_aggregation_condition<type>(system_state);
    else{
        std::cout<<"unknown aggregate condition"<<std::endl;
        exit(EXIT_FAILURE);
    }
    
}

template <typename type>
check_aggregation<type>* iterator_factory<type>::create_check_aggregation(std::string name_type, int lattice, system<type> *system_state, normal_bind<type> *bind_system,
 aggregation_condition<type> *ref_condition, double tolerance)
{

    if ((strcmp(name_type.c_str(), "normal")==0) && (lattice == 1))
        return new check_aggregation_onlattice<type>(system_state, bind_system, ref_condition);
    else if ((strcmp(name_type.c_str(), "normal")==0) && (lattice == 0))
        return new check_aggregation_offlattice<type>(system_state, bind_system, ref_condition, tolerance);
    else{
        std::cout<<"unknown aggregation type"<<std::endl;
        exit(EXIT_FAILURE);
    }

}

template <typename type>
save_config<type>* iterator_factory<type>::create_save_config(std::string name_type, system<type> *ref_sys, simulation_box<type> *ref_box)
{

    if (strcmp(name_type.c_str(), "dlma")==0)
        return new dlma_save_config<type>(ref_sys, ref_box);
    else if (strcmp(name_type.c_str(), "random_site_percolation")==0)
        return new dlma_save_config<type>(ref_sys, ref_box);
    else if (strcmp(name_type.c_str(), "erdos_renyi")==0)
        return new dlma_save_config<type>(ref_sys, ref_box);        
    else{
        std::cout<<"unknown system type"<<std::endl;
        exit(EXIT_FAILURE);
    }


}

template <typename type>
normal_bind<type>* iterator_factory<type>::create_bind_system(std::string name_type, system<type> *system_ptr)
{

    if (strcmp(name_type.c_str(), "normal")==0)
        return new normal_bind<type>(system_ptr);
    else{
        std::cout<<"unknown bind type"<<std::endl;
        exit(EXIT_FAILURE);
    }



}

template <typename type>
system<type>* iterator_factory<type>::create_new_system(std::string name_type, int lattice, char *filename)
{
    if ((strcmp(name_type.c_str(), "dlma")==0) && (lattice == 1))
        return new dlma_system_onlattice<type>(filename);
    else if ((strcmp(name_type.c_str(), "dlma")==0) && (lattice == 0))
        return new dlma_system_offlattice<type>(filename);
    else if ((strcmp(name_type.c_str(), "erdos_renyi")==0))
        return new dlma_system_offlattice<type>(filename);
    else if ((strcmp(name_type.c_str(), "random_site_percolation")==0) && (lattice == 1))
        return new dlma_system_onlattice<type>(filename);
    else{
        std::cout<<"unknown system type"<<std::endl;
        exit(EXIT_FAILURE);
    }        

}

template class iterator_factory<int>;
template class iterator_factory<double>;

}