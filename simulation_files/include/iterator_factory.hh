#include "brownian_movement.hh"
#include "mass_aggregation_condition.hh"
#include "dlma_save_config.hh"
#include "normal_bind.hh"
#include "dlma_system_onlattice.hh"
#include "dlma_system.hh"
#include "brownian_movement_offlattice.hh"
#include "check_aggregation_offlattice.hh"
#include "check_aggregation_onlattice.hh"
#include "dlma_system_offlattice.hh"

#ifndef ITERATOR_FACTORY_H
#define ITERATOR_FACTORY_H

namespace simulation{

template <typename type>
class iterator_factory{

    public:

        virtual ~iterator_factory() = default;
        particle_movement<type>* create_movement(std::string name_type, int dim, int rng_seed, int lattice);
        aggregation_condition<type>* create_aggregation_condition(std::string name_type, system<type> *system_state); 
        check_aggregation<type>* create_check_aggregation(std::string name_type, int lattice, system<type> *system_state, normal_bind<type> *bind_system, aggregation_condition<type> *ref_condition
        , double tolerance);
        save_config<type>* create_save_config(std::string name_type, system<type> *ref_sys, simulation_box<type> *ref_box);
        normal_bind<type>* create_bind_system(std::string name_type, system<type> *system_ptr);
        system<type>* create_new_system(std::string name_type, int lattice, char *filename);



};



}


#endif