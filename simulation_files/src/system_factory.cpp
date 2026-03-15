#include <system_factory.hh>

namespace simulation{

template <typename type>
constituent<type>* system_factory<type> ::create_constituent(int constituent_id, int lattice, int dim, std::string name_type, simulation_box<type> *box){

    if ((strcmp(name_type.c_str(),"particle")==0))
        return new particle<type>(constituent_id, dim, box);

    else if ((strcmp(name_type.c_str(),"cluster")==0))
        return new cluster<type>(constituent_id, dim, box);

    else{
        std::cout<<"unknown constituent"<<std::endl;
        exit(EXIT_FAILURE);
    }
}

template <typename type>
simulation_box<type>* system_factory<type>::create_simulation_box(int lattice, int dim, type *box_lengths, std::vector<boundary_conditions<type>*> system_bc, double tolerance){

    if (lattice == 1){
        return new on_lattice<type>(dim, box_lengths, system_bc);
    }

    else if (lattice == 0){
        return new off_lattice<type>(dim, box_lengths, system_bc, tolerance);
    }

    else{
        std::cout<<"unknown simulation box"<<std::endl;
        exit(EXIT_FAILURE);
    }

}

template <typename type>
boundary_conditions<type>* system_factory<type>::create_boundary_conditions(std::string name_type){

    if (strcmp(name_type.c_str(), "periodic") == 0)
        return new periodic_bc<type>;
    else{
        std::cout<<"unknown boundary condition"<<std::endl;
        exit(EXIT_FAILURE);
    }

}

template class system_factory<int>;
template class system_factory<double>;

}