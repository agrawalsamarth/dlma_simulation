#include "particle.hh"
#include "cluster.hh"
#include "periodic_bc.hh"
#include "on_lattice.hh"
#include "off_lattice.hh"
#include <string.h>

#ifndef SYSTEM_FACTORY_H
#define SYSTEM_FACTORY_H

namespace simulation{

template<typename type>
class system_factory{

    public:

        virtual ~system_factory() = default;

        constituent<type>*               create_constituent(int constituent_id, int lattice, int dim, std::string name_type, simulation_box<type> *box);
        simulation_box<type>*            create_simulation_box(int lattice, int dim, type *box_lengths, std::vector<boundary_conditions<type>*> system_bc, double tolerance);
        boundary_conditions<type>*       create_boundary_conditions(std::string name_type);


};

}


#endif