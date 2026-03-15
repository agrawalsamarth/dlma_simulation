#include "simulation_box.hh"

#ifndef OFF_LATTICE_H
#define OFF_LATTICE_H

namespace simulation{

template <typename type>
class off_lattice: public simulation_box<type>{

    public:

        off_lattice(const int dim, type *box_lengths, std::vector<boundary_conditions<type>*> system_bc, double tolerance);
        ~off_lattice();

        type  get_refill(type x, int axis);
        type  get_L(const int axis);

        void add_particle_to_cell(const int id, type *pos);
        void remove_particle_from_cell(const int id, type *pos);

        std::vector<int> get_neighbour_list(type *pos);

        int  get_periodicity(const int axis);
        type get_periodic_distance(type x, type y, int axis);
        
    private:

        std::vector<boundary_conditions<type>*> box_bc;
        type *L;
        type *halfL;
        int  *L_eff;
        int  *num_grid;
        int  *nx;
        int  *neighbour_x;
        type *inv_deltax;
        std::vector<std::vector<int>> grid;
        int  counter;
        int  D;
        std::vector<int> neighbours;
        int *periodic;
        type r;


};


}

#endif