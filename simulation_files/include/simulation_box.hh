#include <iostream>
#include <vector>

#include "periodic_bc.hh"

#ifndef SIMULATION_BOX_H
#define SIMULATION_BOX_H

namespace simulation{

template<typename type>
class simulation_box{

    public:

        virtual ~simulation_box() = default;

        virtual void add_bc(boundary_conditions<type> *bc) {};
        virtual void set_L(const type L_val, const int axis) {};
        virtual type get_L(const int axis) {};

        virtual type get_refill(type x, int axis) {};
        virtual type get_refill(int old_pos, int new_pos, int axis) {};

        virtual type periodic_distance(type x, int axis) {};

        virtual void add_particle_to_cell(const int id, type *pos) {};
        virtual void remove_particle_from_cell(const int id, type *pos) {};
        virtual int  get_particle_id(type *pos) {};
        virtual int  get_agg_id(type *pos) {};

        virtual std::vector<int> get_neighbour_list(type *ref_pos) {};
        virtual std::vector<int> get_neighbour_list_for_agg(int id, type *ref_pos) {};

        virtual int  get_periodicity(const int axis) {};
        virtual type get_periodic_distance(type x, type y, int axis) {};

        virtual void add_agg_to_cell(const int id, type *pos) {};
        virtual void remove_agg_from_cell(const int id, type *pos) {};
        virtual void clear_cell_field() {};

        

};


}

#endif