#include "constituent.hh"

#ifndef PARTICLE_H
#define PARTICLE_H

namespace simulation{

template <typename type>
class particle: public constituent<type>{

    public:

        particle(const int particle_id, const int dim, simulation_box<type> *system_box);
        ~particle();

        void   move(type *delta_x);
        double get_mass();
        void   set_mass(const double constituent_mass);

        type  pos(const int axis) const; 
        type& pos(const int axis);

        void add_constituent_to_cell();
        void remove_constituent_from_cell();

        void add_agg_to_cell();
        void remove_agg_from_cell();

        void set_aggregate_id(const int id);
        
        int get_id();
        int get_aggregate_id();
        std::vector<int> get_neighbour_list();
        std::vector<int> get_neighbour_list_agg();

        void set_diameter(const double dia);
        double get_diameter();

        void set_original_seed_status(const int seed);
        int  get_original_seed_status();

        void set_current_seed_status(const int seed);
        int  get_current_seed_status();

        void print_neighbours();
        int  get_size();

        
    private:

        int                     D;
        int                     id;
        int                     aggregate_id;
        double                  mass;
        simulation_box<type>   *box;
        type                   *pos_;      
        std::vector<int>        neighbours;
        double                  diameter;
        int                     original_seed_status;
        int                     current_seed_status;  


};


}

#endif