#include <dlma_system.hh>

#ifndef DLMA_SYSTEM_ONLATTICE_H
#define DLMA_SYSTEM_ONLATTICE_H

namespace simulation{

template <typename type>
class dlma_system_onlattice: public dlma_system<type>{

    public:

        dlma_system_onlattice(char *params_name);
        ~dlma_system_onlattice();
        void initialize_system();
        bool check_viability(constituent<type> *c_1, type *dr);
        void move_aggregate(int i, type *dr);
        void print_grid();

    private:

        void initialize_system_for_percolation();
        void build_site_vector();
        void remove_neighbours(int i);
        void build_L_arrs();
        void remove_site_from_vector(int index);
        void exclude_neighbours(int index);
        type* get_array_from_index(int index);
        int  get_index_from_array(type *arr);
        void build_site_vector_for_rsp();

        
        std::vector<int> available_sites;

        int *L_eff;
        int  L_total;
        int  div;
        int  counter;
        type *temp_arr;
        int *temp_pos;
        type *temp_pos_arr;
        type *neighbour_pos;
        int *axis_arr;
        int  neighbour_index;




};



}

#endif