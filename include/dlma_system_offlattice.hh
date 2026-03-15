#include <dlma_system.hh>
#include <misc.hh>

#ifndef DLMA_SYSTEM_OFFLATTICE_H
#define DLMA_SYSTEM_OFFLATTICE_H

namespace simulation{

template <typename type>
class dlma_system_offlattice: public dlma_system<type>{

    public:

        dlma_system_offlattice(char *params_name);
        ~dlma_system_offlattice();
        void initialize_system();
        void move_aggregate(int i, type *dr);
        void add_attachment(const int i, const int j);
        void build_attachment_list();
    
    private:

        void calc_rij();
        //bool check_rij();
        std::vector<coll_deets> build_collision_list(int i, double alpha, type *dr);
        std::vector<int> available_sites_temp;
        coll_deets calc_quad_eqn(type *dr, double alpha, constituent<type> *ref_particle, constituent<type> *nb_particle);
        type fix_overlap(const int i, type *dr);
        void build_site_vector();
        void init_erdos_renyi();
        void build_temp_L_arrs();
        constituent<type> *image;
        std::vector<int>   actual_list;
        std::vector<int>   image_list;
        std::vector<int>   neighbours;
        int iters = 0;
        double *rij;
        std::vector<coll_deets> bonds;
        std::vector<coll_deets> hs;
        //void index_from_array_map();
        void  get_array_from_index_map(int i, int *for_pop_arr);
        int   get_index_from_array_map(int* temp_pos_arr);
        int   pop_counter;
        int   pop_i;
        int  *pop_arr;
        int  *pos_pop_temp;
        int  *neigh_pop_temp;
        int  *pos_arr_for_particle; 
        int   pop_div;
        int   pop_rem;
        int   L_disc;
        double *delta_L_temp;
        int    *L_eff_temp;
        int  *nx_temp;





};



}

#endif