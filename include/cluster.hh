#include "constituent.hh"

#ifndef CLUSTER_H
#define CLUSTER_H

namespace simulation{

template <typename type>
class cluster: public constituent<type>{

    public:

        cluster(const int cluster_id, const int dim, simulation_box<type> *system_box);
        //~cluster();

        void move(type *delta_x);

        //type  element_pos(const int i, const int axis) const;

        void add_constituent(constituent<type> *single_element);

        void add_constituent_to_cell();
        void remove_constituent_from_cell();

        void add_agg_to_cell();
        void remove_agg_from_cell();

        //void set_aggregate_id(const int id);

        double get_mass();
        void   calculate_aggregate_mass();

        int   get_id();
        int   get_aggregate_id();
        int   get_size();
        int   get_element_id(const int i);

        constituent<type>* get_element(const int i);
        std::vector<int> get_neighbour_list(const int i);
        std::vector<int>   get_neighbour_list_agg(const int i);

        int get_element_aggregate_id(const int element_id);

        void set_current_seed_status(const int seed);


    private:

        std::vector<constituent<type>*> elements;
        int id;
        int aggregate_id;
        int D;
        simulation_box<type>* box;
        double mass;
        type temp_pos;


};


}

#endif