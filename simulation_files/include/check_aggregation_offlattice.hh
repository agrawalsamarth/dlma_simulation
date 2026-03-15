#include "check_aggregation.hh"

#ifndef CHECK_AGGREGATION_OFFLATTICE_H
#define CHECK_AGGREGATION_OFFLATTICE_H

namespace simulation{

template <typename type>
class check_aggregation_offlattice: public check_aggregation<type>{

    public:

        check_aggregation_offlattice(system<type> *system_state, normal_bind<type> *bind_system, aggregation_condition<type> *ref_condition, double tolerance);
        //~check_aggregation_offlattice();
        void check_for_aggregation(constituent<type> *c_1);

    private:

        system<type> *sys_state;
        normal_bind<type> *bind_sys;
        aggregation_condition<type> *condition;
        std::vector<int> neighbours;
        constituent<type>* temp;

        int neighbour_id;
        int neighbour_cluster_id;
        int cluster_id;

        type distance;

        constituent<type>* particle_1;
        constituent<type>* particle_2;

        type temp_r;
        type r2;

        double agg_tolerance;



};


}

#endif