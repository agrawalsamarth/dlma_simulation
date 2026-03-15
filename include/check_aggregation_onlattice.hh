#include "check_aggregation.hh"
#include <chrono>

#ifndef CHECK_AGGREGATION_ONLATTICE_H
#define CHECK_AGGREGATION_ONLATTICE_H

namespace simulation{

template <typename type>
class check_aggregation_onlattice: public check_aggregation<type>{

    public:

        check_aggregation_onlattice(system<type> *system_state, normal_bind<type> *bind_system, aggregation_condition<type> *ref_condition);
        void check_for_aggregation(constituent<type> *c_1);
        void display_compute_times();

    private:

        system<type> *sys_state;
        normal_bind<type> *bind_sys;
        aggregation_condition<type> *condition;
        std::vector<int> neighbours;
        constituent<type>* temp;

        int neighbour_id;
        int neighbour_cluster_id;
        int cluster_id;

        constituent<type>* particle_1;
        constituent<type>* particle_2;

        double time_1 = 0.;
        double time_2 = 0.;
        double time_3 = 0.;
        double time_4 = 0.;

        std::chrono::steady_clock::time_point cp_1;
        std::chrono::steady_clock::time_point cp_2;
        std::chrono::steady_clock::time_point cp_3;
        std::chrono::steady_clock::time_point cp_4;
        std::chrono::steady_clock::time_point cp_5;






};


}

#endif