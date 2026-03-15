#include <aggregation_condition.hh>

#ifndef MASS_AGGREGATION_CONDITION_H
#define MASS_AGGREGATION_CONDITION_H

namespace simulation{

template<typename type>
class mass_aggregation_condition: public aggregation_condition<type>{

    public:

        mass_aggregation_condition(system<type> *ref_sys);
        bool agg_condition(constituent<type> * c_1, constituent<type> * c_2);
        void show_out();

    private:

        system<type> *sys_state;
        int id_1;
        int id_2;
        int agg_1;
        int agg_2;
        constituent<type> *cluster_1;
        constituent<type> *cluster_2;


};



}


#endif