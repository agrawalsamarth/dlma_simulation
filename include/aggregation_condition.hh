#include <dlma_system.hh>

#ifndef AGGREGATION_CONDITION_H
#define AGGREGATION_CONDITION_H

namespace simulation{

template <typename type>
class aggregation_condition{

    public:

        virtual void show_out() {};
        virtual bool agg_condition(constituent<type> *c_1, constituent<type> *c_2) {};


};


}

#endif