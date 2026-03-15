#include "normal_bind.hh"
#include "aggregation_condition.hh"

#ifndef CHECK_AGGREGATION_H
#define CHECK_AGGREGATION_H

namespace simulation{

template <typename type>
class check_aggregation{

    public:

        virtual ~check_aggregation() = default;

        virtual void check_for_aggregation(constituent<type> *c_1) {};
        virtual void display_compute_times() {};

};


}

#endif