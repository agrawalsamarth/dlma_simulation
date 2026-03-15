#include "boundary_conditions.hh"

#ifndef PERIODIC_H
#define PERIODIC_H

namespace simulation{

template <typename type>
class periodic_bc: public boundary_conditions<type>{

    public:
    
        type refill(type x, type L);


};


}

#endif