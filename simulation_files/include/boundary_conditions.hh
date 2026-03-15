#include <iostream>

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

namespace simulation{

template <typename type>
class boundary_conditions{

    public:

        virtual ~boundary_conditions() = default;

        virtual type refill(type x, type L) {};
        virtual type refill(type old_pos, type new_pos, type L) {};
};

}

#endif