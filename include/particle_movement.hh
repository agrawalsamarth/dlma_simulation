#include <iostream>
#include <random>

#ifndef PARTICLE_MOVEMENT_H
#define PARTICLE_MOVEMENT_H

namespace simulation{

template <typename type>
class particle_movement{

    public:

        virtual ~particle_movement() = default;
        virtual type*   delta_x() {};
        virtual double  get_rand() {};


};



}

#endif