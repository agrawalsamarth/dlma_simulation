#include <particle_movement.hh>

#ifndef BROWNIAN_MOVEMENT_OFFLATTICE_H
#define BROWNIAN_MOVEMENT_OFFLATTICE_H

namespace simulation{

template <typename type>
class brownian_movement_offlattice: public particle_movement<type>{

    public:

        brownian_movement_offlattice(int dim, int rng_seed);
        ~brownian_movement_offlattice();
        type*   delta_x();
        double get_rand();

    private:

        std::mt19937 generator;
        std::uniform_real_distribution<double> dis;
        int  D;
        type *dr;
        int  axis_rand;
        int  sign_rand;
        double temp;
        double r2;



};


}


#endif