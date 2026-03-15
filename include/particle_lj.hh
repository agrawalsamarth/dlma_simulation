#include <particle.hh>

#ifndef PARTICLE_LJ_H
#define PARTICLE_LJ_H

namespace simulation
{

class particle_lj: public particle<int>{

    public:

        particle_lj();

        int  vel(const int axis) const; 
        int& vel(const int axis);

    private:

        double *vel_; 




};


}


#endif