#include <brownian_movement_offlattice.hh>

namespace simulation{

template <typename type>
brownian_movement_offlattice<type>::brownian_movement_offlattice(int dim, int rng_seed){

    D  = dim;
    dr = (type*)malloc(sizeof(type) * D);

    for (int axis = 0; axis < D; axis++)
        dr[axis] = 0;


    generator.seed(rng_seed);
    dis.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
}

template <typename type>
type* brownian_movement_offlattice<type>::delta_x()
{

    r2 = 0.;

    for (int axis = 0; axis < D; axis++){
        temp     = -1. + 2. * dis(generator);
        dr[axis] = temp;
        r2      += temp * temp;
    }

    r2 = sqrt(r2);

    for (int axis = 0; axis < D; axis++)
        dr[axis] = dr[axis]/(1. * r2);

    return dr;

}

template <typename type>
double brownian_movement_offlattice<type>::get_rand()
{ return dis(generator); }

template <typename type>
brownian_movement_offlattice<type>::~brownian_movement_offlattice()
{ free(dr); }

template class brownian_movement_offlattice<int>;
template class brownian_movement_offlattice<double>;

}