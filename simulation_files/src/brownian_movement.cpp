#include <brownian_movement.hh>

namespace simulation{

template <typename type>
brownian_movement<type>::brownian_movement(int dim, int rng_seed){

    D  = dim;
    dr = (type*)malloc(sizeof(type) * D);

    for (int axis = 0; axis < D; axis++)
        dr[axis] = 0;


    generator.seed(rng_seed);
    dis.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
}

template <typename type>
type* brownian_movement<type>::delta_x()
{

    temp      = dis(generator);
    sign_rand = (temp > 0.5) - (temp <= 0.5); 
    axis_rand = (int)(dis(generator) * D);

    for (int axis = 0; axis < D; axis++){
        dr[axis] = (sign_rand * (axis == axis_rand));
        //std::cout<<"dr = "<<dr[axis]<<std::endl;
    }

    return dr;

}

template <typename type>
double brownian_movement<type>::get_rand()
{ return dis(generator); }

template <typename type>
brownian_movement<type>::~brownian_movement()
{free(dr);}

template class brownian_movement<int>;
template class brownian_movement<double>;

}