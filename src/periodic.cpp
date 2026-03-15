#include <periodic_bc.hh>

namespace simulation{

template<typename type>
type periodic_bc<type>::refill(type x, type L){
    return (x + L * (x < 0) - L * (x >= L));
}

template class periodic_bc<int>;
template class periodic_bc<double>;

}