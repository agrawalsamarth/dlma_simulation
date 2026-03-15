#include <dlma_system.hh>
#include <dlma_system_onlattice.hh>
#include <dlma_system_offlattice.hh>

#ifndef NORMAL_BIND_H
#define NORMAL_BIND_H

namespace simulation{

template<typename type>
class normal_bind{

    public:
        normal_bind(system<type> *system_ptr);
        constituent<type>* bind_aggregates(constituent<type> *c_1, constituent<type> *c_2);

    private:

        system_factory<type> factory;
        system<type>        *sys_state;

};

}

#endif