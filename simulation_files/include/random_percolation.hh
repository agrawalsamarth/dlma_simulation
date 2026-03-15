#include <dlma_system.hh>

#ifndef RANDOM_PERCOLATION_HH
#define RANDOM_PERCOLATION_HH

namespace simulation{

template <typename type>
class dlma_system: public random_percolation<type> {

    public:

        random_percolation(char *params_name);
        void initialize_system();
        void read_params_parser(char *params_name);

};



}

#endif