#include "iterator_factory.hh"

#ifndef SYSTEM_IETRATOR_H
#define SYSTEM_IETRATOR_H

namespace simulation
{

template <typename type>
class system_iterator{

    public:

        virtual ~system_iterator() {};

        virtual void run_system() {};
        virtual void iteration_step() {};
        virtual void save_config_file() {};
        virtual void save_config_file(char *filename) {};
        virtual void create_movie_files(char *filename) {};
        virtual void run_system_for_percolation() {};
        virtual void run_system_for_erdos_renyi() {};


};



}

#endif