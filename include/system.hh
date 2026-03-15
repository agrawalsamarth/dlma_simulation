#include "system_factory.hh"
#include <cstring>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>
#include <map>
#include <algorithm>
#include <split_string.hh>

#ifndef SYSTEM_H
#define SYSTEM_H

namespace simulation{

template <typename type>
class system{

    public:

        //virtual std::vector<std::string> split_string_by_delimiter(const std::string& s, char delimiter) {};
        virtual ~system() = default;
        virtual void read_params_parser(char *params_name) {};
        virtual void initialize_system() {};

        virtual void calculate_propensity()  {};

        virtual void add_aggregate(constituent<type> *new_aggregate) {};
        virtual void remove_aggregate(const int id) {};
        virtual constituent<type>* get_aggregate(const int id) {};

        virtual void build_id_map() {};

        virtual int get_latest_cluster_id() {};
        virtual int get_latest_cluster_id_without_increment() {};
        virtual int get_id_map(int c_id) {};

        virtual int get_lattice() {};
        virtual int get_dim() {};
        virtual int get_max_attachments() {};
        virtual int get_N() {};
        virtual double get_phi() {};
        virtual double get_alpha() {};
        virtual std::vector<int> get_attachment_vector(const int i) {};

        virtual simulation_box<type>* get_box() {};

        virtual constituent<type>* get_constituent(const int i) {};
        virtual constituent<type>* get_particle_by_index(const int i) {};
        virtual void print_id_map() {};
        virtual void print_grid() {};
        virtual bool check_viability(constituent<type> *c_1, type *dr) {};
        virtual void move_aggregate(const int i, type *dr) {};
        virtual void add_attachment(constituent<type> *c_1) {};
        virtual void add_attachment(const int i, const int j) {};
        virtual void print_attachments() {};

        virtual int total_aggregates() {};
        virtual double get_seedmass() {};
        virtual constituent<type>* get_particle_by_id(const int id) {};

        virtual type get_interparticle_distance(constituent<type> *p_1, constituent<type> *p_2) {};
        virtual int choose_aggregate() {};

        virtual void build_attachment_list() {};
        virtual void build_idx_map_for_agg() {};

};



}

#endif