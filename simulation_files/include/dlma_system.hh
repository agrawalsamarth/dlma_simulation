#include <system.hh>

#ifndef DLMA_SYSTEM_H
#define DLMA_SYSTEM_H

namespace simulation{

template <typename type>
class dlma_system: public system<type>{

    public:

        //virtual std::vector<std::string> split_string_by_delimiter(const std::string& s, char delimiter);
        virtual ~dlma_system() = default;
        virtual void read_params_parser(char *params_name);
        virtual void initialize_system() {};

        virtual void calculate_propensity();

        virtual void add_aggregate(constituent<type> *new_aggregate);
        virtual void remove_aggregate(const int id);
        virtual constituent<type>* get_aggregate(const int id);

        virtual void build_id_map();

        virtual int get_latest_cluster_id();
        virtual int get_latest_cluster_id_without_increment();
        virtual int get_id_map(int c_id);

        virtual int get_lattice();
        virtual int get_dim();
        virtual int get_max_attachments();
        virtual int get_N();
        virtual void build_idx_map_for_agg();
        virtual double get_phi();
        virtual double get_alpha();
        virtual std::vector<int> get_attachment_vector(const int i);

        virtual simulation_box<type>* get_box();

        virtual constituent<type>* get_constituent(const int i);
        virtual void print_id_map();
        virtual void print_grid() {};
        virtual bool check_viability(constituent<type> *c_1, type *dr) {};
        virtual void move_aggregate(const int i, type *dr) {};
        virtual void add_attachment(constituent<type> *c_1);
        virtual void add_attachment(const int i, const int j) {};
        virtual void print_attachments();

        virtual int total_aggregates();
        virtual double get_seedmass();
        virtual constituent<type>* get_particle_by_id(const int id);
        virtual type get_interparticle_distance(constituent<type> *p_1, constituent<type> *p_2);
        virtual int choose_aggregate();

        virtual void build_attachment_list() {};
        virtual constituent<type>* get_particle_by_index(const int i);
        

    protected:

        std::mt19937 generator;
        std::uniform_real_distribution<double> dis;

        system_factory<type> factory;
        simulation_box<type> *box;

        std::vector<constituent<type>*> all_particles;
        std::vector<constituent<type>*> aggregates;

        std::map<int, int> id_map;
        std::map<int, int> agg_id_map;
        std::vector<std::vector<int>> attachments;

        std::vector<double> propensity;
        double total_propensity;

        int    N;
        bool   N_flag = false;
        int    N_s;
        bool   N_s_flag = false;
        double seed_pct;
        bool   seed_pct_flag = false;
        int    D;
        int    D_flag = false;
        int    lattice;
        bool   lattice_flag = false;
        double phi;
        bool   phi_flag = false;
        double alpha = 0.5;
        bool   alpha_flag = false;
        double seed_mass;
        bool   seed_mass_flag = false;
        int    rng_seed = 0;
        bool   rng_seed_flag = false;
        type   *L;
        type   *halfL;
        bool   *L_flag;
        bool    L_flag_and;
        bool    L_flag_or;
        int    *periodic;
        bool   *periodic_flag;
        bool    periodic_flag_and;
        int     latest_cluster_id=0;
        double  tolerance;
        bool    tolerance_flag=false;
        int     distance_metric_rgg=0;

        bool    is_viable;
        type    temp_r;
        type    r2;
        std::string system_type = "none";

        double get_manhattan_distance(constituent<type> *p_1, constituent<type> *p_2);

    
    private:

        void check_for_dlma_params();
        void check_for_percolation_params();
        int agg_id;
        void print_agg_map();
        void check_for_erdos_renyi_params();
        

};



}

#endif