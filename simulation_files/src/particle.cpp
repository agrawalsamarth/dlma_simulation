#include <particle.hh>

namespace simulation{

template<typename type>
particle<type>::particle(const int particle_id, const int dim, simulation_box<type> *system_box){
    
    id  = particle_id;
    D   = dim;
    box = system_box;

    pos_ = (type*)malloc(sizeof(type)*D);

    for (int axis = 0; axis < D; axis++)
        pos_[axis] = 0;


}

template<typename type>
particle<type>::~particle()
{ free(pos_);}

template<typename type>
void particle<type>::set_mass(const double constituent_mass)
{ mass = constituent_mass;}

template<typename type>
double particle<type>::get_mass()
{ return mass;}

template<typename type>
type particle<type>::pos(const int axis) const
{ return pos_[axis];}

template<typename type>
type& particle<type>::pos(const int axis)
{ return pos_[axis];}

template<typename type>
void particle<type>::move(type *delta_x)
{

    for (int axis = 0; axis < D; axis++)
        pos_[axis] = box->get_refill(pos_[axis]+delta_x[axis], axis);

}

template<typename type>
void particle<type>::add_constituent_to_cell()
{ box->add_particle_to_cell(id, pos_); }

template<typename type>
void particle<type>::remove_constituent_from_cell()
{ box->remove_particle_from_cell(id, pos_); }

template<typename type>
void particle<type>::add_agg_to_cell()
{ box->add_agg_to_cell(aggregate_id, pos_); }

template<typename type>
void particle<type>::remove_agg_from_cell()
{ box->remove_agg_from_cell(aggregate_id, pos_); }

template<typename type>
void particle<type>::set_aggregate_id(const int cluster_id)
{ aggregate_id = cluster_id;}

template<typename type>
int particle<type>::get_id()
{ return id;}

template<typename type>
int particle<type>::get_aggregate_id()
{ return aggregate_id;}

template<typename type>
std::vector<int> particle<type>::get_neighbour_list()
{   return box->get_neighbour_list(pos_);}

template<typename type>
std::vector<int> particle<type>::get_neighbour_list_agg()
{   return box->get_neighbour_list_for_agg(aggregate_id, pos_);}

template<typename type>
void particle<type>::set_diameter(const double dia)
{ diameter = dia;}

template<typename type>
double particle<type>::get_diameter()
{return diameter;}

template<typename type>
void particle<type>::set_original_seed_status(const int seed)
{ original_seed_status = seed;}

template<typename type>
int  particle<type>::get_original_seed_status()
{ return original_seed_status;}

template<typename type>
void particle<type>::set_current_seed_status(const int seed)
{ current_seed_status = seed;}

template<typename type>
int  particle<type>::get_current_seed_status()
{ return current_seed_status;}

template<typename type>
int particle<type>::get_size()
{ return 1;}

template<typename type>
void particle<type>::print_neighbours()
{

    neighbours = box->get_neighbour_list(pos_);

    for (int i = 0; i < neighbours.size(); i++)
        std::cout<<neighbours[i]<<"\t";


}

template class particle<int>;
template class particle<double>;

}