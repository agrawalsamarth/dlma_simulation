#include<cluster.hh>

namespace simulation{

template <typename type>
cluster<type>::cluster(const int cluster_id, const int dim, simulation_box<type> *system_box)
{
    id  = cluster_id;
    D   = dim;
    box = system_box;
}

template <typename type>
void cluster<type>::add_constituent(constituent<type> *single_element){
    elements.push_back(single_element);
    single_element->set_aggregate_id(id);
}

template <typename type>
void cluster<type>::move(type *delta_x){

    for (int i = 0; i < elements.size(); i++){
        elements[i]->move(delta_x);
    }

}

/*template <typename type>
type cluster<type>::element_pos(const int i, const int axis) const
{ return elements[i]->pos(axis);  }*/

template<typename type>
void cluster<type>::add_constituent_to_cell()
{ 

    for (int i = 0; i < elements.size(); i++)
        elements[i]->add_constituent_to_cell();    

}

template<typename type>
void cluster<type>::remove_constituent_from_cell()
{ 

    for (int i = 0; i < elements.size(); i++)    
        elements[i]->remove_constituent_from_cell(); 
    
}

template<typename type>
void cluster<type>::add_agg_to_cell()
{

    for (int i = 0; i < elements.size(); i++)
        elements[i]->add_agg_to_cell();

}

template<typename type>
void cluster<type>::remove_agg_from_cell()
{

    for (int i = 0; i < elements.size(); i++)
        elements[i]->remove_agg_from_cell();

}

template<typename type>
void cluster<type>::calculate_aggregate_mass()
{

    mass = 0.;

    for (int i = 0; i < elements.size(); i++)
        mass += elements[i]->get_mass();


}

template<typename type>
int cluster<type>::get_id()
{ return id;}

template<typename type>
int cluster<type>::get_aggregate_id()
{ return aggregate_id;}

template<typename type>
double cluster<type>::get_mass()
{return mass;}

template<typename type>
int cluster<type>::get_size()
{return elements.size();}

template<typename type>
constituent<type>* cluster<type>::get_element(const int i)
{ return elements[i];}

template<typename type>
int cluster<type>::get_element_aggregate_id(const int element_id)
{ return elements[element_id]->get_aggregate_id();}

template<typename type>
std::vector<int> cluster<type>::get_neighbour_list(const int i)
{  return elements[i]->get_neighbour_list();}

template<typename type>
std::vector<int> cluster<type>::get_neighbour_list_agg(const int i)
{  return elements[i]->get_neighbour_list_agg();}

template<typename type>
int cluster<type>::get_element_id(const int i)
{ return elements[i]->get_id();}

template<typename type>
void cluster<type>::set_current_seed_status(const int seed)
{
    for (int i = 0; i < elements.size(); i++)
        elements[i]->set_current_seed_status(seed);

}

/*template<typename type>
cluster<type>::~cluster()
{
    for (auto temp_delete_p : elements)
        delete temp_delete_p;
}*/

template class cluster<int>;
template class cluster<double>;

}