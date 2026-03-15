#include <on_lattice.hh>

namespace simulation{

template <typename type>
on_lattice<type>::on_lattice(const int dim, type *box_lengths, std::vector<boundary_conditions<type>*> system_bc){

    D = dim;
    L = (type*)malloc(sizeof(int) * D);

    temp_pos      = (type*)malloc(sizeof(type) * D);
    neighbour_pos = (type*)malloc(sizeof(type) * D);

    for (int axis = 0; axis < D; axis++)
        L[axis] = box_lengths[axis];

    L_eff = (int*)malloc(sizeof(int) * D);

    for (int axis = 0; axis < D; axis++)
        L_eff[axis] = 1;

    for (int axis = 0; axis < D; axis++)
    {
        for (int itr = axis+1; itr < D; itr++)
            L_eff[axis] *= L[itr];
    }

    int L_total = 1;

    for (int axis = 0; axis < D; axis++)
        L_total *= L[axis];

    grid = (int*)malloc(sizeof(int) * L_total);
    for (int i = 0; i < L_total; i++)
        grid[i] = -1;

    /*agg_grid = (int*)malloc(sizeof(int) * L_total);
    for (int i = 0; i < L_total; i++)
        agg_grid[i] = -1;*/

    for (int axis = 0; axis < D; axis++)
        box_bc.push_back(system_bc[axis]);

    periodic = (int*)malloc(sizeof(int) * D);

    for (int axis = 0; axis < D; axis++)
        periodic[axis] = 1;
}

template <typename type>
type on_lattice<type>::get_L(const int axis)
{ return L[axis]; }

template <typename type>
type on_lattice<type>::get_refill(type x, int axis){

    return (box_bc[axis]->refill(x, L[axis]));

}

template <typename type>
void on_lattice<type>::add_particle_to_cell(const int id, type *pos){

    counter = 0;

    for (int axis = 0; axis < D; axis++)
            counter += pos[axis] * L_eff[axis];

    grid[counter] = id;

}

template <typename type>
void on_lattice<type>::remove_particle_from_cell(const int id, type *pos){

    counter = 0;

    for (int axis = 0; axis < D; axis++)
            counter += pos[axis] * L_eff[axis];

    grid[counter] = -1;

}

template<typename type>
void on_lattice<type>::add_agg_to_cell(const int id, type *pos){

    counter = 0;

    for (int axis = 0; axis < D; axis++)
            counter += pos[axis] * L_eff[axis];

    agg_grid[counter] = id;

}

template <typename type>
void on_lattice<type>::remove_agg_from_cell(const int id, type *pos){

    counter = 0;

    for (int axis = 0; axis < D; axis++)
            counter += pos[axis] * L_eff[axis];

    agg_grid[counter] = -1;

}

template <typename type>
int on_lattice<type>::get_particle_id(type *pos){

    counter = 0;

    for (int axis = 0; axis < D; axis++)
            counter += pos[axis] * L_eff[axis];

    return grid[counter];


}

template <typename type>
int on_lattice<type>::get_agg_id(type *pos) {

    counter = 0;

    for (int axis = 0; axis < D; axis++)
            counter += pos[axis] * L_eff[axis];

    return agg_grid[counter];

}

template <typename type>
std::vector<int> on_lattice<type>::get_neighbour_list(type *pos){

    neighbours.clear();

    for (int axis = 0; axis < D; axis++){
        for (int i = 0; i < 2; i++){

            for (int axis_2 = 0; axis_2 < D; axis_2++)
                temp_pos[axis_2] = (axis_2 == axis) * (-1 + 2*i);

            for (int axis_2 = 0; axis_2 < D; axis_2++){
                neighbour_pos[axis_2]  = (pos[axis_2] + temp_pos[axis_2]);
                neighbour_pos[axis_2] += (L[axis_2] * (neighbour_pos[axis_2]==-1)) - (L[axis_2] * (neighbour_pos[axis_2]==L[axis_2]));
            }

            agg_index = get_particle_id(neighbour_pos);

            //if (agg_index != -1)
                neighbours.push_back(agg_index);            

        }
    }

    return neighbours;

}

template <typename type>
std::vector<int> on_lattice<type>::get_neighbour_list_for_agg(int id, type *pos){

    neighbours.clear();

    for (int axis = 0; axis < D; axis++){
        for (int i = 0; i < 2; i++){

            for (int axis_2 = 0; axis_2 < D; axis_2++)
                temp_pos[axis_2] = (axis_2 == axis) * (-1 + 2*i);

            for (int axis_2 = 0; axis_2 < D; axis_2++){
                neighbour_pos[axis_2]  = (pos[axis_2] + temp_pos[axis_2]);
                neighbour_pos[axis_2] += (L[axis_2] * (neighbour_pos[axis_2]==-1)) - (L[axis_2] * (neighbour_pos[axis_2]==L[axis_2]));
            }

            agg_index = get_agg_id(neighbour_pos);

            if ((agg_index != -1) && (agg_index != id))
                neighbours.push_back(agg_index);

        }
    }

    return neighbours;


}

template <typename type>
void on_lattice<type>::clear_cell_field()
{

    int L_total = 1;

    for (int axis = 0; axis < D; axis++)
        L_total *= L[axis];

    for (int i = 0; i < L_total; i++){
        grid[i]     = -1;
    }

}

template <typename type>
int on_lattice<type>::get_periodicity(const int axis)
{return periodic[axis];}

template <typename type>
on_lattice<type>::~on_lattice()
{
    free(neighbour_pos);
    free(L_eff);
    free(grid);
    free(L);
    free(temp_pos);
    free(periodic);
    
    for (int i = 0; i < box_bc.size(); i++)
        delete box_bc[i];

}

template class on_lattice<int>;
template class on_lattice<double>;


}