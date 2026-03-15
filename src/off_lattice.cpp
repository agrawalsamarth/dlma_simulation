#include <off_lattice.hh>

namespace simulation{

template <typename type>
off_lattice<type>::off_lattice(const int dim, type *box_lengths, std::vector<boundary_conditions<type>*> system_bc, double tolerance){

    D = dim;
    L = (type*)malloc(sizeof(type) * D);
    halfL = (type*)malloc(sizeof(type) * D);
    inv_deltax = (type*)malloc(sizeof(type)*D);
    num_grid = (int*)malloc(sizeof(int)*D);
    nx = (int*)malloc(sizeof(int) * D);
    neighbour_x = (int*)malloc(sizeof(int) * D);

    for (int axis = 0; axis < D; axis++){
        L[axis]     = box_lengths[axis];
        halfL[axis] = 0.5 * L[axis];
    } 

    type delta_x[D];

    for (int axis = 0; axis < D; axis++){
        delta_x[axis] = (1. + tolerance);
        inv_deltax[axis] = 1./delta_x[axis];
    }

    for (int axis = 0; axis < D; axis++){
        num_grid[axis] = (int)(L[axis] * inv_deltax[axis]);
        nx[axis] = num_grid[axis];
    }

    for (int axis = 0; axis < D; axis++){
        delta_x[axis]    = L[axis]/(1. * nx[axis]);
        inv_deltax[axis] = 1./delta_x[axis];
    }

    L_eff = (int*)malloc(sizeof(int) * D);

    for (int axis = 0; axis < D; axis++)
        L_eff[axis] = 1;

    for (int axis = 0; axis < D; axis++)
    {
        for (int itr = axis+1; itr < D; itr++)
            L_eff[axis] *= nx[itr];
    }

    int L_total = 1;

    for (int axis = 0; axis < D; axis++)
        L_total *= nx[axis];

    grid.resize(L_total);

    for (int axis = 0; axis < D; axis++)
        box_bc.push_back(system_bc[axis]);

    periodic = (int*)malloc(sizeof(int) * D);

    for (int axis = 0; axis < D; axis++)
        periodic[axis] = 1;
 
}

template <typename type>
off_lattice<type>::~off_lattice()
{

    free(L);
    free(halfL);
    free(inv_deltax);
    free(periodic);
    free(num_grid);
    free(nx);
    free(neighbour_x);
    free(L_eff);

    for (int i = 0; i < box_bc.size(); i++)
        delete box_bc[i];
}

template <typename type>
type off_lattice<type>::get_L(const int axis)
{  return L[axis]; }

template <typename type>
type off_lattice<type>::get_refill(type x, int axis){
    return (box_bc[axis]->refill(x, L[axis]));
}

template <typename type>
void off_lattice<type>::add_particle_to_cell(const int id, type *pos){

    counter = 0;

    for (int axis = 0; axis < D; axis++)
        nx[axis] = (int)(pos[axis]*inv_deltax[axis]);

    for (int axis = 0; axis < D; axis++)
        counter += (nx[axis]*L_eff[axis]);

    grid[counter].push_back(id);

}

template <typename type>
void off_lattice<type>::remove_particle_from_cell(const int id, type *pos){

    counter = 0;

    for (int axis = 0; axis < D; axis++)
        nx[axis] = (int)(pos[axis]*inv_deltax[axis]);

    for (int axis = 0; axis < D; axis++)
        counter += (nx[axis]*L_eff[axis]);

    for (int i = 0; i < grid[counter].size(); i++){

        if (grid[counter][i] == id)
            grid[counter].erase(grid[counter].begin()+i);

    }    

}


template <typename type>
std::vector<int> off_lattice<type>::get_neighbour_list(type *pos){

    neighbours.clear();

    for (int axis = 0; axis < D; axis++)
        nx[axis] = (int)(pos[axis]*inv_deltax[axis]);


    switch (D){


        case 2:

            for (int i = -1; i <= 1; i++){
                for (int j = -1; j <=1; j++){

                    neighbour_x[0]  = (nx[0] + i);
                    neighbour_x[1]  = (nx[1] + j);

                    neighbour_x[0] += periodic[0] * (num_grid[0] * (neighbour_x[0] == -1) - num_grid[0] * (neighbour_x[0] == num_grid[0]));
                    neighbour_x[1] += periodic[1] * (num_grid[1] * (neighbour_x[1] == -1) - num_grid[1] * (neighbour_x[1] == num_grid[1]));

                    counter = 0;

                    for (int axis = 0; axis < D; axis++)
                        counter += (neighbour_x[axis]*L_eff[axis]);

                    for (int n = 0; n < grid[counter].size(); n++)
                        neighbours.push_back(grid[counter][n]);

                }
            }

            break;

        case 3:

            for (int i = -1; i <= 1; i++){
                for (int j = -1; j <=1; j++){
                    for (int k = -1; k <= 1; k++) {

                        neighbour_x[0]  = (nx[0] + i);
                        neighbour_x[1]  = (nx[1] + j);
                        neighbour_x[2]  = (nx[2] + k);


                        neighbour_x[0] += periodic[0] * (num_grid[0] * (neighbour_x[0] == -1) - num_grid[0] * (neighbour_x[0] == num_grid[0]));
                        neighbour_x[1] += periodic[1] * (num_grid[1] * (neighbour_x[1] == -1) - num_grid[1] * (neighbour_x[1] == num_grid[1]));
                        neighbour_x[2] += periodic[2] * (num_grid[2] * (neighbour_x[2] == -1) - num_grid[2] * (neighbour_x[2] == num_grid[2]));

                        counter = 0;

                        for (int axis = 0; axis < D; axis++)
                            counter += (neighbour_x[axis]*L_eff[axis]);

                        for (int n = 0; n < grid[counter].size(); n++)
                            neighbours.push_back(grid[counter][n]);

                    }
                }
            }

            break;

    }


    return neighbours;

}

template <typename type>
int off_lattice<type>::get_periodicity(const int axis)
{return periodic[axis];}

template <typename type>
type off_lattice<type>::get_periodic_distance(type x, type y, int axis)
{
    r = x - y;
    r = r + periodic[axis] * (L[axis]*(r < -halfL[axis]) - L[axis]*(r > halfL[axis]));

    return r;
} 

template class off_lattice<int>;
template class off_lattice<double>;


}