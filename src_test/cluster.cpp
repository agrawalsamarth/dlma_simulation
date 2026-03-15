#include<cluster.hh>

namespace simulation{

template <typename type>
void cluster<type>::add_constituent(constituent<type> *single_element){
    elements.push_back(single_element);
}

template <typename type>
void cluster<type>::move(type *delta_x){

    for (int i = 0; i < elements.size(); i++){
        elements[i]->move(delta_x);
    }

}

template <typename type>
type cluster<type>::element_pos(const int i, const int axis) const 
{ return elements[i]->pos(axis);  }
    
template <typename type>
type& cluster<type>::element_pos(const int i, const int axis)
{ return elements[i]->pos(axis);  }


template class cluster<int>;
template class cluster<double>;

}