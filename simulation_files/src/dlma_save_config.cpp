#include <dlma_save_config.hh>

namespace simulation{

template <typename type>
dlma_save_config<type>::dlma_save_config(system<type> *ref_sys, simulation_box<type> *ref_box)
{
    sys_state = ref_sys;
    box       = ref_box;
}

template <typename type>
void dlma_save_config<type>::convert_data_to_string()
{

    std::string temp_line;
    constituent<type> *temp;
    std::vector<int>  attachments;

    int max_attachments = sys_state->get_max_attachments();
    int lattice = sys_state->get_lattice();
    int N = sys_state->get_N();
    int D = sys_state->get_dim();

    print_data.resize(N);

    int num_atts;

    for (int i = 0; i < N; i++)
    {

        temp_line = "";

        //temp        = sys_state->get_particle_by_id(i);
        temp        = sys_state->get_particle_by_index(i);
        attachments = sys_state->get_attachment_vector(i); 
        
        temp_line += std::to_string(i)+",";

        for (int axis = 0; axis < D; axis++){
            if (lattice == 1){
                temp_line += std::to_string(1. * (temp->pos(axis)) + 0.5)+",";
            }
            else
                temp_line += std::to_string(1. * (temp->pos(axis)))+",";
        }

        //for (int axis = 0; axis < D; axis++)
                //fprintf(f, "%lf,", 1. * (temp->pos(axis)) + 0.5);

        temp_line += std::to_string(temp->get_original_seed_status())+",";
        temp_line += std::to_string(temp->get_current_seed_status())+",";
        temp_line += std::to_string(temp->get_diameter())+",";

        num_atts = attachments.size();

        if (max_attachments > 0)
            temp_line += std::to_string(num_atts)+",";
        else
            temp_line += std::to_string(num_atts);

        for (int j = 0; j < max_attachments; j++){

            if (num_atts == max_attachments){

                if (j != (num_atts-1))
                    temp_line += std::to_string(attachments[j])+",";
                else
                    temp_line += std::to_string(attachments[j]);

            }

            else {

                if (j < num_atts)
                    temp_line += std::to_string(attachments[j])+",";
                else if ((j >= num_atts) && (j != (max_attachments-1)))
                    temp_line += "NaN,";
                else
                    temp_line += "NaN";

            }

        }

        //temp_line += "\n";

        print_data[i] = temp_line;
        

    }



}

template <typename type>
void dlma_save_config<type>::save_configuration(char *filename)
{
    /*FILE *f;
    f= fopen(filename,"w");*/

    //std::cout<<"1"<<std::endl;

    std::ofstream f;
    f.open(filename);

    int num_atts;
    int max_attachments = sys_state->get_max_attachments();
    int lattice = sys_state->get_lattice();
    int N = sys_state->get_N();
    int D = sys_state->get_dim();
    
    double phi   = sys_state->get_phi();
    double alpha = sys_state->get_alpha();

    constituent<type> *temp;
    std::vector<int>  attachments;

    int headers = 10+(3*D);

    /*fprintf(f, "headers=%d\n", headers);
    fprintf(f, "lattice=%d\n", lattice);
    fprintf(f, "N=%d\n", N);
    fprintf(f, "D=%d\n", D);
    fprintf(f, "maxAttachments=%d\n", max_attachments);
    fprintf(f, "folded=%d\n", 1);
    fprintf(f, "phi=%lf\n", phi);
    fprintf(f, "alpha=%lf\n", alpha);*/

    f << "headers="<<headers<<"\n";
    f << "lattice="<<lattice<<"\n";
    f << "N=" <<N << "\n";
    f << "D=" << D << "\n";
    f << "maxAttachments=" << max_attachments << "\n";
    f << "folded=" << 1 << "\n";
    f << "phi=" << phi << "\n";
    f << "alpha=" << alpha << "\n";



    for (int axis = 0; axis < D; axis++){
        f << "x"<<axis<<"_lo="<<0.0<<"\n";
        f << "x"<<axis<<"_hi="<<1. * box->get_L(axis)<<"\n";
    }

    for (int axis = 0; axis < D; axis++){
        f << "x"<<axis<<"_periodic="<<box->get_periodicity(axis)<<"\n";
    }

    f << "seedMass="<< sys_state->get_seedmass()<<"\n";
    f << "columns=" << (5+D+max_attachments) <<"\n";


    f << "id,";

    for (int axis = 0; axis < D; axis++)
        f << "x"<<axis<<",";
    
    f << "assignedSeedStatus,currentSeedStatus,diameter,";
    //fprintf(f, "currentSeedStatus,");
    //fprintf(f, "diameter,");
    if (max_attachments > 0)
        f << "attachments,";
    else
        f << "attachments\n";
        
    for (int att = 1; att <= max_attachments; att++){
        if (att == max_attachments)
            f << "att_"<<att<<"\n";
        else
            f << "att_"<<att<<",";
    }

    //std::cout<<"2"<<std::endl;

    convert_data_to_string();

    for (int i = 0; i < N; i++){

        //fprintf(f, "%s\n",print_data[i].c_str());
        f << print_data[i] << "\n";

    }

    //std::cout<<"3"<<std::endl;

    f.close();

}

template <typename type>
void dlma_save_config<type>::save_configuration()
{
    int num_atts;
    int max_attachments = sys_state->get_max_attachments();
    int lattice = sys_state->get_lattice();
    int N = sys_state->get_N();
    int D = sys_state->get_dim();
    
    double phi   = sys_state->get_phi();
    double alpha = sys_state->get_alpha();

    constituent<type> *temp;
    std::vector<int>  attachments;

    int headers = 10+(3*D);

    printf("headers=%d\n", headers);
    printf("lattice=%d\n", lattice);
    printf("N=%d\n", N);
    printf("D=%d\n", D);
    printf("maxAttachments=%d\n", max_attachments);
    printf("folded=%d\n", 1);
    printf("phi=%lf\n", phi);
    printf("alpha=%lf\n", alpha);

    for (int axis = 0; axis < D; axis++){
        printf("x%d_lo=%lf\n",axis,0.0);
        printf("x%d_hi=%lf\n",axis,1. * box->get_L(axis));
    }

    for (int axis = 0; axis < D; axis++){
        printf("x%d_periodic=%d\n",axis,box->get_periodicity(axis));
    }

    printf("seedMass=%lf\n", sys_state->get_seedmass());
    printf("columns=%d\n", (5+D+max_attachments));


    printf("id,");

    for (int axis = 0; axis < D; axis++)
        printf("x%d,", axis);
    
    printf("assignedSeedStatus,");
    printf("currentSeedStatus,");
    printf("diameter,");
    if (max_attachments > 0)
        printf("attachments,");
    else
        printf("attachments\n");
        
    for (int att = 1; att <= max_attachments; att++){
        if (att == max_attachments)
            printf("att_%d\n", att);
        else
            printf("att_%d,", att);
    }

    for (int i = 0; i < N; i++){

        temp        = sys_state->get_particle_by_id(i);
        attachments = sys_state->get_attachment_vector(i); 

        printf("%d,", i);

        for (int axis = 0; axis < D; axis++){
            if (lattice == 1)
                printf("%lf,", 1. * (temp->pos(axis)) + 0.5);
            else
                printf("%lf,", 1. * (temp->pos(axis)));
        }

        printf("%d,", temp->get_original_seed_status());
        printf("%d,", temp->get_current_seed_status());
        printf("%lf,",temp->get_diameter());

        num_atts = attachments.size();

        if (max_attachments > 0)
            printf("%d,",num_atts);
        else
            printf("%d",num_atts);

        

        for (int j = 0; j < max_attachments; j++){

            if (num_atts == max_attachments){

                if (j != (num_atts-1))
                    printf("%d,", attachments[j]);
                else
                    printf("%d", attachments[j]);

            }

            else {

                if (j < num_atts)
                    printf("%d,", attachments[j]);
                else if ((j >= num_atts) && (j != (max_attachments-1)))
                    printf("NaN,");
                else
                    printf("NaN");

            }

        }

        printf("\n");



    }


}

/*template<typename type>
dlma_save_config<type>::~dlma_save_config()
{
    print_data.resize(0);
    print_data.shrink_to_fit();
}*/

template class dlma_save_config<int>;
template class dlma_save_config<double>;


}