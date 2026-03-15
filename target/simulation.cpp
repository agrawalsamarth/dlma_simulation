#include <dlma_iterator.hh>
#include <split_string.hh>
//#include <chrono>

using namespace simulation;

int main(int argc, char *argv[])
{
    if (argc < 2){
        std::cout<<"please provide params file"<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::ifstream parser(argv[1], std::ifstream::in);

    if (parser.fail()){
        std::cout<<"either file does not exist or does not have permissions"<<std::endl;
        exit(EXIT_FAILURE);
    }

    std::string str;
    std::vector<std::string> results;

    std::string sys_name;
    int lattice_val;

    while (getline(parser, str)){

        results = split_string_by_delimiter(str, '=');

        if (results[0] == "system"){
            sys_name = results[1];
        }

        if (results[0] == "lattice"){
            lattice_val = stoi(results[1]);
        }

    }

    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    if ((sys_name == "dlma") && (lattice_val == 1)){
        system_iterator<int> *test = new dlma_iterator<int>(argv[1]);;
        test->run_system();

        if (argc < 3)
            test->save_config_file();
        else
            test->save_config_file(argv[2]);

        delete test;
        
    }

    else if((sys_name == "dlma") && (lattice_val == 0)){
        system_iterator<double> *test = new dlma_iterator<double>(argv[1]);;
        test->run_system();
        //test->create_movie_files(argv[2]);

        if (argc < 3)
            test->save_config_file();
        else
            test->save_config_file(argv[2]);

        delete test;
    }

    else if((sys_name == "random_site_percolation") && (lattice_val == 1)){
        system_iterator<int> *test = new dlma_iterator<int>(argv[1]);
        test->run_system_for_percolation();

        if (argc < 3)
            test->save_config_file();
        else
            test->save_config_file(argv[2]);

        delete test;
    }

    else if(sys_name == "erdos_renyi"){
        system_iterator<double> *test = new dlma_iterator<double>(argv[1]);
        test->run_system_for_erdos_renyi();

        if (argc < 3)
            test->save_config_file();

        else
            test->save_config_file(argv[2]);

        delete test;
    }

    else{
        std::cout<<"please check system and lattice values"<<std::endl;
        exit(EXIT_FAILURE);
    }

    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    return 0;

}
