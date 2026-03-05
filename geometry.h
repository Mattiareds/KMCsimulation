#include <vector>
#include <fstream>

#ifndef GEOMETRY_H
#define GEOMETRY_H

class metropolis;
class coordinates;

class geometry{
    private:
    coordinates* st = nullptr;
    metropolis* mt = nullptr;
    std::vector<std::vector<double>> planes;
    std::vector<std::vector<int>> info_plane_sites; //how many and in which planes is each site
    std::vector<int> type_of_plane; //0=100c 1=111c 2=two111 3=one111one100 4=two100 5=two111one100 6=two111one100 7=100b 8=111b
    std::vector<std::vector<int>> table_of_nn;
    void make_planes();
    void site_characterisation();
    void initialize_nn();

    public:
    
    void set(coordinates& structure , metropolis& algorithm){ 
        st = &structure; 
        mt = &algorithm; 
    }


    
    bool test_piano(int i,int j);
    void test_border(geometry& c);
    void initialize_external_nn(coordinates& core);
    void pv_substitution();

    //initializer
    void starter();

    // Getter and setter for all
    const std::vector<int>& get_table_of_nn(int idx) const { return table_of_nn.at(idx); }
    const int table_of_nn_size() const { return table_of_nn.size(); } //getter alla dimensione di table_of_nn
    const int get_TOPl(int idx) const { return type_of_plane.at(idx); }
    const int get_plane(int i, int j) { return info_plane_sites[i][j]; }


};

#endif