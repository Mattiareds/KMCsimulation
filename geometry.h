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
    double** planes;
    std::vector<std::vector<int>> info_plane_sites; //how many and in which planes is each site
    std::vector<int> type_of_plane; //0=100c 1=111c 2=two111 3=one111one100 4=two100 5=two111one100 6=two111one100 7=100b 8=111b
    std::vector<std::vector<int>> table_of_nn;
    std::vector<int> upper_sites;
    std::vector<std::vector<int>> edge_map;
    void make_planes();
    void site_characterisation();
    void initialize_nn();
    // bool is_same_plane(int site_i, int site_f);

    public:
    
    void set(coordinates& structure , metropolis& algorithm){ 
        st = &structure; 
        mt = &algorithm; 
    }
    
    bool test_piano(int i,int j);
    void test_border(geometry& c);
    void initialize_external_nn(coordinates& upper);
    void pv_substitution();

    //initializer
    void starter();
    void starter_high();

    //setter
    void set_table_of_nn(int i, int j, const int nn) {
        table_of_nn.at(i).at(j) = nn;
    }

    void append_nn(int i, const int nn){
        table_of_nn.at(i).push_back(nn);
        edge_map.at(i).push_back(nn);
    }

    void set_TOPl(int site, int plane){
        type_of_plane.at(site) = plane;
    }

    void reset_edge_map(int idx, int jdx){
        edge_map.at(idx).at(jdx) = -1;
    }

    void set_plane(int idx, int plane){
        info_plane_sites[idx][1] = plane;
    }

    // Getter and setter for all
    const std::vector<int>& get_table_of_nn(int idx) const { return table_of_nn.at(idx); }
    const std::vector<int>& get_upper_sites() const { return upper_sites;}
    const std::vector<int>& get_edge_map(int idx) const { return edge_map.at(idx);}
    const int& get_edge_site(int idx, int jdx) const { return edge_map.at(idx).at(jdx);}
    const int table_of_nn_size() const { return table_of_nn.size(); } //getter alla dimensione di table_of_nn
    const int get_TOPl(int idx) const { return type_of_plane.at(idx); }
    const int get_plane(int i, int j) { return info_plane_sites[i][j]; }


};

#endif