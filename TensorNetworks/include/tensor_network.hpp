#ifndef TENSOR_NETWORK_HPP
#define TENSOR_NETWORK_HPP

#include "tensor.hpp"
#include <map>
#include <set>

namespace TensorNet{
struct link{
    link(){};
    int t1_id, t2_id; //id тензоров, t1_id <= t2_id
    int n; //число связей тензоров
    std::vector<int> glue1, glue2; //соединения индексов
};
}

using TensorNet::link;

void find_correspondence(const std::vector<int> & glue1, int n1, std::map<int, int> & table1, const std::vector<int> & glue2, int n2, std::map<int, int> & table2){
    std::vector<int> v1, v2;
    for(int i = 0; i < n1; i++){
        v1.push_back(i);
    }
    for(int i = 0; i < n2; i++){
        v2.push_back(i);
    }
    int n = n1;
    int g1 = glue1.size();
    int sep = n - g1;
    std::vector<int> g = glue_vectors(glue1, v1, glue2, v2);
    for(int i = 0; i < sep; i++){
        table1[g[i]] = i;
    }
    for(int i = sep, sz_ = sep+v2.size()-g1; i < sz_; i++){
        table2[g[i]] = i;
    }
}

template <typename T>
class TensorNetwork{
    std::map<int, Tensor<T>> tensors;
    std::map<std::set<int>, TensorNet::link> links;
    std::map<int, std::set<int>> connections;
    int last_id = 0;
    Tensor<T> result;
public:
    // в link в links_ glue1 должен соответствовать тензору уже имеющемуся в сети, а glue2 - новому
    // если тензоры в link совпадают - не важно
    void add_tensor(const Tensor<T> & t, std::map<int, TensorNet::link> & links_){
        tensors[last_id] = t;
        last_id++;
        for(std::pair <int, TensorNet::link> p : links_){
            TensorNet::link l = p.second;
            l.t2_id = last_id-1; l.t1_id = p.first;
            std::set<int> & s1 = connections[l.t1_id]; std::set<int> & s2 = connections[l.t2_id];
            s1.insert(l.t2_id); s2.insert(l.t1_id);
            std::set<int> s; s.insert(l.t1_id); s.insert(l.t2_id);
            links[s] = l;
        }
    }

    void replace_tensor(const Tensor<T> & t, int ID){
        tensors[ID] = t;
    }

    void one_step(int t1_id, int t2_id){
        if(t1_id > t2_id){
            int temp = t1_id; t1_id = t2_id; t2_id = temp;
        }
        std::set<int> s; s.insert(t1_id); s.insert(t2_id);
        TensorNet::link & l = links[s];
        Tensor<T> & t1 = tensors[t1_id]; Tensor<T> & t2 = tensors[t2_id];
        Tensor<T> conv_tensor = tensor_convolution(t1, l.glue1, t2, l.glue2);

        std::map<int, int> table1, table2;
        find_correspondence(l.glue1, t1.get_D(), table1, l.glue2, t2.get_D(), table2);

        if(t1_id != t2_id){
            tensors.erase(t1_id); tensors.erase(t2_id);
        }else{
            tensors.erase(t1_id);
        }
        
        tensors[last_id++] = conv_tensor;
        const std::set<int> s1 = connections[t1_id]; const std::set<int> s2 = connections[t2_id];
        for(int i : s1){
            if(i != t2_id){
                std::set<int> s; s.insert(i); s.insert(t1_id);
                TensorNet::link & l = links[s];
                if(i == l.t2_id){
                    l.glue1.swap(l.glue2);
                    l.t1_id = l.t2_id;
                }
                l.t2_id = last_id-1;
                std::vector<int> new_glue;
                for(int j : l.glue2){
                    new_glue.push_back(table1[j]); 
                }
                l.glue2 = new_glue;
                std::set<int> s_; s_.insert(l.t1_id); s_.insert(l.t2_id);
                links[s_] = l;  
                links.erase(s);
            }
        }

        for(int i : s1){
            connections[i].erase(t1_id);
            if(i != t2_id){
                connections[i].insert(last_id-1);
                connections[last_id-1].insert(i);
            }
        }
        std::set<int> empty;
        connections[t1_id] = empty;

        if(t1_id != t2_id){
            for(int i : s2){
                if(i != t1_id){
                    std::set<int> s; s.insert(i); s.insert(t2_id);
                    TensorNet::link & l = links[s];
                    int ind1,ind2;
                    ind1 = l.t1_id;
                    if(i == l.t2_id){
                        l.glue1.swap(l.glue2);
                        ind1 = l.t2_id;
                    }
                    ind2 = last_id-1;
                    std::set<int> s_; s_.insert(ind1); s_.insert(ind2);
                    std::map<std::set<int>, TensorNet::link>::iterator it = links.find(s_);
                    bool b = (it!=links.end());

                    std::vector<int> new_glue;
                    for(int j : l.glue2){
                        new_glue.push_back(table2[j]); 
                    }

                    if(b){
                        TensorNet::link & l_ = links[s_];
                        int n_ = l.glue1.size();
                        for(int j = 0; j < n_; j++){
                            l_.glue1.push_back(l.glue1[j]);
                            l_.glue2.push_back(new_glue[j]);
                        }
                        l_.n += n_;
                        links.erase(s);
                    }else{
                        l.glue2 = new_glue;
                    }

                }
            }

            links.erase(s);

            for(int i : s2){
                connections[i].erase(t2_id);
                if(i != t1_id){
                    connections[i].insert(last_id-1);
                    connections[last_id-1].insert(i);
                }
            }
            connections[t2_id] = empty;
        }

        revision();
    }

    void revision(){
        std::vector<int> for_erase;
        for(std::pair<int, std::set<int>> p : connections){
            if(p.second.size() == 0){
                for_erase.push_back(p.first);
            }
        }
        for(int i : for_erase){
            connections.erase(i);
        }
    }

    Tensor<T> & get_tensor(int ID){
        return tensors[ID];
    }

    void convolve_Network(){
        while(connections.size() > 1){
            int id1 = -1; int id2 = -1;
            for(std::pair<int, std::set<int>> p : connections){
                for(int j : p.second){
                    id2 = j;
                    break;
                }
                if(id2 != -1){
                    id1 = p.first;
                }
            }
            one_step(id1, id2);
        }
        result = tensors[last_id-1];
    }

    Tensor<T> get_result(){
        return result;
    }

    template <typename other_Type> TensorNetwork<T>& operator = (const TensorNetwork<other_Type> &t){
        tensors = t.tensors;
        links = t.links;
        connections = t.connections;
        last_id = t.last_id;
        result = t.result;
        return *this;
    }
};

#endif //TENSOR_NETWORK_HPP