#ifndef HEIS_MAG_HPP
#define HEIS_MAG_HPP

#include <aiwlib/vec>
#include "COP.hpp"
#include "spec_func.hpp"
#include "tensor_network.hpp"
#include <ctime>
#include <aiwlib/sphere>

#define sphere

template <typename T>
class heis_mag{
private:
    int n = 0;
    std::vector<aiw::Vec<3, T>> p;
    std::map<std::set<int>, T> lambda;
    int l_max = 2;
    int N_teta = 22;
    int N_fi = 22;
    std::vector<std::vector<int>> pows;
    TensorNetwork<std::complex<T>> TN;
    int ind_sz;
    std::map<int, std::set<int>> connections;
    COP cop;
    std::vector<std::map<std::set<int>, int>> order;
    std::map<std::set<int>, std::vector<T>> infeld;
    Tensor<T> result;
    int rank;
    int sz_sph;
public:
    void set_lambda(int i, int j, T lambda_){
        std::set<int> s; s.insert(i); s.insert(j);
        lambda[s] = lambda_;
        std::set<int> & s1 = connections[i]; std::set<int> & s2 = connections[j];
        s1.insert(j); s2.insert(i);
        SpecFunc sf;
        std::vector<T> v(l_max+1);
        for(int i = 0; i <= l_max; i++){
            v[i] = sf.Infeld_func(lambda_, i+0.5, 0.001);
        } 
        infeld[s] = v;
    }

    void set_p(int i, const aiw::Vec<3, T> & v_){
        p[i] = v_;
    }

    std::vector<int> get_tensor_marks(int ID){
        return (TN.get_tensor(ID)).get_marks();
    }

    void init(int n_, int rank_, int sz_sphere_){
        n = n_;
        rank = rank_;
        sz_sph = sz_sphere_;
        p.resize(n);
        pows.resize(n);
        ind_sz = l_max*(l_max+2)+1;
        cop.init_Legendre_Pol(l_max);
        cop.calc_Legendre_Pol_derivatives();
        order.resize(n);
    }

    void set_pow(int node, const std::vector<int> & pow_){
        pows[node] = pow_;
    }

    void calc_lm(int ind, int & l, int & m){
        int count = 1; l = 0; m = 0;
        for(int i = 0; i < l_max+1; i++){
            if(ind <= count-1){
                m = l-(count-1-ind);
                break;
            }else{
                count += 2*(l+1)+1;
                l++;
            }
        }
    }

    Tensor<std::complex<T>> calc_tensor(int i){
        std::map<std::set<int>, int> & ord = order[i];
        Tensor<std::complex<T>> res;
        // T d_teta = M_PI/N_teta; T d_fi = 2*M_PI/N_fi;
        std::set<int> & c = connections[i];
        int k = c.size(); int pp = pows[i].size();
        std::vector<int> tens_ids;
        for(int j : c){
            tens_ids.push_back(j);
        }
        std::sort(tens_ids.begin(), tens_ids.end());
        std::vector<int> sz_(k+pp);
        for(int i = 0; i < k; i++){
            sz_[i] = ind_sz;
        }
        for(int i = 0; i < pp; i++){
            sz_[k+i] = 3;
        }
        res.init(sz_);

        #ifdef no_sphere
        for(int i_teta = 0; i_teta < N_teta; i_teta++){
            T teta = i_teta*d_teta;
            for(int i_fi = 0; i_fi < N_fi; i_fi++){
                T fi = i_fi*d_fi;
                aiw::Vec<3, T> M; M[0] = sin(teta)*cos(fi); M[1] = sin(teta)*sin(fi); M[2] = cos(teta);
                std::complex<T> f = sin(teta)*exp(M*p[i]);
        #endif

        #ifdef sphere
            for(int i_sph=0; i_sph<sz_sph; i_sph++){
		        const aiw::Vec<3, T> &M = aiw::sph_vert(i_sph, rank); T dM = aiw::sph_vert_area(i_sph, rank);
                T teta = acos(M[2]); T fi = 0;
                if(std::abs(teta) < 1e-8){
                    fi = 0;
                }else{
                    T mc = M[0]/sin(teta);
                    T ms = M[1]/sin(teta);
                    if(abs(mc) > 1){
                        if(mc > 0){
                            fi = 0;
                        }else{
                            fi = M_PI;
                        }
                    }else{
                        fi = acos(mc);
                    }
                    if(ms < 0){
                        fi = 2*M_PI-fi;
                    }
                }
                std::complex<T> f = exp(M*p[i]);
        #endif
                typename Tensor<std::complex<T>>::iterator it = res.begin();
                typename Tensor<std::complex<T>>::iterator it_end = res.end();
                
                while(it != it_end){
                    std::vector<int> v = it.get_pos();
                    std::complex<T> r = f;

                    for(int ind = 0; ind < k; ind++){
                        int ind_i = tens_ids[ind];
                        std::set<int> s; s.insert(i); s.insert(ind_i);
                        #ifdef no_sphere
                        if(i_teta == 0 && i_fi == 0){
                        #endif
                        #ifdef sphere
                        if(i_sph == 0){
                        #endif  
                            ord[s] = ind;
                        }
                        T lmbd = lambda[s];
                        bool b = (ind_i > i);
                        int l, m; calc_lm(v[ind], l, m);
                        r *= b? cop.Y(l, m, teta, fi)*infeld[s][l]*pow(2*M_PI, 1.5)/sqrt(lmbd):conj(cop.Y(l, m, teta, fi));
                    }
                    for(int ind = 0; ind < pp; ind++){
                        int pow_ = pows[i][ind];
                        int curr_ind = v[k+ind];
                        r *= pow(M[curr_ind], pow_);
                    }
                    #ifdef no_sphere
                    res[v] += r*d_fi*d_teta;
                    #endif
                    #ifdef sphere
                    res[v] += r*dM;
                    #endif
                    it++;
                }
            }

        #ifdef no_sphere
        }
        #endif

        return res;
    }

    void build_TN(){
        std::vector<Tensor<std::complex<T>>> tensors;
        for(int node = 0; node < n; node++){
            tensors.push_back(calc_tensor(node));
        }

        for(int j = 0; j < n; j++){
            // std::map<std::set<int>, int> & ord_j = order[j];
            // std::map<int, TensorNet::link> links;
            // for(int i = 0; i < j; i++){
            //     std::set<int> s; s.insert(i); s.insert(j);
            //     std::map<std::set<int>, int> & ord_i = order[i];
            //     int indj = ord_j[s]; int indi = ord_i[s];
            //     TensorNet::link & l = links[i]; 
            //     l.glue1 = std::vector<int>{indi}; 
            //     l.glue2 = std::vector<int>{indj}; 
            //     l.n=1;
            // }
            // TN.add_tensor(tensors[j], links);
            include_tensor(j, tensors[j]);
        }
    }

    void include_tensor(int j, const Tensor<std::complex<T>> & t, bool replace = false){
        std::map<std::set<int>, int> & ord_j = order[j];
        std::map<int, TensorNet::link> links;
        for(int i = 0; i < j; i++){
            std::set<int> s; s.insert(i); s.insert(j);
            std::map<std::set<int>, int> & ord_i = order[i];
            int indj = ord_j[s]; int indi = ord_i[s];
            TensorNet::link & l = links[i]; 
            l.glue1 = std::vector<int>{indi}; 
            l.glue2 = std::vector<int>{indj}; 
            l.n=1;
        }
        if(replace){
            TN.replace_tensor(t, j);
        }else{
            TN.add_tensor(t, links); 
        }  
    }

    void calc(){
        TN.convolve_Network();
        Tensor<std::complex<T>> res = TN.get_result();
        result.init(res.get_size());
        std::vector<int> marks_ = res.get_marks();
        for(int i = 0; i < res.get_D(); i++){
            result.set_mark(i, marks_[i]);
        }
        int sz = res.get_sz();
        for(int i = 0; i < sz; i++){
            result[i] = std::real(res[i]);
        }
    }

    void add_index(int ID, int num, int mark = -1){
        for(int i = 0; i < num; i++){
            pows[ID].push_back(1);
        }
        std::vector<int> marks_ = get_tensor_marks(ID);
        Tensor<std::complex<T>> t = calc_tensor(ID);
        int k = connections[ID].size()-1;
        int addition = pows[ID].size();
        for(int i = 0, sz_ = marks_.size(); i < sz_; i++){
            t.set_mark(i, marks_[i]);
        }
        if(mark != -1){
            t.set_mark(k+addition, mark);
        }
        include_tensor(ID, t, true);
    }

    Tensor<T> get_result(){
        return result;
    }

    template <typename other_Type> heis_mag<T>& operator = (heis_mag<other_Type> & h){
        n = h.n;
        p = h.p;
        lambda = h.lambda;
        l_max = h.l_max;
        N_teta = h.N_teta;
        N_fi = h.N_fi;
        pows = h.pows;
        TN = h.TN;
        ind_sz = h.ind_sz;
        connections = h.connections;
        cop = h.cop;
        order = h.order;
        result = h.result;
        infeld = h.infeld;
        rank = h.rank;
        sz_sph = h.sz_sph;
        return *this;
    }
};

template<typename T>
struct moments3{
    int rank, sz_sphere;

    heis_mag<T> Z_t, m1_t, m2_t, m3_t, m1m2_t, m1m3_t, m2m3_t, 
    m1m2m3_t, m1m1m2_t, m1m1m3_t, m2m2m1_t, m2m2m3_t, m3m3m1_t, m3m3m2_t;

    T Z, eta_12, eta_13, eta_23;
    aiw::Vec<3, T> m1eta_12, m1eta_13, m1eta_23, m2eta_12, m2eta_13, m2eta_23, m3eta_12, m3eta_13, m3eta_23, m1, m2, m3;

    T lmbd12, lmbd13, lmbd23;
    aiw::Vec<3, T> p1, p2, p3;

    double t_init = 0, t_calc = 0;
    void calc(){
        clock_t start1 = clock();
        Z_t.init(3, rank, sz_sphere);
        Z_t.set_lambda(0, 1, lmbd12);
        Z_t.set_lambda(1, 2, lmbd23);
        Z_t.set_lambda(2, 0, lmbd13);
        Z_t.set_p(0, p1);
        Z_t.set_p(1, p2);
        Z_t.set_p(2, p3);

        Z_t.build_TN(); 
        m1_t = Z_t;
        m2_t = Z_t; 
        m3_t = Z_t; 
        
        m1_t.add_index(0, 1, 0);
        m2_t.add_index(1, 1, 1);
        m3_t.add_index(2, 1, 2);
        m1m2_t = m1_t;
        m1m3_t = m1_t;
        m2m3_t = m2_t;
        m1m2_t.add_index(1, 1, 1);
        m1m3_t.add_index(2, 1, 2);
        m2m3_t.add_index(2, 1, 2);

        m1m1m2_t = m1m2_t;
        m1m1m2_t.add_index(0, 1, 0);

        m1m1m3_t = m1m3_t;
        m1m1m3_t.add_index(0, 1, 0);

        m2m2m1_t = m1m2_t;
        m2m2m1_t.add_index(1, 1, 1);

        m2m2m3_t = m2m3_t;
        m2m2m3_t.add_index(1, 1, 1);

        m3m3m1_t = m1m3_t;
        m3m3m1_t.add_index(2, 1, 2);

        m3m3m2_t = m2m3_t;
        m3m3m2_t.add_index(2, 1, 2);

        m1m2m3_t = m2m3_t;
        m1m2m3_t.add_index(0, 1, 0);
        clock_t end1 = clock();
        t_init += double(end1 - start1) / CLOCKS_PER_SEC;

        Z_t.calc(); Z = (Z_t.get_result())[0];

        m1_t.calc();
        m2_t.calc();
        m3_t.calc();
        m1m2_t.calc();
        m1m3_t.calc();
        m2m3_t.calc();
        m1m2m3_t.calc();
        m1m1m2_t.calc();
        m1m1m3_t.calc();
        m2m2m1_t.calc();
        m2m2m3_t.calc();
        m3m3m1_t.calc();
        m3m3m2_t.calc();

        Tensor<T> m1m2 = m1m2_t.get_result(); 
        Tensor<T> m1m3 = m1m3_t.get_result(); 
        Tensor<T> m2m3 = m2m3_t.get_result(); 
        Tensor<T> m1m2m3 = m1m2m3_t.get_result();
        Tensor<T> m1m1m2 = m1m1m2_t.get_result();
        Tensor<T> m1m1m3 = m1m1m3_t.get_result();
        Tensor<T> m2m2m1 = m2m2m1_t.get_result();
        Tensor<T> m2m2m3 = m2m2m3_t.get_result();
        Tensor<T> m3m3m1 = m3m3m1_t.get_result();
        Tensor<T> m3m3m2 = m3m3m2_t.get_result();

        int ind1 = m1m2m3.mark_to_ind(0); int ind2 = m1m2m3.mark_to_ind(1); int ind3 = m1m2m3.mark_to_ind(2);

        Tensor<T> m1eta_23_conv = tensor_convolution(m1m2m3, std::vector<int>{ind2}, std::vector<int>{ind3});
        m1eta_23[0] = m1eta_23_conv[0]; m1eta_23[1] = m1eta_23_conv[1]; m1eta_23[2] = m1eta_23_conv[2];
        m1eta_23 /= Z;

        Tensor<T> m2eta_13_conv = tensor_convolution(m1m2m3, std::vector<int>{ind1}, std::vector<int>{ind3});
        m2eta_13[0] = m2eta_13_conv[0]; m2eta_13[1] = m2eta_13_conv[1]; m2eta_13[2] = m2eta_13_conv[2];
        m2eta_13 /= Z;

        Tensor<T> m3eta_12_conv = tensor_convolution(m1m2m3, std::vector<int>{ind1}, std::vector<int>{ind2});
        m3eta_12[0] = m3eta_12_conv[0]; m3eta_12[1] = m3eta_12_conv[1]; m3eta_12[2] = m3eta_12_conv[2];
        m3eta_12 /= Z;

        ind1 = m1m1m2.mark_to_ind(0); ind2 = m1m1m2.mark_to_ind(1);

        Tensor<T> m1eta_12_conv = tensor_convolution(m1m1m2, std::vector<int>{ind1}, std::vector<int>{ind2});
        m1eta_12[0] = m1eta_12_conv[0]; m1eta_12[1] = m1eta_12_conv[1]; m1eta_12[2] = m1eta_12_conv[2];
        m1eta_12 /= Z;

        ind1 = m1m1m3.mark_to_ind(0); ind3 = m1m1m3.mark_to_ind(2);

        Tensor<T> m1eta_13_conv = tensor_convolution(m1m1m3, std::vector<int>{ind1}, std::vector<int>{ind3});
        m1eta_13[0] = m1eta_13_conv[0]; m1eta_13[1] = m1eta_13_conv[1]; m1eta_13[2] = m1eta_13_conv[2];
        m1eta_13 /= Z;

        ind1 = m2m2m1.mark_to_ind(0); ind2 = m2m2m1.mark_to_ind(1); 

        Tensor<T> m2eta_12_conv = tensor_convolution(m2m2m1, std::vector<int>{ind1}, std::vector<int>{ind2});
        m2eta_12[0] = m2eta_12_conv[0]; m2eta_12[1] = m2eta_12_conv[1]; m2eta_12[2] = m2eta_12_conv[2];
        m2eta_12 /= Z;

        ind3 = m2m2m3.mark_to_ind(2); ind2 = m2m2m3.mark_to_ind(1);

        Tensor<T> m2eta_23_conv = tensor_convolution(m2m2m3, std::vector<int>{ind2}, std::vector<int>{ind3});
        m2eta_23[0] = m2eta_23_conv[0]; m2eta_23[1] = m2eta_23_conv[1]; m2eta_23[2] = m2eta_23_conv[2];
        m2eta_23 /= Z;

        ind1 = m3m3m1.mark_to_ind(0); ind3 = m3m3m1.mark_to_ind(2); 

        Tensor<T> m3eta_13_conv = tensor_convolution(m3m3m1, std::vector<int>{ind1}, std::vector<int>{ind3});
        m3eta_13[0] = m3eta_13_conv[0]; m3eta_13[1] = m3eta_13_conv[1]; m3eta_13[2] = m3eta_13_conv[2];
        m3eta_13 /= Z;

        ind2 = m3m3m2.mark_to_ind(1); ind3 = m3m3m2.mark_to_ind(2); 

        Tensor<T> m3eta_23_conv = tensor_convolution(m3m3m2, std::vector<int>{ind2}, std::vector<int>{ind3});
        m3eta_23[0] = m3eta_23_conv[0]; m3eta_23[1] = m3eta_23_conv[1]; m3eta_23[2] = m3eta_23_conv[2];
        m3eta_23 /= Z;

        eta_12 = (m1m2[0]+m1m2[4]+m1m2[8])/Z;
        eta_13 = (m1m3[0]+m1m3[4]+m1m3[8])/Z;
        eta_23 = (m2m3[0]+m2m3[4]+m2m3[8])/Z;

        Tensor<T> m1_res = m1_t.get_result(); 
        m1[0] = m1_res[0]/Z; m1[1] = m1_res[1]/Z; m1[2] = m1_res[2]/Z;

        Tensor<T> m2_res = m2_t.get_result(); 
        m2[0] = m2_res[0]/Z; m2[1] = m2_res[1]/Z; m2[2] = m2_res[2]/Z;

        Tensor<T> m3_res = m3_t.get_result(); 
        m3[0] = m3_res[0]/Z; m3[1] = m3_res[1]/Z; m3[2] = m3_res[2]/Z;

        clock_t end2 = clock();
        t_calc += double(end1 - end2) / CLOCKS_PER_SEC;
    }
};


#endif //HEIS_MAG_HPP