#ifndef HEIS_MAG_HPP
#define HEIS_MAG_HPP

#include <aiwlib/vec>
#include "COP.hpp"
#include "spec_func.hpp"
#include "tensor_network.hpp"

template <typename T>
class heis_mag{
private:
    int n = 0;
    std::vector<aiw::Vec<3, T>> p;
    std::map<std::set<int>, T> lambda;
    int l_max = 2;
    int N_teta = 20;
    int N_fi = 20;
    std::vector<std::vector<int>> pows;
    TensorNetwork<std::complex<T>> TN;
    int ind_sz;
    std::map<int, std::set<int>> connections;
    COP cop;
    std::vector<std::map<std::set<int>, int>> order;
    std::map<std::set<int>, std::vector<T>> infeld;
    Tensor<T> result;
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

    void init(int n_){
        n = n_;
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

    std::complex<T> calc_integral(int i, const std::vector<int> & v){
        std::map<std::set<int>, int> & ord = order[i];
        std::complex<T> res{0, 0};
        T d_teta = M_PI/N_teta; T d_fi = 2*M_PI/N_fi;
        std::set<int> & c = connections[i];
        int k = c.size(); int pp = pows[i].size();
        std::vector<int> tens_ids;
        for(int j : c){
            tens_ids.push_back(j);
        }
        std::sort(tens_ids.begin(), tens_ids.end());
        // SpecFunc sf;

        for(int i_teta = 0; i_teta < N_teta; i_teta++){
            T teta = i_teta*d_teta;
            for(int i_fi = 0; i_fi < N_fi; i_fi++){
                T fi = i_fi*d_fi;
                aiw::Vec<3, T> M; M[0] = sin(teta)*cos(fi); M[1] = sin(teta)*sin(fi); M[2] = cos(teta);
                std::complex<T> contribution = sin(teta)*exp(M*p[i]);
                for(int ind = 0; ind < k; ind++){
                    int ind_i = tens_ids[ind];
                    std::set<int> s; s.insert(i); s.insert(ind_i);
                    if(i_teta == 0 && i_fi == 0){
                        ord[s] = ind;
                    }
                    T lmbd = lambda[s];
                    bool b = (ind_i > i);
                    int l, m; calc_lm(v[ind], l, m);
                    contribution *= b? cop.Y(l, m, teta, fi)*infeld[s][l]*pow(2*M_PI, 1.5)/sqrt(lmbd):conj(cop.Y(l, m, teta, fi));
                }
                for(int ind = 0; ind < pp; ind++){
                    int pow_ = pows[i][ind];
                    int curr_ind = v[k+ind];
                    contribution *= pow(M[curr_ind], pow_);
                }
                res += contribution;
            }
        }

        res *= d_fi*d_teta;
        return res;
    }

    void build_TN(){
        std::vector<Tensor<std::complex<T>>> tensors(n);
        for(int node = 0; node < n; node++){
            Tensor<std::complex<T>> & t = tensors[node];
            int k = connections[node].size(); int pp = pows[node].size();
            std::vector<int> sz_(k+pp);
            for(int i = 0; i < k; i++){
                sz_[i] = ind_sz;
            }
            for(int i = 0; i < pp; i++){
                sz_[k+i] = 3;
            }
            t.init(sz_);
            typename Tensor<std::complex<T>>::iterator it = t.begin();
            typename Tensor<std::complex<T>>::iterator it_end = t.end();
            while(it != it_end){
                const std::vector<int> & v = it.get_pos();
                t[v] = calc_integral(node, v);
                it++;
            }
        }

        for(int j = 0; j < n; j++){
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
            TN.add_tensor(tensors[j], links);
        }
    }

    void calc(){
        TN.convolve_Network();
        Tensor<std::complex<T>> res = TN.get_result();
        result.init(res.get_size());
        int sz = res.get_sz();
        for(int i = 0; i < sz; i++){
            result[i] = std::real(res[i]);
        }
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
        return *this;
    }
};

template<typename T>
struct moments3{
    heis_mag<T> Z_t, m1_t, m2_t, m3_t, m1m2_t, m1m3_t, m2m3_t, 
    m1m2m3_1, m1m2m3_2, m1m2m3_3, m1m1m2, m1m1m3, m2m2m1, m2m2m3, m3m3m1, m3m3m2;

    T Z, eta_12, eta_13, eta_23;
    aiw::Vec<3, T> m1eta_12, m1eta_13, m1eta_23, m2eta_12, m2eta_13, m2eta_23, m3eta_12, m3eta_13, m3eta_23, m1, m2, m3;

    T lmbd12, lmbd13, lmbd23;
    aiw::Vec<3, T> p1, p2, p3;

    void calc(){
        Z_t.init(3);
        Z_t.set_lambda(0, 1, lmbd12);
        Z_t.set_lambda(1, 2, lmbd23);
        Z_t.set_lambda(2, 0, lmbd13);
        Z_t.set_p(0, p1);
        Z_t.set_p(1, p2);
        Z_t.set_p(2, p3);
        m1_t = Z_t; m2_t = Z_t; m3_t = Z_t; m1m2_t = Z_t; m1m3_t = Z_t; m2m3_t = Z_t;
        m1m2m3_1 = Z_t; m1m2m3_2 = Z_t; m1m2m3_3 = Z_t; m1m1m2 = Z_t; m1m1m3 = Z_t; m2m2m1 = Z_t; m2m2m3 = Z_t; m3m3m1 = Z_t; m3m3m2 = Z_t;
        Z_t.build_TN(); Z_t.calc(); Z = (Z_t.get_result())[0];

        std::vector<int> pow1{1}, pow2{1,1};
        m1m2_t.set_pow(0, pow1); m1m2_t.set_pow(1, pow1);
        m1m2_t.build_TN(); m1m2_t.calc(); Tensor<T> m1m2 = m1m2_t.get_result(); 
        eta_12 = (m1m2[0]+m1m2[4]+m1m2[8])/Z;

        m1_t.set_pow(0, pow1);
        m1_t.build_TN(); m1_t.calc(); Tensor<T> m1_res = m1_t.get_result(); m1[0] = m1_res[0]/Z; m1[1] = m1_res[1]/Z; m1[2] = m1_res[2]/Z;
    }
};


#endif //HEIS_MAG_HPP