#include "../include/heis_mag.hpp"
#include <ctime>
#include <aiwlib/sphere>

template<typename T>
struct mom_m{
    int rank, sz_sphere;
    heis_mag<T> Z_t, m1_t;

    T lmbd12, lmbd13, lmbd23;
    aiw::Vec<3, T> p1, p2, p3;

    T Z, m;
    aiw::Vec<3,T> m1;
    double t_calc = 0;
    int l_max;

    void calc(){
        clock_t start = clock();
        Z_t.l_max = l_max;
        Z_t.init(3, rank, sz_sphere);
        Z_t.set_lambda(0, 1, lmbd12);
        Z_t.set_lambda(1, 2, lmbd23);
        Z_t.set_lambda(2, 0, lmbd13);
        Z_t.set_p(0, p1);
        Z_t.set_p(1, p2);
        Z_t.set_p(2, p3);

        Z_t.build_TN(); 
        m1_t = Z_t;
        m1_t.add_index(0, 1, 0);
        Z_t.calc(); Z = (Z_t.get_result())[0];
        m1_t.calc();

        Tensor<T> m_res = m1_t.get_result(); 
        m1[0] = m_res[0]/Z; m1[1] = m_res[1]/Z; m1[2] = m_res[2]/Z;
        m = sqrt(m1[0]*m1[0]+m1[1]*m1[1]+m1[2]*m1[2]);
        clock_t end = clock();
        t_calc += double(end - start) / CLOCKS_PER_SEC;
    }
};

mom_m<double> calc_directly(aiw::Vec<3, double> p1, aiw::Vec<3, double> p2, aiw::Vec<3, double> p3, double lmbd12, double lmbd13, double lmbd23, int sz_sphere, int rank){
    mom_m<double> moms;

    double Z = 0;
    aiw::Vec<3, double> m1;
    double m;

    for(int i=0; i<sz_sphere; i++){
		const aiw::Vec<3> &m_i = aiw::sph_vert(i, rank); double dm_i = aiw::sph_vert_area(i, rank);
        for(int j=0; j<sz_sphere; j++){
		    const aiw::Vec<3> &m_j = aiw::sph_vert(j, rank); double dm_j = aiw::sph_vert_area(j, rank);
            for(int k=0; k<sz_sphere; k++){
		        const aiw::Vec<3> &m_k = aiw::sph_vert(k, rank); double dm_k = aiw::sph_vert_area(k, rank);
                double f = exp(p1*m_i+p2*m_j+p3*m_k+lmbd12*m_i*m_j+lmbd13*m_i*m_k+lmbd23*m_j*m_k);
                double dm = dm_i*dm_j*dm_k;
                Z += f*dm;
                m1 += f*m_i*dm;
            }
        }
        std::cout << i << "\n";
    }

    m1 /= Z;
    m = sqrt(m1[0]*m1[0]+m1[1]*m1[1]+m1[2]*m1[2]);
    moms.m = m;
    moms.Z = Z;
    return moms;
}

int main(){
    aiw::Vec<3, double> p1{1, 2, 3}, p2{4., 0, 6}, p3{7, 8, 0};
    p1 *= 1./p1.abs(); p2 *= 1./p2.abs(); p3 *= 1./p3.abs();
    double lmbd_12{1}, lmbd_13{2}, lmbd_23{3};

    int rank_begin = 0; int rank_end = 2;
    int l_min = 0; int l_max = 3;
    aiw::sph_init_table(rank_end+1);	
    int sz_sph = aiw::sph_vertex_num(rank_end+1);
    mom_m<double> exact = calc_directly(p1, p2, p3, lmbd_12, lmbd_13, lmbd_23, sz_sph, rank_end+1);
    double exact_m = exact.m;

    std::cout << "s l delta_1 t_1 delta_2 t_2\n";
    for(int r = rank_begin; r <= rank_end; r++){
        aiw::sph_init_table(r);	
        
        int sz_sph = aiw::sph_vertex_num(r);
        clock_t start = clock();
        mom_m<double> mom_direct = calc_directly(p1, p2, p3, lmbd_12, lmbd_13, lmbd_23, sz_sph, r);
        clock_t end = clock();
        double t2 = double(end - start) / CLOCKS_PER_SEC;
        

        for(int l = l_min; l <= l_max; l++){
            mom_m<double> mom;
            mom.lmbd12 = lmbd_12; mom.lmbd13 = lmbd_13; mom.lmbd23 = lmbd_23;
            mom.p1 = p1; mom.p2 = p2; mom.p3 = p3;
            mom.rank = r;
            mom.sz_sphere = sz_sph;
            mom.l_max = l;
            mom.calc();
            std::cout << sz_sph << " " << l << " " << mom.m-exact_m << " " << mom.t_calc  << " " << mom_direct.m-exact_m << " " << t2 << "\n";
        }
    }

    return 0;
}