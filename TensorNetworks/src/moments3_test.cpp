#include "../include/heis_mag.hpp"
#include <ctime>
#include <aiwlib/sphere>

#define calc_m1
#define calc_eta_23

const int rank = 1;

moments3<double> calc_directly(aiw::Vec<3, double> p1, aiw::Vec<3, double> p2, aiw::Vec<3, double> p3, double lmbd12, double lmbd13, double lmbd23, int sz_sphere){
    moments3<double> moms;

    double Z = 0;
    double eta12 = 0;
    double eta13 = 0;
    double eta23 = 0;
    aiw::Vec<3, double> m1, m2, m3, m1eta23, m2eta13, m3eta12, m1eta12, m1eta13, m2eta12, m2eta23, m3eta13, m3eta23;

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
                m2 += f*m_j*dm;
                m3 += f*m_k*dm;
                double et23 = m_j*m_k;
                double et12 = m_j*m_i;
                double et13 = m_i*m_k;

                eta12 += f*et12*dm;
                eta13 += f*et13*dm;
                eta23 += f*et23*dm;

                m1eta23 += f*dm*et23*m_i;
                m2eta13 += f*dm*et13*m_j;
                m3eta12 += f*dm*et12*m_k;

                m1eta12 += f*dm*et12*m_i;
                m1eta13 += f*dm*et13*m_i;

                m2eta12 += f*dm*et12*m_j;
                m2eta23 += f*dm*et23*m_j;

                m3eta13 += f*dm*et13*m_k;
                m3eta23 += f*dm*et23*m_k;
            }
        }
    }

    moms.eta_12 = eta12/Z;
    moms.eta_13 = eta13/Z;
    moms.eta_23 = eta23/Z;
    moms.m1 = m1/Z;
    moms.m2 = m2/Z;
    moms.m3 = m3/Z;
    moms.m1eta_23 = m1eta23/Z;
    moms.m2eta_13 = m2eta13/Z;
    moms.m3eta_12 = m3eta12/Z;
    moms.m1eta_12 = m1eta12/Z;
    moms.m1eta_13 = m1eta13/Z;
    moms.m2eta_12 = m2eta12/Z;
    moms.m2eta_23 = m2eta23/Z;
    moms.m3eta_13 = m3eta13/Z;
    moms.m3eta_23 = m3eta23/Z;
    
    moms.Z = Z;
    return moms;
}

int main(){
    aiw::sph_init_table(rank);	
    int sz_sph = aiw::sph_vertex_num(rank);
    
    aiw::Vec<3, double> p1{1, 2, 3}, p2{4., 0, 6}, p3{7, 8, 0};
    p1 *= 1./p1.abs(); p2 *= 1./p2.abs(); p3 *= 1./p3.abs();
    double lmbd_12{1}, lmbd_13{2}, lmbd_23{3};

    double T1{1}, T2{10};
    int N = 10;
    double dT = (T2-T1)/N;
    std::cout << "#:T x x_t y y_t z z_t eta eta_t Z Z_t time time_t\n";
    double t_init =0, t_calc = 0;
    for(int i = 0; i < N; i++){
        double T = T1 + dT*i;
        moments3<double> moms, moms_t;
        clock_t start1 = clock();
        moms.lmbd12 = lmbd_12/T; moms.lmbd13 = lmbd_13/T; moms.lmbd23 = lmbd_23/T;
        moms.p1 = p1/T; moms.p2 = p2/T; moms.p3 = p3/T;
        moms.rank = rank;
        moms.sz_sphere = sz_sph;
        moms.calc();
        clock_t end1 = clock();
        double t1 = double(end1 - start1) / CLOCKS_PER_SEC;
        clock_t start2 = clock();
        moms_t = calc_directly(moms.p1, moms.p2, moms.p3, moms.lmbd12, moms.lmbd13, moms.lmbd23, sz_sph);
        clock_t end2 = clock();
        double t2 = double(end2 - start2) / CLOCKS_PER_SEC;
        std::cout << T << " ";

        #ifdef calc_m1eta_23
        std::cout << moms.m1eta_23[0] << " " << moms_t.m1eta_23[0] << " " << moms.m1eta_23[1] << " " << moms_t.m1eta_23[1]
        << " " << moms.m1eta_23[2] << " " << moms_t.m1eta_23[2];
        #endif

        #ifdef calc_m2eta_13
        std::cout << moms.m2eta_13[0] << " " << moms_t.m2eta_13[0] << " " << moms.m2eta_13[1] << " " << moms_t.m2eta_13[1]
        << " " << moms.m2eta_13[2] << " " << moms_t.m2eta_13[2];
        #endif

        #ifdef calc_m3eta_12
        std::cout << moms.m3eta_12[0] << " " << moms_t.m3eta_12[0] << " " << moms.m3eta_12[1] << " " << moms_t.m3eta_12[1]
        << " " << moms.m3eta_12[2] << " " << moms_t.m3eta_12[2];
        #endif

        #ifdef calc_m1eta_12
        std::cout << moms.m1eta_12[0] << " " << moms_t.m1eta_12[0] << " " << moms.m1eta_12[1] << " " << moms_t.m1eta_12[1]
        << " " << moms.m1eta_12[2] << " " << moms_t.m1eta_12[2];
        #endif

        #ifdef calc_m1eta_13
        std::cout << moms.m1eta_13[0] << " " << moms_t.m1eta_13[0] << " " << moms.m1eta_13[1] << " " << moms_t.m1eta_13[1]
        << " " << moms.m1eta_13[2] << " " << moms_t.m1eta_13[2];
        #endif

        #ifdef calc_m2eta_12
        std::cout << moms.m2eta_12[0] << " " << moms_t.m2eta_12[0] << " " << moms.m2eta_12[1] << " " << moms_t.m2eta_12[1]
        << " " << moms.m2eta_12[2] << " " << moms_t.m2eta_12[2];
        #endif

        #ifdef calc_m2eta_23
        std::cout << moms.m2eta_23[0] << " " << moms_t.m2eta_23[0] << " " << moms.m2eta_23[1] << " " << moms_t.m2eta_23[1]
        << " " << moms.m2eta_23[2] << " " << moms_t.m2eta_23[2];
        #endif

        #ifdef calc_m3eta_13
        std::cout << moms.m3eta_13[0] << " " << moms_t.m3eta_13[0] << " " << moms.m3eta_13[1] << " " << moms_t.m3eta_13[1]
        << " " << moms.m3eta_13[2] << " " << moms_t.m3eta_13[2];
        #endif

        #ifdef calc_m3eta_23
        std::cout << moms.m3eta_23[0] << " " << moms_t.m3eta_23[0] << " " << moms.m3eta_23[1] << " " << moms_t.m3eta_23[1]
        << " " << moms.m3eta_23[2] << " " << moms_t.m3eta_23[2];
        #endif

        #ifdef calc_m1
        std::cout << moms.m1[0] << " " << moms_t.m1[0] << " " << moms.m1[1] << " " << moms_t.m1[1]
        << " " << moms.m1[2] << " " << moms_t.m1[2];
        #endif

        #ifdef calc_m2
        std::cout << moms.m2[0] << " " << moms_t.m2[0] << " " << moms.m2[1] << " " << moms_t.m2[1]
        << " " << moms.m2[2] << " " << moms_t.m2[2];
        #endif

        #ifdef calc_m3
        std::cout << moms.m3[0] << " " << moms_t.m3[0] << " " << moms.m3[1] << " " << moms_t.m3[1]
        << " " << moms.m3[2] << " " << moms_t.m3[2];
        #endif

        #ifdef calc_eta_12
        std::cout << " " << moms.eta_12 << " " << moms_t.eta_12;
        #endif

        #ifdef calc_eta_13
        std::cout << " " << moms.eta_13 << " " << moms_t.eta_13;
        #endif

        #ifdef calc_eta_23
        std::cout << " " << moms.eta_23 << " " << moms_t.eta_23;
        #endif

        std::cout << " " << moms.Z << " " << moms_t.Z << " " << t1 << " " << t2 << "\n";
        
        
        t_init += moms.t_init; t_calc += moms.t_calc;
    }
    std::cout << t_calc << " " << t_init << "\n";
    return 0;
}