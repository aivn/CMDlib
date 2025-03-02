#include "../include/heis_mag.hpp"
#include <ctime>
#include <aiwlib/sphere>

template<typename T>
struct mom_m{
    int rank, sz_sphere;
    heis_mag<T> Z_t, m1_t;

    T lmbd12, lmbd13, lmbd23, lmbd14, lmbd15, lmbd16, lmbd24, lmbd25, lmbd26, lmbd34, lmbd35, lmbd36, lmbd45, lmbd46, lmbd56;
    aiw::Vec<3, T> p1, p2, p3, p4, p5, p6;

    T Z, m;
    aiw::Vec<3,T> m1;
    double t_calc = 0;
    int l_max;

    void calc(){
        clock_t start = clock();
        Z_t.l_max = l_max;
        Z_t.init(6, rank, sz_sphere);
        Z_t.set_lambda(0, 1, lmbd12);
        Z_t.set_lambda(0, 2, lmbd13);
        Z_t.set_lambda(0, 3, lmbd14);
        Z_t.set_lambda(0, 4, lmbd15);
        Z_t.set_lambda(0, 5, lmbd16);
        Z_t.set_lambda(1, 2, lmbd23);
        Z_t.set_lambda(1, 3, lmbd24);
        Z_t.set_lambda(1, 4, lmbd25);
        Z_t.set_lambda(1, 5, lmbd26);
        Z_t.set_lambda(2, 3, lmbd34);
        Z_t.set_lambda(2, 4, lmbd35);
        Z_t.set_lambda(2, 5, lmbd36);
        Z_t.set_lambda(3, 4, lmbd45);
        Z_t.set_lambda(3, 5, lmbd46);
        Z_t.set_lambda(4, 6, lmbd56);
        Z_t.set_p(0, p1);
        Z_t.set_p(1, p2);
        Z_t.set_p(2, p3);
        Z_t.set_p(3, p4);
        Z_t.set_p(4, p5);
        Z_t.set_p(5, p6);

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

template<typename T>
struct mom_m_exact{
    int rank, sz_sphere;
    // heis_mag<T> Z_t, m1_t;

    T lmbd12, lmbd13, lmbd23, lmbd14, lmbd15, lmbd16, lmbd24, lmbd25, lmbd26, lmbd34, lmbd35, lmbd36, lmbd45, lmbd46, lmbd56;
    aiw::Vec<3, T> p1, p2, p3, p4, p5, p6;

    T Z, m;
    aiw::Vec<3,T> m1;
    double t_calc = 0;
    int l_max;

    void calc(){
        clock_t start = clock();
        
        for(int i1=0; i1<sz_sphere; i1++){
            const aiw::Vec<3> &m_1 = aiw::sph_vert(i1, rank); double dm_1 = aiw::sph_vert_area(i1, rank);
            for(int i2=0; i2<sz_sphere; i2++){
                const aiw::Vec<3> &m_2 = aiw::sph_vert(i2, rank); double dm_2 = aiw::sph_vert_area(i2, rank);
                for(int i3=0; i3<sz_sphere; i3++){
                    const aiw::Vec<3> &m_3 = aiw::sph_vert(i3, rank); double dm_3 = aiw::sph_vert_area(i3, rank);
                    for(int i4=0; i4<sz_sphere; i4++){
                        const aiw::Vec<3> &m_4 = aiw::sph_vert(i4, rank); double dm_4 = aiw::sph_vert_area(i4, rank);
                        for(int i5=0; i5<sz_sphere; i5++){
                            const aiw::Vec<3> &m_5 = aiw::sph_vert(i5, rank); double dm_5 = aiw::sph_vert_area(i5, rank);
                            for(int i6=0; i6<sz_sphere; i6++){
                                const aiw::Vec<3> &m_6 = aiw::sph_vert(i6, rank); double dm_6 = aiw::sph_vert_area(i6, rank);
                                double w = p1*m_1+p2*m_2+p3*m_3+p4*m_4+p5*m_5+p6*m_6;
                                w += lmbd12*m_1*m_2+lmbd13*m_1*m_3+lmbd14*m_1*m_4+lmbd15*m_1*m_5+lmbd23*m_2*m_3+lmbd24*m_2*m_4+lmbd25*m_2*m_5+lmbd26*m_2*m_6;
                                w += lmbd34*m_3*m_4+lmbd35*m_3*m_5+lmbd3*m_3*m_6+lmbd45*m_4*m_5+lmbd46*m_4*m_6+lmbd56*m_5*m_6;
                                double f = exp(w);
                                double dm = dm_1*dm_2*dm_3*dm_4*dm_5*dm_6;
                                Z += f*dm;
                                m1 += f*m_1*dm;
                            }
                        }
                    }
                }
            }
        }
    
        m1 /= Z;
        m = sqrt(m1[0]*m1[0]+m1[1]*m1[1]+m1[2]*m1[2]);
        moms.m = m;
        moms.Z = Z;

        clock_t end = clock();
        t_calc += double(end - start) / CLOCKS_PER_SEC;
    }
};

int main(){
    aiw::Vec<3, double> p1{1, 2, 3}, p2{4., 0, 6}, p3{7, 8, 0}, p4{7, 2, 3}, p5{24., 0, 6}, p6{7, 88, 3};
    p1 *= 1./p1.abs(); p2 *= 1./p2.abs(); p3 *= 1./p3.abs();
    p4 *= 1./p4.abs(); p5 *= 1./p5.abs(); p6 *= 1./p6.abs();
    double lmbd12{1}, lmbd13{2}, lmbd14{3},  lmbd15{1}, lmbd16{2}; 
    double lmbd23{1}, lmbd24{2}, lmbd25{3},  lmbd26{2};
    double lmbd34{1}, lmbd35{2}, lmbd36{3};
    double lmbd45{1}, lmbd46{2};
    double lmbd56{3};

    int rank_begin = 0; int rank_end = 0;
    int l_min = 0; int l_max = 3;
    aiw::sph_init_table(rank_end+1);	
    int sz_sph = aiw::sph_vertex_num(rank_end+1);
    mom_m_exact<double> exact;
    double exact_m = exact.m;

    std::cout << "s l delta_1 t_1 delta_2 t_2\n";
    for(int r = rank_begin; r <= rank_end; r++){
        aiw::sph_init_table(r);	
        
        int sz_sph = aiw::sph_vertex_num(r);
        
        mom_m_exact<double> mom_direct;
        mom_direct.rank = r;
        mom_direct.sz_sphere = sz_sph;
        mom_direct.p1 = p1;
        mom_direct.p2 = p2;
        mom_direct.p3 = p3;
        mom_direct.p4 = p4;
        mom_direct.p5 = p5;
        mom_direct.p6 = p6;
        mom_direct.lmbd12 = lmbd12;
        mom_direct.lmbd13 = lmbd13;
        mom_direct.lmbd14 = lmbd14;
        mom_direct.lmbd15 = lmbd15;
        mom_direct.lmbd16 = lmbd16;
        mom_direct.lmbd23 = lmbd23;
        mom_direct.lmbd24 = lmbd24;
        mom_direct.lmbd25 = lmbd25;
        mom_direct.lmbd26 = lmbd26;
        mom_direct.lmbd34 = lmbd34;
        mom_direct.lmbd35 = lmbd35;
        mom_direct.lmbd36 = lmbd36;
        mom_direct.lmbd45 = lmbd45;
        mom_direct.lmbd46 = lmbd46;
        mom_direct.lmbd56 = lmbd56;

        for(int l = l_min; l <= l_max; l++){
            mom_m<double> mom;
            mom.l_max = l;
            mom_direct.rank = r;
            mom.sz_sphere = sz_sph;
            mom.p1 = p1;
            mom.p2 = p2;
            mom.p3 = p3;
            mom.p4 = p4;
            mom.p5 = p5;
            mom.p6 = p6;
            mom.lmbd12 = lmbd12;
            mom.lmbd13 = lmbd13;
            mom.lmbd14 = lmbd14;
            mom.lmbd15 = lmbd15;
            mom.lmbd16 = lmbd16;
            mom.lmbd23 = lmbd23;
            mom.lmbd24 = lmbd24;
            mom.lmbd25 = lmbd25;
            mom.lmbd26 = lmbd26;
            mom.lmbd34 = lmbd34;
            mom.lmbd35 = lmbd35;
            mom.lmbd36 = lmbd36;
            mom.lmbd45 = lmbd45;
            mom.lmbd46 = lmbd46;
            mom.lmbd56 = lmbd56;
            mom.calc();
            std::cout << sz_sph << " " << l << " " << mom.m-exact_m << " " << mom.t_calc  << " " << mom_direct.m-exact_m << " " << mom.t_calc << "\n";
        }
    }

    return 0;
}