#include "../include/heis_mag.hpp"
#include <ctime>

moments3<double> calc_directly(aiw::Vec<3, double> p1, aiw::Vec<3, double> p2, aiw::Vec<3, double> p3, double lmbd12, double lmbd13, double lmbd23){
    moments3<double> moms;
    int N_fi = 20;
    int N_teta = 20;
    double teta_max = M_PI;
    double teta_min = 0.;
    double fi_max = 2*M_PI;
    double fi_min = 0.;
    double d_fi = (fi_max-fi_min)/N_fi;
    double d_teta = (teta_max-teta_min)/N_teta;

    double Z = 0;
    double eta12 = 0;
    aiw::Vec<3, double> m1;

    for(int i1 = 0; i1 < N_fi; i1++){
        double fi1 = fi_min+i1*d_fi;
        for(int j1 = 0; j1 < N_teta; j1++){
            double teta1 = teta_min+j1*d_teta;
            aiw::Vec<3, double> M1{sin(teta1)*cos(fi1), sin(teta1)*sin(fi1), cos(teta1)};
            for(int i2 = 0; i2 < N_fi; i2++){
                double fi2 = fi_min+i2*d_fi;
                for(int j2 = 0; j2 < N_teta; j2++){
                    double teta2 = teta_min+j2*d_teta;
                    aiw::Vec<3, double> M2{sin(teta2)*cos(fi2), sin(teta2)*sin(fi2), cos(teta2)};
                    for(int i3 = 0; i3 < N_fi; i3++){
                        double fi3 = fi_min+i3*d_fi;
                        for(int j3 = 0; j3 < N_teta; j3++){
                            double teta3 = teta_min+j3*d_teta;
                            aiw::Vec<3, double> M3{sin(teta3)*cos(fi3), sin(teta3)*sin(fi3), cos(teta3)};
                            double f = sin(teta1)*sin(teta2)*sin(teta3)*exp(p1*M1+p2*M2+p3*M3+lmbd12*M1*M2+lmbd13*M1*M3+lmbd23*M2*M3);
                            Z += f;
                            eta12 += M1*M2*f;
                            m1 += M1*f;
                        }
                    }
                }
            }
        }
    }
    double k = pow(d_fi*d_teta, 3);
    Z *= k;
    m1 = m1*k/Z;
    eta12 = eta12*k/Z;

    moms.eta_12 = eta12;
    moms.m1 = m1;
    moms.Z = Z;
    return moms;
}

int main(){
    aiw::Vec<3, double> p1{1, 2, 3}, p2{4., 0, 6}, p3{7, 8, 0};
    p1 *= 1./p1.abs(); p2 *= 1./p2.abs(); p3 *= 1./p3.abs();
    double lmbd_12{1}, lmbd_13{2}, lmbd_23{3};

    double T1{1}, T2{10};
    int N = 10;
    double dT = (T2-T1)/N;
    std::cout << "#:T m1x m1x_t m1y m1y_t m1z m1z_t eta12 eta12_t Z Z_t time time_t\n";
    for(int i = 0; i < N; i++){
        double T = T1 + dT*i;
        moments3<double> moms, moms_t;
        clock_t start1 = clock();
        moms.lmbd12 = lmbd_12/T; moms.lmbd13 = lmbd_13/T; moms.lmbd23 = lmbd_23/T;
        moms.p1 = p1/T; moms.p2 = p2/T; moms.p3 = p3/T;
        moms.calc();
        clock_t end1 = clock();
        double t1 = double(end1 - start1) / CLOCKS_PER_SEC;
        clock_t start2 = clock();
        moms_t = calc_directly(moms.p1, moms.p2, moms.p3, moms.lmbd12, moms.lmbd13, moms.lmbd23);
        clock_t end2 = clock();
        double t2 = double(end2 - start2) / CLOCKS_PER_SEC;
        std::cout << T << " " << moms.m1[0] << " " << moms_t.m1[0] << " " << moms.m1[1] << " " << moms_t.m1[1]
        << " " << moms.m1[2] << " " << moms_t.m1[2] << " " << moms.eta_12 << " " << moms_t.eta_12 << " " << moms.Z << " " << moms_t.Z << " " << t1 << " " << t2 << "\n";
    }
    
    return 0;
}