#include <bits/stdc++.h>
#include "integral_equations_2plus1_system.h"
#include "mphib_testing.h"
#include "SA_method_2plus1_system.h"

#define PRINT_PARAM(x) std::cout<<#x<<" = "<<x<<std::endl; 

void test_M2()
{
    double pmin = 1e-10;
    double pmax = 1.0; 
    double ppoints = 5; 

    double delp = std::abs(pmax - pmin)/ppoints; 

    std::vector<comp> pvec; 
    for(int i=0; i<ppoints; ++i)
    {
        double p = pmin + i*delp; 

        pvec.push_back(p); 
    }

    double E = 3.1; 
    double a = -2.0; 
    double r = 0.0; 
    double m1 = 1.0;
    double m2 = 0.5; 

    double eps = 0.0; 
    double total_P = 0.0; 

    std::vector<comp> weights; 
    Eigen::MatrixXcd kernmat;

    kernel_2plus1_system_ERE(kernmat, E, pvec, pvec, weights, weights, m1, m2, eps, eps, eps, total_P, a, r, a, r); 

}

void test_q2k_sigk()
{
    double En = 3.1; 
    double pmin = -1.5;
    double pmax = 1.5;
    double ppoints = 20;
    double mi = 1.0;
    double total_P = 0.0; 

    double delp = std::abs(pmax - pmin)/ppoints; 

    std::vector<comp> repvec;
    std::vector<comp> impvec; 

    for(int i=0;i<ppoints; ++i)
    {
        double p = pmin + i*delp; 

        repvec.push_back(p);
        impvec.push_back(ii*p); 
    }

    for(int i=0; i<repvec.size(); ++i)
    {
        comp sig1 = sigma(En, repvec[i], mi, total_P);
        comp sig2 = sigma(En, impvec[i], mi, total_P); 

        std::cout<<repvec[i]<<'\t'<<sig1<<'\t'<<impvec[i]<<'\t'<<sig2<<std::endl; 
    }
}

void test_sigk_q2k()
{
    double En = 3.1;
    double mi = 1.0; 

    double re_sigmin = 2.1;
    double im_sigmin = -1.0;
    double im_sigmax = 1.0;
    double sig_points = 20;
    double del_sig = std::abs(im_sigmin - im_sigmax)/sig_points; 

    for(int i=0; i<sig_points; ++i)
    {
        double im_sig = im_sigmin + i*del_sig; 
        comp sigk = re_sigmin + ii*im_sig; 

        comp p = pmom(En, sigk, mi);

        std::cout<<"sigk:"<<sigk<<'\t'<<"p="<<p<<std::endl;  
    }
}

void M2k_Plotter_KK()
{
    // atmpi = 0.06906
    // atmK = 0.09698
    // scattering_length_1_piK = 4.04;// - 0.2; //total uncertainty 0.05 stat 0.15 systematic 
    // scattering_length_2_KK = 4.07;// - 0.07; //total uncertainty 0.07 stat 
    // atinv = 5.666 GeV 
    // anisotropy xi = 3.444

    //double 
}

void test_kinematic_variables()
{
    double m1 = 1.0;
    double m2 = 1.0;//0.999;

    double mi = m1;
    double mj = m1;
    double mk = m2; 

    double a01 = 2.0; 
    double a02 = 2.0; 

    double En = 2.6;

    comp kmax1 = pmom(En, 0.0, m1);
    comp kmax2 = pmom(En, 0.0, m2);
    comp sigb1 = sigma_b_plus(a01, mj, mk); 
    comp sigb2 = sigma_b_minus(a01, mj, mk); 

    comp qb = qb_i(En, sigb1, mi);

    std::cout<<"kmax1 = "<<kmax1<<std::endl; 
    std::cout<<"kmax2 = "<<kmax2<<std::endl; 
    std::cout<<"sigb1 = "<<sigb1<<std::endl; 
    std::cout<<"sigb2 = "<<sigb2<<std::endl; 
    std::cout<<"qb = "<<qb<<std::endl; 


}

void test_dpqb_building()
{
    double m1 = 1.0;
    double m2 = 0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 100;//200;//250;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "check_dqq_test1.dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 

    double En_points = 3;//100.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //for(int i=0; i<(int)En_points; ++i)
    {
        int i = 1; 
        comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 25; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = 0.01; //eps_for_m2k; 
        dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
        dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

        //std::cout<<"pvec size = "<<pvec_for_m1m2.size()<<std::endl; 
        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        comp gval = gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 

        std::cout<<En.real()<<'\t'<<(En*En).real()<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<std::endl; 
        fout<<std::real(En)<<'\t'<<std::real(En*En)<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<std::endl; 
        
    
    }
    fout.close(); 
}



void test_dpqb_building_1()
{
    double m1 = 1.0;
    double m2 = 1.0;//0.99999999990;//0.999; 
    double a0_1 = 16.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 250;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "mphib_above_phibthreshold_a" + std::to_string((int)a0_1) + ".dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    for(int i=0; i<(int)En_points; ++i)
    {
        comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 15; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = eps_for_m2k; 
        dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
        dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        comp gval = gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb,En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(20)<<En.real()<<'\t'<<(En*En).real()<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<std::endl; 
        fout<<std::setprecision(20)<<std::real(En)<<'\t'<<std::real(En*En)<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<std::endl; 
        
    
    }
    fout.close(); 
}


//We change m2 value by 0.1 starting from 0.9
//and check rhoMphib value 
void test_d11_building_1()
{
    comp ii = {0.0,1.0};
    double m1 = 1.0;
    double m2 = 9.0; //1.0;//0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    double m2_initial = 0.7; 
    double m2_final = 1.0; 
    double del_m2 = 0.1; 
    double m2_points = std::abs(m2_initial - m2_final)/del_m2; 
    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 5000;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.00001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    
    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::ofstream fout1; 
    

    /* file index naming*/
    fout<<"En"<<'\t'
        <<"s"<<'\t'
        <<"Re(rhoMphib)"<<'\t'
        <<"Im(rhoMphib)"<<'\t'
        <<"perc_diff"
        <<std::endl; 

    fout1<<"En"<<'\t'
         <<"s"<<'\t'
         <<"Redqq11"<<'\t'
         <<"Imdqq11"<<'\t'
         <<"Redqq12"<<'\t'
         <<"Imdqq12"<<'\t'
         <<"Redqq21"<<'\t'
         <<"Imdqq21"<<'\t'
         <<"Redqq22"<<'\t'
         <<"Imdqq22"
         <<std::endl; 

    for(int j=0; j<(int)m2_points; ++j)
    {
        m2 = m2_initial + j*del_m2; 

        comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
        comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
        comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
        comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
        comp phib1plus = std::sqrt(sigb1plus) + m1;
        comp phib1minus = std::sqrt(sigb1minus) + m1;
        comp phib2plus = std::sqrt(sigb2plus) + m2;
        comp phib2minus = std::sqrt(sigb2minus) + m2;

        std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
        std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
        std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
        std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

        std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
        std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
        std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
        std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

        comp threeparticle_threshold = (m1 + m1 + m2); 

        //comp En = (phib1 + threeparticle_threshold)/2.0;

        std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

        std::string filename =    "mphib_above_phibthreshold_m2_" + std::to_string(m2) 
                                + "_a1_" + std::to_string(a0_1)
                                + ".dat";
        std::string filename1 =   "dqq_m2_" + std::to_string(m2)
                                + "_a1_" + std::to_string(a0_1) 
                                + ".dat"; 

        fout.open(filename.c_str()); 
        fout1.open(filename1.c_str()); 

        comp En_initial = phib1plus; 
        comp En_final = threeparticle_threshold; 

        double En_points = 500.0; 
        comp delEn = std::abs(En_initial - En_final)/En_points; 
        //abort();
    
    for(int i=0; i<(int)En_points; ++i)
    {
        double img_En = 0.0; 
        comp En = std::real(En_initial + ((comp)i)*delEn) + ii*img_En; 
        //comp qb; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 25; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = eps_for_m2k; 
        dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
        dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        comp dqq12 = dqqmat(0,1); 
        comp dqq21 = dqqmat(1,0); 
        comp dqq22 = dqqmat(1,1); 

        
        comp gval = gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb,En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(6)
                 <<"m2 = "<<m2<<'\t'
                 <<"sigb = "<<sigb1plus<<'\t'
                 <<"En = "<<En<<'\t'
                 <<"s = "<<(En*En)<<'\t'
                 <<"rhoM = "<<rhophib_mphib<<'\t'
                 <<"diff = "<<diff<<'\t'
                 <<"d12 = "<<dqq12<<'\t'
                 <<"d21 = "<<dqq21<<'\t'
                 <<"d22 = "<<dqq22
                 <<std::endl; 
        
        
        fout<<std::setprecision(20)
            <<std::real(En)<<'\t'
            <<std::real(En*En)<<'\t'
            <<rhophib_mphib.real()<<'\t'
            <<rhophib_mphib.imag()<<'\t'
            <<diff
            <<std::endl; 
        
        fout1<<std::setprecision(20)
             <<std::real(En)<<'\t'
             <<std::real(En*En)<<'\t'
             <<std::real(dqq11)<<'\t'
             <<std::imag(dqq11)<<'\t'
             <<std::real(dqq12)<<'\t'
             <<std::imag(dqq12)<<'\t'
             <<std::real(dqq21)<<'\t'
             <<std::imag(dqq21)<<'\t'
             <<std::real(dqq22)<<'\t'
             <<std::imag(dqq22)
             <<std::endl; 
    
    }
    fout.close(); 
    fout1.close(); 
    }
}

//We check the N dependence of our solutions
void test_dpqb_vs_N_building()
{
    double m1 = 1.0;
    double m2 = 0.9;//0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    //int number_of_points = 250;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "mphib_above_phibthreshold_a" + std::to_string((int)a0_1) + "_vs_N_1.dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En = (phib1plus + threeparticle_threshold)/2.0; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    int N_initial = 50; 
    int N_final = 1000; 
    int N_points = 30; 
    int del_N = std::abs(N_initial - N_final)/N_points; 

    for(int i=0; i<(int)N_points; ++i)
    {
        //comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        int number_of_points = N_initial + i*del_N; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 55; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = 0.27;//energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = eps_for_m2k; 
        test_dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points, 'n'); 
        test_dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points, 'n');
        //std::cout<<"eps for m2k = "<<eps_for_m2k<<std::endl; 
        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        double eta_1 = 1.0; 
        comp gval = gfunc_i(eta_1, sigb1plus, m1, m2); //gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb,En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(20)
                 <<En.real()<<'\t'
                 <<(En*En).real()<<'\t'
                 <<rhophib_mphib<<'\t'
                 <<(1.0/mphib).imag()<<'\t'
                 <<rhophib_val<<'\t'
                 <<diff<<'\t'
                 <<number_of_points<<std::endl; 
        fout<<std::setprecision(20)<<std::real(En)<<'\t'<<std::real(En*En)<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<'\t'<<number_of_points<<std::endl; 
        
    
    }
    fout.close(); 
}

//This is not used anymore
//Energy is set to be at one fourth range between 
//phib threshold and three particle threshold 
void test_dpqb_vs_N_building_1()
{
    double m1 = 1.0;
    double m2 = 0.9;//0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    //int number_of_points = 250;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "mphib_above_phibthreshold_a" + std::to_string((int)a0_1) + "_vs_N_en_onefourth_2.dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp En = (En_1 + En_2)/2.0; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //N range = [50, 500] and [500, 5000]
    int N_initial = 500; 
    int N_final = 5000; 
    int N_points = 10; 
    int del_N = std::abs(N_initial - N_final)/N_points; 

    for(int i=0; i<(int)N_points; ++i)
    {
        //comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        int number_of_points = N_initial + i*del_N; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 55; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = 0.0;//eps_for_m2k; 
        dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
        dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        double eta_1 = 1.0; 
        comp gval = gfunc_i(eta_1, sigb1plus, m1, m2);//gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb, En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(20)<<En.real()<<'\t'<<(En*En).real()<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<'\t'<<number_of_points<<std::endl; 
        fout<<std::setprecision(20)<<std::real(En)<<'\t'<<std::real(En*En)<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<'\t'<<number_of_points<<std::endl; 
        
    
    }
    fout.close(); 
}

//Here we calculate d_22 for mphib 
//where the two body sub-channel is 
//a bound-state of m1 and m1 
void test_dpqb_for_m1m1_vs_N_building_1()
{
    double m1 = 1.0;
    double m2 = 1.0;//0.99999999990;//0.999; 
    double a0_1 = -2.0; 
    double a0_2 = -2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    //int number_of_points = 250;  

    double eps_for_m2k = 0.0000001;
    double eps_for_ope = 0.0000001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    


    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "mphib_above_for_m1m1_phibthreshold_a" + std::to_string((int)a0_1) + "_vs_N_en_onefourth_2.dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp En = (En_1 + En_2)/2.0; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //N range = [50, 500] and [500, 5000]
    int N_initial = 10; 
    int N_final = 100; 
    int N_points = 10; 
    int del_N = std::abs(N_initial - N_final)/N_points; 

    //qb vals for i=1 and i=2 case:
    comp qb_val1plus = qb_i(En, sigb1plus, m1);
    comp qb_val1minus = qb_i(En, sigb1minus, m1);
    comp qb_val2plus = qb_i(En, sigb2plus, m2);
    comp qb_val2minus = qb_i(En, sigb2minus, m2);

    std::cout<<"qb 1 + = "<<qb_val1plus<<std::endl;
    std::cout<<"qb 1 - = "<<qb_val1minus<<std::endl;
    std::cout<<"qb 2 + = "<<qb_val2plus<<std::endl;
    std::cout<<"qb 2 - = "<<qb_val2minus<<std::endl;

    //test M2k for large sigk approx //
    comp sigk_test = 50000000000; 
    comp M2k_test1 = M2k_ERE_s2k(0.5, En, sigk_test, 0.0, 2.0, 0.0, m1, m1, m2, eps_for_m2k); 
    comp M2k_test2 = M2k_ERE_large_sigk_approx(0.5, sigk_test, 2.0); 
    std::cout<<"M2k large sigk approx test"<<std::endl; 
    std::cout<<"M2k test 1 = "<<M2k_test1<<std::endl;
    std::cout<<"M2k test 2 = "<<M2k_test2<<std::endl; 
    for(int i=0; i<(int)N_points; ++i)
    {
        //int i=0; 
        //comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        int number_of_points = N_initial + i*del_N; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 55; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = eps_for_m2k; 
        char debug = 'n'; 
        test_dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points, debug); 
        test_dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points, debug);

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/
        
        comp dqq11 = dqqmat(0,0); 
        comp dqq22 = dqqmat(1,1); 
        
        double eta_1 = 0.5; 
        comp gval = gfunc_i(eta_1, sigb1plus, m1, m2);//gfunc(sigb1plus,a0_1);
        
        std::cout<<"gval = "<<gval<<'\t'
                 <<"dqq00 = "<<dqqmat(0,0)<<'\t'
                 <<"dqq01 = "<<dqqmat(0,1)<<'\t'
                 <<"dqq10 = "<<dqqmat(1,0)<<'\t'
                 <<"dqq11 = "<<dqqmat(1,1)<<std::endl;
        
        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        
        comp rhophib_val = rhophib(qb,En); 
        std::cout<<"rhophib = "<<rhophib_val<<std::endl; 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(20)<<En.real()<<'\t'<<(En*En).real()<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<'\t'<<number_of_points<<std::endl; 
        fout<<std::setprecision(20)
            <<std::real(En)<<'\t'
            <<std::real(En*En)<<'\t'
            <<rhophib_mphib.real()<<'\t'
            <<rhophib_mphib.imag()<<'\t'
            <<diff<<'\t'
            <<number_of_points<<std::endl; 
        
    
    }
    fout.close(); 
}


//This creates the datafile for the plot analogous 
//to fig:8 of https://arxiv.org/pdf/2010.09820
void test_delta_rhophib_density()
{
    double m1 = 1.0;
    double m2 = 0.9;//0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    //int number_of_points = 250;  

    double eps_for_m2k = 0.001;
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "deltarhophib_density_a" + std::to_string((int)a0_1) + "_eps_vs_N.dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En = (phib1plus + threeparticle_threshold)/2.0; 


    int N_initial = 100; 
    int N_final = 1000; 
    int N_points = 400; 
    int del_N = std::abs(N_initial - N_final)/N_points; 

    std::vector<double> eps_multiplier(4); 
    eps_multiplier[0] = 1.0e-4; 
    eps_multiplier[1] = 1.0e-3; 
    eps_multiplier[2] = 1.0e-2; 
    eps_multiplier[3] = 1.0e-1; 

    double eps_coeff_ini = 1.0e-10; 
    double eps_coeff_fin = 1.0; 
    double eps_coeff_points = 100.0; 
    double del_eps_coeff = std::abs(eps_coeff_ini - eps_coeff_fin)/eps_coeff_points; 

    std::vector<double> eps_vec((int)eps_coeff_points*eps_multiplier.size()); 

    int eps_counter = 0; 
    for(int i=0; i<eps_multiplier.size(); ++i)
    {
        double eps_multiplier_val = eps_multiplier[i]; 
        for(int j=0; j<(int)eps_coeff_points; ++j)
        {
            double eps_coeff = eps_coeff_ini + j*del_eps_coeff;
            double eps_val = eps_multiplier_val*eps_coeff; 
            eps_vec[eps_counter] = eps_val;
            eps_counter = eps_counter + 1;  
            //std::cout<<std::setprecision(20); 
            //std::cout<<i<<'\t'<<j<<'\t'<<eps_multiplier_val<<'\t'<<eps_coeff<<'\t'<<eps_val<<std::endl; 
        }
    }

    

    for(int i=0; i<(int)N_points; ++i)
    {
        int number_of_points = N_initial + i*del_N; 

    for(int j=0; j<eps_vec.size(); ++j)
    {
        //comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 25; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        double eps_value = eps_vec[j]; 
        eps_for_m2k = eps_value; //energy_dependent_epsilon(eta, En, qb, sigb1plus, kmax, m1, number_of_points);
        eps_for_ope = eps_for_m2k; 
        dpqb_solver_ERE(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points); 
        dqq_interpolator(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2, number_of_points);

        comp dqq11 = dqqmat(0,0); 
        
        comp gval = gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb,En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(20)
                 <<"E="<<En.real()<<'\t'
                 <<"N="<<number_of_points<<'\t'
                 <<"eps="<<eps_value<<'\t'
                 <<"dRho="<<diff<<std::endl; 
        std::cout<<std::setprecision(20)
                 <<"rhopb="<<rhophib_val<<'\t'
                 <<"qb="<<qb<<'\t'
                 <<"mphib_inv="<<1.0/mphib<<std::endl; 
        
        
        fout<<std::setprecision(20)
            <<std::real(En)<<'\t'
            <<std::real(En*En)<<'\t'
            <<rhophib_mphib.real()<<'\t'
            <<rhophib_mphib.imag()<<'\t'
            <<diff<<'\t'
            <<number_of_points<<'\t'
            <<eps_value<<'\t'
            <<qb.real()<<'\t'
            <<kmax.real()<<'\t'
            <<sigb1plus.real()<<std::endl; 
        
    
    }
    }
    fout.close(); 
    
}


void Mphib_degenerate_testing()
{
    double m = 1.0; 
    double a0 = 2.0; 
    double total_P = 0.0; 
    double eps_for_ope = 0.0; 
    double eps_for_cutoff = 0.0; 
    int number_of_points = 3000; 

    comp sigb = sigma_b_plus(a0, m, m); 
    comp phibval = std::sqrt(sigb) + m; 
    double En_initial = std::real(phibval); 
    double En_final = 3.0*m; 
    double En_points = 250; 

    double eta = 25; 

    double del_En = std::abs(En_initial - En_final)/En_points; 

    for(int i=0; i<En_points; ++i)
    {
        double En = En_initial + i*del_En; 
        comp qb = qb_i(En, sigb, m); 
        comp kmax = pmom(En, 0.0, m); 
        double eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb, kmax, m, number_of_points);
        //PRINT_PARAM(eps_for_m2k);
        Mphib_degenerate_mass(En, m, a0, total_P, eps_for_ope, eps_for_m2k, eps_for_cutoff, number_of_points); 

    }
}

void Mphib_degenerate_testing_vs_N()
{
    double m = 1.0; 
    double a0 = 2.0; 
    double total_P = 0.0; 
    double eps_for_ope = 0.0; 
    double eps_for_cutoff = 0.0; 
    int number_of_points = 3000; 

    comp sigb = sigma_b_plus(a0, m, m); 
    comp phibval = std::sqrt(sigb) + m; 
    //std::cout<<sigb<<'\t'<<phibval<<std::endl; 
    double En_initial = std::real(phibval); 
    double En_final = 3.0*m; 
    double En_points = 250; 
    double En = std::abs(En_initial + En_final)/2.0; 

    double eta = 25; 

    double del_En = std::abs(En_initial - En_final)/En_points; 

    int N_initial = 500; 
    int N_final = 5000; 
    int N_points = 10;
    int del_N = std::abs(N_initial - N_final)/N_points; 

    for(int i=0; i<N_points; ++i)
    {
        //double En = En_initial + i*del_En; 
        int number_of_points = N_initial + i*del_N; 
        comp qb = qb_i(En, sigb, m); 
        comp kmax = pmom(En, 0.0, m); 
        double eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb, kmax, m, number_of_points);
        //PRINT_PARAM(eps_for_m2k);
        Mphib_degenerate_mass(En, m, a0, total_P, eps_for_ope, eps_for_m2k, eps_for_cutoff, number_of_points); 

    }
}

void Mphib_degenerate_testing_vs_N_En_onefourth()
{
    double m = 1.0; 
    double a0 = 2.0; 
    double total_P = 0.0; 
    double eps_for_ope = 0.0; 
    double eps_for_cutoff = 0.0; 
    int number_of_points = 3000; 

    comp sigb = sigma_b_plus(a0, m, m); 
    comp phibval = std::sqrt(sigb) + m; 
    //std::cout<<sigb<<'\t'<<phibval<<std::endl; 
    double En_initial = std::real(phibval); 
    double En_final = 3.0*m; 
    double En_points = 250; 
    double En_1 = std::abs(En_initial + En_final)/2.0; 
    double En_2 = En_initial; 

    double En = std::abs(En_1 + En_2)/2.0; 

    double eta = 25; 

    double del_En = std::abs(En_initial - En_final)/En_points; 

    //Two different ranges of N used
    //N = [50, 500] and [500, 5000]
    int N_initial = 50; 
    int N_final = 500; 
    int N_points = 10;
    int del_N = std::abs(N_initial - N_final)/N_points; 

    for(int i=0; i<N_points; ++i)
    {
        //double En = En_initial + i*del_En; 
        int number_of_points = N_initial + i*del_N; 
        comp qb = qb_i(En, sigb, m); 
        comp kmax = pmom(En, 0.0, m); 
        double eps_for_m2k = energy_dependent_epsilon(eta, En, qb, sigb, kmax, m, number_of_points);
        //PRINT_PARAM(eps_for_m2k);
        Mphib_degenerate_mass(En, m, a0, total_P, eps_for_ope, eps_for_m2k, eps_for_cutoff, number_of_points); 

    }
}

void test_M2_sigb_range()
{
    double m1 = 1.0; 
    double m2_initial = 1.0e-10; 
    double m2_final = 1.0; 
    double m2_points = 10; 
    double del_m2 = std::abs(m2_initial - m2_final)/m2_points; 

    double a0_1_initial = 1.0e-10; 
    double a0_1_final = 3; 
    double a0_1_points = 1000; 
    double del_a0_1 = std::abs(a0_1_initial - a0_1_final)/a0_1_points; 

    std::string filename1 = "m2_dependence_on_sigb.dat";
    std::string filename2 = "a0_1_dependence_on_sigb.dat"; 
    std::string filename3 = "g_i_residue.dat";

    std::ofstream fout1, fout2, fout3; 
    
    fout1.open(filename1.c_str()); 
    fout2.open(filename2.c_str()); 
    fout3.open(filename3.c_str()); 

    for(int i=0; i<m2_points; ++i)
    {
        double a0_1 = 2.0; 
        double m2 = m2_initial + i*del_m2; 
        comp sigb1 = sigma_b_plus(a0_1, m1, m2);
        comp sigb2 = sigma_b_minus(a0_1, m1, m2); 
        double twoparticle_threshold = m1 + m2; 

        fout1<<std::setprecision(20); 
        fout1<<m2<<'\t'
             <<sigb1.real()<<'\t'
             <<sigb1.imag()<<'\t'
             <<sigb2.real()<<'\t'
             <<sigb2.imag()<<'\t'
             <<twoparticle_threshold*twoparticle_threshold<<std::endl;
    }
    fout1.close();

    for(int i=0; i<a0_1_points; ++i)
    {
        double a0_1 = a0_1_initial + i*del_a0_1; 
        double m2 = 0.71; 

        comp sigb1 = sigma_b_plus(a0_1, m1, m2); 
        comp sigb2 = sigma_b_minus(a0_1, m1, m2); 
        double twoparticle_threshold = m1 + m2; 

        fout2<<std::setprecision(20); 
        fout2<<a0_1<<'\t'
             <<sigb1.real()<<'\t'
             <<sigb1.imag()<<'\t'
             <<sigb2.real()<<'\t'
             <<sigb2.imag()<<'\t'
             <<twoparticle_threshold*twoparticle_threshold<<std::endl; 
    }
    fout2.close(); 

    
    for(int i=0; i<m2_points; ++i)
    {
        double eta_i = 1.0; 
        double a0_1 = 2.0; 
        double m2 = m2_initial + i*del_m2; 
        comp sigb1 = sigma_b_plus(a0_1, m1, m2); 
        comp sigb2 = sigma_b_minus(a0_1, m1, m2); 
        comp g1 = gfunc_i(eta_i, sigb1, m1, m2);
        comp g2 = gfunc_i(eta_i, sigb2, m1, m2); 

        double twoparticle_threshold = m1 + m2; 

        fout3<<m2<<'\t'
             <<g1.real()<<'\t'
             <<g1.imag()<<'\t'
             <<g2.real()<<'\t'
             <<g2.imag()<<std::endl; 
        
        /*std::cout<<"m2 = "<<m2<<'\t'
                 <<"sigb1 = "<<sigb1<<'\t'
                 <<"g1 = "<<g1<<'\t'
                 <<"sigb2 = "<<sigb2<<'\t'
                 <<"g2 = "<<g2<<std::endl; */
    }
    fout3.close(); 
    
    //compare g with the degenerate version
    double eta_i = 0.5; 
    double m2 = 1.0; 
    double a0_1 = 2.0; 
    comp sigb1 = sigma_b_plus(a0_1, m1, m2); 
    comp sigb2 = sigma_b_minus(a0_1, m1, m2); 
    comp g1 = gfunc_i(eta_i, sigb1, m1, m2);
    comp g2 = gfunc_i(eta_i, sigb2, m1, m2); 
    comp degen_g1 = gfunc(sigb1, a0_1);
    comp degen_g2 = gfunc(sigb2, a0_1); 

    std::cout<<g1<<'\t'<<g2<<'\t'<<degen_g1<<'\t'<<degen_g2<<std::endl; 

}

void test_M2_1()
{
    comp ii = {0.0, 1.0}; 
    double m1 = 1.0; 
    double m2 = 0.9; 
    double a0_1 = 16.0; 
    double a0_2 = 16.0; 

    double kmin = 0.0; 
    double kmax = 2.0; 
    double kpoints = 100;

    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    double En = 3.0; 
    double eps_for_m2k = 0.0;//0.0001; 
    double sigk_points = 5000.0; 
    double sigk_initial = 0.0; 
    double sigk_final_1 = 0.02; 
    double del_sigk_1 = std::abs(sigk_initial - sigk_final_1)/sigk_points; 
    double sigk_final = 5; 
    double del_sigk_2 = std::abs(sigk_final_1 - sigk_final)/sigk_points; 

    double sigk_imag = 0.00001; 

    std::string filename = "M2k_a0_" + to_string(a0_1) + ".dat"; 
    std::ofstream fout; 
    fout.open(filename.c_str()); 

    for(int i=0; i<sigk_points; ++i)
    {
        double sigk = sigk_initial + i*del_sigk_1;
        comp sigk_comp = sigk + ii*sigk_imag;  
        comp m2k_1 = M2k_ERE_s2k(eta_1, En, sigk_comp, 0.0, a0_1, 0.0, m1, m1, m2, eps_for_m2k);
        comp m2k_2 = M2k_ERE_s2k(eta_2, En, sigk_comp, 0.0, a0_2, 0.0, m2, m1, m1, eps_for_m2k);

        std::cout<<std::setprecision(20)<<sigk<<'\t'<<m2k_1<<'\t'<<m2k_2<<std::endl; 
        fout<<std::setprecision(20); 
        fout<<sigk<<'\t'
            <<m2k_1.real()<<'\t'
            <<m2k_1.imag()<<'\t'
            <<m2k_2.real()<<'\t'
            <<m2k_2.imag()<<std::endl;
    }
    for(int i=0; i<sigk_points; ++i)
    {
        double sigk = sigk_final_1 + i*del_sigk_2; 
        comp sigk_comp = sigk + ii*sigk_imag;  
        comp m2k_1 = M2k_ERE_s2k(eta_1, En, sigk_comp, 0.0, a0_1, 0.0, m1, m1, m2, eps_for_m2k);
        comp m2k_2 = M2k_ERE_s2k(eta_2, En, sigk_comp, 0.0, a0_2, 0.0, m2, m1, m1, eps_for_m2k);

        std::cout<<std::setprecision(20)<<sigk<<'\t'<<m2k_1<<'\t'<<m2k_2<<std::endl; 
        fout<<std::setprecision(20); 
        fout<<sigk<<'\t'
            <<m2k_1.real()<<'\t'
            <<m2k_1.imag()<<'\t'
            <<m2k_2.real()<<'\t'
            <<m2k_2.imag()<<std::endl;
    }
    fout.close(); 

}


void test_Gs()
{
    comp ii = {0.0, 1.0}; 
    double m1 = 1.0; 
    double m2 = 1.0; 

    double En = std::sqrt(8.3); 
    double delta = 1.0e-4; 
    double a0_1 = 16.0; 
    double a0_2 = 16.0;
    double eps_for_cutoff = 0.0; 
    double eps_for_ope = 0.0; 

    comp sigb1 = sigma_b_plus(a0_1, m1, m2);

    comp qb = qb_i(En, sigb1, m1);

    double sigp_initial = 3.3; 
    double sigp_final = 4.0; 
    double sigp_points = 200000; 
    double del_sigp = std::abs(sigp_initial - sigp_final)/sigp_points; 

    std::string filename1 = "test_Gs.dat";
    std::ofstream fout; 
    fout.open(filename1.c_str()); 

    for(int i=0; i<sigp_points; ++i)
    {
        comp sigp = sigp_initial + (comp)i*del_sigp + ii*delta;

        comp p = pmom(En, sigp, m1); 

        comp Gs = GS_pk(En, p, qb, m1, m1, m2, eps_for_ope, eps_for_cutoff); 

        std::cout<<sigp.real()<<'\t'
                 <<sigp.imag()<<'\t'
                 <<Gs.real()<<'\t'
                 <<Gs.imag()<<std::endl; 

        fout<<std::setprecision(20)
            <<sigp.real()<<'\t'
            <<sigp.imag()<<'\t'
            <<Gs.real()<<'\t'
            <<Gs.imag()<<std::endl; 
    }
    fout.close(); 
}

//These set works when we set k of (p,k) to qb
//and find the singularities of p
comp beta_x_j(  comp En, 
                comp x, 
                double mj, 
                comp r )//this is the spectator momentum 
{
    comp omgr = omega_func(r, mj);
    comp A = En - omgr; 

    comp res = A*A - x*x*r*r; 

    return res;  
}

comp p_ope_cut_plus(    comp En, 
                        comp x, 
                        comp r, 
                        double mi, 
                        double mj, 
                        double mk, 
                        double eps )
{
    comp ii = {0.0, 1.0}; 
    comp betax = beta_x_j(En, x, mj, r); 
    comp beta0 = beta_x_j(En, 0.0, mj, r); 
    comp beta1 = beta_x_j(En, 1.0, mj, r); 
    comp A = beta1 + ii*eps;
    comp xi = A + (mi*mi - mk*mk);

    comp num = -r*x*xi + std::sqrt(beta0)*std::sqrt(xi*xi - 4.0*mi*mi*betax); 
    comp denom = 2.0*betax; 

    return num/denom; 
}

comp p_ope_cut_minus(   comp En, 
                        comp x, 
                        comp r, 
                        double mi, 
                        double mj, 
                        double mk, 
                        double eps )
{
    comp ii = {0.0, 1.0}; 
    comp betax = beta_x_j(En, x, mj, r); 
    comp beta0 = beta_x_j(En, 0.0, mj, r); 
    comp beta1 = beta_x_j(En, 1.0, mj, r); 
    comp A = beta1 + ii*eps;
    comp xi = A + (mi*mi - mk*mk);

    comp num = -r*x*xi - std::sqrt(beta0)*std::sqrt(xi*xi - 4.0*mi*mi*betax); 
    comp denom = 2.0*betax; 

    return num/denom; 
}

//Here we define the singularity cuts for the case
//when we set p of (p,k) to qb and find the singularities 
//based on k

comp beta_x_i(  comp En, 
                comp x, 
                double mi, 
                comp p  )
{
    comp omgp = omega_func(p, mi);
    comp A = (En - omgp); 
    comp B = x*p;

    return A*A - B*B; 
}

comp r_ope_cut_plus(    comp En,
                        comp x, 
                        comp p, 
                        double mi, 
                        double mj, 
                        double mk, 
                        double eps )
{
    comp ii = {0.0, 1.0}; 
    comp betax = beta_x_i(En, x, mi, p); 
    comp beta0 = beta_x_i(En, 0.0, mi, p); 
    comp beta1 = beta_x_i(En, 1.0, mi, p); 
    comp A = beta1 + ii*eps; 
    comp xi = A + (mj*mj - mk*mk);

    comp num = -p*x*xi + std::sqrt(beta0)*std::sqrt(xi*xi - 4.0*mj*mj*betax);
    comp denom = 2.0*betax; 

    return num/denom; 
}

comp r_ope_cut_minus(   comp En,
                        comp x, 
                        comp p, 
                        double mi, 
                        double mj, 
                        double mk, 
                        double eps )
{
    comp ii = {0.0, 1.0}; 
    comp betax = beta_x_i(En, x, mi, p); 
    comp beta0 = beta_x_i(En, 0.0, mi, p); 
    comp beta1 = beta_x_i(En, 1.0, mi, p); 
    comp A = beta1 + ii*eps; 
    comp xi = A + (mj*mj - mk*mk);

    comp num = -p*x*xi - std::sqrt(beta0)*std::sqrt(xi*xi - 4.0*mj*mj*betax);
    comp denom = 2.0*betax; 

    return num/denom; 
}

void test_kernel_singularities()
{
    comp ii = {0.0,1.0};
    double pi = std::acos(-1.0); 

    double m1 = 1.0; 
    double m2 = 0.9;

    comp En = std::sqrt(8.2); 

    double a0_1 = 2.0; 
    double a0_2 = 2.0; 

    double eps_for_cutoff = 0.0; 
    double eps_for_m2k = 0.0; 
    double eps_for_ope = 0.0; 

    comp sigb11 = sigma_b_plus(a0_1, m1, m2); 
    comp sigb11_1 = sigma_b_minus(a0_1, m1, m2); 

    comp sigb12 = sigma_b_plus(a0_1, m2, m1); 
    comp sigb21 = sigma_b_plus(a0_2, m1, m1); 
    comp sigb22 = sigma_b_minus(a0_2, m1, m1); 
    //std::cout<<sigb11<<'\t'<<sigb11_1<<'\t'<<sigb12<<'\t'<<sigb21<<'\t'<<sigb22<<std::endl; 
    //abort(); 

    comp qb11 = qb_i(En, sigb11, m1); 
    comp qb12 = qb_i(En, sigb12, m1); 
    comp qb21 = qb_i(En, sigb21, m2); 

    comp pole_for_M2_11 = qb_i(En, sigb11, m1); 
    comp pole_for_M2_12 = qb_i(En, sigb11_1, m1); 
    comp pole_for_M2_2 = qb_i(En, sigb21, m2); 

    double xinitial = -1.0; 
    double xfinal = 1.0; 
    double xpoints = 1000.0; 
    double delx = std::abs(xinitial - xfinal)/xpoints; 

    double mi, mj, mk; 

    std::string filename1 = "p_ope_cut.dat";
    std::string filename2 = "omega_cut.dat";
    std::ofstream fout, fout1;
    fout.open(filename1.c_str()); 


    for(int i=0; i<xpoints; ++i)
    {
        double x = xinitial + i*delx; 

        mi = m1; 
        mj = m1; 
        mk = m2; 

        comp pcutplus_qb11 = p_ope_cut_plus(En, x, qb11, mi, mj, mk, eps_for_ope);
        comp pcutminus_qb11 = p_ope_cut_minus(En, x, qb11, mi, mj, mk, eps_for_ope);

        mi = m1; 
        mj = m2; 
        mk = m1; 

        comp pcutplus_qb12 = p_ope_cut_plus(En, x, qb12, mi, mj, mk, eps_for_ope);
        comp pcutminus_qb12 = p_ope_cut_minus(En, x, qb12, mi, mj, mk, eps_for_ope);

        mi = m2; 
        mj = m1; 
        mk = m1; 
        comp pcutplus_qb21 = p_ope_cut_plus(En, x, qb21, mi, mj, mk, eps_for_ope);
        comp pcutminus_qb21 = p_ope_cut_minus(En, x, qb21, mi, mj, mk, eps_for_ope);

        fout<<setprecision(20)
            <<En.real()<<'\t'
            <<En.imag()<<'\t'
            <<pole_for_M2_11.real()<<'\t'
            <<pole_for_M2_11.imag()<<'\t'
            <<pole_for_M2_12.real()<<'\t'
            <<pole_for_M2_12.imag()<<'\t'
            <<pole_for_M2_2.real()<<'\t'
            <<pole_for_M2_2.imag()<<'\t'
            <<pcutplus_qb11.real()<<'\t'
            <<pcutplus_qb11.imag()<<'\t'
            <<pcutminus_qb11.real()<<'\t'
            <<pcutminus_qb11.imag()<<'\t'
            <<pcutplus_qb12.real()<<'\t'
            <<pcutplus_qb12.imag()<<'\t'
            <<pcutminus_qb12.real()<<'\t'
            <<pcutminus_qb12.imag()<<'\t'
            <<pcutplus_qb21.real()<<'\t'
            <<pcutplus_qb21.imag()<<'\t'
            <<pcutminus_qb21.real()<<'\t'
            <<pcutminus_qb21.imag()<<std::endl;

        std::cout<<En<<'\t'
                 <<qb11<<'\t'
                 <<qb12<<'\t'
                 <<qb21<<'\t'
                 <<pcutplus_qb11<<'\t'
                 <<pcutminus_qb11<<'\t'
                 <<pcutplus_qb12<<'\t'
                 <<pcutminus_qb12<<'\t'
                 <<pcutplus_qb21<<'\t'
                 <<pcutminus_qb21<<std::endl;


    }

    fout.close(); 

    


}

/* We copied the old pcut plus for identical particles 
to check the newly built one */

/*
comp omega_comp(    comp p,
                    double m    )
{
    return sqrt(p*p + m*m);
}

comp betaxqs(   double x,
                comp s,
                comp q,
                double m      )
{
    return (omega_comp(q,m) - sqrt(s))*(omega_comp(q,m) - sqrt(s)) - x*x*q*q;
}

comp betaxqs_comp(   comp x,
                comp s,
                comp q,
                double m      )
{
    //return (omega_comp(q,m) - sqrt(s))*(omega_comp(q,m) - sqrt(s)) - x*x*q*q;
    return pow(omega_comp(q,m) - sqrt(s),2.0) - pow(x*q,2.0);

}

comp sqbeta0qs(     comp s,
                    comp q,
                    double m    )
{
    return omega_comp(q,m) - sqrt(s);
}

comp pcut_plus(     double x,
                    comp s,
                    comp q,
                    double m,
                    double eps  )
{
    comp ii = {0.0,1.0};
    comp firstterm = -q*x*betaxqs(1.0,s,q,m);
    comp secondterm = sqbeta0qs(s,q,m);
    comp thirdtermsq = (betaxqs(1.0,s,q,m)+ii*eps)*(betaxqs(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs(x,s,q,m);
    comp thirdterm;
    
    if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = mysqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)>=0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)>=0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) - ii*abs(imag(thirdtermsq)));
    }
    else if(real(thirdtermsq)<0.0 && imag(thirdtermsq)<0.0)
    {
        thirdterm = ii*sqrt(abs(real(thirdtermsq)) + ii*abs(imag(thirdtermsq)));
    }
    comp forthterm = 2.0*betaxqs(x,s,q,m);

    return (firstterm + secondterm*thirdterm)/forthterm;
}

comp pcut_plus_comp(     comp x,
                    comp s,
                    comp q,
                    double m,
                    double eps  )
{
    comp ii = {0.0,1.0};
    if(imag(x)==0.0)
    {
        x = real(x) + ii*1.0e-15;
    }
    comp firstterm = -q*x*betaxqs_comp(1.0,s,q,m);
    //cout<<"firstterm:"<<firstterm<<endl;
    comp secondterm = sqbeta0qs(s,q,m);
    //cout<<"secondterm:"<<secondterm<<endl;
    comp thirdtermsq = (betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m);
    comp thirdterm;
    if(imag(thirdtermsq)>0.0)
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
    }
    else 
    {
        if(real(thirdtermsq)<=0.0)
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
        else 
        {
            thirdterm = sqrt((betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m));
        }
    }
    //cout<<"p_thirdtermsq:"<<(betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps) - 4.0*m*m*betaxqs_comp(x,s,q,m)<<endl;
    //cout<<"thirdterm"<<thirdterm<<endl;
    //cout<<"3.1part:"<<(betaxqs_comp(1.0,s,q,m)+ii*eps)*(betaxqs_comp(1.0,s,q,m)+ii*eps)<<endl;
    //cout<<"3.2part:"<<4.0*m*m*betaxqs_comp(x,s,q,m)<<endl;
    comp forthterm = 2.0*betaxqs_comp(x,s,q,m);
    //cout<<"forthterm:"<<forthterm<<endl;
    //cout<<"res:"<<(firstterm + secondterm*thirdterm)/forthterm<<endl;
    return (firstterm + secondterm*thirdterm)/forthterm;
}

*/
/*
void compare_pcuts()
{
    comp ii = {0.0,1.0};
    double pi = std::acos(-1.0); 

    double m1 = 1.0; 
    double m2 = 1.0;

    double x = -1.0; 

    comp En = std::sqrt(8.3); 

    double a0_1 = 2.0; 
    double a0_2 = 2.0; 

    double eps_for_cutoff = 0.0; 
    double eps_for_m2k = 0.0; 
    double eps_for_ope = 0.0; 

    comp sigb11 = sigma_b_plus(a0_1, m1, m2); 
    comp sigb12 = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2 = sigma_b_plus(a0_1, m1, m1); 

    comp qb11 = qb_i(En, sigb11, m1); 
    comp qb12 = qb_i(En, sigb12, m1); 
    comp qb2 = qb_i(En, sigb2, m2); 

    double xinitial = -1.0; 
    double xfinal = 1.0; 
    double xpoints = 1000.0; 
    double delx = std::abs(xinitial - xfinal)/xpoints; 

    double mi, mj, mk; 

    
    mi = m1; 
    mj = m1; 
    mk = m2;

    comp pcutplus_qb11 = p_ope_cut_plus(En, x, qb11, mi, mj, mk, eps_for_ope);
    comp pcutminus_qb11 = p_ope_cut_minus(En, x, qb11, mi, mj, mk, eps_for_ope);

    comp pcutplus_qb12 = p_ope_cut_plus(En, x, qb12, mi, mj, mk, eps_for_ope);
    comp pcutminus_qb12 = p_ope_cut_minus(En, x, qb12, mi, mj, mk, eps_for_ope);

    mi = m2; 
    mj = m1; 
    mk = m1; 
    comp pcutplus_qb2 = p_ope_cut_plus(En, x, qb2, mi, mj, mk, eps_for_ope);
    comp pcutminus_qb2 = p_ope_cut_minus(En, x, qb2, mi, mj, mk, eps_for_ope);

    comp pcut_old = pcut_plus_comp(  x, En*En, qb11, m1, eps_for_ope);
    std::cout<<pcutplus_qb11<<'\t'<<pcut_old<<std::endl; 

}
*/


void test_Gs_surface()
{
    comp ii = {0.0, 1.0}; 
    double pi = std::acos(-1.0); 
    double m1 = 1.0; 
    double m2 = 0.9; 

    comp En = std::sqrt(8.2); 
    double total_P = 0.0; 

    double delta = 1.0e-4; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0;
    double eps_for_cutoff = 0.0; 
    double eps_for_ope = 0.0; 
    double eps_for_m2k = 0.0; 

    comp sigb1 = sigma_b_plus(a0_1, m1, m2);
    comp sigb2 = sigma_b_plus(a0_2, m1, m1); 


    comp qb = qb_i(En, sigb1, m1);
    comp qb2 = qb_i(En, sigb2, m2); 

    double sigp_initial = 3.3; 
    double sigp_final = 4.0; 
    double sigp_points = 200000; 
    double del_sigp = std::abs(sigp_initial - sigp_final)/sigp_points; 

    double preal_initial = -1.5; 
    double preal_final = 1.51; 
    double pimag_initial = -1.5; 
    double pimag_final = 1.51; 

    double p_points = 253; 
    double del_preal = std::abs(preal_initial - preal_final)/p_points; 
    double del_pimag = std::abs(pimag_initial - pimag_final)/p_points;

    double mi, mj, mk; 
    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    std::string filename1 = "test_Gs_surface.dat";
    std::ofstream fout; 
    fout.open(filename1.c_str()); 

    for(int i=0; i<p_points; ++i)
    {
        for(int j=0; j<p_points; ++j)
        {
            comp preal = preal_initial + (comp)i*del_preal; 
            comp pimag = pimag_initial + (comp)j*del_pimag; 

            comp p = preal + ii*pimag; 

            //(i,j) = 1 1 ; k = 2
            mi = m1; 
            mj = m1; 
            mk = m2; 
            
            comp Gs11 = GS_pk(En, qb, p, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
            comp M2k_11 = M2k_ERE(eta_1, En, p, total_P, a0_1, 0.0, mj, mk, mi, eps_for_m2k);
            comp omg_11 = omega_func(p, mj);

            comp kern1 = (p*p/(std::pow(2.0*pi,2.0)*omg_11))*Gs11*M2k_11; 

            //(i,j) = 1 2 ; k = 1
            mi = m1; 
            mj = m2; 
            mk = m1; 
            
            comp Gs12 = GS_pk(En, qb, p, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
            comp M2k_12 = M2k_ERE(eta_2, En, p, total_P, a0_2, 0.0, mj, mk, mi, eps_for_m2k);
            comp omg_12 = omega_func(p, mj);

            comp kern2 = (p*p/(std::pow(2.0*pi,2.0)*omg_12))*Gs12*M2k_12; 

            //(i,j) = 2 1 ; k = 1
            mi = m2; 
            mj = m1; 
            mk = m1; 
            
            comp Gs21 = GS_pk(En, qb2, p, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
            comp M2k_21 = M2k_ERE(eta_1, En, p, total_P, a0_1, 0.0, mj, mk, mi, eps_for_m2k);
            comp omg_21 = omega_func(p, mj);

            comp kern3 = (p*p/(std::pow(2.0*pi,2.0)*omg_21))*Gs21*M2k_21; 


            fout<<std::setprecision(20)
                <<En.real()<<'\t'
                <<En.imag()<<'\t'
                <<preal.real()<<'\t'
                <<pimag.real()<<'\t'
                <<qb.real()<<'\t'
                <<qb.imag()<<'\t'
                <<kern1.real()<<'\t'
                <<kern1.imag()<<'\t'
                <<kern2.real()<<'\t'
                <<kern2.imag()<<'\t'
                <<kern3.real()<<'\t'
                <<kern3.imag()<<std::endl; 

            std::cout<<std::setprecision(20)
                     <<preal<<'\t'
                     <<pimag<<'\t'
                     <<qb.real()<<'\t'
                     <<qb.imag()<<'\t'
                     <<kern1.real()<<'\t'
                     <<kern1.imag()<<'\t'
                     <<kern2.real()<<'\t'
                     <<kern2.imag()<<'\t'
                     <<kern3.real()<<'\t'
                     <<kern3.imag()<<std::endl; 

        }
         
    }
    fout.close(); 
}

void test_Gs_surface_1()
{
    comp ii = {0.0, 1.0}; 
    double pi = std::acos(-1.0); 
    double m1 = 1.0; 
    double m2 = 0.9; 

    comp En = std::sqrt(8.2); 
    double total_P = 0.0; 

    double delta = 1.0e-4; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0;
    double eps_for_cutoff = 0.0; 
    double eps_for_ope = 0.0; 
    double eps_for_m2k = 0.0; 

    comp sigb1 = sigma_b_plus(a0_1, m1, m2);
    comp sigb1_deep = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2 = sigma_b_plus(a0_2, m1, m1); 


    comp qb1 = qb_i(En, sigb1, m1);
    comp qb1_deep = qb_i(En, sigb1_deep, m1); 
    comp qb2 = qb_i(En, sigb2, m2); 

    double sigp_initial = 3.3; 
    double sigp_final = 4.0; 
    double sigp_points = 200000; 
    double del_sigp = std::abs(sigp_initial - sigp_final)/sigp_points; 

    double preal_initial = -1.5;//1.24;//-1.5;//0.5;//-1.5;//-0.5;//-1.5; 
    double preal_final = 1.51;//1.28;//1.51;//1.41;//1.51;//0.51;//1.51; 
    double pimag_initial = -1.5;//-0.02;//-1.5;//-0.5;//-1.5;//-0.6;//-1.5; 
    double pimag_final = 1.51;//0.021;//1.51;//0.51;//1.51;//-0.51;//1.51; 

    double p_points = 253; 
    double del_preal = std::abs(preal_initial - preal_final)/p_points; 
    double del_pimag = std::abs(pimag_initial - pimag_final)/p_points;

    double mi, mj, mk; 
    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    std::string filename1 = "test_Gs_surface_1.dat";
    std::ofstream fout; 
    fout.open(filename1.c_str()); 

    for(int i=0; i<p_points; ++i)
    {
        for(int j=0; j<p_points; ++j)
        {
            comp preal = preal_initial + (comp)i*del_preal; 
            comp pimag = pimag_initial + (comp)j*del_pimag; 

            comp p = preal + ii*pimag; 

            //(i,j) = 2 1 ; k = 1
            mi = m2; 
            mj = m1; 
            mk = m1; 
            
            comp Gs21 = GS_pk(En, qb2, p, mi, mj, mk, eps_for_ope, eps_for_cutoff); 
            comp M2k_21 = M2k_ERE(eta_1, En, p, total_P, a0_1, 0.0, mj, mk, mi, eps_for_m2k);
            comp omg_21 = omega_func(p, mj);

            comp kern3 = (p*p/(std::pow(2.0*pi,2.0)*omg_21))*Gs21*M2k_21; 


            fout<<std::setprecision(20)
                <<preal.real()<<'\t'
                <<pimag.real()<<'\t'
                <<kern3.real()<<'\t'
                <<kern3.imag()<<std::endl; 

            std::cout<<std::setprecision(20)
                     <<preal<<'\t'
                     <<pimag<<'\t'
                     <<kern3.real()<<'\t'
                     <<kern3.imag()<<std::endl; 

        }
         
    }
    fout.close(); 

    std::ofstream fout1; 
    std::string filename2 = "kernel_singularities_1.dat";
    fout1.open(filename2.c_str()); 

    double xinitial = -1.0; 
    double xfinal = 1.0; 
    double xpoints = 1000.0; 
    double delx = std::abs(xinitial - xfinal)/xpoints;

    for(int i=0; i<xpoints; ++i)
    {
        double x = xinitial + i*delx; 

        mi = m2; 
        mj = m1; 
        mk = m1; 

        comp pcutplus_val = r_ope_cut_plus(En, x, qb2, mi, mj, mk, eps_for_ope); 
        comp pcutminus_val = r_ope_cut_minus(En, x, qb2, mi, mj, mk, eps_for_ope); 

        fout1<<std::setprecision(20)<<'\t'
             <<pcutplus_val.real()<<'\t'
             <<pcutplus_val.imag()<<'\t'
             <<pcutminus_val.real()<<'\t'
             <<pcutminus_val.imag()<<'\t'
             <<qb1.real()<<'\t'
             <<qb1.imag()<<'\t'
             <<qb1_deep.real()<<'\t'
             <<qb1_deep.imag()<<'\t'
             <<qb2.real()<<'\t'
             <<qb2.imag()<<std::endl; 
    }
    fout1.close(); 


}


void compare_Bmats()
{
    double m1 = 1.0; 
    double m2 = 0.9; 

    double a0_1 = -2.0; 
    double a0_2 = -2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    double total_P = 0.0; 
    double r = 0.0; 

    double eps_for_m2k = 0.001; 
    double eps_for_ope = 0.001; 
    double eps_for_cutoff = 0.0; 

    double number_of_points = 500.0; 

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    


    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp En = (En_1 + En_2)/2.0; 

    //qb vals for i=1 and i=2 case:
    comp qb_val1plus = qb_i(En, sigb1plus, m1);
    comp qb_val1minus = qb_i(En, sigb1minus, m1);
    comp qb_val2plus = qb_i(En, sigb2plus, m2);
    comp qb_val2minus = qb_i(En, sigb2minus, m2);

    std::cout<<"qb 1 + = "<<qb_val1plus<<std::endl;
    std::cout<<"qb 1 - = "<<qb_val1minus<<std::endl;
    std::cout<<"qb 2 + = "<<qb_val2plus<<std::endl;
    std::cout<<"qb 2 - = "<<qb_val2minus<<std::endl;

    comp kmax_for_m1 = pmom(En, 0.0, m1); 
    comp kmax_for_m2 = pmom(En, 0.0, m2); 
    comp epsilon_for_kvec = 1.0e-5; 

    std::cout<<"kmax_for_m1 = "<<kmax_for_m1<<std::endl; 
    std::cout<<"kmax_for_m2 = "<<kmax_for_m2<<std::endl; 

    std::vector<comp> pvec_for_m1m2;
    std::vector<comp> weights_for_pvec_for_m1m2; 
    std::vector<comp> kvec_for_m1m1; 
    std::vector<comp> weights_for_kvec_for_m1m1; 

    flavor_based_momentum_vector(pvec_for_m1m2, weights_for_pvec_for_m1m2, En, m1, number_of_points);
    flavor_based_momentum_vector(kvec_for_m1m1, weights_for_kvec_for_m1m1, En, m2, number_of_points);

    std::cout<<"pvec for m1m2 = "<<std::endl; 
    vec_printer(pvec_for_m1m2); 
    std::cout<<"kvec for m1m1 = "<<std::endl; 
    vec_printer(kvec_for_m1m1); 
    std::cout<<"=============================="<<std::endl;

    int size1 = pvec_for_m1m2.size(); 
    int size2 = kvec_for_m1m1.size(); 
    int total_size = size1 + size2; 
    Eigen::MatrixXcd B_mat1(total_size, total_size); 
    Eigen::MatrixXcd B_mat2(total_size, total_size); 

    test_Bmat_2plus1_system_ERE( B_mat1, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2 );
    
    test_Bmat_GMtilde_multiplied_2plus1_system_ERE(B_mat2, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, a0_2); 
    

    std::cout<<"Bmat1 determinant = "<<B_mat1.determinant()<<std::endl;
    std::cout<<"Bmat2 determinant = "<<B_mat2.determinant()<<std::endl; 

    for(int i=0; i<total_size; ++i)
    {
        for(int j=0; j<total_size; ++j)
        {
            comp A = B_mat1(i,j); 
            comp B = B_mat2(i,j); 

            double diff = std::abs(A - B);

            std::cout<<std::setprecision(10); 
            std::cout<<"i:"<<i<<'\t'
                     <<"j:"<<j<<'\t'
                     <<"Bmat1="<<A<<'\t'
                     <<"Bmat2="<<B<<'\t'
                     <<"diff="<<diff<<std::endl; 
        }
    }


}



void plot_single_integral_equation_components()
{
    double m1 = 1.0; 
    double m2 = 0.9; 

    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 

    double total_P = 0.0; 
    double r = 0.0; 

    

    double number_of_points = 1000.0; 
    double number_of_points_1 = 1000.0; 

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;

    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    


    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp En = std::sqrt(8.2);//(En_1 + En_2)/2.0; 
    

    

    //qb vals for i=1 and i=2 case:
    comp qb_val1plus = qb_i(En, sigb1plus, m1);
    comp qb_val1minus = qb_i(En, sigb1minus, m1);
    comp qb_val2plus = qb_i(En, sigb2plus, m2);
    comp qb_val2minus = qb_i(En, sigb2minus, m2);

    std::cout<<"qb 1 + = "<<qb_val1plus<<std::endl;
    std::cout<<"qb 1 - = "<<qb_val1minus<<std::endl;
    std::cout<<"qb 2 + = "<<qb_val2plus<<std::endl;
    std::cout<<"qb 2 - = "<<qb_val2minus<<std::endl;

    comp kmax_for_m1 = pmom(En, 0.0, m1); 
    comp kmax_for_m2 = pmom(En, 0.0, m2); 
    comp epsilon_for_kvec = 1.0e-5; 

    std::cout<<"kmax_for_m1 = "<<kmax_for_m1<<std::endl; 
    std::cout<<"kmax_for_m2 = "<<kmax_for_m2<<std::endl; 

    std::vector<comp> pvec_for_m1m2;
    std::vector<comp> weights_for_pvec_for_m1m2; 
    std::vector<comp> kvec_for_m1m1; 
    std::vector<comp> weights_for_kvec_for_m1m1; 

    double eta_for_eps = 25.0; 
    double eps_for_m2k = energy_dependent_epsilon(eta_for_eps, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points ); 
    std::cout<<"eps_for_m2k = "<<eps_for_m2k<<std::endl; 
    double eps_for_ope = eps_for_m2k; 
    double eps_for_cutoff = 0.0; 

    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    flavor_based_momentum_vector(pvec_for_m1m2, weights_for_pvec_for_m1m2, En, m1, number_of_points);
    flavor_based_momentum_vector(kvec_for_m1m1, weights_for_kvec_for_m1m1, En, m2, number_of_points_1);

    /*
    std::cout<<"pvec for m1m2 = "<<std::endl; 
    vec_printer(pvec_for_m1m2); 
    std::cout<<"kvec for m1m1 = "<<std::endl; 
    vec_printer(kvec_for_m1m1); 
    std::cout<<"=============================="<<std::endl;
    */

    for(int i=0; i<pvec_for_m1m2.size(); ++i)
    {
        std::cout<<"p="<<pvec_for_m1m2[i]<<'\t'
                 <<"w="<<weights_for_pvec_for_m1m2[i]<<std::endl;
    }
    
    int size1 = pvec_for_m1m2.size(); 
    int size2 = kvec_for_m1m1.size(); 

    int qb1plus_ind, qb1minus_ind, qb2plus_ind; 

    for(int i=0; i<size1-1; ++i)
    {
        double pval1 = pvec_for_m1m2[i].real(); 
        double pval2 = pvec_for_m1m2[i+1].real(); 
        double kval1 = kvec_for_m1m1[i].real(); 
        double kval2 = kvec_for_m1m1[i+1].real(); 

        if(qb_val1plus.real()>=pval1 && qb_val1plus.real()<=pval2)
        {
            qb1plus_ind = i; 
        }

        if(qb_val1minus.real()>=pval1 && qb_val1minus.real()<=pval2)
        {
            qb1minus_ind = i; 
        }

        if(qb_val2plus.real()>=kval1 && qb_val2plus.real()<=kval2)
        {
            qb2plus_ind = i; 
        }

    }

    std::cout<<"q1+ ind = "<<qb1plus_ind<<'\t'
             <<"qb = "<<qb_val1plus<<'\t'
             <<"p_i = "<<pvec_for_m1m2[qb1plus_ind]<<'\t'
             <<"p_i+1 = "<<pvec_for_m1m2[qb1plus_ind+1]<<std::endl; 

    /*
    std::cout<<"q1- ind = "<<qb1minus_ind<<'\t'
             <<"qb = "<<qb_val1minus<<'\t'
             <<"p_i = "<<pvec_for_m1m2[qb1minus_ind]<<'\t'
             <<"p_i+1 = "<<pvec_for_m1m2[qb1minus_ind+1]<<std::endl; 

    
    std::cout<<"q2+ ind = "<<qb2plus_ind<<'\t'
             <<"qb = "<<qb_val2plus<<'\t'
             <<"p_i = "<<kvec_for_m1m1[qb2plus_ind]<<'\t'
             <<"p_i+1 = "<<kvec_for_m1m1[qb2plus_ind+1]<<std::endl; 
 
    */
    int total_size = size1 + size2; 
    Eigen::MatrixXcd B_mat1(total_size, total_size); 
    Eigen::MatrixXcd B_mat2(total_size, total_size); 

    
    //test_Bmat_2plus1_system_ERE( B_mat1, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, a0_2, r0_2 );
    test_Bmat_2plus1_system_ERE_1( B_mat1, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, r0_1, eta_1, a0_2, r0_2, eta_2 );
    
    //test_Bmat_GMtilde_multiplied_2plus1_system_ERE(B_mat2, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, m1, m2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, a0_2); 
    

    std::cout<<"Bmat1 determinant = "<<B_mat1.determinant()<<std::endl;
    //std::cout<<"Bmat2 determinant = "<<B_mat2.determinant()<<std::endl; 
    /*std::cout<<"Bmat1 = "<<std::endl; 
    for(int i=0; i<size1; ++i)
    {
        for(int j=0; j<size1; ++j)
        {
            std::cout<<i<<'\t'
                     <<j<<'\t'
                     <<B_mat1(i,j)<<std::endl; 
        }
    }*/

    //Print and compare the two Bmats we have
    char bmat_compare = 'n'; 
    if(bmat_compare=='y')
    {
    for(int i=0; i<total_size; ++i)
    {
        for(int j=0; j<total_size; ++j)
        {
            comp A = B_mat1(i,j); 
            comp B = B_mat2(i,j); 

            double diff = std::abs(A - B);

            std::cout<<std::setprecision(10); 
            std::cout<<"i:"<<i<<'\t'
                     <<"j:"<<j<<'\t'
                     <<"Bmat1="<<A<<'\t'
                     <<"Bmat2="<<B<<'\t'
                     <<"diff="<<diff<<std::endl; 
        }
    }
    }

    

    std::cout<<B_mat1.rows()<<std::endl; 
    std::cout<<B_mat1.cols()<<std::endl; 


    //Print all values of Bmat for individual i vals 
    //and the surface plot of Bmat 
    char bmat_plots_plus_surface = 'n'; 
    if(bmat_plots_plus_surface=='y')
    {
    std::ofstream fout1; 
    std::string filename;

    std::string filename1 = "B_mat_surface.dat"; 
    std::ofstream fout2; 

    fout2.open(filename1.c_str()); 
    
    
    for(int i=0; i<B_mat1.rows(); ++i)
    {
        int ival = i; 
        filename = "bmat_file_i_"  + std::to_string(i) + ".dat"; 
        fout1.open(filename.c_str()); 

        for(int j=0; j<B_mat1.cols(); ++j)
        {
            comp mom_p; 
            comp mom_k; 
            comp weight_p; 
            comp weight_k; 

            if(i<number_of_points && j<number_of_points)
            {
                mom_p = pvec_for_m1m2[j]; 
                mom_k = pvec_for_m1m2[j]; 
                weight_p = weights_for_pvec_for_m1m2[j]; 
                weight_k = weights_for_pvec_for_m1m2[j]; 
            }
            else if(i<number_of_points && j>=number_of_points)
            {
                int ind = j - number_of_points; 
                mom_p = pvec_for_m1m2[ind]; 
                mom_k = kvec_for_m1m1[ind]; 
                weight_p = weights_for_pvec_for_m1m2[ind]; 
                weight_k = weights_for_kvec_for_m1m1[ind]; 
            }
            else if(i>=number_of_points && j<number_of_points)
            {
                mom_p = kvec_for_m1m1[j]; 
                mom_k = pvec_for_m1m2[j]; 
                weight_p = weights_for_kvec_for_m1m1[j]; 
                weight_k = weights_for_pvec_for_m1m2[j]; 
            }
            else
            {
                int ind = j - number_of_points;
                mom_p = kvec_for_m1m1[ind]; 
                mom_k = kvec_for_m1m1[ind]; 
                weight_p = weights_for_kvec_for_m1m1[ind]; 
                weight_k = weights_for_kvec_for_m1m1[ind]; 
            }

            comp Bmat_result = B_mat1(i,j); 
            std::cout<<"j="<<j<<'\t'
                     <<"p="<<mom_p<<'\t'
                     <<"wt="<<weight_p<<'\t'
                     <<"k="<<mom_k<<'\t'
                     <<"wt="<<weight_k<<'\t'
                     <<"Bmat="<<B_mat1(i,j)<<std::endl; 
            fout1<<j<<'\t'
                 <<mom_p.real()<<'\t'
                 <<mom_p.imag()<<'\t'
                 <<weight_p.real()<<'\t'
                 <<weight_p.imag()<<'\t'
                 <<mom_k.real()<<'\t'
                 <<mom_k.imag()<<'\t'
                 <<weight_k.real()<<'\t'
                 <<weight_k.imag()<<'\t'
                 <<B_mat1(i,j).real()<<'\t'
                 <<B_mat1(i,j).imag()<<std::endl; 
            
            fout2<<i<<'\t'
                 <<j<<'\t'
                 <<mom_p.real()<<'\t'
                 <<mom_p.imag()<<'\t'
                 <<weight_p.real()<<'\t'
                 <<weight_p.imag()<<'\t'
                 <<mom_k.real()<<'\t'
                 <<mom_k.imag()<<'\t'
                 <<weight_k.real()<<'\t'
                 <<weight_k.imag()<<'\t'
                 <<B_mat1(i,j).real()<<'\t'
                 <<B_mat1(i,j).imag()<<std::endl; 
        }
        fout1.close(); 
    }
    fout2.close(); 
    }

    Eigen::MatrixXcd negGSpqb(total_size, 2); 
    negative_GSpqb_mat(negGSpqb, En, pvec_for_m1m2, kvec_for_m1m1, weights_for_pvec_for_m1m2, weights_for_kvec_for_m1m1, qb_val1plus, m1, m2, eps_for_ope, eps_for_cutoff, total_P); 


    char gmat_plot = 'n';
    if(gmat_plot=='y')
    {

        std::ofstream fout3; 
        std::string filename3 = "gmat_plot.dat";

        fout3.open(filename3.c_str()); 


        int gmat_row = negGSpqb.rows(); 
        int gmat_col = negGSpqb.cols(); 

        std::cout<<"gmat = "<<std::endl; 
        for(int i=0; i<gmat_row; ++i)
        {
            for(int j=0; j<gmat_col; ++j)
            {
                if(j==0)
                {
                    std::cout<<i<<'\t'
                             <<j<<'\t'
                             <<negGSpqb(i,j)<<'\t';
                    fout3<<i<<'\t'
                         <<negGSpqb(i,j).real()<<'\t'
                         <<negGSpqb(i,j).imag()<<'\t';
                }
                else 
                {
                    std::cout<<j<<'\t'
                             <<negGSpqb(i,j)<<std::endl;
                    fout3<<negGSpqb(i,j).real()<<'\t'
                         <<negGSpqb(i,j).imag()<<std::endl;
                }
                 
            }
        }
        fout3.close(); 

    }

    Eigen::MatrixXcd dmat; 
    double relerr; 
    LinearSolver_2(B_mat1, dmat, negGSpqb, relerr); 
    std::cout<<"relative error = "<<relerr<<std::endl; 


    //Here we plot the dmat for a single 
    //kinematic energy, En. We plot dmat 
    //against the index values of the matrix
    //to check for any arbitrary singularities 
    //arising in the system. 
    char dmat_plot = 'n';
    if(dmat_plot=='y')
    {

        std::ofstream fout4; 
        std::string filename4 = "dmat.dat"; 

        fout4.open(filename4.c_str()); 

        int dmat_row = dmat.rows(); 
        int dmat_col = dmat.cols(); 

        for(int i=0; i<dmat_row; ++i)
        {
            for(int j=0; j<dmat_col; ++j)
            {
                if(j==0)
                {
                    std::cout<<i<<'\t'
                             <<j<<'\t'
                             <<dmat(i,j)<<'\t';
                    fout4<<i<<'\t'
                         <<dmat(i,j).real()<<'\t'
                         <<dmat(i,j).imag()<<'\t';
                }
                else 
                {
                    std::cout<<j<<'\t'
                             <<dmat(i,j)<<std::endl;
                    
                    fout4<<dmat(i,j).real()<<'\t'
                         <<dmat(i,j).imag()<<std::endl; 
                }
            }
        }
        fout4.close(); 
    }

    /* Build constituents for degenerate case from 
    the 2+1 system, solve the eqn to check results */

    char degen_test = 'n';

    if(degen_test=='y')
    {
        Eigen::MatrixXcd degen_Bmat(size1, size1); 
        Eigen::VectorXcd degen_Gvec(size1); 
        Eigen::VectorXcd degen_dvec(size1); 

        for(int i=0; i<size1; ++i)
        {
            for(int j=0; j<size1; ++j)
            {
                degen_Bmat(i,j) = B_mat1(i,j); 
            }
        }

        for(int i=0; i<size1; ++i)
        {
            degen_Gvec(i) = negGSpqb(i,0);
        }

        LinearSolver_3(degen_Bmat, degen_dvec, degen_Gvec, relerr);

        std::cout<<"degen dvec = "<<std::endl; 

        for(int i=0; i<size1; ++i)
        {
            std::cout<<i<<'\t'<<degen_dvec(i)<<std::endl; 
        }

        comp degen_dqq = {0.0,0.0}; 
        test_degen_dqq_interpolator(degen_dvec, En, m1, pvec_for_m1m2, weights_for_pvec_for_m1m2, qb_val1plus, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, number_of_points, degen_dqq);

        comp gval = gfunc_i(eta_1, sigb1plus, m1, m1); 
        comp gsq = gval*gval; 
        comp mphib = gsq*degen_dqq; 
        comp rhopb = rhophib(qb_val1plus, En); 

        double error = std::abs(((1.0/mphib).imag() + rhopb)/rhopb )*100.0; 

        std::cout<<"g = "<<gval<<std::endl; 
        std::cout<<"Im(1/mphib) = "<<(1.0/mphib).imag()<<std::endl; 
        std::cout<<"rho phib = "<<rhopb<<std::endl; 
        std::cout<<"error = "<<error<<std::endl; 
    }
    /*---------------------------------------------*/

    
    comp dmat11_qq_val = dmat(qb1plus_ind, 0); 
    comp g1_val = gfunc_i(eta_1, sigb1plus, m1, m2); 
    comp mphib_11_test_val = g1_val*g1_val*dmat11_qq_val; 

    comp rhophib_11 = rhophib(qb_val1plus, En); 

    double diff = std::abs(((1.0/mphib_11_test_val).imag() + rhophib_11)/rhophib_11)*100.0;
    
    std::cout<<"g = "<<g1_val<<std::endl; 
    std::cout<<"mphib = "<<mphib_11_test_val<<std::endl; 
    std::cout<<"Im(1/mphib) = "<<(1.0/mphib_11_test_val).imag()<<std::endl;
    std::cout<<"rhophib = "<<rhophib_11<<std::endl;
    std::cout<<"delrho = "<<diff<<"%"<<std::endl;
    
    Eigen::MatrixXcd dmat_qq_intpltd;
    test_dqq_interpolator_1(dmat_qq_intpltd, dmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb_val1plus, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_1, 0.0, eta_1, a0_2, 0.0, eta_2, number_of_points, 'n');

    comp dmat11_qq_val_intpltd = dmat_qq_intpltd(0,0); 

    g1_val = gfunc_i(eta_1, sigb1plus, m1, m2); 
    mphib_11_test_val = g1_val*g1_val*dmat11_qq_val_intpltd; 

    diff = std::abs(((1.0/mphib_11_test_val).imag() + rhophib_11)/rhophib_11)*100.0;
    
    std::cout<<"p - q = "<<pvec_for_m1m2[qb1plus_ind] - qb_val1plus<<std::endl; 
    std::cout<<"g = "<<g1_val<<std::endl; 
    std::cout<<"mphib = "<<mphib_11_test_val<<std::endl; 
    std::cout<<"Im(1/mphib) = "<<(1.0/mphib_11_test_val).imag()<<std::endl;
    std::cout<<"rhophib = "<<rhophib_11<<std::endl;
    std::cout<<"delrho = "<<diff<<"%"<<std::endl;
    
    //Changed multiple things, commented several
    //check the whole code, when we set a>0, and 
    //eta = 1, we get something closer to previous 
    //results
    //check with old code if they generate the same
    //numbers or not
    
   

}

void test_dpqb_vs_N_building_2()
{
    double m1 = 1.0;
    double m2 = 0.9;//0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 
    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    //int number_of_points = 250;  

    

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;


    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
    std::string filename =  "mphib_BFmethod_vs_N_En_eta_25_3_forpaper.dat";

    fout.open(filename.c_str()); 

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp s = 8.2; 
    comp En = std::sqrt(s);//(En_1 + En_2)/2.0; 

    //qb vals for i=1 and i=2 case:
    comp qb_val1plus = qb_i(En, sigb1plus, m1);
    comp qb_val1minus = qb_i(En, sigb1minus, m1);
    comp qb_val2plus = qb_i(En, sigb2plus, m2);
    comp qb_val2minus = qb_i(En, sigb2minus, m2);

    comp kmax_for_m1 = pmom(En, 0.0, m1); 
    comp kmax_for_m2 = pmom(En, 0.0, m2); 
    comp epsilon_for_kvec = 1.0e-5; 

    std::cout<<"kmax_for_m1 = "<<kmax_for_m1<<std::endl; 
    std::cout<<"kmax_for_m2 = "<<kmax_for_m2<<std::endl; 

    double eta_for_eps = 25; 

    double eps_for_m2k_1 = 0.0;// = energy_dependent_epsilon(eta_for_eps, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points ); 
    double eps_for_m2k_2 = 0.0; 
    double eps_for_ope_1 = eps_for_m2k_1; 
    double eps_for_ope_2 = eps_for_m2k_2; 
    double eps_for_cutoff_1 = 0.0; 
    double eps_for_cutoff_2 = 0.0; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //N range = [50, 500] and [500, 5000]
    int N_initial = 50; 
    int N_final = 500; 
    int N_points = 0; 
    int del_N = 1;//std::abs(N_initial - N_final)/N_points; 
    N_points = std::abs((N_initial - N_final)/del_N);

    for(int i=0; i<(int)N_points; ++i)
    {
        //comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        int number_of_points = N_initial + i*del_N; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = 25; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k_1 = energy_dependent_epsilon(eta, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points);
        eps_for_m2k_2 = eps_for_m2k_1; //energy_dependent_epsilon(eta, En, qb_val2plus, sigb2plus, kmax_for_m2, m2, number_of_points);
        eps_for_ope_1 = eps_for_m2k_1; 
        eps_for_ope_2 = eps_for_m2k_2; 
        test_dpqb_solver_ERE_2(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k_1, eps_for_m2k_2, eps_for_ope_1, eps_for_ope_2, eps_for_cutoff_1, eps_for_cutoff_2, total_P, a0_1, r0_1, eta_1, a0_2, r0_2, eta_2, number_of_points, 'n'); 
        test_dqq_interpolator_2(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k_1, eps_for_m2k_2, eps_for_ope_1, eps_for_ope_2, eps_for_cutoff_1, eps_for_cutoff_2, total_P, a0_1, r0_1, eta_1, a0_2, r0_2, eta_2, number_of_points, 'n');

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        //double eta_1 = 1.0; 
        comp gval = gfunc_i(eta_1, sigb1plus, m1, m2);//gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb, En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(20)<<En.real()<<'\t'<<(En*En).real()<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<'\t'<<number_of_points<<std::endl; 
        //fout<<std::setprecision(20)<<std::real(En)<<'\t'<<std::real(En*En)<<'\t'<<rhophib_mphib.real()<<'\t'<<rhophib_mphib.imag()<<'\t'<<diff<<'\t'<<number_of_points<<std::endl; 
        
        fout<<std::setprecision(20)<<number_of_points<<'\t'
                                   <<diff<<'\t'
                                   <<rhophib_mphib.real()<<'\t'
                                   <<rhophib_mphib.imag()<<std::endl; 
    
    }
    fout.close(); 
}

/* This sets N = 50 and tests the eta dependence */
void test_dpqb_vs_N_50_eta_dependence()
{
    double m1 = 1.0;
    double m2 = 0.9;//0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 
    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    //int number_of_points = 250;  

    

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;


    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
   

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp s = 8.2; 
    comp En = std::sqrt(s);//(En_1 + En_2)/2.0; 

    //qb vals for i=1 and i=2 case:
    comp qb_val1plus = qb_i(En, sigb1plus, m1);
    comp qb_val1minus = qb_i(En, sigb1minus, m1);
    comp qb_val2plus = qb_i(En, sigb2plus, m2);
    comp qb_val2minus = qb_i(En, sigb2minus, m2);

    comp kmax_for_m1 = pmom(En, 0.0, m1); 
    comp kmax_for_m2 = pmom(En, 0.0, m2); 
    comp epsilon_for_kvec = 1.0e-5; 

    std::cout<<"kmax_for_m1 = "<<kmax_for_m1<<std::endl; 
    std::cout<<"kmax_for_m2 = "<<kmax_for_m2<<std::endl; 

    double eta_for_eps = 25; 

    double eps_for_m2k_1 = 0.0;// = energy_dependent_epsilon(eta_for_eps, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points ); 
    double eps_for_m2k_2 = 0.0; 
    double eps_for_ope_1 = eps_for_m2k_1; 
    double eps_for_ope_2 = eps_for_m2k_2; 
    double eps_for_cutoff_1 = 0.0; 
    double eps_for_cutoff_2 = 0.0; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //N range = [50, 500] and [500, 5000]
    int N_initial = 50; 
    int N_final = 2000; 
    int N_points = 0; 
    int del_N = 50;//std::abs(N_initial - N_final)/N_points; 
    N_points = std::abs((N_initial - N_final)/del_N);

    double eta_initial = 2.0; 
    double eta_final = 50.0;
    double eta_points = 1000.0; 
    double del_eta = std::abs(eta_final - eta_initial)/eta_points;  

    int number_of_points = 250; //N_initial + i*del_N; 

    std::string filename =  "mphib_above_phibthreshold_a" + std::to_string((int)a0_1) + "_vs_eta_for_N_" + std::to_string((int)number_of_points) + ".dat";

    fout.open(filename.c_str()); 
    for(int i=0; i<(int)eta_points; ++i)
    {
        //comp En = En_initial + ((comp)i)*delEn; 
        //comp qb; 
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double eta = eta_initial + i*del_eta; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k_1 = energy_dependent_epsilon(eta, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points);
        eps_for_m2k_2 = energy_dependent_epsilon(eta, En, qb_val2plus, sigb2plus, kmax_for_m2, m2, number_of_points);
        eps_for_ope_1 = eps_for_m2k_1; 
        eps_for_ope_2 = eps_for_m2k_2; 
        test_dpqb_solver_ERE_2(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k_1, eps_for_m2k_2, eps_for_ope_1, eps_for_ope_2, eps_for_cutoff_1, eps_for_cutoff_2, total_P, a0_1, r0_1, eta_1, a0_2, r0_2, eta_2, number_of_points, 'n'); 
        test_dqq_interpolator_2(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k_1, eps_for_m2k_2, eps_for_ope_1, eps_for_ope_2, eps_for_cutoff_1, eps_for_cutoff_2, total_P, a0_1, r0_1, eta_1, a0_2, r0_2, eta_2, number_of_points, 'n');

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        //double eta_1 = 1.0; 
        comp gval = gfunc_i(eta_1, sigb1plus, m1, m2);//gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb, En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        
        std::cout<<std::setprecision(5)
                 <<"En="<<En.real()<<'\t'
                 <<"s="<<(En*En).real()<<'\t'
                 <<"rhoM="<<rhophib_mphib<<'\t'
                 <<"err="<<diff<<'\t'
                 <<"N="<<number_of_points<<'\t'
                 <<"eta="<<eta<<std::endl; 
        fout<<std::setprecision(20)
            <<std::real(En)<<'\t'
            <<std::real(En*En)<<'\t'
            <<rhophib_mphib.real()<<'\t'
            <<rhophib_mphib.imag()<<'\t'
            <<diff<<'\t'
            <<eta<<std::endl; 
        
    
    }
    fout.close(); 
}


/* This is the new iteration of the delta_rhophib
density making plot, we use omp here to accelerate 
the calculation */
void test_delta_rhophib_density_with_omp()
{
    double m1 = 1.0;
    double m2 = 0.9;//0.99999999990;//0.999; 
    double a0_1 = 2.0; 
    double a0_2 = 2.0; 
    double r0_1 = 0.0; 
    double r0_2 = 0.0; 
    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    //int number_of_points = 250;  

    

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;


    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
   

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp s = 8.2; 
    comp En = std::sqrt(s);//(En_1 + En_2)/2.0; 

    //qb vals for i=1 and i=2 case:
    comp qb_val1plus = qb_i(En, sigb1plus, m1);
    comp qb_val1minus = qb_i(En, sigb1minus, m1);
    comp qb_val2plus = qb_i(En, sigb2plus, m2);
    comp qb_val2minus = qb_i(En, sigb2minus, m2);

    comp kmax_for_m1 = pmom(En, 0.0, m1); 
    comp kmax_for_m2 = pmom(En, 0.0, m2); 
    comp epsilon_for_kvec = 1.0e-5; 

    std::cout<<"kmax_for_m1 = "<<kmax_for_m1<<std::endl; 
    std::cout<<"kmax_for_m2 = "<<kmax_for_m2<<std::endl; 

    double eta_for_eps = 25; 

    double eps_for_m2k_1 = 0.0;// = energy_dependent_epsilon(eta_for_eps, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points ); 
    double eps_for_m2k_2 = 0.0; 
    double eps_for_ope_1 = eps_for_m2k_1; 
    double eps_for_ope_2 = eps_for_m2k_2; 
    double eps_for_cutoff_1 = 0.0; 
    double eps_for_cutoff_2 = 0.0; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //N range = [50, 500] and [500, 5000]
    int N_initial = 10; 
    int N_final = 1010; 
    int N_points = 500; 
    int del_N = std::abs(N_initial - N_final)/N_points; 
    //N_points = std::abs((N_initial - N_final)/del_N);

    double eta_initial = 2.0; 
    double eta_final = 50.0;
    double eta_points = 1000.0; 
    double del_eta = std::abs(eta_final - eta_initial)/eta_points;  

    double log_eps_initial = log(pow(10,-4.0));
    double log_eps_final = log(pow(10,-1.0)); 
    double log_eps_points = 500.0; 
    double del_log_eps = std::abs(log_eps_initial - log_eps_final)/log_eps_points; 


    std::vector<double> N_vec; 
    std::vector<double> eps_vec; 

    for(int i=0; i<N_points; ++i)
    {
        double N = N_initial + i*del_N; 
        for(int j=0; j<log_eps_points; ++j)
        {
            double log_eps = log_eps_initial + j*del_log_eps; 
            double eps = std::exp(log_eps); 
            N_vec.push_back(N); 
            eps_vec.push_back(eps); 
            std::cout<<"N="<<N<<'\t'
                     <<"log_eps="<<log_eps<<'\t'
                     <<"eps="<<eps<<std::endl; 
        }
    }
    //abort(); 
    std::string filename =  "delta_rhophib_density_hardcutoff.dat";
    fout.open(filename.c_str()); 

    //std::vector<std::vector<comp> > 
    int ind = 0; 
    int run_size = N_vec.size(); 

    std::vector<double> diff_vec(run_size); 
    int run_counter = 0; 
    int run_counter2 = 1; 

    /*-------------Input for Python plot script ---------*/
    /*the plot script will use this input to draw individual epsilon
    lines based on single eta values on the density plot */
    std::ofstream fout1; 
    std::string filename1 = "eta_epsilon_input_for_" + filename; 
    fout1.open(filename1.c_str()); 

    fout1<<"En"<<'\t'
         <<"qb"<<'\t'
         <<"sigb"<<'\t'
         <<"kmax"<<'\t'
         <<"mi"<<std::endl; 
    fout1<<En.real()<<'\t'
         <<qb_val1plus.real()<<'\t'
         <<sigb1plus.real()<<'\t'
         <<kmax_for_m1.real()<<'\t'
         <<m1<<std::endl; 
    fout1.close();
    /*----------------------------------------------------*/

    #pragma omp parallel for schedule(dynamic)
    for(ind=0; ind<(int)run_size; ++ind)
    {
        Eigen::MatrixXcd dpqbmat; 
        Eigen::MatrixXcd dqqmat; 
        std::vector<comp> pvec_for_m1m2; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        double N = N_vec[ind];
        double eps = eps_vec[ind]; 
        int number_of_points = N; 

        comp qb = qb_i(En, sigb1plus, m1);
        comp kmax = pmom(En, 0.0, m1); 

        eps_for_m2k_1 = eps;//energy_dependent_epsilon(eta, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points);
        eps_for_m2k_2 = eps;//energy_dependent_epsilon(eta, En, qb_val2plus, sigb2plus, kmax_for_m2, m2, number_of_points);
        eps_for_ope_1 = eps_for_m2k_1; 
        eps_for_ope_2 = eps_for_m2k_2; 
        test_dpqb_solver_ERE_2(dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k_1, eps_for_m2k_2, eps_for_ope_1, eps_for_ope_2, eps_for_cutoff_1, eps_for_cutoff_2, total_P, a0_1, r0_1, eta_1, a0_2, r0_2, eta_2, number_of_points, 'n'); 
        test_dqq_interpolator_2(dqqmat, dpqbmat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, qb, eps_for_m2k_1, eps_for_m2k_2, eps_for_ope_1, eps_for_ope_2, eps_for_cutoff_1, eps_for_cutoff_2, total_P, a0_1, r0_1, eta_1, a0_2, r0_2, eta_2, number_of_points, 'n');

        /*for(int i=0; i<pvec_for_m1m2.size(); ++i)
        {
            std::cout<<"pvec i="<<i<<'\t'<<"val = "<<pvec_for_m1m2[i]<<std::endl;
        }*/

        comp dqq11 = dqqmat(0,0); 
        
        comp gval = gfunc_i(eta_1, sigb1plus, m1, m2);//gfunc(sigb1plus,a0_1);

        comp mphib = gval*gval*dqq11; 
        comp rhophib_mphib = rhophib(qb, En)*mphib; 
        comp rhophib_val = rhophib(qb, En); 

        double diff = (double)std::abs((std::imag(1.0/mphib) + std::real(rhophib_val))/(std::real(rhophib_val)))*100;
        diff_vec[ind] = diff; 
       
        double percent_completed = ((double)run_counter+1)/run_size*100.0; 
        int divisor = run_size/10;

        #pragma omp critical 
        if(run_counter2==200)
        {
            std::cout<<"run progress = "<<percent_completed<<"%"<<std::endl;

            run_counter2 = 1; 
            
            std::cout<<"ex: N = "<<number_of_points<<'\t'
                     <<"eps = "<<eps<<'\t'
                     <<"diff = "<<diff<<'\t'
                     <<"run_no = "<<run_counter
                     <<"/"<<run_size<<std::endl; 
        }
        run_counter = run_counter + 1; 
        run_counter2 = run_counter2 + 1; 
    
    }

    for(int i=0; i<run_size; ++i)
    {
        double N_val = N_vec[i]; 
        double eps_val = eps_vec[i]; 
        double diff_val = diff_vec[i]; 

        fout<<N_val<<'\t'
            <<eps_val<<'\t'
            <<diff_val<<std::endl; 
    }
    fout.close(); 
}

void test_SA_method_1()
{
    double m1 = 1.0;
    double m2 = 0.9;//0.99999999990;//0.999; 
    double a0_m1 = 2.0; 
    double a0_m2 = 2.0; 
    double r0_m1 = 0.0; 
    double r0_m2 = 0.0; 
    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 5000;  

    

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_m1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_m1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_m2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_m2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;


    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */
    std::ofstream fout;
   

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp s = 8.2; 
    comp En = std::sqrt(s);//(En_1 + En_2)/2.0; 

    //qb vals for i=1 and i=2 case:
    comp qb_val1plus = qb_i(En, sigb1plus, m1);
    comp qb_val1minus = qb_i(En, sigb1minus, m1);
    comp qb_val2plus = qb_i(En, sigb2plus, m2);
    comp qb_val2minus = qb_i(En, sigb2minus, m2);

    comp kmax_for_m1 = pmom(En, 0.0, m1); 
    comp kmax_for_m2 = pmom(En, 0.0, m2); 
    comp epsilon_for_kvec = 1.0e-5; 

    std::cout<<"kmax_for_m1 = "<<kmax_for_m1<<std::endl; 
    std::cout<<"kmax_for_m2 = "<<kmax_for_m2<<std::endl; 

    double eta_for_eps = 25; 

    double eps_for_m2k_1 = energy_dependent_epsilon(eta_for_eps, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points ); 
    double eps_for_m2k_2 = 0.0; 
    double eps_for_ope_1 = eps_for_m2k_1; 
    double eps_for_ope_2 = eps_for_m2k_2; 
    double eps_for_cutoff_1 = 0.0; 
    double eps_for_cutoff_2 = 0.0; 
    double eps_for_m2k = eps_for_m2k_1; 
    double eps_for_ope = eps_for_m2k; 
    double eps_for_cutoff = 0.0; 
    comp qb_1 = qb_val1plus; 
    comp qb_2 = qb_1;//qb_val2plus; 
    comp sigb1 = sigb1plus;
    comp sigb2 = sigb2plus; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //N range = [50, 500] and [500, 5000]
    int N_initial = 50; 
    int N_final = 2000; 
    int N_points = 0; 
    int del_N = 50;//std::abs(N_initial - N_final)/N_points; 
    N_points = std::abs((N_initial - N_final)/del_N);

    double eta_initial = 2.0; 
    double eta_final = 50.0;
    double eta_points = 1000.0; 
    double del_eta = std::abs(eta_final - eta_initial)/eta_points;  


    std::vector<comp> pvec_for_m1m2;
    std::vector<comp> kvec_for_m1m1; 
    std::vector<comp> weights_for_pvec_for_m1m2; 
    std::vector<comp> weights_for_kvec_for_m1m1; 

    flavor_based_momentum_vector(pvec_for_m1m2, weights_for_pvec_for_m1m2, En, m1, number_of_points);
    flavor_based_momentum_vector(kvec_for_m1m1, weights_for_kvec_for_m1m1, En, m2, number_of_points);

    char debug = 'n'; 
    Eigen::MatrixXcd dqq_mat; 


    test_dS_qbiqbi_SA_Method_2plus1_system_ERE(dqq_mat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, sigb1, sigb2, qb_1, qb_2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, eta_1, a0_m2, r0_m2, eta_2, debug);

    comp g1_val = gfunc_i(eta_1, sigb1plus, m1, m2); 
    comp dqq_val = dqq_mat(0,0); 
    comp mphib = g1_val*g1_val*dqq_val; 
    comp rhophib_11 = rhophib(qb_val1plus, En); 
    double diff = std::abs(((1.0/mphib).imag() + rhophib_11)/rhophib_11)*100.0; 

    std::cout<<"Im(1/mphib) = "<<1.0/mphib<<std::endl; 
    std::cout<<"rhophib = "<<rhophib_11<<std::endl; 
    std::cout<<"diff = "<<diff<<std::endl; 


}

void test_SA_method_vs_N()
{
    double m1 = 1.0;
    double m2 = 0.9;//0.99999999990;//0.999; 
    double a0_m1 = 2.0; 
    double a0_m2 = 2.0; 
    double r0_m1 = 0.0; 
    double r0_m2 = 0.0; 
    double eta_1 = 1.0; 
    double eta_2 = 0.5; 

    //double En = 1.95;
    double total_P = 0.0; 
    double r = 0; 
    int number_of_points = 5000;  

    

    /*-----------------------------------------*/

    comp sigb1plus = sigma_b_plus(a0_m1, m1, m2);
    comp sigb1minus = sigma_b_minus(a0_m1, m1, m2); 
    comp sigb2plus = sigma_b_plus(a0_m2, m1, m1); 
    comp sigb2minus = sigma_b_minus(a0_m2, m1, m1); 
    comp phib1plus = std::sqrt(sigb1plus) + m1;
    comp phib1minus = std::sqrt(sigb1minus) + m1;
    comp phib2plus = std::sqrt(sigb2plus) + m2;
    comp phib2minus = std::sqrt(sigb2minus) + m2;


    std::cout<<"sigb1+ = "<<sigb1plus<<std::endl; 
    std::cout<<"sigb1- = "<<sigb1minus<<std::endl; 
    
    std::cout<<"sigb2+ = "<<sigb2plus<<std::endl; 
    std::cout<<"sigb2- = "<<sigb2minus<<std::endl; 

    std::cout<<"phib+ threshold 1 = "<<phib1plus<<std::endl; 
    std::cout<<"phib- threshold 1 = "<<phib1minus<<std::endl; 
    std::cout<<"phib+ threshold 2 = "<<phib2plus<<std::endl; 
    std::cout<<"phib- threshold 2 = "<<phib2minus<<std::endl; 

    comp threeparticle_threshold = (m1 + m1 + m2); 

    //comp En = (phib1 + threeparticle_threshold)/2.0;

    std::cout<<"threeparticle threshold = "<<threeparticle_threshold<<std::endl; 

    
    //abort();

    
    /* Generating a file here to check */

    comp En_initial = phib1plus; 
    comp En_final = threeparticle_threshold; 
    comp En_1 = (phib1plus + threeparticle_threshold)/2.0; 
    comp En_2 = En_initial; 
    comp s = 8.2; 
    comp En = std::sqrt(s);//(En_1 + En_2)/2.0; 

    //qb vals for i=1 and i=2 case:
    comp qb_val1plus = qb_i(En, sigb1plus, m1);
    comp qb_val1minus = qb_i(En, sigb1minus, m1);
    comp qb_val2plus = qb_i(En, sigb2plus, m2);
    comp qb_val2minus = qb_i(En, sigb2minus, m2);

    comp kmax_for_m1 = pmom(En, 0.0, m1); 
    comp kmax_for_m2 = pmom(En, 0.0, m2); 
    comp epsilon_for_kvec = 1.0e-5; 

    std::cout<<"kmax_for_m1 = "<<kmax_for_m1<<std::endl; 
    std::cout<<"kmax_for_m2 = "<<kmax_for_m2<<std::endl; 

    
    comp qb_1 = qb_val1plus; 
    comp qb_2 = qb_1;//qb_val2plus; 
    comp sigb1 = sigb1plus;
    comp sigb2 = sigb2plus; 

    double En_points = 3000.0; 
    comp delEn = std::abs(En_initial - En_final)/En_points; 

    //N range = [50, 500] and [500, 5000]
    int N_initial = 50; 
    int N_final = 500; 
    int N_points = 0; 
    int del_N = 1;//std::abs(N_initial - N_final)/N_points; 
    N_points = std::abs((N_initial - N_final)/del_N);

    double eta_initial = 2.0; 
    double eta_final = 50.0;
    double eta_points = 1000.0; 
    double del_eta = std::abs(eta_final - eta_initial)/eta_points;  

    std::ofstream fout; 
    std::string filename = "Mphib_SA_method_vs_N_eta_25_1_3_forpaper.dat";
    fout.open(filename.c_str()); 

    for(int i=0; i<N_points; ++i)
    {
        number_of_points = N_initial + i*del_N; 

        double eta_for_eps = 25; 

        double eps_for_m2k_1 = energy_dependent_epsilon(eta_for_eps, En, qb_val1plus, sigb1plus, kmax_for_m1, m1, number_of_points ); 
        double eps_for_m2k_2 = eps_for_m2k_1; 
        double eps_for_ope_1 = eps_for_m2k_1; 
        double eps_for_ope_2 = eps_for_m2k_2; 
        double eps_for_cutoff_1 = 0.0; 
        double eps_for_cutoff_2 = 0.0; 
        double eps_for_m2k = eps_for_m2k_1; 
        double eps_for_ope = eps_for_m2k; 
        double eps_for_cutoff = 0.0; 

        std::vector<comp> pvec_for_m1m2;
        std::vector<comp> kvec_for_m1m1; 
        std::vector<comp> weights_for_pvec_for_m1m2; 
        std::vector<comp> weights_for_kvec_for_m1m1; 

        flavor_based_momentum_vector(pvec_for_m1m2, weights_for_pvec_for_m1m2, En, m1, number_of_points);
        flavor_based_momentum_vector(kvec_for_m1m1, weights_for_kvec_for_m1m1, En, m2, number_of_points);

        //simple_momentum_vector(pvec_for_m1m2, weights_for_pvec_for_m1m2, 0.0, 0.25*kmax_for_m1, number_of_points); 
        //simple_momentum_vector(kvec_for_m1m1, weights_for_kvec_for_m1m1, 0.0, kmax_for_m2, number_of_points); 
        
        char debug = 'n'; 
        Eigen::MatrixXcd dqq_mat; 

    
        test_dS_qbiqbi_SA_Method_2plus1_system_ERE(dqq_mat, En, m1, m2, pvec_for_m1m2, weights_for_pvec_for_m1m2, kvec_for_m1m1, weights_for_kvec_for_m1m1, sigb1, sigb2, qb_1, qb_2, eps_for_m2k, eps_for_ope, eps_for_cutoff, total_P, a0_m1, r0_m1, eta_1, a0_m2, r0_m2, eta_2, debug);

        comp g1_val = gfunc_i(eta_1, sigb1plus, m1, m2); 
        comp dqq_val = dqq_mat(0,0); 
        comp mphib = g1_val*g1_val*dqq_val; 
        comp rhophib_11 = rhophib(qb_val1plus, En); 
        comp rhoM = rhophib_11*mphib; 
        double diff = std::abs(((1.0/mphib).imag() + rhophib_11)/rhophib_11)*100.0; 

        std::cout<<"N = "<<number_of_points<<std::endl; 
        std::cout<<"Im(1/mphib) = "<<1.0/mphib<<std::endl; 
        std::cout<<"rhophib = "<<rhophib_11<<std::endl; 
        std::cout<<"diff = "<<diff<<std::endl; 
        std::cout<<"========================"<<std::endl; 

        fout<<std::setprecision(20)<<number_of_points<<'\t'
            <<diff<<'\t'
            <<rhoM.real()<<'\t'
            <<rhoM.imag()<<std::endl; 
    }
    fout.close(); 
    


}


int main()
{
    //test_M2(); 
    //test_q2k_sigk();
    //test_sigk_q2k();
    //test_kinematic_variables();
    //test_dpqb_building();
    //test_dpqb_building_1();
    //test_d11_building_1();
    //test_dpqb_vs_N_building();

    //Mphib_degenerate_testing();
    //Mphib_degenerate_testing_vs_N();
    //Mphib_degenerate_testing_vs_N_En_onefourth();
    //test_dpqb_vs_N_building_1();
    //test_delta_rhophib_density();
    //test_M2_sigb_range();
    //test_M2_1(); 
    //test_dpqb_for_m1m1_vs_N_building_1();
    
    //test_Gs();
    //test_Gs_surface();
    //test_kernel_singularities();
    //compare_pcuts();
    
    //compare_Bmats(); 
    //plot_single_integral_equation_components();
    //test_dpqb_vs_N_building_2();
    //test_dpqb_vs_N_50_eta_dependence();
    //test_delta_rhophib_density_with_omp();
    
    //SA Method Testing:
    //test_SA_method_1(); 
    test_SA_method_vs_N();

    //Going back to OPE again
    //test_Gs_surface_1();
    return 0;
}