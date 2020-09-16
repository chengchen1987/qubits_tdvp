#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <fstream>
#include <algorithm>
#include <random>

#include "itensor/all.h"

#include "ed_basis.h"

using namespace std;

class Parameters
{
public:
    int N;
    int nop;
    //
    double gamma;
    //
    int N_ini;
    int evo_type;
    double dt;
    int nt;
    char evo_IG;
    char evo_IGsq;
    char evo_ExtremeP;
    // 
    int N_gse; // global subspace expansion steps
};

class Dyn_DataStruc
{
public:
    void Initialize(Parameters& para, const int& _prt_ind, const double& _epsilon);
    void ReleaseSpace(Parameters& para);
    void PrintDynResults(Parameters& para);

    double varepsilon;
    int prt_ind;

    //int t_len;
    //  //double* time_vec;
    //
    double* IG_p_t;
    double* IGsq_p_t;
    //double* EE_p_t; 
    double* ExtremeP_p_t;
    double* mzi_p_t;
};

//
void Vec_fwrite_double(const char* fname, double* data, const int& dsize);

//
void set_target_product_state(itensor::InitState& state, const int &N, const int &ini_state);
auto rand_perm_k_in_n(const int &n, const int &k);
void rand_ini_state(const int &n, const int &k, itensor::InitState& state);
void construct_qubits_mpo(itensor::AutoMPO &ampo, const int &N, double *J1_coeff, double *J2_coeff, double *hii, const double &sign);
void construct_Hdiag(Basis &basis, const l_int&Dim, double*hii, double *H_diag);

//
void compute_evo_obs(Parameters &para, itensor::BasicSiteSet<itensor::SpinHalfSite>& _sites, itensor::MPS& _wf, const int &ini_state, const int &it, Dyn_DataStruc &dyn_data);

double one_site_ob(auto &_sites, auto &_wf, int i1, std::string const &_op1);
double two_site_corr(auto &_sites, auto &_wf, int i1, int i2, std::string const &_op1, std::string const &_op2);
void compute_2site_corr(Parameters &para, auto &_sites, auto &_wf, int b1, int b2, std::string const &_op1, std::string const &_op2, double *Cij);

double Get_IG_from_mzVec(double* _mzi, const int &N, const int &nop, const int &ini_state);
double Get_ExtremeP_from_mzVec(double* _mzi, const int &N);
double Get_IGsq_from_Szij(double* Szij, double *mzi, const int &N, const int &nop, const int &ini_state);

//
l_int Get_TargetFock(const double& Target_E, const l_int&Dim, vector <pair<double, l_int> >& Fock_E_n);


#endif
