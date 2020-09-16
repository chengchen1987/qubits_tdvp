#include "utils.h"

void set_target_product_state(itensor::InitState& state, const int &N, const int &ini_state)
{
    for(int i = 1; i <= N; ++i)
    {
        if ((ini_state>>(i-1))&1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
        //
        std::cout << ((ini_state>>(i-1))&1);
    }
}

// randomly select k index in (0,n), for initial MPS with res[i] sites occupied
auto rand_perm_k_in_n(const int &n, const int &k)
{
    std::vector<int> res(k);
    int j;
    for (j = 0; j < k; ++j) res[j]=j;
    for (; j < n; j++)
    {
        int m = rand() % (j + 1);
        if (m < k) res[m] = j;
    }
    return res;
}

void rand_ini_state(const int &n, const int &k, itensor::InitState& state)
{
    srand(time(0));
    for(int j = 1; j < n; ++j) state.set(j,"Dn");  // set all sites to empty 
    auto res = rand_perm_k_in_n(n, k);
    for(int i = 0; i < k; ++i) state.set(res[i]+1,"Up");
    // print intiail Fock state
    std::vector<int> fock0(n);
    for(int i = 0; i < n; ++i) fock0[i] = 0;
    for(int i = 0; i < k; ++i) fock0[res[i]] = 1;
    std::cout << "Initial Fock MPS:" << std::endl;
    for(int i = 0; i < n; ++i) std::cout << fock0[i];
    std::cout << std::endl;
}

void Dyn_DataStruc::Initialize(Parameters& para, const int& _prt_ind, const double& _epsilon)
{
    prt_ind = _prt_ind;
    varepsilon = _epsilon;
    //
    int t_len = para.nt;
    if ('y' == para.evo_IG)
        IG_p_t = new double[t_len];
    if ('y' == para.evo_IGsq)
        IGsq_p_t = new double[t_len];
    if ('y' == para.evo_ExtremeP)
    {
        ExtremeP_p_t = new double[t_len];
        mzi_p_t = new double[t_len * para.N];
    }
}

void Dyn_DataStruc::ReleaseSpace(Parameters& para)
{
    if ('y' == para.evo_IG)
        delete[]IG_p_t;
    if ('y' == para.evo_IGsq)
        delete[]IGsq_p_t;
    if ('y' == para.evo_ExtremeP) {
        delete[]ExtremeP_p_t;
        delete[]mzi_p_t;
    }
}

void::Dyn_DataStruc::PrintDynResults(Parameters& para)
{
    int t_len = para.nt;
    if ('y' == para.evo_IG)
    {
        char fname[80];
        sprintf(fname, "evo_IG_e%0.2f_ind%d.bin", varepsilon, prt_ind);
        //        format("wf_idelta%d_iphi%d.mps", idelta,iphi)
        Vec_fwrite_double(fname, IG_p_t, t_len);
    }

    if ('y' == para.evo_IGsq)
    {
        char fname[80];
        sprintf(fname, "evo_IGsq_e%0.2f_ind%d.bin", varepsilon, prt_ind);
        Vec_fwrite_double(fname, IGsq_p_t, t_len);
    }

    if ('y' == para.evo_ExtremeP)
    {
        char fname[80];
        sprintf(fname, "evo_ExtremeP_e%0.2f_ind%d.bin", varepsilon, prt_ind);
        Vec_fwrite_double(fname, ExtremeP_p_t, t_len);
        //
        char fname1[80];
        sprintf(fname1, "evo_mzi_e%0.2f_ind%d.bin", varepsilon, prt_ind);
        Vec_fwrite_double(fname1, mzi_p_t, para.N * t_len);
    }
}

void Vec_fwrite_double(const char* fname, double* data, const int& dsize)
{
    FILE* f_out;
    f_out = fopen(fname, "wb");
    fwrite(data, sizeof(double), dsize, f_out);
    fclose(f_out);
}

// compute <wf|op1_{i1}]|wf>
double one_site_ob(auto &_sites, auto &_wf, int i1, std::string const &_op1)
{
    _wf.position(i1);
    auto ket_op = _wf(i1)*op(_sites,_op1,i1);
    auto bra = dag(prime(_wf(i1),"Site"));
    return real(eltC(bra*ket_op));  // since I know it is real
}

// compute <wf|op1_{i1}*op2_{i2]|wf>
double two_site_corr(auto &_sites, auto &_wf, int i1, int i2, std::string const &_op1, std::string const &_op2)
{
    if (1 == (i2-i1))
    {
        _wf.position(i1);
        auto ket_op = _wf(i1)*_wf(i2)*op(_sites,_op1,i1)*op(_sites,_op2,i2);
        auto bra = dag(_wf(i1)*_wf(i2));
        bra.prime("Site");
        return real(eltC(bra*ket_op));
    }
    if (1 < (i2-i1))
    {
        auto op1 = op(_sites,_op1,i1);
        auto op2 = op(_sites,_op2,i2);
        // set initial tensor, make an auxilary copy C
        _wf.position(i1);
        auto C = _wf(i1); 
        C *= op1;
        auto ir = commonIndex(_wf(i1),_wf(i1+1),"Link");
        C *= dag(prime(prime(_wf(i1),"Site"),ir)); 
        for (int bi = i1+1; bi < i2; ++bi)
        {
            C *= _wf(bi);
            C *= dag(prime(_wf(bi),"Link"));
        }
        C *= _wf(i2);
        C *= op2;
        auto il = commonIndex(_wf(i2),_wf(i2-1),"Link");
        C *= dag(prime(prime(_wf(i2),"Site"),il));
        return real(eltC(C));
    }
    else
    {
        return 0.0;
    }
}

void compute_2site_corr(Parameters &para, auto &_sites, auto &_wf, int b1, int b2, std::string const &_op1, std::string const &_op2, double *Cij)
{
    // itensor/mps index starts from 1 
    b1 += 1;
    //    b2 += 1;
    int N = para.N;
    int i1 = b1;
    if (1 > (b2-b1))
    {
        std::cout << "b2 must larger than b1!";
        exit(-10);
    }
    else if (1 == (b2-b1)) 
    {
        Cij[(b1-1)*N + (b2-1)] = two_site_corr(_sites, _wf, b1, b2, "Sz", "Sz");
        Cij[(b1-1)*N + (b2-1)] = Cij[(b2-1)*N + (b1-1)];
    }   
    else 
    {
        _wf.position(i1);
        auto C = _wf(i1); 
        auto op1 = op(_sites,_op1,i1);
        C *= op1;
        auto link_1 = commonIndex(_wf(i1),_wf(i1+1),"Link");
        C *= dag(prime(prime(_wf(i1),"Site"),link_1)); 
        //
        std::cout << "corr2: before i2" << std::endl;
        for (int i2 = i1+1; i2 <= b2; i2++){
            std::cout << "i2 = " << i2 << std::endl;
            // set initial tensor, make an auxilary copy C
            for (int bi = i1+1; bi < i2; ++bi)
            {
                C *= _wf(bi);
                C *= dag(prime(_wf(bi),"Link"));
                auto op2 = op(_sites,_op2,i2);
                auto Caux = C * _wf(i2) * op2;
                std::cout << "Caux = C * _wf(i2) * op2" << std::endl;
                //C *= _wf(i2);
                //C *= op2;
                auto link_2 = commonIndex(_wf(i2),_wf(i2-1),"Link");
                Caux *= dag(prime(prime(_wf(i2),"Site"),link_2));
                Cij[(i1-1)*N + (i2-1)] = real(eltC(Caux));
                std::cout << "Cij: " << Cij[(i1-1)*N + (i2-1)]  << std::endl;
            }
        } 
    }
}

double Get_ExtremeP_from_mzVec(double* _mzi, const int &N)
{
    double* mz_abs = new double[N];
    for (int i = 0; i < N; i++)
        mz_abs[i] = std::abs(_mzi[i]);
    double extremeP = 0.5 - *std::max_element(mz_abs, mz_abs + N);
    delete[] mz_abs;
    return extremeP;
}

double Get_IG_from_mzVec(double* _mzi, const int &N, const int &nop, const int &ini_state)
{
    double aux = 0;
    double cf1 = 1.0/nop;
    double cf0 = -1.0/(N-nop);
    for(int i = 0; i < N; ++i)
    {
        int c01 = ((ini_state >> i) & 1);
        aux += c01 ? cf1*(_mzi[i]+0.5) : cf0*(_mzi[i]+0.5);
    }
    return aux;
}

double Get_IGsq_from_Szij(double* Szij, double *mzi, const int &N, const int &nop, const int &ini_state)
{
    double cf1 = 1.0/nop;
    double cf0 = -1.0/(N-nop);
    // get coefficient matrix 
    double* ci = new double[N];
    for (int i = 0; i < N; i++) 
    {
        int c01 = ((ini_state >> i) & 1);
        ci[i] = c01 ? cf0 : cf1;
    }
    //
    double aux = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) 
        {
            aux += ci[i] * ci[j] * Szij[i*N + j];
            aux += ci[i] * ci[j] * 0.5 * (mzi[i] + mzi[j]);
            aux += ci[i] * ci[j] * 0.25;
        }
    delete []ci;
    return aux;
}

void compute_evo_obs(Parameters &para, itensor::BasicSiteSet<itensor::SpinHalfSite>& _sites, itensor::MPS& _wf, const int &ini_state, const int &it, Dyn_DataStruc &dyn_data)
{
    // Sz_i(t)
    for(int i = 0; i < para.N; ++i)
    {
        dyn_data.mzi_p_t[it*para.N + i] = real(one_site_ob(_sites, _wf, i+1, "Sz"));
    }
    // IG(t)
    dyn_data.IG_p_t[it] = Get_IG_from_mzVec(&dyn_data.mzi_p_t[it*para.N], para.N, para.nop, ini_state);
    // ExtremeP(t)
    dyn_data.ExtremeP_p_t[it] =  Get_ExtremeP_from_mzVec(&dyn_data.mzi_p_t[it*para.N], para.N);  
    // compute QFI = 4(<I^2> - <I>^2)
    int L = para.N;
    // get <Sz_i*Sz_j>
    if ('y' == para.evo_IGsq)
    {
        double *Szij = new double[L*L];
        for (int i = 0; i < L; i++)
        {
            Szij[i*L + i] = 0.25;
            /*            for (int j = i+1; j < L; j++)
                          {
                          Szij[i*L + j] = two_site_corr(_sites, _wf, i+1, j+1, "Sz", "Sz");
                          Szij[j*L + i] = Szij[i*L + j];
                          }*/
            std::cout << "before compute_2site_corr() " << std::endl;
            compute_2site_corr(para, _sites, _wf, i, L, "Sz", "Sz", Szij);  
            for (int j = i+1; j < L; j++)
            {
                Szij[j*L + i] = Szij[i*L + j];
            }
        }
        dyn_data.IGsq_p_t[it] = Get_IGsq_from_Szij(Szij, &dyn_data.mzi_p_t[it*para.N], para.N, para.nop, ini_state);
        delete []Szij;
    }
}

void construct_qubits_mpo(itensor::AutoMPO &ampo, const int &N, double *J1_coeff, double *J2_coeff, double *hii, const double &sign)
{
    for(int i = 1; i <= N-1; ++ i)
    {
        ampo += sign*J1_coeff[i-1],"S+",i,"S-",i+1;
        ampo += sign*J1_coeff[i-1],"S-",i,"S+",i+1;
    }
    for(int i = 1; i <= N-2; ++ i)
    {
        ampo += sign*J2_coeff[i-1],"S+",i,"S-",i+2;
        ampo += sign*J2_coeff[i-1],"S-",i,"S+",i+2;
    }
    for(int i = 1; i <= N; ++ i)
    {
        ampo += sign*hii[i-1],"Sz",i;
    }
}

l_int Get_TargetFock(const double& Target_E, const l_int&Dim, vector <pair<double, l_int> >& Fock_E_n) 
{
    for (l_int p = 0; p < Dim; p++) {
        if (Fock_E_n[p].first > Target_E) {
            return p;
        }
    }
    return -1;
}

void construct_Hdiag(Basis &basis, const l_int&Dim, double*hii, double *H_diag)
{
    for (l_int k = 0; k < Dim; k++)
    {
        double H_aux = 0;
        l_int Num_k  = basis.get_state(k);
        for (int i = 0; i < basis.get_L(); i++) H_aux += hii[i] * (((Num_k >> i) & 1));
        H_diag[k] = H_aux;
    }
}
