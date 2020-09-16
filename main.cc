#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>

#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"
using namespace itensor;

#include "ed_basis.h"
#include "utils.h"

int main(int argc, char* argv[])
{
    std::complex<double> ii(0,1);
    double rescale_coeff = 2 * PI / 1000;

    int time0 = time(0);

    // Ideal qubits, J1-J2 chain or triangular ladder 
    //parameters from file: Hamiltonian
    Parameters Params;
    auto input = InputGroup(argv[1],"input");
    Params.N = input.getInt("N"); 
    Params.nop = input.getInt("nop"); // number of particles/up-spins
    Params.gamma = input.getReal("gamma",0);
    Params.N_ini = input.getInt("N_ini",2);
    Params.evo_type = input.getInt("evo_type",0); // 0/1 for real/imaginary evolution
    Params.dt = input.getReal("dt",1.0);
    Params.nt = input.getInt("nt",100);
    Params.N_gse = input.getInt("N_gse",3);

    Params.evo_IG = 'y';
    Params.evo_IGsq = 'y';
    Params.evo_ExtremeP = 'y';

    //parameters from file: static dmrg 
    auto dmrg_nsweeps = input.getInt("dmrg_nsweeps");
    auto table = InputGroup(input,"dmrg_sweeps");
    auto dmrg_sweeps = Sweeps(dmrg_nsweeps,table);
    println(dmrg_sweeps);

    // read J1/J2 coupling coefficients
    double *J1_coeff = new double [Params.N-1];
    std::ifstream read_J1("J1_coeff.in", std::ios::in);
    std::cout << "J1/NN couplings (GHz):" << std::endl;
    for (int i = 0; i < Params.N-1; i++)
    {
        read_J1 >> J1_coeff[i];
        J1_coeff[i] *= rescale_coeff;
        std::cout << J1_coeff[i] << std::endl;
    }
    read_J1.close();
    //
    double *J2_coeff = new double [Params.N-2];
    std::ifstream read_J2("J2_coeff.in", std::ios::in);
    std::cout << "J2/NNN couplings (GHz):" << std::endl;
    for (int i = 0; i < Params.N-2; i++)
    { 
        read_J2 >> J2_coeff[i];
        J2_coeff[i] *= rescale_coeff; 
        std::cout << J2_coeff[i] << std::endl;
    }
    read_J2.close();
    //
    Params.gamma *= rescale_coeff;
    double x0 = 0.5 * (Params.N - 1);
    double *hii = new double[Params.N];
    std::cout << "onsite potential in code (GHz):" << std::endl;
    for (int i = 0; i < Params.N; i++)
    {
        hii[i] = -(i - x0)*Params.gamma;
        std::cout << hii[i] << std::endl;
    }
    double energy_shift = 0;
    for(int i = 1; i <= Params.N; ++ i) energy_shift += Params.gamma*(i-Params.N/2.0-0.5)*0.5;
    std::cout << "energy shift: " << energy_shift << std::endl;

    //
    auto t_coeff = (Params.evo_type) ? ii : 1.0;
    //std::cout << "\n tau = " << tau << std::endl;

    // Make the Hamiltonian for qubits
    // Make N spin 1/2's
    auto sites = SpinHalf(Params.N);
    // note: \sigma = 2*S
    auto ampo = AutoMPO(sites);    
    int sign = 1;
    construct_qubits_mpo(ampo, Params.N, J1_coeff, J2_coeff, hii, sign);
    auto H = toMPO(ampo);

    // dmrg to get E_min/E_max
    // Set the initial state to be Neel state
    auto dmrg_state0 = InitState(sites);
    rand_ini_state(Params.N, Params.nop, dmrg_state0);
    auto dmrg_psi0 = MPS(dmrg_state0);
    // get E_min
    auto [E_min,wf_min] = dmrg(H,dmrg_psi0,dmrg_sweeps,"Quiet");
    //printfln("\n DMRG: E_min = %.10f",E_min);

    // get E_max
    auto ampo1 = AutoMPO(sites);
    sign = -1; 
    construct_qubits_mpo(ampo1, Params.N, J1_coeff, J2_coeff, hii, sign);
    auto H1 = toMPO(ampo1); 
    auto [E_max,wf_max] = dmrg(H1,dmrg_psi0,dmrg_sweeps,"Quiet");
    E_max = - E_max;
    E_min += energy_shift;
    E_max += energy_shift;
    
    printfln("\n DMRG: E_min = %.10f",E_min);
    printfln("\n DMRG: E_max = %.10f",E_max);
    std::ofstream ofe("E_min_max.dat");
    ofe << E_min << std::endl << E_max;
    ofe.close();

    // Time evolution: initial states information
    double E_band = E_max - E_min;
    double epsilon = 0.5;
    double E_Target = epsilon * E_band + E_min;

    Basis basis(Params.N, Params.nop);
    auto Dim = basis.get_Dim();
    double *H_diag = new double [Dim];
    construct_Hdiag(basis, Dim, hii, H_diag);
    /* 
       FILE* f_out;
       f_out = fopen("FockSpec.bin", "wb");
       fwrite(H_diag, sizeof(double), Dim, f_out);
       fclose(f_out);
       */
    vector <pair<double, l_int> > Fock_E_n;
    for (l_int p = 0; p < Dim; p++)
    {
        pair<double, l_int> aux = make_pair(H_diag[p], p);
        Fock_E_n.push_back(aux);
    }
    std::sort(Fock_E_n.begin(), Fock_E_n.end()); 
    l_int p_Fock_target = Get_TargetFock(E_Target, Dim, Fock_E_n);
    vector <pair<int, l_int> > target_ind_Fock;

    int N_ini = Params.N_ini; 
    if (-1 != p_Fock_target)
    {
        for (int i = -N_ini / 2; i < N_ini / 2; i++)
        {
            l_int p_Fock = p_Fock_target + i;
            if (p_Fock >= 0 && p_Fock < Dim)
            {
                l_int p_ini = Fock_E_n[p_Fock].second;
                pair<int, l_int> aux = make_pair(i, p_ini);
                target_ind_Fock.push_back(aux);
            }
        }
    }

    // write time vec
    double *time_vec = new double[Params.nt];
    //for (int i = 0; i < Params.nt; i++) time_vec[i] = Params.dt * i;   
    
    // try different dt at different times 
    // 0.1 to 1, 1 to 10, then dt = dt
    for (int i = 0; i < 10; i++)
        time_vec[i] = 0.1*i;
    for (int i = 10; i < 20; i++)
        time_vec[i] = 1.0*(i-9);
    for (int i = 20; i < Params.nt; i++)
        time_vec[i] = time_vec[i-1] + Params.dt; 
    //
    Vec_fwrite_double("evo_time_vec.bin", time_vec, Params.nt);

    // compute time evolution
#pragma omp parallel for
    for (size_t i_p = 0; i_p < target_ind_Fock.size(); i_p++)
    {
        std::cout<<"There are "<<omp_get_num_threads()<<" threads"<<std::endl;
        
        int prt_ind = target_ind_Fock[i_p].first;
        l_int p_ind = target_ind_Fock[i_p].second;
        l_int p_state = basis.get_state(p_ind);
        double E_real = H_diag[p_ind];

        char fini[80];
        sprintf(fini, "IniStateInfo_state_epsi_e%0.2f_ind%d.dat", epsilon, prt_ind);
        std::ofstream ofini(fini);
        ofini << p_state << endl;   // this is a binary number
        ofini << std::setprecision(14) << (E_real - E_min)/E_band << endl;
        ofini.close();

        // evo data
        Dyn_DataStruc dyn_data;
        dyn_data.Initialize(Params, prt_ind, epsilon);
        //
        /////////////////////////////////////////////////////////////
        // tdvp 
        auto state = InitState(sites);
        int ini_state = p_state;
        std::cout << "initial state:" <<std::endl;
        set_target_product_state(state, Params.N, ini_state);
        auto psi1 = MPS(state);

        // start TDVP, either one site or two site algorithm can be used by adjusting the "NumCenter" argument
        println("----------------------------------------GSE-TDVP---------------------------------------");

        auto energy = real(innerC(psi1,H,psi1));
        printfln("Initial energy = %.5f", energy);

        auto sweeps = Sweeps(1);
        sweeps.maxdim() = 2000;
        //sweeps.cutoff() = 1E-12;
        sweeps.niter() = 12;

        // get observables at t=0
        compute_evo_obs(Params, sites, psi1, ini_state, 0, dyn_data);
        // compute time evolution 
        for(int it = 1; it < Params.nt; ++it)
        {
            auto tau = t_coeff*(time_vec[it]-time_vec[it-1]);
            std::cout << "time-step: " << it << ", dt: " << time_vec[it]-time_vec[it-1] << ", time: " << time_vec[it] << std::endl;       
            //if(n < 10)
            if(it < Params.N_gse)
            {
                // Global subspace expansion
                std::vector<Real> epsilonK = {1E-12, 1E-12};
                addBasis(psi1,H,epsilonK,{"Cutoff",1E-12,
                        "Method","DensityMatrix",
                        "KrylovOrd",3,
                        "DoNormalize",true,
                        "Quiet",true});
            }

            // TDVP sweep
            energy = tdvp(psi1,H,-tau,sweeps,{"DoNormalize",true,
                    "Quiet",true,
                    "NumCenter",1});

            // observables 
            compute_evo_obs(Params, sites, psi1, ini_state, it, dyn_data);
        }
        printfln("\nEnergy after time evolution = %.10f",energy);
        printfln("Using overlap = %.10f", real(innerC(psi1,H,psi1)) );

        dyn_data.PrintDynResults(Params); 
        dyn_data.ReleaseSpace(Params);

    }

    std::cout << "Total time cost: " << time(0) - time0 << " s" << std::endl;
    return 0;
}
