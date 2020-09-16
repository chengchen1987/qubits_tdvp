#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
##SBATCH -n 1
##SBATCH -c 64
#SBATCH -o slurm-%j.o
#SBATCH -e slurm-%j.err
#SBATCH -J mbl_tdvp
##SBATCH -t 3-05:00:00

export OMP_NUM_THREADS=1
EXECFILE=mbl_tdvp

ORIG_FOLDER=`pwd`
FILES="mbl_tdvp J1_coeff.in J2_coeff.in"

# fixed parameters
L=16
nop=8

N_ini=2
evo_type=1
dt=5.0
nt=311

N_gse=15

# loop for realizations, disorder
r_start=0
r_end=1

Rand_V=0.0
#gnum=18
#g_list=(0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0)
#gnum=9
#g_list=(0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5)
#g_list=(5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0)
gnum=1
#g_list=(0.5)
g_list=(1.0)
#g_list=(2.0)
#g_list=(3.0)
#g_list=(4.0)
#g_list=(5.0)
#g_list=(6.0)
#g_list=(7.0)
#g_list=(8.0)
#g_list=(9.0)
#g_list=(4.5)


for ((r=${r_start}; r<${r_end}; r++))
do
    #for ((iV=0; iV<${Vnum}; iV++))
    for ((ig=0; ig<${gnum}; ig++))
    do
        #Rand_V=${V_list[$iV]}
        gamma=${g_list[$ig]}
        cd $ORIG_FOLDER
        dirname=Qubits_L${L}_n${nop}_g${gamma}_V${Rand_V}_r${r}
        mkdir $dirname
        for i_file in $FILES
        do
            cp $i_file $dirname
        done
        echo "$dirname"
        cd $dirname
# input parameters
cat > input.in << EOF
input 
    {
    N = ${L}
    nop = ${nop}

    gamma=${gamma}

    N_ini = ${N_ini}
    evo_type = ${evo_type}
    dt = ${dt}
    nt = ${nt}

    N_gse = ${N_gse}

    dmrg_nsweeps = 10
    dmrg_sweeps
        {
        maxdim  mindim  cutoff  niter  noise
        50      20      1E-6    4      1E-4
        80      20      1E-8    4      1E-6
        100     10      1E-10   4      1E-8
        150     10      1E-12   4      1E-10
        200     10      1E-12   4      0
        200     10      1E-12   4      0
        500     10      1E-12   4      0
        500     10      1E-12   4      0
        1000     10      1E-12   4      0
        1000     10      1E-12   4      0
        2000     10      1E-12   4      0
        2000     10      1E-12   4      0
        }
}

EOF
            srun ./$EXECFILE input.in > output
            #./$EXECFILE input.in > output
            rm $EXECFILE
        sleep 1
    done
done
