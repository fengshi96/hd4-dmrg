#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <stdio.h>
#include <assert.h>
#include <cassert>
#include <string>
#include <iomanip>      // std::setprecision
#include "itensor/all.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
using namespace itensor;
void Chain1(int N, int Norb, Eigen::VectorXi Nc_, Eigen::MatrixXi& N1neigh_, Eigen::VectorXi indx_);

int main(int argc, char* argv[]) {
    if(argc < 2)  {
        printfln("Usage: %s input_file",argv[0]);
        return 0;
    }
    auto input = InputGroup(argv[1],"input");
    auto orbitals = input.getInt("Orbitals");
    auto LLX = input.getInt("LLX");
    auto totalSweeps = input.getInt("totalSweeps");
    auto N = LLX;
    auto Norb = LLX*orbitals;

    double JExch = input.getReal("JExch");
    double Lambda = input.getReal("Lambda");
    double BSzPinning = input.getReal("BSzPinning"); 

    // -- setting up sweeping paramaters --
    auto sw_table = InputGroup(input,"table_name");
    auto sweeps = Sweeps(totalSweeps,sw_table);
    println(sweeps);
    Args args;
    args.add("Verbose",true);

    Eigen::MatrixXi N1neigh_;
    Eigen::VectorXi indx_,Nc_;
    Chain1(LLX, orbitals, Nc_, N1neigh_, indx_);

    Eigen::MatrixXd Connx(Norb,Norb), Conny(Norb,Norb), Connz(Norb,Norb);
    Eigen::MatrixXd Sorx(Norb,Norb), Sory(Norb,Norb), Sorz(Norb,Norb);
    Connx.setZero();

    // --- Kitaev Model --- // Neighbors for each site
    auto sites = SpinOne(Norb); //make a chain of N spin 1/2's
    auto ampo = AutoMPO(sites);
    for(int i=0;i<N;i++){ 	// ith site
        int ia = i*orbitals + 0;
        int ib = i*orbitals + 1;

        int j = N1neigh_(i,0);
        int ja = j*orbitals + 0;
        int jb = j*orbitals + 1;

        assert(ia<Norb); assert(ib<Norb);
        assert(ja<Norb); assert(jb<Norb);
        assert(j != -1); // --- assertions not working for some reason!

	// === Pinning ====
	if(i==0) ampo += BSzPinning,"Sz",ia+1; 

	// === Spin - Heisenberg term === 
        double exc = -2.0*JExch/4.0;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1;
        if(i<j) Connx(ia,ja) = exc;

        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1;
        if(i<j) Conny(ia,ja) = exc;

        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1;
        if(i<j) Connz(ia,ja) = exc;
	// if(i<j) std::cout << ia+1 << "-" << ja+1 << std::endl;

        // === spin orbit term ===
        exc = Lambda/2.0;
        ampo += exc,"Sx",ia+1,"Sx",ib+1;
        Sorx(ia,ib) = exc;

        ampo += exc,"Sy",ia+1,"Sy",ib+1;
        Sory(ia,ib) = exc;

        ampo += exc,"Sz",ia+1,"Sz",ib+1;
        Sorz(ia,ib) = exc;
        // std::cout << ia+1 << "-" << ib+1 << std::endl;

        // === Si.Sj Li.Lj ===
        exc = JExch/4.0;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1,"Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1,"Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1,"Sz",ib+1,"Sz",jb+1;

        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1,"Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1,"Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1,"Sz",ib+1,"Sz",jb+1;

        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1,"Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1,"Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1,"Sz",ib+1,"Sz",jb+1;
	// if(i<j) std::cout << ia+1 << "-" << ja+1 << "-" << ib+1 << "-" << jb+1 << std::endl;

        // === Si.Sj Li.Lj Li.Lj ===
        exc = JExch/4.0;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sx",ib+1,"Sx",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sx",ib+1,"Sx",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sx",ib+1,"Sx",jb+1, "Sz",ib+1,"Sz",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sy",ib+1,"Sy",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sy",ib+1,"Sy",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sy",ib+1,"Sy",jb+1, "Sz",ib+1,"Sz",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sz",ib+1,"Sz",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sz",ib+1,"Sz",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sx",ia+1,"Sx",ja+1, "Sz",ib+1,"Sz",jb+1, "Sz",ib+1,"Sz",jb+1;

        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sx",ib+1,"Sx",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sx",ib+1,"Sx",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sx",ib+1,"Sx",jb+1, "Sz",ib+1,"Sz",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sy",ib+1,"Sy",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sy",ib+1,"Sy",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sy",ib+1,"Sy",jb+1, "Sz",ib+1,"Sz",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sz",ib+1,"Sz",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sz",ib+1,"Sz",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sy",ia+1,"Sy",ja+1, "Sz",ib+1,"Sz",jb+1, "Sz",ib+1,"Sz",jb+1;

        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sx",ib+1,"Sx",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sx",ib+1,"Sx",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sx",ib+1,"Sx",jb+1, "Sz",ib+1,"Sz",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sy",ib+1,"Sy",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sy",ib+1,"Sy",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sy",ib+1,"Sy",jb+1, "Sz",ib+1,"Sz",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sz",ib+1,"Sz",jb+1, "Sx",ib+1,"Sx",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sz",ib+1,"Sz",jb+1, "Sy",ib+1,"Sy",jb+1;
        if(i<j) ampo += exc,"Sz",ia+1,"Sz",ja+1, "Sz",ib+1,"Sz",jb+1, "Sz",ib+1,"Sz",jb+1;
	// if(i<j) std::cout << ia+1 << "-" << ja+1 << "-" << ib+1 << "-" << jb+1 << "-" << ib+1 << "-" << jb+1 << std::endl;
    }
    // ----------------------------------
    // exit(1);


    auto H = MPO(ampo);

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i) {
        if(i%2 == 1)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
    }
    auto psi = MPS(state);

    printfln("Initial energy = %.5f", overlap(psi,H,psi) );

    // Begin the DMRG calculation
    auto energy = dmrg(psi,H,sweeps,args);

    // Print the final energy reported by DMRG
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing overlap = %.10f", overlap(psi,H,psi) );

    writeToFile(std::string("sites.txt"),sites); //file name will be sites_100
    writeToFile(std::string("psi.txt"),psi);     //file name will be psi_100






    // =========================================================
    // =================== Observables =========================
    // =========================================================
    std::cout << " Starting calculations of Observables " << std::endl;
    std::ofstream outfile;
    outfile.open("Observables.dat");

    Cplx SzTot=0, LzTot=0, JzTot=0;
    println("\nj Sx Sy Sz Lx Ly Lz (Sz+Lz) = ");
    for(int j=0; j < N; ++j) {
        //re-gauge psi to get ready to measure at position j
        int ja = j*orbitals + 0 + 1;
        int jb = j*orbitals + 1 + 1;

        // ==== lower orbital - spin index ======
        psi.position(ja);
        ITensor ket = psi.A(ja);
        ITensor bra = dag(prime(ket,Site));

        ITensor Sxjop = sites.op("Sx",ja); //*sites.op("Sx",j);
        ITensor Syjop = sites.op("Sy",ja); //*sites.op("Sy",j);
        ITensor Szjop = sites.op("Sz",ja); //*sites.op("Sz",j);

        //take an inner product
        auto sxj = (bra*Sxjop*ket).cplx();
        auto syj = (bra*Syjop*ket).cplx();
        auto szj = (bra*Szjop*ket).cplx();
        //printfln("%d %.12f %.12f %.12f",j,sxj,syj,szj);


        // ==== upper orbital - orbital index ======
        psi.position(jb);
        ket = psi.A(jb);
        bra = dag(prime(ket,Site));

        ITensor Lxjop = sites.op("Sx",jb); //*sites.op("Sx",j);
        ITensor Lyjop = sites.op("Sy",jb); //*sites.op("Sy",j);
        ITensor Lzjop = sites.op("Sz",jb); //*sites.op("Sz",j);

        //take an inner product
        auto Lxj = (bra*Lxjop*ket).cplx();
        auto Lyj = (bra*Lyjop*ket).cplx();
        auto Lzj = (bra*Lzjop*ket).cplx();


                if((j>10) && (j<=N-10)) {
                        SzTot += Cplx(szj); //.cplx();
                        LzTot += Cplx(Lzj);
                        JzTot += Cplx(szj)+Cplx(Lzj);
                }

        printfln("%d %.12f %.12f %.12f %.12f %.12f %.12f %.12f",j,sxj,syj,szj,Lxj,Lyj,Lzj,szj+Lzj);
        outfile << j << "  " << sxj << "  " << syj << "  " << szj << "  "
                << "  " << Lxj << "  " << Lyj << "  " << Lzj << "  " << szj+Lzj << "  \n";

    }
    outfile << std::endl;
    outfile.close();


    printfln("\nTotal Sz, Lz, Jz = %.12f %.12f %.12f",SzTot/(N-20),LzTot/(N-20),JzTot/(N-20));




    return 0;
    }













// ---------------- Honeycomb ---------------------------
void Chain1(int N, int Norb, Eigen::VectorXi Nc_,
                Eigen::MatrixXi& N1neigh_, Eigen::VectorXi indx_) {

    std::cout << "creating chain: " << std::endl;
    int Nsite_ = N;
    bool IsPeriodicX=false;

    // Site labeling
    indx_.resize(Nsite_);
    Nc_.resize(Nsite_);

    for(int i=0; i<Nsite_; i++){
        indx_[i] = i;
        Nc_[i] = i;
    }
    std::cout << Nc_ << std::endl;

    N1neigh_.resize(Nsite_,1);
    for(int i=0;i<Nsite_;i++){ 	// ith site
        if (i < Nsite_-1) N1neigh_(i,0) = i+1;

        if(IsPeriodicX==true && i==Nsite_-1) {
            N1neigh_(i,0) = 0;
        } else if (IsPeriodicX==false && i==Nsite_-1){
            N1neigh_(i,0) = -1;
        }
    }

    std::cout << " 1st Nearest neighbors " << std::endl;
    std::cout << N1neigh_ << std::endl;
    std::cout << std::endl;

} // end function


