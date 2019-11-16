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
    for(int i=1;i<N;i++){ 	// ith site
	// === Pinning ====
	if(i==1) ampo += BSzPinning,"Sz",i+1; 

	// === Spin - Heisenberg term === 
        double exc;
        // === Li.Lj ===
        exc = JExch;
        ampo += exc,"Sx",i,"Sx",i+1;
        ampo += exc,"Sy",i,"Sy",i+1;
        ampo += exc,"Sz",i,"Sz",i+1;

        // === Li.Lj Li.Lj ===
/*        exc = JExch/4.0;
        ampo += exc,"Sx",i,"Sx",i+1, "Sx",i,"Sx",i+1;
        ampo += exc,"Sx",i,"Sx",i+1, "Sy",i,"Sy",i+1;
        ampo += exc,"Sx",i,"Sx",i+1, "Sz",i,"Sz",i+1;
        ampo += exc,"Sy",i,"Sy",i+1, "Sx",i,"Sx",i+1;
        ampo += exc,"Sy",i,"Sy",i+1, "Sy",i,"Sy",i+1;
        ampo += exc,"Sy",i,"Sy",i+1, "Sz",i,"Sz",i+1;
        ampo += exc,"Sz",i,"Sz",i+1, "Sx",i,"Sx",i+1;
        ampo += exc,"Sz",i,"Sz",i+1, "Sy",i,"Sy",i+1;
        ampo += exc,"Sz",i,"Sz",i+1, "Sz",i,"Sz",i+1;
*/
	// if(i<j) std::cout << ia+1 << "-" << ja+1 << "-" << ib+1 << "-" << jb+1 << "-" << ib+1 << "-" << jb+1 << std::endl;
    }
    // ----------------------------------
    // exit(1);

    /*
    std::cout << "Connx = \n" << Connx << " \n "
              << "Conny = \n" << Conny << " \n "
              << "Connz = \n" << Connz << " \n ";

    std::cout << "Sorx = \n" << Sorx << " \n "
              << "Sory = \n" << Sory << " \n "
              << "Sorz = \n" << Sorz << " \n ";
    */
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






   /* // =========================================================
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





//     // Sum total S.S to check thata it's
//     // equal to ground state energy
//     outfile.open("Observables.dat", std::ios_base::app);
//     Real totalSdS = 0.;
//     println("\nj S_{i}.S_{i+1} = ");
//     for(int b = 1; b < N; ++b) {

//         psi.position(b);
//         ITensor bondket = psi.A(b)*psi.A(b+1);
//         ITensor bondbra = dag(prime(bondket,Site));

//         ITensor pmop = 0.5*sites.op("S+",b)*sites.op("S-",b+1);
//         ITensor mpop = 0.5*sites.op("S-",b)*sites.op("S+",b+1);
//         ITensor zzop = sites.op("Sz",b)*sites.op("Sz",b+1);

//         auto zz = (bondbra*zzop*bondket).cplx();
//         auto xx = (bondbra*pmop*bondket).cplx();
//         auto yy = (bondbra*mpop*bondket).cplx();

//         printfln("%d %.12f %.12f %.12f",b,xx.real(), yy.real(), zz.real());
//         totalSdS += xx.real()+yy.real()+zz.real();
//     }

//     printfln("\nSum of S.S = %.12f",totalSdS);
//     printfln("Ground state energy from DMRG = %.12f",energy);



//     Eigen::MatrixXd SxSx(N,N), SySy(N,N), SzSz(N,N);
//     SxSx.setZero(); SySy.setZero(); SzSz.setZero();
//     // http://itensor.org/docs.cgi?page=formulas/correlator_mps
//     // Given an MPS or IQMPS called "psi",
//     // constructed from a SiteSet "sites"
//     for(int i=1; i<=N; i++) {
//         auto op_i = sites.op("Sx",i);
//         //below we will assume j > i
//         psi.position(i); //'gauge' the MPS to site i

//         for(int j=i+1; j<=N;j++) {

//             if(j>N) continue;
//             std::cout << i << "-" << j << std::endl;

////             auto op_j = sites.op("Sx",j);

////             //psi.Anc(1) *= psi.A(0); //Uncomment if doing iDMRG calculation
////             //index linking i to i+1:
////             auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
////             auto C = psi.A(i)*op_i*dag(prime(psi.A(i),Site,ir));
////             for(int k = i+1; k < j; ++k) {
////                 C *= psi.A(k);
////                 C *= dag(prime(psi.A(k),Link));
////             }

////             C *= psi.A(j);
////             C *= op_j;

////             auto jl = commonIndex(psi.A(j),psi.A(j-1),Link); //index linking j to j-1:
////             C *= dag(prime(psi.A(j),jl,Site));

////             auto result = C.real();                            // or C.cplx() if expecting complex
////             SxSx(i-1,j-1) = result;                            // std::cout << i << "-" << j << " ==> " << result << std::endl;

//         }
//     }

//     std::cout << std::setprecision(5);
//     std::cout << "SxSx = \n" << SxSx << " \n" << std::endl;

*/


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


