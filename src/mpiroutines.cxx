/*! \file mpiroutines.cxx
 *  \brief this file contains routines used with MPI compilation.

    MPI routines generally pertain to domain decomposition or to specific MPI tasks that determine what needs to be broadcast between various threads.

 */

#ifdef USEMPI

#include <cassert>
#include <tuple>
#include <random>

//-- For MPI

#include "stf.h"
#include "logging.h"
#include "timer.h"

#ifdef SWIFTINTERFACE
#include "swiftinterface.h"
using namespace Swift;
#endif

/// function that sends information between threads 
static std::tuple<std::vector<Int_t>, std::vector<float>>
exchange_indices_and_props(
    const std::vector<Int_t> &indices, const std::vector<float> &props,
    std::size_t props_per_index, int rank, int tag, MPI_Comm &mpi_comm)
{
    Int_t num_indices = indices.size();
    Int_t num_props = num_indices * props_per_index;
    assert(num_indices <= std::numeric_limits<std::int32_t>::max());
    assert(props.size() == num_props);

    // Send/recv number of indices to allocate reception buffers
    Int_t num_indices_recv;
    MPI_Status status;
    MPI_Sendrecv(
        &num_indices, 1, MPI_Int_t, rank, tag * 2,
        &num_indices_recv, 1, MPI_Int_t, rank, tag * 2,
        mpi_comm, &status);

    // Send/recv actual indices and properties
    std::vector<Int_t> indices_recv(num_indices_recv);
    Int_t num_props_recv = num_indices_recv * props_per_index;
    std::vector<float> props_recv(num_props_recv);
    MPI_Sendrecv(
        indices.data(), num_indices, MPI_Int_t, rank, tag * 3,
        indices_recv.data(), num_indices_recv, MPI_Int_t, rank, tag * 3,
        mpi_comm, &status);
    MPI_Sendrecv(
        props.data(), num_props, MPI_FLOAT, rank, tag * 4,
        props_recv.data(), num_props_recv, MPI_FLOAT, rank, tag * 4,
        mpi_comm, &status);
    return {std::move(indices_recv), std::move(props_recv)};
}

/// @def Routines to assist point-to-point communication 
//@{

/// @brief Generate a vector list of mpi ranks that need to communicate
/// @param send_info information containing number of items to be sent from one task to another
/// @return vector of pair of mpi process ranks
inline vector<tuple<int, int>> MPIGenerateCommPairs(Int_t *send_info)
{
    vector<tuple<int, int>> commpair;
    for(auto task1 = 0; task1 < NProcs; task1++)
    {
        for(auto task2 = task1+1; task2 < NProcs; task2++)
        {
                if (send_info[task1 * NProcs + task2] == 0 && send_info[task2 * NProcs + task1] == 0) continue;
            commpair.push_back(make_tuple(task1, task2));
        }
    }
    // to ensure that ThisTask = 0 doesn't dominate the first set of communication pairs,
    // randomize pairs
    unsigned seed = 4322;
    shuffle(commpair.begin(), commpair.end(), std::default_random_engine(seed));
    return commpair;
}
/// @brief Set which task is sending and receiving
/// @param task1 one of the tasks in the pair
/// @param task2 second task in the pair 
/// @return 
inline tuple<int, int> MPISetSendRecvTask(int task1, int task2){
    auto sendTask = task1, recvTask = task2;
    if (ThisTask == task2)
    {
        recvTask = task1;
        sendTask = task2;
    }
    return make_tuple(sendTask, recvTask);
}


/// @brief Initialize the number of chunks and chunksize for a communication
/// @param nsend number of times to send
/// @param nrecv number of items to receive
/// @param maxchunksize max number to send in one go
/// @return tuple of number of send/receives, current send chunk size, current receive chunk size, send off access, receive offset access
inline tuple<int, int, int, Int_t, Int_t> MPIInitialzeCommChunks(Int_t nsend, Int_t nrecv, Int_t maxchunksize)
{
    //send info in loops to minimize memory footprint
    int cursendchunksize, currecvchunksize, nsendchunks, nrecvchunks, numsendrecv;
    Int_t sendoffset, recvoffset;
    cursendchunksize = currecvchunksize = maxchunksize;
    nsendchunks=ceil(static_cast<Double_t>(nsend)/static_cast<Double_t>(maxchunksize));
    nrecvchunks=ceil(static_cast<Double_t>(nrecv)/static_cast<Double_t>(maxchunksize));
    if (cursendchunksize > nsend) {
        nsendchunks = 1;
        cursendchunksize = nsend;
    }
    if (currecvchunksize > nrecv) {
        nrecvchunks = 1;
        currecvchunksize = nrecv;
    }
    numsendrecv = std::max(nsendchunks, nrecvchunks);
    sendoffset = recvoffset = 0;
    return make_tuple(numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset);
}

/// @brief Initialize the number of chunks and chunksize for a communication
/// @param nsend number of times to send
/// @param nrecv number of items to receive
/// @param cursendchunksize current number of items sent
/// @param currecvchunksize current number of items received
/// @param sendoffset current send offset
/// @param recvoffset current receive offset
inline void MPIUpdateCommChunks(Int_t nsend, Int_t nrecv, 
int &cursendchunksize, int &currecvchunksize, 
Int_t &sendoffset, Int_t &recvoffset)
{
    sendoffset += cursendchunksize;
    recvoffset += currecvchunksize;
    cursendchunksize = std::min(static_cast<Int_t>(cursendchunksize), nsend - sendoffset);
    currecvchunksize = std::min(static_cast<Int_t>(currecvchunksize), nrecv - recvoffset);
}

//@}

/// @name Domain decomposition routines and io routines to place particles correctly in local mpi data space
//@{

/// @brief Using bisection distance mpi_dxsplit, determine mpi decomposition. Here domains are constructured in data units
/// @param opt Options structure containing runtime arguments
void MPIInitialDomainDecomposition(Options &opt)
{
    if (opt.impiusemesh) {
        MPIInitialDomainDecompositionWithMesh(opt);
        return;
    }
    Int_t i,j,k;
    int Nsplit,isplit;
    Double_t diffsplit;
    int b,a;

    if (ThisTask==0) {
        //first split need not be simply having the dimension but determine
        //number of splits to have Nprocs=a*2^b, where a and b are integers
        //initial integers
        b=floor(log((float)NProcs)/log(2.0))-1;
        a=floor(NProcs/pow(2,b));
        diffsplit=(double)NProcs/(double)a/(double)pow(2,b);
        while (diffsplit!=1) {
            b--;
            a=floor(NProcs/pow(2,b));
            diffsplit=(double)NProcs/(double)a/(double)pow(2,b);
        }
        Nsplit=b+1;
        mpi_ideltax[0]=0;mpi_ideltax[1]=1;mpi_ideltax[2]=2;
        isplit=0;
        for (j=0;j<3;j++) mpi_nxsplit[j]=0;
        for (j=0;j<Nsplit;j++) {
            mpi_nxsplit[mpi_ideltax[isplit++]]++;
            if (isplit==3) isplit=0;
        }
        for (j=0;j<3;j++) mpi_nxsplit[j]=pow(2.0,mpi_nxsplit[j]);
        //and adjust first dimension
        mpi_nxsplit[0]=mpi_nxsplit[0]/2*a;

        //for all the cells along the boundary of axis with the third split axis (smallest variance)
        //set the domain limits to the sims limits
        int ix=mpi_ideltax[0],iy=mpi_ideltax[1],iz=mpi_ideltax[2];
        int mpitasknum;
        for (j=0;j<mpi_nxsplit[iy];j++) {
            for (i=0;i<mpi_nxsplit[ix];i++) {
                mpitasknum=i+j*mpi_nxsplit[ix]+0*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iz][0]=mpi_xlim[iz][0];
                mpitasknum=i+j*mpi_nxsplit[ix]+(mpi_nxsplit[iz]-1)*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iz][1]=mpi_xlim[iz][1];
            }
        }
        //here for domains along second axis
        for (k=0;k<mpi_nxsplit[iz];k++) {
            for (i=0;i<mpi_nxsplit[ix];i++) {
                mpitasknum=i+0*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iy][0]=mpi_xlim[iy][0];
                mpitasknum=i+(mpi_nxsplit[iy]-1)*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[iy][1]=mpi_xlim[iy][1];
            }
        }
        //finally along axis with largest variance
        for (k=0;k<mpi_nxsplit[iz];k++) {
            for (j=0;j<mpi_nxsplit[iy];j++) {
                mpitasknum=0+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[ix][0]=mpi_xlim[ix][0];
                mpitasknum=(mpi_nxsplit[ix]-1)+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                mpi_domain[mpitasknum].bnd[ix][1]=mpi_xlim[ix][1];
            }
        }
        //here use the three different histograms to define the boundary
        Double_t bndval[3];
        for (i=0;i<mpi_nxsplit[ix];i++) {
            bndval[0]=(mpi_xlim[ix][1]-mpi_xlim[ix][0])*(Double_t)(i+1)/(Double_t)mpi_nxsplit[ix];
            if(i<mpi_nxsplit[ix]-1) {
            for (j=0;j<mpi_nxsplit[iy];j++) {
                for (k=0;k<mpi_nxsplit[iz];k++) {
                    //define upper limit
                    mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                    mpi_domain[mpitasknum].bnd[ix][1]=bndval[0];
                    //define lower limit
                    mpitasknum=(i+1)+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                    mpi_domain[mpitasknum].bnd[ix][0]=bndval[0];
                }
            }
            }
            //now for secondary splitting
            if (mpi_nxsplit[iy]>1)
            for (j=0;j<mpi_nxsplit[iy];j++) {
                bndval[1]=(mpi_xlim[iy][1]-mpi_xlim[iy][0])*(Double_t)(j+1)/(Double_t)mpi_nxsplit[iy];
                if(j<mpi_nxsplit[iy]-1) {
                for (k=0;k<mpi_nxsplit[iz];k++) {
                    mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                    mpi_domain[mpitasknum].bnd[iy][1]=bndval[1];
                    mpitasknum=i+(j+1)*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                    mpi_domain[mpitasknum].bnd[iy][0]=bndval[1];
                }
                }
                if (mpi_nxsplit[iz]>1)
                for (k=0;k<mpi_nxsplit[iz];k++) {
                    bndval[2]=(mpi_xlim[iz][1]-mpi_xlim[iz][0])*(Double_t)(k+1)/(Double_t)mpi_nxsplit[iz];
                    if (k<mpi_nxsplit[iz]-1){
                    mpitasknum=i+j*mpi_nxsplit[ix]+k*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                    mpi_domain[mpitasknum].bnd[iz][1]=bndval[2];
                    mpitasknum=i+j*mpi_nxsplit[ix]+(k+1)*(mpi_nxsplit[ix]*mpi_nxsplit[iy]);
                    mpi_domain[mpitasknum].bnd[iz][0]=bndval[2];
                    }
                }
            }
        }
        LOG(info) << "Initial MPI Domains are:";
        for (int j=0;j<NProcs;j++) {
            std::ostringstream os;
            os << " ThisTask= " << j << " :: ";
            for (int k=0;k<3;k++)
                os << k << " " << mpi_domain[j].bnd[k][0] << " " << mpi_domain[j].bnd[k][1] << " | ";
            LOG(info) << os.str();
        }
    }
    //broadcast data
    MPI_Bcast(mpi_domain, NProcs*sizeof(MPI_Domain), MPI_BYTE, 0, MPI_COMM_WORLD);
}

/// @brief Using mesh and space-filling Z curve, determine mpi decomposition. Here domains are constructured in data units
/// @param opt Options structure containing runtime arguments
void MPIInitialDomainDecompositionWithMesh(Options &opt){
    if (ThisTask==0) {
        // each processor takes subsection of volume, where the grow in number of cells per proc 
        // is set such that the cells per dim grows as (log(N)/log(2))^(y). Currently y=1
        // unless the numcellsperdim has been set, in which case that is used. 
        if (opt.numcellsperdim == 0) {
            unsigned int NProcfac = std::ceil(log(static_cast<double>(NProcs))/log(2.0));
            opt.numcellsperdim = opt.minnumcellperdim*std::max(NProcfac, 1u);
        }
        unsigned int n3 = opt.numcells = opt.numcellsperdim*opt.numcellsperdim*opt.numcellsperdim;
        double idelta = 1.0/(double)opt.numcellsperdim;
        for (auto i=0; i<3; i++) {
            opt.spacedimension[i] = (mpi_xlim[i][1] - mpi_xlim[i][0]);
            opt.cellwidth[i] = (opt.spacedimension[i] * idelta);
            opt.icellwidth[i] = 1.0/opt.cellwidth[i];
        }

        //now order according to Z-curve or Morton curve
        //first fill curve
        vector<bitset<16>> zcurvevalue(3);
        struct zcurvestruct{
            unsigned int coord[3];
            unsigned long long index;
            bitset<48> zcurvevalue;
        };
        vector<zcurvestruct> zcurve(n3);
        unsigned long long index;
        for (auto ix=0;ix<opt.numcellsperdim;ix++) {
            for (auto iy=0;iy<opt.numcellsperdim;iy++) {
                for (auto iz=0;iz<opt.numcellsperdim;iz++) {
                    index = ix*opt.numcellsperdim*opt.numcellsperdim + iy*opt.numcellsperdim + iz;
                    zcurve[index].coord[0] = ix;
                    zcurve[index].coord[1] = iy;
                    zcurve[index].coord[2] = iz;
                    zcurve[index].index = index;
                    zcurvevalue[0] = bitset<16>(ix);
                    zcurvevalue[1] = bitset<16>(iy);
                    zcurvevalue[2] = bitset<16>(iz);
                    for (auto j=0; j<16;j++)
                    {
                        for (auto i = 0; i<3; i++) {
                            zcurve[index].zcurvevalue[j*3+i] = zcurvevalue[i][j];
                        }
                    }
                }
            }
        }
        //then sort index array based on the zcurvearray value
        sort(zcurve.begin(), zcurve.end(), [](const zcurvestruct &a, const zcurvestruct &b){
            return a.zcurvevalue.to_ullong() < b.zcurvevalue.to_ullong();
        });
        //finally assign cells to tasks
        opt.cellnodeids.resize(n3);
        opt.cellnodeorder.resize(n3);
        opt.cellloc = new cell_loc[n3];
        int nsub = max((int)floor(n3/(double)NProcs), 1);
        int itask = 0, count = 0;
        vector<int> numcellspertask(NProcs,0);
        for (auto i=0;i<n3;i++)
        {
            if (count == nsub) {
                count = 0;
                itask++;
            }
            if (itask == NProcs) itask -= 1;
            opt.cellnodeids[zcurve[i].index] = itask;
            opt.cellnodeorder[i] = zcurve[i].index;
            numcellspertask[itask]++;
            count++;
        }
        LOG(info) << "Z-curve Mesh MPI decomposition:";
        LOG(info) << " Mesh has resolution of " << opt.numcellsperdim << " per spatial dim";
        LOG(info) << " with each mesh spanning (" << opt.cellwidth[0] << ", " << opt.cellwidth[1] << ", " << opt.cellwidth[2] << ")";
        LOG(info) << "MPI tasks :";
        for (auto i=0; i<NProcs; i++) {
            LOG(info) << " Task "<< i << " has " << numcellspertask[i] / double(n3) << " of the volume";
        }
    }
    //broadcast data
    MPI_Bcast(&opt.numcells, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(&opt.numcellsperdim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(opt.spacedimension, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(opt.cellwidth, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(opt.icellwidth, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (ThisTask != 0) {
        opt.cellnodeids.resize(opt.numcells);
        opt.cellnodeorder.resize(opt.numcells);
    }
    opt.cellnodenumparts.resize(opt.numcells,0);
    MPI_Bcast(opt.cellnodeids.data(), opt.numcells, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(opt.cellnodeorder.data(), opt.numcells, MPI_INTEGER, 0, MPI_COMM_WORLD);
}

/// @brief Find load imbalance ((max-min)/expected average)
/// @param opt Options structure containing runtime arguments
inline double MPILoadBalanceWithMesh(Options &opt) {
    //calculate imbalance based on min and max in mpi domains
    vector<Int_t> mpinumparts(NProcs, 0);
    for (auto i=0;i<opt.numcells;i++)
    {
        auto itask = opt.cellnodeids[i];
        mpinumparts[itask] += opt.cellnodenumparts[i];
    }
    double minval, maxval, ave, std, sum;
    minval = maxval = mpinumparts[0];
    ave = std = sum = 0;
    for (auto &x:mpinumparts) {
        if (minval > x) minval = x;
        if (maxval < x) maxval = x;
        ave += x;
        std += x*x;
    }
    ave /= (double)NProcs;
    std /= (double)NProcs;
    std = sqrt(std - ave*ave);
    return (maxval-minval)/ave;
}

/// @brief Using mesh and space-filling Z curve redo mpi decomposition to improve load balance
/// @param opt Options structure containing runtime arguments
bool MPIRepartitionDomainDecompositionWithMesh(Options &opt){
    Int_t *buff = new Int_t[opt.numcells];
    for (auto i=0;i<opt.numcells;i++) buff[i]=0;
    MPI_Allreduce(opt.cellnodenumparts.data(), buff, opt.numcells, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
    for (auto i=0;i<opt.numcells;i++) opt.cellnodenumparts[i]=buff[i];
    delete[] buff;
    double optimalave = 0; for (auto i=0;i<opt.numcells;i++) optimalave += opt.cellnodenumparts[i];
    optimalave /= (double)NProcs;
    auto loadimbalance = MPILoadBalanceWithMesh(opt);
    LOG_RANK0(info) << "MPI imbalance of " << loadimbalance;
    if (loadimbalance > opt.mpimeshimbalancelimit) {
        LOG_RANK0(info) << "Imbalance too large, adjusting MPI domains ...";
        int itask = 0;
        Int_t numparts = 0 ;
        vector<int> numcellspertask(NProcs,0);
        vector<int> mpinumparts(NProcs,0);
        for (auto i=0;i<opt.numcells;i++)
        {
            auto index = opt.cellnodeorder[i];
            numcellspertask[itask]++;
            opt.cellnodeids[index] = itask;
            numparts += opt.cellnodenumparts[index];
            if (numparts > optimalave && itask < NProcs-1) {
                mpinumparts[itask] = numparts;
                itask++;
                numparts = 0;
            }
        }
        mpinumparts[NProcs-1] = numparts;
        if (ThisTask == 0) {
            for (auto x:mpinumparts) if (x == 0) {
                LOG(error) << "MPI Process has zero particles associated with it, likely due to too many mpi tasks requested or too coarse a mesh used.";
                LOG(error) << "Current number of tasks: " << NProcs;
                LOG(error) << "Current mesh resolution " << opt.numcellsperdim << "^3";
                Int_t sum = 0;
                for (auto x:opt.cellnodenumparts) {
                    sum +=x;
                }
                LOG(error) << "Total number of particles loaded " << sum;
                LOG(error) << "Suggested number of particles per mpi processes is > 1e7";
                LOG(error) << "Suggested number of mpi processes using 1e7 is " << (int)ceil(sum / 1e7);
                LOG(error) << "Increase mesh resolution or reduce MPI Processes ";
                MPI_Abort(MPI_COMM_WORLD,8);
            }
            LOG(info) << "Now have MPI imbalance of " << MPILoadBalanceWithMesh(opt);
            LOG(info) << "MPI tasks:";
            for (auto i=0; i<NProcs; i++) {
                LOG(info) << " Task " << i << " has " << numcellspertask[i] / double(opt.numcells) << " of the volume";
            }
        }
        for (auto &x:opt.cellnodenumparts) x=0;
        Nlocal = mpinumparts[ThisTask];
        //only need to reread file if baryon search active to determine number of
        //baryons in local mpi domain. Otherwise, simple enough to update Nlocal
        return (opt.iBaryonSearch>0);
    }
    return false;
}

void MPINumInDomain(Options &opt)
{
    
    //when reading number in domain, use all available threads to read all available files
    //first set number of read threads to either total number of mpi process or files, which ever is smaller
    //store old number of read threads
    if (NProcs==1) {
    Nlocal=Ntotal;
    Nmemlocal=Nlocal;
    return;
    }
    int nsnapread=opt.nsnapread;
    //opt.nsnapread=min(NProcs,opt.num_files);
    if(opt.inputtype==IOTIPSY) MPINumInDomainTipsy(opt);
    else if (opt.inputtype==IOGADGET) MPINumInDomainGadget(opt);
    else if (opt.inputtype==IORAMSES) MPINumInDomainRAMSES(opt);
#ifdef USEHDF
    else if (opt.inputtype==IOHDF) MPINumInDomainHDF(opt);
#endif
    if (ThisTask == 0) {
        if (Ntotal/1e7 < NProcs) {
            LOG(warning) << "Suggested number of particles per mpi processes is roughly > 1e7";
            LOG(warning) << "Number of MPI tasks greater than this suggested number";
            LOG(warning) << "May result in poor performance";
        }
    }
    //if using mesh, check load imbalance and also repartition
    if (opt.impiusemesh) {
        if (MPIRepartitionDomainDecompositionWithMesh(opt)){
            if(opt.inputtype==IOTIPSY) MPINumInDomainTipsy(opt);
            else if (opt.inputtype==IOGADGET) MPINumInDomainGadget(opt);
            else if (opt.inputtype==IORAMSES) MPINumInDomainRAMSES(opt);
#ifdef USEHDF
            else if (opt.inputtype==IOHDF) MPINumInDomainHDF(opt);
#endif
        }
    }
    opt.nsnapread=nsnapread;
    //adjust the memory allocated to allow some buffer room.
    Nmemlocal=Nlocal*(1.0+opt.mpipartfac);
    if (opt.iBaryonSearch) Nmemlocalbaryon=Nlocalbaryon[0]*(1.0+opt.mpipartfac);

}

void MPIDomainExtent(Options &opt)
{
    if(opt.inputtype==IOTIPSY) MPIDomainExtentTipsy(opt);
    else if (opt.inputtype==IOGADGET) MPIDomainExtentGadget(opt);
    else if (opt.inputtype==IORAMSES) MPIDomainExtentRAMSES(opt);
#ifdef USEHDF
    else if (opt.inputtype==IOHDF) MPIDomainExtentHDF(opt);
#endif
}

void MPIDomainDecomposition(Options &opt)
{
    MPIInitialDomainDecomposition(opt);
    if(opt.inputtype==IOTIPSY) MPIDomainDecompositionTipsy(opt);
    else if (opt.inputtype==IOGADGET) MPIDomainDecompositionGadget(opt);
    else if (opt.inputtype==IORAMSES) MPIDomainDecompositionRAMSES(opt);
#ifdef USEHDF
    else if (opt.inputtype==IOHDF) MPIDomainDecompositionHDF(opt);
#endif
}

/// @brief Adjust mpi domain decomposition to set it to internal runtime units
/// @param opt Options structure containing runtime arguments
void MPIAdjustDomain(Options &opt){
    if (opt.impiusemesh) {
        MPIAdjustDomainWithMesh(opt);
        return;
    }
    Double_t aadjust, lscale;
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    lscale=opt.lengthinputconversion*aadjust;
    if (opt.inputcontainslittleh) lscale /=opt.h;
    for (int j=0;j<NProcs;j++) for (int k=0;k<3;k++) {mpi_domain[j].bnd[k][0]*=lscale;mpi_domain[j].bnd[k][1]*=lscale;}
}

/// @brief Adjust mpi domain decomposition to set it to internal runtime units
/// @param opt Options structure containing runtime arguments
void MPIAdjustDomainWithMesh(Options &opt)
{
    //once data loaded update cell widths
    Double_t aadjust, lscale;
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    lscale=opt.lengthinputconversion*aadjust;
    if (opt.inputcontainslittleh) lscale /=opt.h;
    for (auto i=0;i<3;i++) {
        opt.spacedimension[i] *= lscale;
        opt.cellwidth[i] *= lscale;
        opt.icellwidth[i] /= lscale;
    }
}

//@}

/// @brief given a position and a mpi thread domain information, determine which mpi process a particle is assigned to
/// @param opt Options structure containing runtime arguments
/// @param x x position
/// @param y y position
/// @param z z position
/// @return mpi process rank in MPI_COMM_WORLD
int MPIGetParticlesProcessor(Options &opt, Double_t x, Double_t y, Double_t z){
    if (NProcs==1) return 0;
    if (opt.impiusemesh) {
        unsigned int ix, iy, iz;
        unsigned long long index;
        ix=floor(x*opt.icellwidth[0]);
        iy=floor(y*opt.icellwidth[1]);
        iz=floor(z*opt.icellwidth[2]);
        index = ix*opt.numcellsperdim*opt.numcellsperdim+iy*opt.numcellsperdim+iz;
        opt.cellnodenumparts[index]++;
        if (index >= 0 && index < opt.numcells) return opt.cellnodeids[index];
    }
    else {
        for (int j=0;j<NProcs;j++){
            if( (mpi_domain[j].bnd[0][0]<=x) && (mpi_domain[j].bnd[0][1]>=x)&&
                (mpi_domain[j].bnd[1][0]<=y) && (mpi_domain[j].bnd[1][1]>=y)&&
                (mpi_domain[j].bnd[2][0]<=z) && (mpi_domain[j].bnd[2][1]>=z) )
                return j;
        }
    }
    LOG(error) << "Particle outside the mpi domains of every process (" << x << "," << y << "," << z << ")";
    MPI_Abort(MPI_COMM_WORLD,9);
    return -1;
}

/// @def Routines related to managing extra properties of baryon particles 
//@{
void MPIStripExportParticleOfExtraInfo(Options &opt, Int_t n, Particle *Part)
{
#if defined(GASON) || defined(STARON) || defined(BHON) || defined(EXTRADMON)
    Int_t numextrafields = 0;
#ifdef GASON
    numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields > 0)
    {
        for (auto i=0;i<n;i++)
        {
            if (Part[i].HasHydroProperties()) Part[i].SetHydroProperties();
        }
    }
#endif
#ifdef STARON
    numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields > 0)
    {
        for (auto i=0;i<n;i++)
        {
            if (Part[i].HasStarProperties()) Part[i].SetStarProperties();
        }
    }
#endif
#ifdef BHON
    numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields > 0)
    {
        for (auto i=0;i<n;i++)
        {
            if (Part[i].HasBHProperties()) Part[i].SetBHProperties();
        }
    }
#endif
#ifdef EXTRADMON
    numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields > 0)
    {
        for (auto i=0;i<n;i++)
        {
            if (Part[i].HasExtraDMProperties()) Part[i].SetExtraDMProperties();
        }
    }
#endif
#endif
}

void MPIFillBuffWithHydroInfo(Options &opt, Int_t nlocalbuff, Particle *Part, vector<Int_t> &indices, vector<float> &propbuff, bool resetbuff)
{
#ifdef GASON
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    string field;
    indices.clear();
    propbuff.clear();

    numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasHydroProperties()) indices.push_back(i);
    num = indices.size();

    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetInternalProperties(field);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetChemistry(field);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetChemistryProduction(field);
        }
        if (resetbuff) Part[index].SetHydroProperties();
    }
#endif
}

void MPIFillBuffWithStarInfo(Options &opt, Int_t nlocalbuff, Particle *Part, vector<Int_t> &indices, vector<float> &propbuff, bool resetbuff)
{
#ifdef STARON
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    string field;
    indices.clear();
    propbuff.clear();

    numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasStarProperties()) indices.push_back(i);
    num = indices.size();

    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetInternalProperties(field);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetChemistry(field);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetChemistryProduction(field);
        }
        if (resetbuff) Part[index].SetStarProperties();
    }
#endif
}

void MPIFillBuffWithBHInfo(Options &opt, Int_t nlocalbuff, Particle *Part, vector<Int_t> &indices, vector<float> &propbuff, bool resetbuff)
{
#ifdef BHON
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    string field;
    indices.clear();
    propbuff.clear();

    numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasBHProperties()) indices.push_back(i);
    num = indices.size();

    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetInternalProperties(field);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetChemistry(field);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetChemistryProduction(field);
        }
        if (resetbuff) Part[index].SetBHProperties();
    }
#endif
}

void MPIFillBuffWithExtraDMInfo(Options &opt, Int_t nlocalbuff, Particle *Part, vector<Int_t> &indices, vector<float> &propbuff, bool resetbuff)
{
#ifdef EXTRADMON
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    string field;
    indices.clear();
    propbuff.clear();

    numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasExtraDMProperties()) indices.push_back(i);
    num = indices.size();

    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetExtraDMProperties().GetExtraProperties(field);
        }
        if (resetbuff) Part[index].SetExtraDMProperties();
    }
#endif
}

void MPIFillFOFBuffWithHydroInfo(Options &opt, Int_t numexport, fofid_in *FoFGroupData, Particle *&Part, vector<Int_t> &indices, vector<float> &propbuff, bool iforexport)
{
#ifdef GASON
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    string field;
    indices.clear();
    propbuff.clear();

    numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<numexport;i++)
    {
        if (!FoFGroupData[i].p.HasHydroProperties()) continue;
        indices.push_back(i);
    }
    num = indices.size();
    if (num == 0) return;

    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetHydroProperties().GetInternalProperties(field);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetHydroProperties().GetChemistry(field);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetHydroProperties().GetChemistryProduction(field);
        }
        FoFGroupData[index].p.SetHydroProperties();
        if (iforexport) Part[FoFGroupData[index].Index].SetHydroProperties();
    }
#endif
}

void MPIFillFOFBuffWithStarInfo(Options &opt, Int_t numexport, fofid_in *FoFGroupData, Particle *&Part, vector<Int_t> &indices, vector<float> &propbuff, bool iforexport)
{
#ifdef STARON
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    string field;
    indices.clear();
    propbuff.clear();

    numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<numexport;i++)
    {
        if (!FoFGroupData[i].p.HasStarProperties()) continue;
        indices.push_back(i);
    }
    num = indices.size();
    if (num == 0) return;

    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetStarProperties().GetInternalProperties(field);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetStarProperties().GetChemistry(field);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetStarProperties().GetChemistryProduction(field);
        }
        FoFGroupData[index].p.SetStarProperties();
        if (iforexport) Part[FoFGroupData[index].Index].SetStarProperties();
    }
#endif
}

void MPIFillFOFBuffWithBHInfo(Options &opt, Int_t numexport, fofid_in *FoFGroupData, Particle *&Part, vector<Int_t> &indices, vector<float> &propbuff, bool iforexport)
{
#ifdef BHON
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    string field;
    indices.clear();
    propbuff.clear();

    numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<numexport;i++)
    {
        if (!FoFGroupData[i].p.HasBHProperties()) continue;
        indices.push_back(i);
    }
    num = indices.size();
    if (num == 0) return;

    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetBHProperties().GetInternalProperties(field);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetBHProperties().GetChemistry(field);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetBHProperties().GetChemistryProduction(field);
        }
        FoFGroupData[index].p.SetBHProperties();
        if (iforexport) Part[FoFGroupData[index].Index].SetBHProperties();
    }
#endif
}

void MPIFillFOFBuffWithExtraDMInfo(Options &opt, Int_t numexport, fofid_in *FoFGroupData, Particle *&Part, vector<Int_t> &indices, vector<float> &propbuff, bool iforexport)
{
#ifdef EXTRADMON
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    string field;
    indices.clear();
    propbuff.clear();

    numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<numexport;i++)
    {
        if (!FoFGroupData[i].p.HasExtraDMProperties()) continue;
        indices.push_back(i);
    }
    num = indices.size();
    if (num == 0) return;

    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = FoFGroupData[index].p.GetExtraDMProperties().GetExtraProperties(field);
        }
        FoFGroupData[index].p.SetExtraDMProperties();
        if (iforexport) Part[FoFGroupData[index].Index].SetExtraDMProperties();
    }
#endif
}
//@}

/// @def Routines related to sending information between threads when reading in data. 
//@{
void MPISendParticleInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID)
{
    vector<Int_t> indices_gas, indices_star, indices_bh, indices_extradm;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff_gas, propbuff_star, propbuff_bh, propbuff_extradm;
    string field;

#ifdef GASON
    numextrafields = opt.gas_internalprop_names.size() + opt.gas_chem_names.size() + opt.gas_chemproduction_names.size();
    if (numextrafields > 0)
    {
        for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasHydroProperties()) indices_gas.push_back(i);
        num = indices_gas.size();
        if (num>0) {
            propbuff_gas.resize(numextrafields*num);
            for (auto i=0;i<num;i++)
            {
                index = indices_gas[i];
                offset = 0;
                for (auto iextra=0;iextra<opt.gas_internalprop_names.size();iextra++)
                {
                    field = opt.gas_internalprop_names[iextra];
                    propbuff_gas[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetInternalProperties(field);
                }
                offset += opt.gas_internalprop_names.size();
                for (auto iextra=0;iextra<opt.gas_chem_names.size();iextra++)
                {
                    field = opt.gas_chem_names[iextra];
                    propbuff_gas[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetChemistry(field);
                }
                offset += opt.gas_chem_names.size();
                for (auto iextra=0;iextra<opt.gas_chemproduction_names.size();iextra++)
                {
                    field = opt.gas_chemproduction_names[iextra];
                    propbuff_gas[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetChemistryProduction(field);
                }
                Part[index].SetHydroProperties();
            }
        }
    }
#endif
#ifdef STARON
    numextrafields = opt.star_internalprop_names.size() + opt.star_chem_names.size() + opt.star_chemproduction_names.size();
    if (numextrafields > 0)
    {
        for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasStarProperties()) indices_star.push_back(i);
        num = indices_star.size();
        if (num > 0)
        {
            propbuff_star.resize(numextrafields*num);
            for (auto i=0;i<num;i++)
            {
                index = indices_star[i];
                offset = 0;
                for (auto iextra=0;iextra<opt.star_internalprop_names.size();iextra++)
                {
                    field = opt.star_internalprop_names[iextra];
                    propbuff_star[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetInternalProperties(field);
                }
                offset += opt.star_internalprop_names.size();
                for (auto iextra=0;iextra<opt.star_chem_names.size();iextra++)
                {
                    field = opt.star_chem_names[iextra];
                    propbuff_star[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetChemistry(field);
                }
                offset += opt.star_chem_names.size();
                for (auto iextra=0;iextra<opt.star_chemproduction_names.size();iextra++)
                {
                    field = opt.star_chemproduction_names[iextra];
                    propbuff_star[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetChemistryProduction(field);
                }
                Part[index].SetStarProperties();
            }
        }
    }
#endif
#ifdef BHON
    numextrafields = opt.bh_internalprop_names.size() + opt.bh_chem_names.size() + opt.bh_chemproduction_names.size();
    if (numextrafields > 0)
    {
        for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasBHProperties()) indices_bh.push_back(i);
        num = indices_bh.size();
        if (num > 0)
        {
            propbuff_bh.resize(numextrafields*num);
            for (auto i=0;i<num;i++)
            {
                index = indices_bh[i];
                offset = 0;
                for (auto iextra=0;iextra<opt.bh_internalprop_names.size();iextra++)
                {
                    field = opt.bh_internalprop_names[iextra];
                    propbuff_bh[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetInternalProperties(field);
                }
                offset += opt.bh_internalprop_names.size();
                for (auto iextra=0;iextra<opt.bh_chem_names.size();iextra++)
                {
                    field = opt.bh_chem_names[iextra];
                    propbuff_bh[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetChemistry(field);
                }
                offset += opt.bh_chem_names.size();
                for (auto iextra=0;iextra<opt.bh_chemproduction_names.size();iextra++)
                {
                    field = opt.bh_chemproduction_names[iextra];
                    propbuff_bh[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetChemistryProduction(field);
                }
                Part[index].SetBHProperties();
            }
        }
    }
#endif
#ifdef EXTRADMON
    numextrafields = opt.extra_dm_internalprop_names.size();
    if (numextrafields > 0)
    {
        for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasExtraDMProperties()) indices_extradm.push_back(i);
        num = indices_extradm.size();
        if (num>0) {
            propbuff_extradm.resize(numextrafields*num);
            for (auto i=0;i<num;i++)
            {
                index = indices_extradm[i];
                offset = 0;
                for (auto iextra=0;iextra<opt.extra_dm_internalprop_names.size();iextra++)
                {
                    field = opt.extra_dm_internalprop_names[iextra];
                    propbuff_extradm[i*numextrafields + iextra + offset] = Part[index].GetExtraDMProperties().GetExtraProperties(field);
                }
                Part[index].SetExtraDMProperties();
            }
        }
    }
#endif
    MPI_Ssend(Part, sizeof(Particle)*nlocalbuff, MPI_BYTE, taskID, taskID, MPI_COMM_WORLD);
#ifdef GASON
    numextrafields = opt.gas_internalprop_names.size() + opt.gas_chem_names.size() + opt.gas_chemproduction_names.size();
    if (numextrafields > 0)
    {
        num = indices_gas.size();
        MPI_Send(&num,sizeof(Int_t),MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
        if (num > 0)
        {
            MPI_Send(indices_gas.data(),sizeof(Int_t)*num,MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
            MPI_Send(propbuff_gas.data(),sizeof(float)*num*numextrafields,MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
        }
    }
#endif
#ifdef STARON
    numextrafields = opt.star_internalprop_names.size() + opt.star_chem_names.size() + opt.star_chemproduction_names.size();
    if (numextrafields > 0)
    {
        num = indices_star.size();
        MPI_Send(&num,sizeof(Int_t),MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
        if (num > 0)
        {
            MPI_Send(indices_star.data(),sizeof(Int_t)*num,MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
            MPI_Send(propbuff_star.data(),sizeof(float)*num*numextrafields,MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
        }
    }
#endif
#ifdef BHON
    numextrafields = opt.bh_internalprop_names.size() + opt.bh_chem_names.size() + opt.bh_chemproduction_names.size();
    if (numextrafields > 0)
    {
        num = indices_bh.size();
        MPI_Send(&num,sizeof(Int_t),MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
        if (num > 0)
        {
            MPI_Send(indices_bh.data(),sizeof(Int_t)*num,MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
            MPI_Send(propbuff_bh.data(),sizeof(float)*num*numextrafields,MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
        }
    }
#endif
#ifdef EXTRADMON
    numextrafields = opt.extra_dm_internalprop_names.size();
    if (numextrafields > 0)
    {
        num = indices_gas.size();
        MPI_Send(&num,sizeof(Int_t),MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
        if (num > 0)
        {
            MPI_Send(indices_extradm.data(),sizeof(Int_t)*num,MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
            MPI_Send(propbuff_extradm.data(),sizeof(float)*num*numextrafields,MPI_BYTE,taskID,taskID,MPI_COMM_WORLD);
        }
    }
#endif
}

void MPISendHydroInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID)
{
#ifdef GASON
    vector<Int_t> indices;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;

    numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasHydroProperties()) indices.push_back(i);
    num = indices.size();
    MPI_Send(&num,1,MPI_Int_t,taskID,taskID,MPI_COMM_WORLD);
    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetInternalProperties(field);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetChemistry(field);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetChemistryProduction(field);
        }
    }
    MPI_Send(indices.data(),num,MPI_Int_t,taskID,taskID,MPI_COMM_WORLD);
    MPI_Send(propbuff.data(),num*numextrafields,MPI_FLOAT,taskID,taskID,MPI_COMM_WORLD);
#endif
}

void MPISendStarInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID)
{
#ifdef STARON
    vector<Int_t> indices;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;

    numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size()  +opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasStarProperties()) indices.push_back(i);
    num = indices.size();
    MPI_Send(&num,1,MPI_Int_t,taskID,taskID,MPI_COMM_WORLD);
    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetInternalProperties(field);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetChemistry(field);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetChemistryProduction(field);
        }
    }
    MPI_Send(indices.data(),num,MPI_Int_t,taskID,taskID,MPI_COMM_WORLD);
    MPI_Send(propbuff.data(),num*numextrafields,MPI_FLOAT,taskID,taskID,MPI_COMM_WORLD);
#endif
}

void MPISendBHInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID)
{
#ifdef BHON
    vector<Int_t> indices;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;

    numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasBHProperties()) indices.push_back(i);
    num = indices.size();
    MPI_Send(&num,1,MPI_Int_t,taskID,taskID,MPI_COMM_WORLD);
    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetInternalProperties(field);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetChemistry(field);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetChemistryProduction(field);
        }
    }
    MPI_Send(indices.data(),num,MPI_Int_t,taskID,taskID,MPI_COMM_WORLD);
    MPI_Send(propbuff.data(),num*numextrafields,MPI_FLOAT,taskID,taskID,MPI_COMM_WORLD);
#endif
}

void MPISendExtraDMInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int taskID)
{
#ifdef EXTRADMON
    vector<Int_t> indices;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;

    numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasExtraDMProperties()) indices.push_back(i);
    num = indices.size();
    MPI_Send(&num,1,MPI_Int_t,taskID,taskID,MPI_COMM_WORLD);
    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetExtraDMProperties().GetExtraProperties(field);
        }
    }
    MPI_Send(indices.data(),num,MPI_Int_t,taskID,taskID,MPI_COMM_WORLD);
    MPI_Send(propbuff.data(),num*numextrafields,MPI_FLOAT,taskID,taskID,MPI_COMM_WORLD);
#endif
}

void MPIISendHydroInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int dst, int tag, MPI_Request &rqst)
{
#ifdef GASON
    vector<Int_t> indices;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;

    numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasHydroProperties()) indices.push_back(i);
    num = indices.size();
    MPI_Isend(&num, 1, MPI_Int_t, dst, tag, MPI_COMM_WORLD, &rqst);
    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetInternalProperties(field);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetChemistry(field);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetHydroProperties().GetChemistryProduction(field);
        }
    }
    MPI_Isend(indices.data(), num, MPI_Int_t, dst, tag*2, MPI_COMM_WORLD, &rqst);
    MPI_Isend(propbuff.data(), num*numextrafields, MPI_FLOAT, dst, tag*3, MPI_COMM_WORLD, &rqst);
#endif
}

void MPIISendStarInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int dst, int tag, MPI_Request &rqst)
{
#ifdef STARON
    vector<Int_t> indices;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;

    numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasStarProperties()) indices.push_back(i);
    num = indices.size();
    MPI_Isend(&num, 1, MPI_Int_t, dst, tag, MPI_COMM_WORLD, &rqst);
    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetInternalProperties(field);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetChemistry(field);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetStarProperties().GetChemistryProduction(field);
        }
    }
    MPI_Isend(indices.data(), num, MPI_Int_t, dst, tag*2, MPI_COMM_WORLD, &rqst);
    MPI_Isend(propbuff.data(), num*numextrafields, MPI_FLOAT, dst, tag*3, MPI_COMM_WORLD, &rqst);
#endif
}

void MPIISendBHInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int dst, int tag, MPI_Request &rqst)
{
#ifdef BHON
    vector<Int_t> indices;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;

    numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasBHProperties()) indices.push_back(i);
    num = indices.size();
    MPI_Isend(&num, 1, MPI_Int_t, dst, tag, MPI_COMM_WORLD, &rqst);
    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetInternalProperties(field);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetChemistry(field);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetBHProperties().GetChemistryProduction(field);
        }
    }
    MPI_Isend(indices.data(), num, MPI_Int_t, dst, tag*2, MPI_COMM_WORLD, &rqst);
    MPI_Isend(propbuff.data(), num*numextrafields, MPI_FLOAT, dst, tag*3, MPI_COMM_WORLD, &rqst);
#endif
}

void MPIISendExtraDMInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int dst, int tag, MPI_Request &rqst)
{
#ifdef EXTRADMON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num = 0, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;

    numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;
    for (auto i=0;i<nlocalbuff;i++) if (Part[i].HasExtraDMProperties()) indices.push_back(i);
    num = indices.size();
    MPI_Isend(&num, 1, MPI_Int_t, dst, tag, MPI_COMM_WORLD, &rqst);
    if (num == 0) return;
    propbuff.resize(numextrafields*num);
    for (auto i=0;i<num;i++)
    {
        index = indices[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            propbuff[i*numextrafields + iextra + offset] = Part[index].GetExtraDMProperties().GetExtraProperties(field);
        }
    }
    MPI_Isend(indices.data(), num, MPI_Int_t, dst, tag*2, MPI_COMM_WORLD, &rqst);
    MPI_Isend(propbuff.data(), num*numextrafields, MPI_FLOAT, dst, tag*3, MPI_COMM_WORLD, &rqst);
#endif
}

void MPIReceiveHydroInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int readtaskID)
{
#ifdef GASON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;
    HydroProperties x;
    numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    MPI_Recv(&num, 1, MPI_Int_t, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    if (num == 0) return;
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocalbuff;i++) Part[i].NullHydroProperties();
    indices.resize(num);
    propbuff.resize(numextrafields*num);
    MPI_Recv(indices.data(), num, MPI_Int_t, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    MPI_Recv(propbuff.data(), num*numextrafields, MPI_FLOAT, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    for (auto i=0;i<num;i++)
    {
        index=indices[i];
        Part[index].SetHydroProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetInternalProperties(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetChemistry(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetChemistryProduction(field,propbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPIReceiveStarInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int readtaskID)
{
#ifdef STARON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;
    StarProperties x;
    numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    MPI_Recv(&num, 1, MPI_Int_t, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    if (num == 0) return;
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocalbuff;i++) Part[i].NullStarProperties();
    indices.resize(num);
    propbuff.resize(numextrafields*num);
    MPI_Recv(indices.data(), num, MPI_Int_t, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    MPI_Recv(propbuff.data(), num*numextrafields, MPI_FLOAT, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    for (auto i=0;i<num;i++)
    {
        index=indices[i];
        Part[index].SetStarProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            Part[index].GetStarProperties().SetInternalProperties(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            Part[index].GetStarProperties().SetChemistry(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            Part[index].GetStarProperties().SetChemistryProduction(field,propbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPIReceiveBHInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int readtaskID)
{
#ifdef BHON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;
    BHProperties x;
    numextrafields = opt.bh_internalprop_unique_input_names.size() +opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    MPI_Recv(&num, 1, MPI_Int_t, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    if (num == 0) return;
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocalbuff;i++) Part[i].NullBHProperties();
    indices.resize(num);
    propbuff.resize(numextrafields*num);
    MPI_Recv(indices.data(), num, MPI_Int_t, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    MPI_Recv(propbuff.data(), num*numextrafields, MPI_FLOAT, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    for (auto i=0;i<num;i++)
    {
        index=indices[i];
        Part[index].SetBHProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            Part[index].GetBHProperties().SetInternalProperties(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            Part[index].GetBHProperties().SetChemistry(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            Part[index].GetBHProperties().SetChemistryProduction(field,propbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPIReceiveExtraDMInfoFromReadThreads(Options &opt, Int_t nlocalbuff, Particle *Part, int readtaskID)
{
#ifdef EXTRADMON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;
    ExtraDMProperties x;
    numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;
    MPI_Recv(&num, 1, MPI_Int_t, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    if (num == 0) return;
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocalbuff;i++) Part[i].NullExtraDMProperties();
    indices.resize(num);
    propbuff.resize(numextrafields*num);
    MPI_Recv(indices.data(), num, MPI_Int_t, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    MPI_Recv(propbuff.data(), num*numextrafields, MPI_FLOAT, readtaskID, ThisTask, MPI_COMM_WORLD, &status);
    for (auto i=0;i<num;i++)
    {
        index=indices[i];
        Part[index].SetExtraDMProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            Part[index].GetExtraDMProperties().SetExtraProperties(field,propbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPIReceiveParticlesFromReadThreads(Options &opt, Particle *&Pbuf, Particle *Part, int *&readtaskID, int *&irecv, int *&mpi_irecvflag, Int_t *&Nlocalthreadbuf, MPI_Request *&mpi_request, Particle *&Pbaryons)
{
    int irecvflag;
    Int_t i,k,Nlocaltotalbuf;
    MPI_Status status;

    //for all threads not reading snapshots, simply receive particles as necessary from all threads involved with reading the data
    //first determine which threads are going to send information to this thread.
    for (i=0;i<opt.nsnapread;i++) if (irecv[i]) {
        mpi_irecvflag[i]=0;
        MPI_Irecv(&Nlocalthreadbuf[i], 1, MPI_Int_t, readtaskID[i], ThisTask+NProcs, MPI_COMM_WORLD, &mpi_request[i]);
    }
    Nlocaltotalbuf=0;
    //non-blocking receives for the number of particles one expects to receive
    do {
        irecvflag=0;
        for (i=0;i<opt.nsnapread;i++) if (irecv[i]) {
            if (mpi_irecvflag[i]==0) {
                //test if a request has been sent for a Recv call by one of the read threads
                MPI_Test(&mpi_request[i], &mpi_irecvflag[i], &status);
                if (mpi_irecvflag[i]) {
                    if (Nlocalthreadbuf[i]>0) {
                        MPI_Recv(&Part[Nlocal],sizeof(Particle)*Nlocalthreadbuf[i],MPI_BYTE,readtaskID[i],ThisTask, MPI_COMM_WORLD,&status);
                        MPIReceiveHydroInfoFromReadThreads(opt, Nlocalthreadbuf[i], &Part[Nlocal], readtaskID[i]);
                        MPIReceiveStarInfoFromReadThreads(opt, Nlocalthreadbuf[i], &Part[Nlocal], readtaskID[i]);
                        MPIReceiveBHInfoFromReadThreads(opt, Nlocalthreadbuf[i], &Part[Nlocal], readtaskID[i]);
                        MPIReceiveExtraDMInfoFromReadThreads(opt, Nlocalthreadbuf[i], &Part[Nlocal], readtaskID[i]);
                        Nlocal+=Nlocalthreadbuf[i];
                        Nlocaltotalbuf+=Nlocalthreadbuf[i];
                        mpi_irecvflag[i]=0;
                        MPI_Irecv(&Nlocalthreadbuf[i], 1, MPI_Int_t, readtaskID[i], ThisTask+NProcs, MPI_COMM_WORLD, &mpi_request[i]);
                    }
                    else {
                        irecv[i]=0;
                    }
                }
            }
        }
        for (i=0;i<opt.nsnapread;i++) irecvflag+=irecv[i];
    } while(irecvflag>0);
    //now that data is local, must adjust data iff a separate baryon search is required.
    if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
        for (i=0;i<Nlocal;i++) {
            k=Part[i].GetType();
            if (!(k==GASTYPE||k==STARTYPE||k==BHTYPE)) Part[i].SetID(0);
            else {
                Nlocalbaryon[0]++;
                if  (k==GASTYPE) {Part[i].SetID(1);Nlocalbaryon[1]++;}
                else if  (k==STARTYPE) {Part[i].SetID(2);Nlocalbaryon[2]++;}
                else if  (k==BHTYPE) {Part[i].SetID(3);Nlocalbaryon[3]++;}
            }
        }
        //sorted so that dark matter particles first, baryons after
        sort(Part, Part + Nlocal, IDCompareVec);
        Nlocal-=Nlocalbaryon[0];
        //index type separated
        for (i=0;i<Nlocal;i++) Part[i].SetID(i);
        for (i=0;i<Nlocalbaryon[0];i++) Part[i+Nlocal].SetID(i+Nlocal);
    }
}

void MPISendReceiveHydroInfoBetweenThreads(Options &opt, Int_t nlocalbuff, Particle *Pbuf, Int_t nlocal, Particle *Part, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef GASON
    auto numextrafields = opt.gas_internalprop_unique_input_names.size()  + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicessend, indicesrecv;
    Int_t numsend, numrecv, index, offset = 0;
    vector<float> propsendbuff, proprecvbuff;
    string field;
    HydroProperties x;

    //first determine what needs to be sent.
    for (auto i=0;i<nlocalbuff;i++) if (Pbuf[i].HasHydroProperties()) indicessend.push_back(i);
    numsend = indicessend.size();
    if (numsend >0) {
        propsendbuff.resize(numextrafields*numsend);
        for (auto i=0;i<numsend;i++)
        {
            index = indicessend[i];
            offset = 0;
            for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
            {
                field = opt.gas_internalprop_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetHydroProperties().GetInternalProperties(field);
            }
            offset += opt.gas_internalprop_unique_input_names.size();
            for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
            {
                field = opt.gas_chem_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetHydroProperties().GetChemistry(field);
            }
            offset += opt.gas_chem_unique_input_names.size();
            for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
            {
                field = opt.gas_chemproduction_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetHydroProperties().GetChemistryProduction(field);
            }
        }
    }
    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    numrecv = indicesrecv.size();

    if (numrecv == 0) return;
    //and then update the local information
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocal;i++) Part[i].NullHydroProperties();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        Part[index].SetHydroProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPISendReceiveStarInfoBetweenThreads(Options &opt, Int_t nlocalbuff, Particle *Pbuf, Int_t nlocal, Particle *Part, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef STARON

    auto numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicessend, indicesrecv;
    Int_t numsend, numrecv, index, offset = 0;
    vector<float> propsendbuff, proprecvbuff;
    string field;
    StarProperties x;

    //first determine what needs to be sent.
    for (auto i=0;i<nlocalbuff;i++) if (Pbuf[i].HasStarProperties()) indicessend.push_back(i);
    numsend = indicessend.size();
    if (numsend >0) {
        propsendbuff.resize(numextrafields*numsend);
        for (auto i=0;i<numsend;i++)
        {
            index = indicessend[i];
            offset = 0;
            for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
            {
                field = opt.star_internalprop_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetStarProperties().GetInternalProperties(field);
            }
            offset += opt.star_internalprop_unique_input_names.size();
            for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
            {
                field = opt.star_chem_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetStarProperties().GetChemistry(field);
            }
            offset += opt.star_chem_unique_input_names.size();
            for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
            {
                field = opt.star_chemproduction_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetStarProperties().GetChemistryProduction(field);
            }
        }
    }

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    numrecv = indicesrecv.size();
    if (numrecv == 0) return;

    //and then update the local information
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocal;i++) Part[i].NullStarProperties();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        Part[index].SetStarProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            Part[index].GetStarProperties().SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            Part[index].GetStarProperties().SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            Part[index].GetStarProperties().SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPISendReceiveBHInfoBetweenThreads(Options &opt, Int_t nlocalbuff, Particle *Pbuf, Int_t nlocal, Particle *Part, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef BHON

    auto numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicessend, indicesrecv;
    Int_t numsend, numrecv, index, offset = 0;
    vector<float> propsendbuff, proprecvbuff;
    string field;
    BHProperties x;

    //first determine what needs to be sent.
    for (auto i=0;i<nlocalbuff;i++) if (Pbuf[i].HasBHProperties()) indicessend.push_back(i);
    numsend = indicessend.size();
    if (numsend >0) {
        propsendbuff.resize(numextrafields*numsend);
        for (auto i=0;i<numsend;i++)
        {
            index = indicessend[i];
            offset = 0;
            for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
            {
                field = opt.bh_internalprop_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetBHProperties().GetInternalProperties(field);
            }
            offset += opt.bh_internalprop_unique_input_names.size();
            for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
            {
                field = opt.bh_chem_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetBHProperties().GetChemistry(field);
            }
            offset += opt.bh_chem_unique_input_names.size();
            for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
            {
                field = opt.bh_chemproduction_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetBHProperties().GetChemistryProduction(field);
            }
        }
    }

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    numrecv = indicesrecv.size();
    if (numrecv == 0) return;
    //and then update the local information
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocal;i++) Part[i].NullBHProperties();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        Part[index].SetBHProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            Part[index].GetBHProperties().SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            Part[index].GetBHProperties().SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            Part[index].GetBHProperties().SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPISendReceiveExtraDMInfoBetweenThreads(Options &opt, Int_t nlocalbuff, Particle *Pbuf, Int_t nlocal, Particle *Part, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef EXTRADMON

    auto numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicessend, indicesrecv;
    Int_t numsend, numrecv, index, offset = 0;
    vector<float> propsendbuff, proprecvbuff;
    string field;
    ExtraDMProperties x;

    //first determine what needs to be sent.
    for (auto i=0;i<nlocalbuff;i++) if (Pbuf[i].HasExtraDMProperties()) indicessend.push_back(i);
    numsend = indicessend.size();
    if (numsend >0) {
        propsendbuff.resize(numextrafields*numsend);
        for (auto i=0;i<numsend;i++)
        {
            index = indicessend[i];
            offset = 0;
            for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
            {
                field = opt.extra_dm_internalprop_unique_input_names[iextra];
                propsendbuff[i*numextrafields + iextra + offset] = Pbuf[index].GetExtraDMProperties().GetExtraProperties(field);
            }
        }
    }

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    numrecv = indicesrecv.size();
    if (numrecv == 0) return;
    //and then update the local information
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocal;i++) Part[i].NullExtraDMProperties();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        Part[index].SetExtraDMProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            Part[index].GetExtraDMProperties().SetExtraProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPISendReceiveBuffWithHydroInfoBetweenThreads(Options &opt, Particle *PartLocal, vector<Int_t> &indicessend,  vector<float> &propsendbuff, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef GASON

    auto numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicesrecv;
    vector<float> proprecvbuff;
    Int_t index, offset;
    string field;
    HydroProperties x;

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    auto numrecv = indicesrecv.size();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            x.SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            x.SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            x.SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        PartLocal[index].SetHydroProperties(x);
    }
    indicessend.clear();
    propsendbuff.clear();
#endif
}

void MPISendReceiveBuffWithStarInfoBetweenThreads(Options &opt, Particle *PartLocal, vector<Int_t> &indicessend,  vector<float> &propsendbuff, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef STARON
    auto numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicesrecv;
    vector<float> proprecvbuff;
    Int_t index, offset;
    string field;
    StarProperties x;

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    auto numrecv = indicesrecv.size();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            x.SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            x.SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            x.SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        PartLocal[index].SetStarProperties(x);
    }
    indicessend.clear();
    propsendbuff.clear();
#endif
}

void MPISendReceiveBuffWithBHInfoBetweenThreads(Options &opt, Particle *PartLocal, vector<Int_t> &indicessend,  vector<float> &propsendbuff, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef BHON
    auto numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicesrecv;
    vector<float> proprecvbuff;
    Int_t index, offset;
    string field;
    BHProperties x;

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    auto numrecv = indicesrecv.size();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            x.SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            x.SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            x.SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        PartLocal[index].SetBHProperties(x);
    }
    indicessend.clear();
    propsendbuff.clear();
#endif
}

void MPISendReceiveBuffWithExtraDMInfoBetweenThreads(Options &opt, Particle *PartLocal, vector<Int_t> &indicessend,  vector<float> &propsendbuff, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef EXTRADMON

    auto numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicesrecv;
    vector<float> proprecvbuff;
    Int_t index, offset;
    string field;
    ExtraDMProperties x;

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    auto numrecv = indicesrecv.size();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            x.SetExtraProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        PartLocal[index].SetExtraDMProperties(x);
    }
    indicessend.clear();
    propsendbuff.clear();
#endif
}

/// @brief send between read threads of input particle data
void MPISendParticlesBetweenReadThreads(Options &opt, Particle *&Pbuf, Particle *Part, Int_t *&nreadoffset, int *&ireadtask, int *&readtaskID, Particle *&Pbaryons, Int_t *&mpi_nsend_baryon)
{
    MPI_Comm mpi_comm_read = MPI_COMM_WORLD;
    if (ireadtask[ThisTask]>=0) {
        LOG(debug)<< "preparing to send to other reading tasks";
        //split the communication into small buffers
        int icycle=0,ibuf;
        //maximum send size
        int maxchunksize=2147483648/opt.nsnapread/sizeof(Particle);
        int nsend,nrecv,nsendchunks,nrecvchunks,numsendrecv;
        int sendTask,recvTask;
        int sendoffset,recvoffset;
        int isendrecv;
        int cursendchunksize,currecvchunksize;
        MPI_Status status;
        for (ibuf=0;ibuf< opt.nsnapread; ibuf++){
            //if there are an even number of read tasks, then communicate such that 0 communcates with N-1, 1<->N-2, etc
            //and moves on to next communication 0<->N-2, 1<->N-3, etc with the communication in chunks
            //first base on read thread position
            sendTask=ireadtask[ThisTask];
            ///map so that 0 <->N-1, 1 <->N-2, etc to start moving to
            recvTask=abs(opt.nsnapread-1-ibuf-sendTask);
            //if have cycled passed zero, then need to adjust recvTask
            if (icycle==1) recvTask=opt.nsnapread-recvTask;

            //now adjust to the actually task ID in the MPI_COMM_WORLD
            sendTask=ThisTask;
            recvTask=readtaskID[recvTask];
            //if ibuf>0 and now at recvTask=0, then next time, cycle
            if (ibuf>0 && recvTask==0) icycle=1;
            //if sendtask!=recvtask, and information needs to be sent, send information
            if (sendTask!=recvTask && (mpi_nsend[ThisTask * NProcs + recvTask] > 0 || mpi_nsend[recvTask * NProcs + ThisTask] > 0)) {
                nsend=mpi_nsend[ThisTask * NProcs + recvTask];
                nrecv=mpi_nsend[recvTask * NProcs + ThisTask];
                //calculate how many send/recvs are needed
                nsendchunks=ceil((double)nsend/(double)maxchunksize);
                nrecvchunks=ceil((double)nrecv/(double)maxchunksize);
                numsendrecv=max(nsendchunks,nrecvchunks);
                //initialize the offset in the particle array
                sendoffset=0;
                recvoffset=0;
                isendrecv=1;
                do
                {
                    //determine amount to be sent
                    cursendchunksize=min(maxchunksize,nsend-sendoffset);
                    currecvchunksize=min(maxchunksize,nrecv-recvoffset);
                    //blocking point-to-point send and receive. Here must determine the appropriate offset point in the local export buffer
                    //for sending data and also the local appropriate offset in the local the receive buffer for information sent from the local receiving buffer
                    MPI_Sendrecv(&Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset],sizeof(Particle)*cursendchunksize, MPI_BYTE, recvTask, TAG_IO_A+isendrecv,
                        &Part[Nlocal],sizeof(Particle)*currecvchunksize, MPI_BYTE, recvTask, TAG_IO_A+isendrecv,
                                MPI_COMM_WORLD, &status);
                    MPISendReceiveHydroInfoBetweenThreads(opt, cursendchunksize,  &Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset], currecvchunksize, &Part[Nlocal], recvTask, TAG_IO_A+isendrecv, mpi_comm_read);
                    MPISendReceiveStarInfoBetweenThreads(opt, cursendchunksize,  &Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset], currecvchunksize, &Part[Nlocal], recvTask, TAG_IO_A+isendrecv, mpi_comm_read);
                    MPISendReceiveBHInfoBetweenThreads(opt, cursendchunksize,  &Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset], currecvchunksize, &Part[Nlocal], recvTask, TAG_IO_A+isendrecv, mpi_comm_read);
                    MPISendReceiveExtraDMInfoBetweenThreads(opt, cursendchunksize,  &Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset], currecvchunksize, &Part[Nlocal], recvTask, TAG_IO_A+isendrecv, mpi_comm_read);
                    Nlocal+=currecvchunksize;
                    sendoffset+=cursendchunksize;
                    recvoffset+=currecvchunksize;
                    isendrecv++;
                } while (isendrecv<=numsendrecv);
            }
            //if separate baryon search, send baryons too
            if (opt.iBaryonSearch && opt.partsearchtype!=PSTALL) {
                nsend=mpi_nsend_baryon[ThisTask * NProcs + recvTask];
                nrecv=mpi_nsend_baryon[recvTask * NProcs + ThisTask];
                //calculate how many send/recvs are needed
                nsendchunks=ceil((double)nsend/(double)maxchunksize);
                nrecvchunks=ceil((double)nrecv/(double)maxchunksize);
                numsendrecv=max(nsendchunks,nrecvchunks);
                //initialize the offset in the particle array
                sendoffset=0;
                recvoffset=0;
                isendrecv=1;
                do
                {
                    //determine amount to be sent
                    cursendchunksize=min(maxchunksize,nsend-sendoffset);
                    currecvchunksize=min(maxchunksize,nrecv-recvoffset);
                    //blocking point-to-point send and receive. Here must determine the appropriate offset point in the local export buffer
                    //for sending data and also the local appropriate offset in the local the receive buffer for information sent from the local receiving buffer
                    MPI_Sendrecv(&Pbuf[nreadoffset[ireadtask[recvTask]]+mpi_nsend[ThisTask * NProcs + recvTask]+sendoffset],sizeof(Particle)*cursendchunksize, MPI_BYTE, recvTask, TAG_IO_B+isendrecv,
                        &Pbaryons[Nlocalbaryon[0]],sizeof(Particle)*currecvchunksize, MPI_BYTE, recvTask, TAG_IO_B+isendrecv,
                                MPI_COMM_WORLD, &status);
                    MPISendReceiveHydroInfoBetweenThreads(opt, cursendchunksize,  &Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset], currecvchunksize, &Pbaryons[Nlocalbaryon[0]], recvTask, TAG_IO_B+isendrecv, mpi_comm_read);
                    MPISendReceiveStarInfoBetweenThreads(opt, cursendchunksize,  &Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset], currecvchunksize, &Pbaryons[Nlocalbaryon[0]], recvTask, TAG_IO_B+isendrecv, mpi_comm_read);
                    MPISendReceiveBHInfoBetweenThreads(opt, cursendchunksize,  &Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset], currecvchunksize, &Pbaryons[Nlocalbaryon[0]], recvTask, TAG_IO_B+isendrecv, mpi_comm_read);
                    MPISendReceiveExtraDMInfoBetweenThreads(opt, cursendchunksize,  &Pbuf[nreadoffset[ireadtask[recvTask]]+sendoffset], currecvchunksize, &Pbaryons[Nlocalbaryon[0]], recvTask, TAG_IO_B+isendrecv, mpi_comm_read);
                    Nlocalbaryon[0]+=currecvchunksize;
                    sendoffset+=cursendchunksize;
                    recvoffset+=currecvchunksize;
                    isendrecv++;
                } while (isendrecv<=numsendrecv);
            }
        }
    }
}

/// @brief send between read threads of input particle data to vector of buffers
void MPISendParticlesBetweenReadThreads(Options &opt, 
    vector<Particle> *&Preadbuf, Particle *Part, int *&ireadtask, int *&readtaskID, Particle *&Pbaryons, 
    MPI_Comm &mpi_comm_read, Int_t *&mpi_nsend_readthread, Int_t *&mpi_nsend_readthread_baryon)
{
    if (ireadtask[ThisTask]>=0) {
        //split the communication into small buffers
        int icycle=0,ibuf;
        //maximum send size
        int maxchunksize=2147483648/opt.nsnapread/sizeof(Particle);
        int nsend,nrecv,nsendchunks,nrecvchunks,numsendrecv;
        int sendTask,recvTask;
        int sendoffset,recvoffset;
        int isendrecv;
        int cursendchunksize,currecvchunksize;
        MPI_Status status;
        for (ibuf=0;ibuf< opt.nsnapread; ibuf++){
            //if there are an even number of read tasks, then communicate such that 0 communcates with N-1, 1<->N-2, etc
            //and moves on to next communication 0<->N-2, 1<->N-3, etc with the communication in chunks
            //first base on read thread position
            sendTask=ireadtask[ThisTask];
            ///map so that 0 <->N-1, 1 <->N-2, etc to start moving to
            recvTask=abs(opt.nsnapread-1-ibuf-sendTask);
            //if have cycled passed zero, then need to adjust recvTask
            if (icycle==1) recvTask=opt.nsnapread-recvTask;

            //if ibuf>0 and now at recvTask=0, then next time, cycle
            if (ibuf>0 && recvTask==0) icycle=1;
            //if sendtask!=recvtask, and information needs to be sent, send information
            if (sendTask!=recvTask && (mpi_nsend_readthread[sendTask * opt.nsnapread + recvTask] > 0 || mpi_nsend_readthread[recvTask * opt.nsnapread + sendTask] > 0)) {
                nsend=mpi_nsend_readthread[sendTask * opt.nsnapread + recvTask];
                nrecv=mpi_nsend_readthread[recvTask * opt.nsnapread + sendTask];
                //calculate how many send/recvs are needed
                nsendchunks=ceil((double)nsend/(double)maxchunksize);
                nrecvchunks=ceil((double)nrecv/(double)maxchunksize);
                numsendrecv=max(nsendchunks,nrecvchunks);
                //initialize the offset in the particle array
                sendoffset=0;
                recvoffset=0;
                isendrecv=1;
                LOG(trace) <<"sending/receving to/from "<<recvTask<<" [nsend,nrecv] = "<<nsend<<", "<<nrecv<<" in "<<numsendrecv<<" loops";
                do
                {
                    //determine amount to be sent
                    cursendchunksize=min(maxchunksize,nsend-sendoffset);
                    currecvchunksize=min(maxchunksize,nrecv-recvoffset);
                    //blocking point-to-point send and receive. Here must determine the appropriate offset point in the local export buffer
                    //for sending data and also the local appropriate offset in the local the receive buffer for information sent from the local receiving buffer
                    MPI_Sendrecv(&Preadbuf[recvTask][sendoffset],sizeof(Particle)*cursendchunksize, MPI_BYTE, recvTask, TAG_IO_A+isendrecv,
                        &(Part[Nlocal]),sizeof(Particle)*currecvchunksize, MPI_BYTE, recvTask, TAG_IO_A+isendrecv,
                                mpi_comm_read, &status);
                    MPISendReceiveHydroInfoBetweenThreads(opt, cursendchunksize,  &Preadbuf[recvTask][sendoffset], currecvchunksize, &Part[Nlocal], recvTask, TAG_IO_A+isendrecv, mpi_comm_read);
                    MPISendReceiveStarInfoBetweenThreads(opt, cursendchunksize,  &Preadbuf[recvTask][sendoffset], currecvchunksize, &Part[Nlocal], recvTask, TAG_IO_A+isendrecv, mpi_comm_read);
                    MPISendReceiveBHInfoBetweenThreads(opt, cursendchunksize,  &Preadbuf[recvTask][sendoffset], currecvchunksize, &Part[Nlocal], recvTask, TAG_IO_A+isendrecv, mpi_comm_read);
                    MPISendReceiveExtraDMInfoBetweenThreads(opt, cursendchunksize,  &Preadbuf[recvTask][sendoffset], currecvchunksize, &Part[Nlocal], recvTask, TAG_IO_A+isendrecv, mpi_comm_read);
                    Nlocal+=currecvchunksize;
                    sendoffset+=cursendchunksize;
                    recvoffset+=currecvchunksize;
                    isendrecv++;
                } while (isendrecv<=numsendrecv);
            }
            //if separate baryon search, send baryons too
            if (opt.iBaryonSearch && opt.partsearchtype!=PSTALL) {
                nsend=mpi_nsend_readthread_baryon[sendTask * opt.nsnapread + recvTask];
                nrecv=mpi_nsend_readthread_baryon[recvTask * opt.nsnapread + sendTask];
                //calculate how many send/recvs are needed
                nsendchunks=ceil((double)nsend/(double)maxchunksize);
                nrecvchunks=ceil((double)nrecv/(double)maxchunksize);
                numsendrecv=max(nsendchunks,nrecvchunks);
                //initialize the offset in the particle array
                sendoffset=0;
                recvoffset=0;
                isendrecv=1;
                do
                {
                    //determine amount to be sent
                    cursendchunksize=min(maxchunksize,nsend-sendoffset);
                    currecvchunksize=min(maxchunksize,nrecv-recvoffset);
                    //blocking point-to-point send and receive. Here must determine the appropriate offset point in the local export buffer
                    //for sending data and also the local appropriate offset in the local the receive buffer for information sent from the local receiving buffer
                    MPI_Sendrecv(&Preadbuf[recvTask][mpi_nsend_readthread[sendTask * opt.nsnapread + recvTask]+sendoffset],sizeof(Particle)*cursendchunksize, MPI_BYTE, recvTask, TAG_IO_B+isendrecv,
                        &Pbaryons[Nlocalbaryon[0]],sizeof(Particle)*currecvchunksize, MPI_BYTE, recvTask, TAG_IO_B+isendrecv,
                                mpi_comm_read, &status);
                    MPISendReceiveHydroInfoBetweenThreads(opt, cursendchunksize,  &Preadbuf[recvTask][mpi_nsend_readthread[sendTask * opt.nsnapread + recvTask]+sendoffset], currecvchunksize, &Pbaryons[Nlocalbaryon[0]], recvTask, TAG_IO_B+isendrecv, mpi_comm_read);
                    MPISendReceiveStarInfoBetweenThreads(opt, cursendchunksize,  &Preadbuf[recvTask][mpi_nsend_readthread[sendTask * opt.nsnapread + recvTask]+sendoffset], currecvchunksize, &Pbaryons[Nlocalbaryon[0]], recvTask, TAG_IO_B+isendrecv, mpi_comm_read);
                    MPISendReceiveBHInfoBetweenThreads(opt, cursendchunksize,  &Preadbuf[recvTask][mpi_nsend_readthread[sendTask * opt.nsnapread + recvTask]+sendoffset], currecvchunksize, &Pbaryons[Nlocalbaryon[0]], recvTask, TAG_IO_B+isendrecv, mpi_comm_read);
                    MPISendReceiveExtraDMInfoBetweenThreads(opt, cursendchunksize,  &Preadbuf[recvTask][mpi_nsend_readthread[sendTask * opt.nsnapread + recvTask]+sendoffset], currecvchunksize, &Pbaryons[Nlocalbaryon[0]], recvTask, TAG_IO_B+isendrecv, mpi_comm_read);
                    Nlocalbaryon[0]+=currecvchunksize;
                    sendoffset+=cursendchunksize;
                    recvoffset+=currecvchunksize;
                    isendrecv++;
                } while (isendrecv<=numsendrecv);
            }
        }
    }
}

//adds a particle read from an input file to the appropriate buffers
void MPIAddParticletoAppropriateBuffer(Options &opt, const int &ibuf, Int_t ibufindex, int *&ireadtask, const Int_t &BufSize, Int_t *&Nbuf, Particle *&Pbuf, Int_t &numpart, Particle *Part, Int_t *&Nreadbuf, vector<Particle>*&Preadbuf){
    if (ibuf==ThisTask) {
        Nbuf[ibuf]--;
        Part[numpart++]=Pbuf[ibufindex];
    }
    else {
        if(Nbuf[ibuf]==BufSize&&ireadtask[ibuf]<0) {
            MPI_Send(&Nbuf[ibuf], 1, MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
            MPI_Send(&Pbuf[ibuf*BufSize],sizeof(Particle)*Nbuf[ibuf],MPI_BYTE,ibuf,ibuf,MPI_COMM_WORLD);
            MPISendHydroInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
            MPISendStarInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
            MPISendBHInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
            MPISendExtraDMInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
            Nbuf[ibuf]=0;
        }
        else if (ireadtask[ibuf]>=0) {
            if (ibuf!=ThisTask) {
                if (Nreadbuf[ireadtask[ibuf]]==Preadbuf[ireadtask[ibuf]].size()) Preadbuf[ireadtask[ibuf]].resize(Preadbuf[ireadtask[ibuf]].size()+BufSize);
                Preadbuf[ireadtask[ibuf]][Nreadbuf[ireadtask[ibuf]]]=Pbuf[ibufindex];
                Nreadbuf[ireadtask[ibuf]]++;
                Nbuf[ibuf]=0;
            }
        }
    }
}

//@}


/// \name routines which check to see if some search region overlaps with local mpi domain
//@{
///search if some region is in the local mpi domain
int MPIInDomain(Double_t xsearch[3][2], Double_t bnd[3][2]){
    Double_t xsearchp[3][2];
    if (NProcs==1) return 1;
    if (!((bnd[0][1] < xsearch[0][0]) || (bnd[0][0] > xsearch[0][1]) ||
        (bnd[1][1] < xsearch[1][0]) || (bnd[1][0] > xsearch[1][1]) ||
        (bnd[2][1] < xsearch[2][0]) || (bnd[2][0] > xsearch[2][1])))
            return 1;
    else {
        if (mpi_period==0) return 0;
        else {
            for (int j=0;j<3;j++) {xsearchp[j][0]=xsearch[j][0];xsearchp[j][1]=xsearch[j][1];}
            for (int j=0;j<3;j++) {
                if (!((bnd[j][1] < xsearch[j][0]+mpi_period) || (bnd[j][0] > xsearch[j][1]+mpi_period))) {xsearchp[j][0]+=mpi_period;xsearchp[j][1]+=mpi_period;}
                else if (!((bnd[j][1] < xsearch[j][0]-mpi_period) || (bnd[j][0] > xsearch[j][1]-mpi_period))) {xsearchp[j][0]-=mpi_period;xsearchp[j][1]-=mpi_period;}
            }
            if (!((bnd[0][1] < xsearchp[0][0]) || (bnd[0][0] > xsearchp[0][1]) ||
            (bnd[1][1] < xsearchp[1][0]) || (bnd[1][0] > xsearchp[1][1]) ||
            (bnd[2][1] < xsearchp[2][0]) || (bnd[2][0] > xsearchp[2][1])))
                return 1;
            else return 0;
        }
    }
}

///\todo clean up memory allocation in these functions, no need to keep allocating xsearch,xsearchp,numoverlap,etc
/// Determine if a particle needs to be exported to another mpi domain based on a physical search radius
int MPISearchForOverlap(Particle &Part, Double_t &rdist){
    Double_t xsearch[3][2];
    for (auto k=0;k<3;k++) {xsearch[k][0]=Part.GetPosition(k)-rdist;xsearch[k][1]=Part.GetPosition(k)+rdist;}
    return MPISearchForOverlap(xsearch);
}

int MPISearchForOverlap(Coordinate &x, Double_t &rdist){
    Double_t xsearch[3][2];
    for (auto k=0;k<3;k++) {xsearch[k][0]=x[k]-rdist;xsearch[k][1]=x[k]+rdist;}
    return MPISearchForOverlap(xsearch);
}

int MPISearchForOverlap(Double_t xsearch[3][2]){
    Double_t xsearchp[7][3][2];//used to store periodic reflections
    int numoverlap=0,numreflecs=0,ireflec[3],numreflecchoice=0;
    int j,k;

    for (j=0;j<NProcs;j++) {
        if (j!=ThisTask) {
            //determine if search region is not outside of this processors domain
            if(!((mpi_domain[j].bnd[0][1] < xsearch[0][0]) || (mpi_domain[j].bnd[0][0] > xsearch[0][1]) ||
                (mpi_domain[j].bnd[1][1] < xsearch[1][0]) || (mpi_domain[j].bnd[1][0] > xsearch[1][1]) ||
                (mpi_domain[j].bnd[2][1] < xsearch[2][0]) || (mpi_domain[j].bnd[2][0] > xsearch[2][1])))
                numoverlap++;
        }
    }
    if (mpi_period!=0) {
        for (k=0;k<3;k++) if (xsearch[k][0]<0||xsearch[k][1]>mpi_period) ireflec[numreflecs++]=k;
        if (numreflecs==1)numreflecchoice=1;
        else if (numreflecs==2) numreflecchoice=3;
        else if (numreflecs==3) numreflecchoice=7;
        for (j=0;j<numreflecchoice;j++) for (k=0;k<3;k++) {xsearchp[j][k][0]=xsearch[k][0];xsearchp[j][k][1]=xsearch[k][1];}
        if (numreflecs==1) {
            if (xsearch[ireflec[0]][0]<0) {
                    xsearchp[0][ireflec[0]][0]=xsearch[ireflec[0]][0]+mpi_period;
                    xsearchp[0][ireflec[0]][1]=xsearch[ireflec[0]][1]+mpi_period;
            }
            else if (xsearch[ireflec[0]][1]>mpi_period) {
                    xsearchp[0][ireflec[0]][0]=xsearch[ireflec[0]][0]-mpi_period;
                    xsearchp[0][ireflec[0]][1]=xsearch[ireflec[0]][1]-mpi_period;
            }
        }
        else if (numreflecs==2) {
            k=0;j=0;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            k++;j++;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            j++;k=0;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            k++;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
        }
        else if (numreflecs==3) {
            j=0;k=0;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            j++;k++;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            j++;k++;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            j++;k=0;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            k=1;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            j++;k=0;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            k=2;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            j++;k=1;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            k=2;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            j++;k=0;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            k++;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
            k++;
            if (xsearch[ireflec[k]][0]<0) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]+mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]+mpi_period;
            }
            else if (xsearch[ireflec[k]][1]>mpi_period) {
                    xsearchp[j][ireflec[k]][0]=xsearch[ireflec[k]][0]-mpi_period;
                    xsearchp[j][ireflec[k]][1]=xsearch[ireflec[k]][1]-mpi_period;
            }
        }
        for (j=0;j<NProcs;j++) for (k=0;k<numreflecchoice;k++) {
            if (j!=ThisTask) {
                if(!((mpi_domain[j].bnd[0][1] < xsearchp[k][0][0]) || (mpi_domain[j].bnd[0][0] > xsearchp[k][0][1]) ||
                (mpi_domain[j].bnd[1][1] < xsearchp[k][1][0]) || (mpi_domain[j].bnd[1][0] > xsearchp[k][1][1]) ||
                (mpi_domain[j].bnd[2][1] < xsearchp[k][2][0]) || (mpi_domain[j].bnd[2][0] > xsearchp[k][2][1])))
                numoverlap++;
            }
        }
    }
    return numoverlap;
}

///\todo clean up memory allocation in these functions, no need to keep allocating xsearch,xsearchp,numoverlap,etc
/// Determine if a particle needs to be exported to another mpi domain based on a physical search radius

int MPISearchForOverlapUsingMesh(Options &opt, Particle &Part, Double_t &rdist){
    Double_t xsearch[3][2];
    for (auto k=0;k<3;k++) {xsearch[k][0]=Part.GetPosition(k)-rdist;xsearch[k][1]=Part.GetPosition(k)+rdist;}
    return MPISearchForOverlapUsingMesh(opt, xsearch);
}

int MPISearchForOverlapUsingMesh(Options &opt, Coordinate &x, Double_t &rdist){
    Double_t xsearch[3][2];
    for (auto k=0;k<3;k++) {xsearch[k][0]=x[k]-rdist;xsearch[k][1]=x[k]+rdist;}
    return MPISearchForOverlapUsingMesh(opt, xsearch);
}

int MPISearchForOverlapUsingMesh(Options &opt, Double_t xsearch[3][2]){
    int numoverlap=0;
    /// Store whether an MPI domain has already been sent to
    vector<int>sent_mpi_domain(NProcs);
    for(int i=0; i<NProcs; i++) sent_mpi_domain[i] = 0;

    /*vector<int> celllist=MPIGetCellListInSearchUsingMesh(opt,xsearch);
    for (auto j:celllist) {
        const int cellnodeID = opt.cellnodeids[j];
        // Only check if particles have overlap with neighbouring cells that are on another MPI domain and have not already been sent to
        if (sent_mpi_domain[cellnodeID] == 1) continue;
        numoverlap++;
        sent_mpi_domain[cellnodeID]++;
    }
    */
    vector<int> cellnodeidlist = MPIGetCellNodeIDListInSearchUsingMesh(opt,xsearch);
    for (auto cellnodeID:cellnodeidlist) {
        /// Only check if particles have overlap with neighbouring cells that are on another MPI domain and have not already been sent to
        if (sent_mpi_domain[cellnodeID] == 1) continue;
        numoverlap++;
        sent_mpi_domain[cellnodeID]++;
    }
    return numoverlap;
}
//@}

/// \name Routines involved in reading input data
//@{
/// \brief Distribute the mpi processes that read the input files so as to spread the read threads evenly throughout the MPI_COMM_WORLD
void MPIDistributeReadTasks(Options&opt, int *&ireadtask, int*&readtaskID){
    //initialize
    if (opt.nsnapread>NProcs) opt.nsnapread=NProcs;
#ifndef USEPARALLELHDF
    //if not using parallel hdf5, allow only one task per file
    if (opt.num_files<opt.nsnapread) opt.nsnapread=opt.num_files;
#else
    //if parallel hdf5 but not reading hdf then again, max one task per file
    if (opt.inputtype!=IOHDF) if (opt.num_files<opt.nsnapread) opt.nsnapread=opt.num_files;
#endif
    for (int i=0;i<NProcs;i++) ireadtask[i]=-1;
    int spacing=max(1,static_cast<int>(static_cast<float>(NProcs)/static_cast<float>(opt.nsnapread)));
    for (int i=0;i<opt.nsnapread;i++) {ireadtask[i*spacing]=i;readtaskID[i]=i*spacing;}
}
/// \brief Set what tasks read what files
int MPISetFilesRead(Options&opt, std::vector<int> &ireadfile, int *&ireadtask){
    //to determine which files the thread should read
    int nread, niread, nfread;
    for (int i=0;i<opt.num_files;i++) ireadfile[i]=0;
#ifndef USEPARALLELHDF
    nread=static_cast<int>(static_cast<float>(opt.num_files)/static_cast<float>(opt.nsnapread));
    niread=ireadtask[ThisTask]*nread;
    nfread=(ireadtask[ThisTask]+1)*nread;
    if (ireadtask[ThisTask]==opt.nsnapread-1) nfread=opt.num_files;
    for (int i=niread;i<nfread;i++) ireadfile[i]=1;
#else
    //for parallel hdf, multiple tasks can be set to read the same file
    //but if nfiles >= nsnapread, proceed as always.
    if (opt.num_files>=opt.nsnapread) {
        int isel = 0;
        vector<int> readID(opt.nsnapread);
        for (auto i=0;i<NProcs;i++) if (ireadtask[i] > -1 ) readID[isel++]=i;
        isel = 0;
        nread = 0;
        niread = -1; 
        for (auto i=0;i<opt.num_files;i++) {
            if (ThisTask == readID[isel]) {
                ireadfile[i] = 1;
                nread++;
                if (niread == -1) niread = i;
                nfread = i;
            }
            isel++;
            if (isel >=opt.nsnapread) isel=0;
        }
    }
    else {
        int ntaskread = ceil(static_cast<float>(opt.nsnapread)/static_cast<float>(opt.num_files));
        int ifile = floor(static_cast<float>(ireadtask[ThisTask])/static_cast<float>(ntaskread));
        ireadfile[ifile] = 1;
        niread = ifile;
    }
#endif
    return niread;
}
//@}

/// \name MPI file write related routines
//@{
/// @brief Initialize the write communicators
void MPIInitWriteComm(){
    ThisWriteTask = ThisTask;
    ThisWriteComm = ThisTask;
    NProcsWrite = NProcs;
    NWriteComms = NProcs;
    mpi_comm_write = MPI_COMM_WORLD;
}
/// @brief Define the write communicators (what tasks below to what communicators)
void MPIBuildWriteComm(Options &opt){
#ifdef USEPARALLELHDF
    if (opt.mpinprocswritesize > 1) {
        ThisWriteComm = (int)(floor(ThisTask/(float)opt.mpinprocswritesize));
        NWriteComms = (int)(ceil(NProcs/(float)opt.mpinprocswritesize));
        MPI_Comm_split(MPI_COMM_WORLD, ThisWriteComm, ThisTask, &mpi_comm_write);
        MPI_Comm_rank(mpi_comm_write, &ThisWriteTask);
        MPI_Comm_size(mpi_comm_write, &NProcsWrite);
    }
    else{
        ThisWriteTask = ThisTask;
        ThisWriteComm = ThisTask;
        NProcsWrite = NProcs;
        NWriteComms = NProcs;
        mpi_comm_write = MPI_COMM_WORLD;
    }
#endif
}
/// @brief Free any communicators involved in writing data 
void MPIFreeWriteComm(){
    if (mpi_comm_write != MPI_COMM_WORLD) MPI_Comm_free(&mpi_comm_write);
    mpi_comm_write = MPI_COMM_WORLD;
    ThisWriteTask = ThisTask;
    ThisWriteComm = ThisTask;
    NProcsWrite = NProcs;
    NWriteComms = NProcs;
}
//@}

/// \name Routines involved in moving particles between tasks
//@{

void MPIReceiveHydroInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sourceTaskID, int tag)
{
#ifdef GASON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;
    HydroProperties x;
    numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    MPI_Recv(&num, 1, MPI_Int_t, sourceTaskID, tag, MPI_COMM_WORLD,&status);
    if (num == 0) return;
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocalbuff;i++) Part[i].NullHydroProperties();
    indices.resize(num);
    propbuff.resize(numextrafields*num);
    MPI_Recv(indices.data(), num, MPI_Int_t, sourceTaskID, tag*2, MPI_COMM_WORLD, &status);
    MPI_Recv(propbuff.data(), propbuff.size(), MPI_FLOAT, sourceTaskID, tag*3, MPI_COMM_WORLD,&status);
    for (auto i=0;i<num;i++)
    {
        index=indices[i];
        Part[index].SetHydroProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetInternalProperties(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetChemistry(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            Part[index].GetHydroProperties().SetChemistryProduction(field,propbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPIReceiveStarInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sourceTaskID, int tag)
{
#ifdef STARON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;
    StarProperties x;
    numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    MPI_Recv(&num, 1, MPI_Int_t, sourceTaskID, tag, MPI_COMM_WORLD,&status);
    if (num == 0) return;
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocalbuff;i++) Part[i].NullStarProperties();
    indices.resize(num);
    propbuff.resize(numextrafields*num);
    MPI_Recv(indices.data(), num, MPI_Int_t, sourceTaskID, tag*2, MPI_COMM_WORLD, &status);
    MPI_Recv(propbuff.data(), num*numextrafields, MPI_FLOAT, sourceTaskID, tag*3, MPI_COMM_WORLD,&status);
    for (auto i=0;i<num;i++)
    {
        index=indices[i];
        Part[index].SetStarProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            Part[index].GetStarProperties().SetInternalProperties(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            Part[index].GetStarProperties().SetChemistry(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            Part[index].GetStarProperties().SetChemistryProduction(field,propbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPIReceiveBHInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sourceTaskID, int tag)
{
#ifdef BHON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;
    BHProperties x;
    numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;
    MPI_Recv(&num, 1, MPI_Int_t, sourceTaskID, tag, MPI_COMM_WORLD,&status);
    if (num == 0) return;
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocalbuff;i++) Part[i].NullBHProperties();
    indices.resize(num);
    propbuff.resize(numextrafields*num);
    MPI_Recv(indices.data(), num, MPI_Int_t, sourceTaskID, tag*2, MPI_COMM_WORLD, &status);
    MPI_Recv(propbuff.data(), num*numextrafields, MPI_FLOAT, sourceTaskID, tag*3, MPI_COMM_WORLD,&status);
    for (auto i=0;i<num;i++)
    {
        index=indices[i];
        Part[index].SetBHProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            Part[index].GetBHProperties().SetInternalProperties(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            Part[index].GetBHProperties().SetChemistry(field,propbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            Part[index].GetBHProperties().SetChemistryProduction(field,propbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPIReceiveExtraDMInfo(Options &opt, Int_t nlocalbuff, Particle *Part, int sourceTaskID, int tag)
{
#ifdef EXTRADMON
    MPI_Status status;
    vector<Int_t> indices;
    Int_t num, numextrafields = 0, index, offset = 0;
    vector<float> propbuff;
    string field;
    ExtraDMProperties x;
    numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;
    MPI_Recv(&num, 1, MPI_Int_t, sourceTaskID, tag, MPI_COMM_WORLD,&status);
    if (num == 0) return;
    //explicitly NULLing copied information which was done with a BYTE copy
    //The unique pointers will have meaningless info so NULL them (by relasing ownership)
    //and then setting the released pointer to null via in built function.
    for (auto i=0;i<nlocalbuff;i++) Part[i].NullExtraDMProperties();
    indices.resize(num);
    propbuff.resize(numextrafields*num);
    MPI_Recv(indices.data(), num, MPI_Int_t, sourceTaskID, tag*2, MPI_COMM_WORLD, &status);
    MPI_Recv(propbuff.data(), propbuff.size(), MPI_FLOAT, sourceTaskID, tag*3, MPI_COMM_WORLD,&status);
    for (auto i=0;i<num;i++)
    {
        index=indices[i];
        Part[index].SetExtraDMProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            Part[index].GetExtraDMProperties().SetExtraProperties(field,propbuff[i*numextrafields+iextra+offset]);
        }
    }
#endif
}

void MPISendReceiveFOFHydroInfoBetweenThreads(Options &opt, fofid_in *FoFGroupDataLocal, vector<Int_t> &indicessend,  vector<float> &propsendbuff, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef GASON

    auto numextrafields = opt.gas_internalprop_unique_input_names.size() + opt.gas_chem_unique_input_names.size() + opt.gas_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicesrecv;
    vector<float> proprecvbuff;
    Int_t index, offset;
    string field;
    HydroProperties x;
    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    auto numrecv = indicesrecv.size();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        FoFGroupDataLocal[index].p.SetHydroProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.gas_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.gas_internalprop_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetHydroProperties().SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chem_unique_input_names.size();iextra++)
        {
            field = opt.gas_chem_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetHydroProperties().SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.gas_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.gas_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.gas_chemproduction_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetHydroProperties().SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
    }
    indicessend.clear();
    propsendbuff.clear();
#endif
}

void MPISendReceiveFOFStarInfoBetweenThreads(Options &opt, fofid_in *FoFGroupDataLocal, vector<Int_t> &indicessend,  vector<float> &propsendbuff, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef STARON

    auto numextrafields = opt.star_internalprop_unique_input_names.size() + opt.star_chem_unique_input_names.size() + opt.star_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicesrecv;
    vector<float> proprecvbuff;
    Int_t index, offset;
    string field;
    StarProperties x;

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    auto numrecv = indicesrecv.size();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        FoFGroupDataLocal[index].p.SetStarProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.star_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.star_internalprop_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetStarProperties().SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chem_unique_input_names.size();iextra++)
        {
            field = opt.star_chem_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetStarProperties().SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.star_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.star_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.star_chemproduction_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetStarProperties().SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
    }
    indicessend.clear();
    propsendbuff.clear();
#endif
}

void MPISendReceiveFOFBHInfoBetweenThreads(Options &opt, fofid_in *FoFGroupDataLocal, vector<Int_t> &indicessend,  vector<float> &propsendbuff, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef BHON

    auto numextrafields = opt.bh_internalprop_unique_input_names.size() + opt.bh_chem_unique_input_names.size() + opt.bh_chemproduction_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicesrecv;
    vector<float> proprecvbuff;
    Int_t index, offset;
    string field;
    BHProperties x;

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    auto numrecv = indicesrecv.size();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        FoFGroupDataLocal[index].p.SetBHProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.bh_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.bh_internalprop_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetBHProperties().SetInternalProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_internalprop_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chem_unique_input_names.size();iextra++)
        {
            field = opt.bh_chem_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetBHProperties().SetChemistry(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
        offset += opt.bh_chem_unique_input_names.size();
        for (auto iextra=0;iextra<opt.bh_chemproduction_unique_input_names.size();iextra++)
        {
            field = opt.bh_chemproduction_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetBHProperties().SetChemistryProduction(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
    }
    indicessend.clear();
    propsendbuff.clear();
#endif
}

void MPISendReceiveFOFExtraDMInfoBetweenThreads(Options &opt, fofid_in *FoFGroupDataLocal, vector<Int_t> &indicessend,  vector<float> &propsendbuff, int recvTask, int tag, MPI_Comm &mpi_comm)
{
#ifdef EXTRADMON

    auto numextrafields = opt.extra_dm_internalprop_unique_input_names.size();
    if (numextrafields == 0) return;

    vector<Int_t> indicesrecv;
    vector<float> proprecvbuff;
    Int_t index, offset;
    string field;
    ExtraDMProperties x;

    std::tie(indicesrecv, proprecvbuff) =
        exchange_indices_and_props(indicessend, propsendbuff, numextrafields,
            recvTask, tag, mpi_comm);
    auto numrecv = indicesrecv.size();
    for (auto i=0;i<numrecv;i++)
    {
        index=indicesrecv[i];
        FoFGroupDataLocal[index].p.SetExtraDMProperties(x);
        offset = 0;
        for (auto iextra=0;iextra<opt.extra_dm_internalprop_unique_input_names.size();iextra++)
        {
            field = opt.extra_dm_internalprop_unique_input_names[iextra];
            FoFGroupDataLocal[index].p.GetExtraDMProperties().SetExtraProperties(field,proprecvbuff[i*numextrafields+iextra+offset]);
        }
    }
    indicessend.clear();
    propsendbuff.clear();
#endif
}

void MPIGetExportNum(const Int_t nbodies, Particle *Part, Double_t rdist){
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs];
    Double_t xsearch[3][2];

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=0;i<nbodies;i++) {
        for (int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist;xsearch[k][1]=Part[i].GetPosition(k)+rdist;}
        for (j=0;j<NProcs;j++) {
            if (j!=ThisTask) {
                //determine if search region is not outside of this processors domain
                if(MPIInDomain(xsearch,mpi_domain[j].bnd))
                {
                    nexport++;
                    nsend_local[j]++;
                }
            }
        }
    }
    NExport=nexport;//*(1.0+MPIExportFac);
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];
    NImport = nimport;
}

void MPIGetExportNumUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t rdist){
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs];
    Double_t xsearch[3][2];
    //siminfo *s = &opt.swiftsiminfo;
    //Options *s = &opt;

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;

    LOG(info) << "Finding number of particles to export to other MPI domains...";

    vector<int>sent_mpi_domain(NProcs);

    for (i=0;i<nbodies;i++) {
        for(int k=0; k<NProcs; k++) sent_mpi_domain[k] = 0;
        for(int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist;xsearch[k][1]=Part[i].GetPosition(k)+rdist;}
        vector<int> cellnodeidlist=MPIGetCellNodeIDListInSearchUsingMesh(opt,xsearch);
        for (auto cellnodeID:cellnodeidlist) {
            /// Only check if particles have overlap with neighbouring cells that are on another MPI domain and have not already been sent to
            if (sent_mpi_domain[cellnodeID] == 1) continue;
            nexport++;
            nsend_local[cellnodeID]++;
            sent_mpi_domain[cellnodeID]++;
        }
    }
    NExport=nexport;//*(1.0+MPIExportFac);
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];
    NImport = nimport;
}

/// @brief Determine which particles have a spatial linking length 
/// such that linking overlaps the domain of another processor 
/// store the necessary information to send that data and then send that information
void MPIBuildParticleExportList(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Double_t rdist){
    Int_t i, j,nexport=0,nimport=0;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    Double_t xsearch[3][2];
    int sendTask, recvTask;
    int maxchunksize=LOCAL_MAX_MSGSIZE/NProcs/sizeof(Particle);
    Int_t nsend, nrecv;
    int nsendchunks,nrecvchunks,numsendrecv;
    Int_t sendoffset,recvoffset;
    int isendrecv;
    int cursendchunksize,currecvchunksize;
    MPI_Status status;
    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=0;i<nbodies;i++) {
        for (int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist;xsearch[k][1]=Part[i].GetPosition(k)+rdist;}
        for (j=0;j<NProcs;j++) {
            if (j!=ThisTask) {
                //determine if search region is not outside of this processors domain
                if(MPIInDomain(xsearch,mpi_domain[j].bnd))
                {
                    FoFDataIn[nexport].Index = i;
                    FoFDataIn[nexport].Task = j;
                    FoFDataIn[nexport].iGroup = pfof[Part[i].GetID()];//set group id
                    FoFDataIn[nexport].iGroupTask = ThisTask;//and the task of the group
                    FoFDataIn[nexport].iLen = Len[i];
                    nexport++;
                    nsend_local[j]++;
                }
            }
        }
    }
    if (nexport>0) {
        //sort the export data such that all particles to be passed to thread j are together in ascending thread number
        //qsort(FoFDataIn, nexport, sizeof(struct fofdata_in), fof_export_cmp);
        std::sort(FoFDataIn, FoFDataIn + nexport, fof_export_cmp_vec);
        for (i=0;i<nexport;i++) {
            PartDataIn[i] = Part[FoFDataIn[i].Index];
#ifdef GASON
            PartDataIn[i].SetHydroProperties();
#endif
#ifdef STARON
            PartDataIn[i].SetStarProperties();
#endif
#ifdef BHON
            PartDataIn[i].SetBHProperties();
#endif
#ifdef EXTRADMON
            PartDataIn[i].SetExtraDMProperties();
#endif
        }
    }
    //then store the offset in the export particle data for the jth Task in order to send data.
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    NImport=0;for (j=0;j<NProcs;j++)NImport+=mpi_nsend[ThisTask+j*NProcs];
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            //blocking point-to-point send and receive. Here must determine the appropriate offset point in the local export buffer
            //for sending data and also the local appropriate offset in the local the receive buffer for information sent from the local receiving buffer
            //first send FOF data and then particle data
            MPI_Sendrecv(&FoFDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct fofdata_in), MPI_BYTE,
                recvTask, TAG_FOF_A,
                &FoFDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct fofdata_in),
                MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, &status);
            MPI_Sendrecv(&PartDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(Particle), MPI_BYTE,
                recvTask, TAG_FOF_B,
                &PartDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(Particle),
                MPI_BYTE, recvTask, TAG_FOF_B, MPI_COMM_WORLD, &status);
            MPISendReceiveHydroInfoBetweenThreads(opt, 
                cursendchunksize, &PartDataIn[noffset[recvTask]+sendoffset], 
                currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], 
                recvTask, TAG_FOF_B_HYDRO, mpi_comm);
            MPISendReceiveStarInfoBetweenThreads(opt, 
                cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], 
                currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], 
                recvTask, TAG_FOF_B_STAR, mpi_comm);
            MPISendReceiveBHInfoBetweenThreads(opt, 
                cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], 
                currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], 
                recvTask, TAG_FOF_B_BH, mpi_comm);
            MPISendReceiveExtraDMInfoBetweenThreads(opt, 
                cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], 
                currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], 
                recvTask, TAG_FOF_B_EXTRA_DM, mpi_comm);

            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }

}

/// @brief Similar to @ref MPIBuildParticleExportList but uses mesh to determine when mpi's to search
void MPIBuildParticleExportListUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Double_t rdist){
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    Double_t xsearch[3][2];
    Int_t sendTask,recvTask;
    int maxchunksize=LOCAL_MAX_MSGSIZE/NProcs/sizeof(Particle);
    MPI_Status status;
    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    vector<int>sent_mpi_domain(NProcs);

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;

    LOG(info) << "Now building exported particle list for FOF search ";
    for (i=0;i<nbodies;i++) {
        for(int k=0; k<NProcs; k++) sent_mpi_domain[k] = 0;
        for(int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist;xsearch[k][1]=Part[i].GetPosition(k)+rdist;}
        vector<int> cellnodeidlist=MPIGetCellNodeIDListInSearchUsingMesh(opt,xsearch);
        for (auto cellnodeID:cellnodeidlist) {
            /// Only check if particles have overlap with neighbouring cells that are on another MPI domain and have not already been sent to
            if (sent_mpi_domain[cellnodeID] == 1) continue;
            FoFDataIn[nexport].Index = i;
            FoFDataIn[nexport].Task = cellnodeID;
            FoFDataIn[nexport].iGroup = pfof[Part[i].GetID()];//set group id
            FoFDataIn[nexport].iGroupTask = ThisTask;//and the task of the group
            FoFDataIn[nexport].iLen = Len[i];
            nexport++;
            nsend_local[cellnodeID]++;
            sent_mpi_domain[cellnodeID]++;
        }
    }

    if (nexport>0) {
        //sort the export data such that all particles to be passed to thread j are together in ascending thread number
        // qsort(FoFDataIn, nexport, sizeof(struct fofdata_in), fof_export_cmp);
        std::sort(FoFDataIn, FoFDataIn + nexport, fof_export_cmp_vec);
        for (i=0;i<nexport;i++) {
            PartDataIn[i] = Part[FoFDataIn[i].Index];
#ifdef GASON
            PartDataIn[i].SetHydroProperties();
#endif
#ifdef STARON
            PartDataIn[i].SetStarProperties();
#endif
#ifdef BHON
            PartDataIn[i].SetBHProperties();
#endif
#ifdef EXTRADMON
            PartDataIn[i].SetExtraDMProperties();
#endif
        }
    }
    //then store the offset in the export particle data for the jth Task in order to send data.
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    NImport=0;for (j=0;j<NProcs;j++)NImport+=mpi_nsend[ThisTask+j*NProcs];
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            //blocking point-to-point send and receive. Here must determine the appropriate offset point in the local export buffer
            //for sending data and also the local appropriate offset in the local the receive buffer for information sent from the local receiving buffer
            //first send FOF data and then particle data
            MPI_Sendrecv(&FoFDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct fofdata_in), MPI_BYTE,
                recvTask, TAG_FOF_A,
                &FoFDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct fofdata_in),
                MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, &status);
            MPI_Sendrecv(&PartDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(Particle), MPI_BYTE,
                recvTask, TAG_FOF_B,
                &PartDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(Particle),
                MPI_BYTE, recvTask, TAG_FOF_B, MPI_COMM_WORLD, &status);
            MPISendReceiveHydroInfoBetweenThreads(opt, 
                cursendchunksize, &PartDataIn[noffset[recvTask]+sendoffset], 
                currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], 
                recvTask, TAG_FOF_B_HYDRO, mpi_comm);
            MPISendReceiveStarInfoBetweenThreads(opt, 
                cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], 
                currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], 
                recvTask, TAG_FOF_B_STAR, mpi_comm);
            MPISendReceiveBHInfoBetweenThreads(opt, 
                cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], 
                currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], 
                recvTask, TAG_FOF_B_BH, mpi_comm);
            MPISendReceiveExtraDMInfoBetweenThreads(opt, 
                cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], 
                currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], 
                recvTask, TAG_FOF_B_EXTRA_DM, mpi_comm);
                
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }

}

/// @brief like @ref MPIGetExportNum but number based on NN search, useful for reducing mem at the expense of cpu cycles
void MPIGetNNExportNum(const Int_t nbodies, Particle *Part, Double_t *rdist){
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs];
    Double_t xsearch[3][2];

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=0;i<nbodies;i++)
    {
#ifdef STRUCDEN
        if (Part[i].GetType()<=0) continue;
#endif
        if (rdist[i] == 0) continue;
        for (int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist[i];xsearch[k][1]=Part[i].GetPosition(k)+rdist[i];}
        for (j=0;j<NProcs;j++) {
            if (j!=ThisTask) {
                //determine if search region is not outside of this processors domain
                if(MPIInDomain(xsearch,mpi_domain[j].bnd))
                {
                    nexport++;
                    nsend_local[j]++;
                }
            }
        }
    }
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    NExport=nexport;
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];
    NImport=nimport;
}

/// @brief like @ref MPIGetExportNum but number based on NN search, useful for reducing mem at the expense of cpu cycles
void MPIGetNNExportNumUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t *rdist){
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs];
    Double_t xsearch[3][2];
    vector<int>sent_mpi_domain(NProcs);

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=0;i<nbodies;i++)
    {
#ifdef STRUCDEN
        if (Part[i].GetType() <= 0) continue;
#endif
        if (rdist[i] == 0) continue;
        for(int k=0; k<NProcs; k++) sent_mpi_domain[k] = 0;
        for (int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist[i];xsearch[k][1]=Part[i].GetPosition(k)+rdist[i];}
        vector<int> cellnodeidlist=MPIGetCellNodeIDListInSearchUsingMesh(opt,xsearch);
        for (auto cellnodeID:cellnodeidlist) {
            /// Only check if particles have overlap with neighbouring cells that are on another MPI domain and have not already been sent to
            if (sent_mpi_domain[cellnodeID] == 1) continue;
            nexport++;
            nsend_local[cellnodeID]++;
            sent_mpi_domain[cellnodeID]++;
        }
    }
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    NExport=nexport;
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];
    NImport=nimport;
}

/// @brief like @ref MPIBuildParticleExportList but each particle has a different distance stored in rdist, used to find nearest neighbours
void MPIBuildParticleNNExportList(const Int_t nbodies, Particle *Part, Double_t *rdist){
    Int_t i, j,nexport=0,nimport=0;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    Double_t xsearch[3][2];
    MPI_Status status;
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(nndata_in);

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=0;i<nbodies;i++)
    {
#ifdef STRUCDEN
        if (Part[i].GetType()<=0) continue;
#endif
        if (rdist[i] == 0) continue;
        for (int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist[i];xsearch[k][1]=Part[i].GetPosition(k)+rdist[i];}
        for (j=0;j<NProcs;j++) {
            if (j!=ThisTask) {
                //determine if search region is not outside of this processors domain
                if(MPIInDomain(xsearch,mpi_domain[j].bnd))
                {
                    //NNDataIn[nexport].Index=i;
                    NNDataIn[nexport].ToTask=j;
                    NNDataIn[nexport].FromTask=ThisTask;
                    NNDataIn[nexport].R2=rdist[i]*rdist[i];
                    //NNDataIn[nexport].V2=vdist2[i];
                    for (int k=0;k<3;k++) {
                        NNDataIn[nexport].Pos[k]=Part[i].GetPosition(k);
                        NNDataIn[nexport].Vel[k]=Part[i].GetVelocity(k);
                    }
                    nexport++;
                    nsend_local[j]++;
                }
            }
        }
    }
    //sort the export data such that all particles to be passed to thread j are together in ascending thread number
    //if (nexport>0) qsort(NNDataIn, nexport, sizeof(struct nndata_in), nn_export_cmp);
    if (nexport>0) {
	       std::sort(NNDataIn, NNDataIn + nexport, [](const nndata_in &a, const nndata_in &b) {
            return a.ToTask < b.ToTask;
        });
     }

    //then store the offset in the export particle data for the jth Task in order to send data.
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            MPI_Sendrecv(
                &NNDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct nndata_in), MPI_BYTE,
                recvTask, TAG_NN_A+ichunk,
                &NNDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct nndata_in), MPI_BYTE, 
                recvTask, TAG_NN_A+ichunk, 
                MPI_COMM_WORLD, &status);
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }

}

/// @brief like @ref MPIBuildParticleNNExportList but use mesh to determine overlap
void MPIBuildParticleNNExportListUsingMesh(Options &opt, const Int_t nbodies, Particle *Part, Double_t *rdist){
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    Double_t xsearch[3][2];
    Int_t sendTask,recvTask;
    MPI_Status status;
    vector<int>sent_mpi_domain(NProcs);
    int maxchunksize=2147483648/NProcs/sizeof(nndata_in);

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=0;i<nbodies;i++)
    {
#ifdef STRUCDEN
        if (Part[i].GetType()<=0) continue;
#endif
        if (rdist[i] == 0) continue;
        for (int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist[i];xsearch[k][1]=Part[i].GetPosition(k)+rdist[i];}

        /// Store whether an MPI domain has already been sent to
        for(int k=0; k<NProcs; k++) sent_mpi_domain[k] = 0;
        for (int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist[i];xsearch[k][1]=Part[i].GetPosition(k)+rdist[i];}
        vector<int> cellnodeidlist=MPIGetCellNodeIDListInSearchUsingMesh(opt,xsearch);
        for (auto cellnodeID:cellnodeidlist) {
            /// Only check if particles have overlap with neighbouring cells that are on another MPI domain and have not already been sent to
            if (sent_mpi_domain[cellnodeID] == 1) continue;
            //NNDataIn[nexport].Index=i;
            NNDataIn[nexport].ToTask=cellnodeID;
            NNDataIn[nexport].FromTask=ThisTask;
            NNDataIn[nexport].R2=rdist[i]*rdist[i];
            //NNDataIn[nexport].V2=vdist2[i];
            for (int k=0;k<3;k++) {
                NNDataIn[nexport].Pos[k]=Part[i].GetPosition(k);
                NNDataIn[nexport].Vel[k]=Part[i].GetVelocity(k);
            }
            nexport++;
            nsend_local[cellnodeID]++;
            sent_mpi_domain[cellnodeID]++;
        }
    }
    //sort the export data such that all particles to be passed to thread j are together in ascending thread number
    if (nexport>0) {
	       std::sort(NNDataIn, NNDataIn + nexport, [](const nndata_in &a, const nndata_in &b) {
            return a.ToTask < b.ToTask;
        });
     }

    //then store the offset in the export particle data for the jth Task in order to send data.
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            MPI_Sendrecv(
                &NNDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct nndata_in), MPI_BYTE,
                recvTask, TAG_NN_A+ichunk,
                &NNDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct nndata_in), MPI_BYTE, 
                recvTask, TAG_NN_A+ichunk, 
                MPI_COMM_WORLD, &status);
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }
}

/// @brief Mirror to @ref MPIGetNNExportNum, use exported particles, run ball search 
/// to find number of all local particles that need to be
/// imported back to exported particle's thread so that a proper NN search can be made.
void MPIGetNNImportNum(const Int_t nbodies, KDTree *tree, Particle *Part, int iallflag){
    Int_t i, j, nexport=0;
    Int_t nsend_local[NProcs],nbuffer[NProcs];
    Int_t oldnsend[NProcs*NProcs];
    bool *iflagged = new bool[nbodies];
    vector<Int_t> taggedindex;
    for(j=0;j<NProcs;j++)
    {
        nbuffer[j]=0;
        for (int k=0;k<j;k++)nbuffer[j]+=mpi_nsend[ThisTask+k*NProcs];//offset on "receiver" end
    }
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (j=0;j<NProcs;j++) {
        if (j==ThisTask) continue;
        if (mpi_nsend[ThisTask+j*NProcs]==0) continue;
        for (i=0;i<nbodies;i++) iflagged[i]=false;
        //search local list and tag all local particles that need to be exported back (or imported) to the exported particles thread
        for (i=nbuffer[j];i<nbuffer[j]+mpi_nsend[ThisTask+j*NProcs];i++) {
            taggedindex = tree->SearchBallPosTagged(NNDataGet[i].Pos, NNDataGet[i].R2);
            if (taggedindex.size()==0) continue;
            for (auto &index:taggedindex) {
                if (iflagged[index]) continue;
                #ifdef STRUCDEN
                if (iallflag==0 && Part[index].GetType()<0) continue;
                #endif
                nexport++;
                nsend_local[j]++;
            }
            for (auto &index:taggedindex) iflagged[index]=true;
        }
    }
    delete[] iflagged;
    //must store old mpi nsend for accessing NNDataGet properly.
    for (j=0;j<NProcs;j++) for (int k=0;k<NProcs;k++) oldnsend[k+j*NProcs]=mpi_nsend[k+j*NProcs];
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    NImport=0;
    for (j=0;j<NProcs;j++)NImport+=mpi_nsend[ThisTask+j*NProcs];
    NExport=nexport;
    for (j=0;j<NProcs;j++) for (int k=0;k<NProcs;k++) mpi_nsend[k+j*NProcs]=oldnsend[k+j*NProcs];
}

/// @brief Mirror to @ref MPIBuildParticleNNExportList, use exported particles, 
/// run ball search to find all local particles that need to be
/// imported back to exported particle's thread so that a proper NN search can be made.
/// Is also used for calculating spherical overdensity quantities, where iSOcalc = true
Int_t MPIBuildParticleNNImportList(Options &opt, const Int_t nbodies, KDTree *tree, Particle *Part, int iallflag, bool iSOcalc){
    Int_t i, j, nexport=0,ncount;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    bool *iflagged = new bool[nbodies];
    vector<Int_t> taggedindex;
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(Particle);
    MPI_Status status;
    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    for(j=0;j<NProcs;j++)
    {
        nbuffer[j]=0;
        for (int k=0;k<j;k++)nbuffer[j]+=mpi_nsend[ThisTask+k*NProcs];//offset on "receiver" end
    }
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (j=0;j<NProcs;j++) {
        if (j==ThisTask) continue;
        if (mpi_nsend[ThisTask+j*NProcs]==0) continue;
        for (i=0;i<nbodies;i++) iflagged[i]=false;
        //search local list and tag all local particles that need to be exported back (or imported) to the exported particles thread
        for (i=nbuffer[j];i<nbuffer[j]+mpi_nsend[ThisTask+j*NProcs];i++) {
            taggedindex = tree->SearchBallPosTagged(NNDataGet[i].Pos, NNDataGet[i].R2);
            if (taggedindex.size()==0) continue;
            for (auto &index:taggedindex) {
                if (iflagged[index]) continue;
                iflagged[index] = true;
#ifdef STRUCDEN
                if (iallflag==0 && Part[index].GetType()<0) continue;
#endif
                PartDataIn[nexport]=Part[index];
                nexport++;
                nsend_local[j]++;
            }
        }
    }
    delete[] iflagged;
    //sort the export data such that all particles to be passed to thread j are together in ascending thread number
    //qsort(NNDataReturn, nexport, sizeof(struct nndata_in), nn_export_cmp);

    //now if there is extra information, strip off the data from the particles to be sent to store in a separate buffer
    //and store it explicitly into a buffer. Here are the buffers
    vector<Int_t> indices_gas_send;
    vector<float> propbuff_gas_send;
    vector<Int_t> indices_star_send;
    vector<float> propbuff_star_send;
    vector<Int_t> indices_bh_send;
    vector<float> propbuff_bh_send;
    vector<Int_t> indices_extra_dm_send;
    vector<float> propbuff_extra_dm_send;

    //if no information stored in the extra fields will be used, then just remove it before sending
    //now if this is not called for SO calculations, assume that particles do not need to be exported
    //with extra information
    if (!iSOcalc) {
        MPIStripExportParticleOfExtraInfo(opt, nexport, PartDataIn);
    }

    //then store the offset in the export particle data for the jth Task in order to send data.
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
#if defined(GASON) || defined(STARON) || defined(BHON) || defined(EXTRADMON)
            if (iSOcalc) {
                MPIFillBuffWithHydroInfo(opt, cursendchunksize, &PartDataIn[noffset[recvTask]+sendoffset], indices_gas_send, propbuff_gas_send, true);
                MPIFillBuffWithStarInfo(opt, cursendchunksize, &PartDataIn[noffset[recvTask]+sendoffset], indices_star_send, propbuff_star_send, true);
                MPIFillBuffWithBHInfo(opt, cursendchunksize, &PartDataIn[noffset[recvTask]+sendoffset], indices_bh_send, propbuff_bh_send, true);
                MPIFillBuffWithExtraDMInfo(opt, cursendchunksize, &PartDataIn[noffset[recvTask]+sendoffset], indices_extra_dm_send, propbuff_extra_dm_send, true);
            }
#endif
            MPI_Sendrecv(&PartDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(Particle), MPI_BYTE,
                recvTask, TAG_NN_B+ichunk,
                &PartDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(Particle), MPI_BYTE, 
                recvTask, TAG_NN_B+ichunk, MPI_COMM_WORLD, &status);
#if defined(GASON) || defined(STARON) || defined(BHON) || defined(EXTRADMON)
            if (iSOcalc) {
                MPISendReceiveBuffWithHydroInfoBetweenThreads(opt, &PartDataGet[nbuffer[recvTask]+recvoffset], indices_gas_send, propbuff_gas_send, recvTask, TAG_NN_B+ichunk, mpi_comm);
                MPISendReceiveBuffWithStarInfoBetweenThreads(opt, &PartDataGet[nbuffer[recvTask]+recvoffset], indices_star_send, propbuff_star_send, recvTask, TAG_NN_B+ichunk, mpi_comm);
                MPISendReceiveBuffWithBHInfoBetweenThreads(opt, &PartDataGet[nbuffer[recvTask]+recvoffset], indices_bh_send, propbuff_bh_send, recvTask, TAG_NN_B+ichunk, mpi_comm);
                MPISendReceiveBuffWithExtraDMInfoBetweenThreads(opt, &PartDataGet[nbuffer[recvTask]+recvoffset], indices_extra_dm_send, propbuff_extra_dm_send, recvTask, TAG_NN_B+ichunk, mpi_comm);
            }
#endif
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }
    ncount=0;for (int k=0;k<NProcs;k++)ncount+=mpi_nsend[ThisTask+k*NProcs];
    return ncount;
}

/// @brief similar @ref MPIGetExportNum but number based on halo properties to see if any
vector<bool> MPIGetHaloSearchExportNum(const Int_t ngroup, PropData *&pdata, vector<Double_t> &rdist)
{
    Int_t i,j,nexport=0,nimport=0;
    Int_t nsend_local[NProcs];
    Double_t xsearch[3][2];
    vector<bool> halooverlap(ngroup+1);


    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=1;i<=ngroup;i++)
    {
        halooverlap[i]=false;
        for (int k=0;k<3;k++) {xsearch[k][0]=pdata[i].gcm[k]-rdist[i];xsearch[k][1]=pdata[i].gcm[k]+rdist[i];}
        for (j=0;j<NProcs;j++) {
            if (j!=ThisTask) {
                //determine if search region is not outside of this processors domain
                if(MPIInDomain(xsearch,mpi_domain[j].bnd))
                {
                    nexport++;
                    nsend_local[j]++;
                    halooverlap[i]=true;
                }
            }
        }
    }
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];
    NImport = nimport;
    NExport=nexport;
    return halooverlap;
}

/// @brief similar @ref MPIGetHaloSearchExportNum but using mesh mpi decomposition
vector<bool> MPIGetHaloSearchExportNumUsingMesh(Options &opt, const Int_t ngroup, PropData *&pdata, vector<Double_t> &rdist)
{
    Int_t nexport=0,nimport=0;
    Int_t nsend_local[NProcs];
    Double_t xsearch[3][2];
    vector<bool> halooverlap(ngroup+1);
    vector<int>sent_mpi_domain(NProcs);

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (auto j=0;j<NProcs;j++) nsend_local[j]=0;
    for (auto i=1;i<=ngroup;i++)
    {
        for (int k=0; k<NProcs; k++) sent_mpi_domain[k] = 0;
        for (int k=0;k<3;k++) {xsearch[k][0]=pdata[i].gcm[k]-rdist[i];xsearch[k][1]=pdata[i].gcm[k]+rdist[i];}
        vector<int> cellnodeidlist=MPIGetCellNodeIDListInSearchUsingMesh(opt,xsearch);
        for (auto cellnodeID:cellnodeidlist) {
            /// Only check if particles have overlap with neighbouring cells that are on another MPI domain and have not already been sent to
            if (sent_mpi_domain[cellnodeID] == 1) continue;
            nexport++;
            nsend_local[cellnodeID]++;
            halooverlap[i]=true;
            sent_mpi_domain[cellnodeID]++;
        }
    }

    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (auto j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];
    NImport = nimport;
    NExport=nexport;
    return halooverlap;
}

/// @brief like @ref MPIBuildParticleExportList but each particle has a different distance stored in rdist used to find nearest neighbours
void MPIBuildHaloSearchExportList(const Int_t ngroup, PropData *&pdata, vector<Double_t> &rdist, vector<bool> &halooverlap)
{
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    Double_t xsearch[3][2];
    MPI_Status status;
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(nndata_in);

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=1;i<=ngroup;i++)
    {
        if (halooverlap[i]==false) continue;
        for (int k=0;k<3;k++) {xsearch[k][0]=pdata[i].gcm[k]-rdist[i];xsearch[k][1]=pdata[i].gcm[k]+rdist[i];}
        for (j=0;j<NProcs;j++) {
            if (j!=ThisTask) {
                //determine if search region is not outside of this processors domain
                if(MPIInDomain(xsearch,mpi_domain[j].bnd))
                {
                    //NNDataIn[nexport].Index=i;
                    NNDataIn[nexport].ToTask=j;
                    NNDataIn[nexport].FromTask=ThisTask;
                    NNDataIn[nexport].R2=rdist[i]*rdist[i];
                    //NNDataIn[nexport].V2=vdist2[i];
                    for (int k=0;k<3;k++) {
                        NNDataIn[nexport].Pos[k]=pdata[i].gcm[k];
                    }
                    nexport++;
                    nsend_local[j]++;
                }
            }
        }
    }
    //sort the export data such that all particles to be passed to thread j are together in ascending thread number
    // if (nexport>0) qsort(NNDataIn, nexport, sizeof(struct nndata_in), nn_export_cmp);
    if (nexport>0) std::sort(NNDataIn, NNDataIn + nexport, nn_export_cmp_vec);

    //then store the offset in the export data for the jth Task in order to send data.
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of items to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            MPI_Sendrecv(&NNDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct nndata_in), MPI_BYTE,
                recvTask, TAG_NN_A+ichunk,
                &NNDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct nndata_in), MPI_BYTE, 
                recvTask, TAG_NN_A+ichunk, 
                MPI_COMM_WORLD, &status);
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }
}

/// @brief similar \ref MPIBuildHaloSearchExportList but for mesh mpi decomposition
void MPIBuildHaloSearchExportListUsingMesh(Options &opt, const Int_t ngroup, PropData *&pdata, vector<Double_t> &rdist, vector<bool> &halooverlap)
{
    Int_t nexport=0,nimport=0;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    Double_t xsearch[3][2];
    MPI_Status status;
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(nndata_in);
    vector<int>sent_mpi_domain(NProcs);

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (auto j=0;j<NProcs;j++) nsend_local[j]=0;
    for (auto i=1;i<=ngroup;i++)
    {
        if (halooverlap[i]==false) continue;
        for (int k=0; k<NProcs; k++) sent_mpi_domain[k] = 0;
        for (int k=0;k<3;k++) {xsearch[k][0]=pdata[i].gcm[k]-rdist[i];xsearch[k][1]=pdata[i].gcm[k]+rdist[i];}
        vector<int> cellnodeidlist=MPIGetCellNodeIDListInSearchUsingMesh(opt,xsearch);
        for (auto cellnodeID:cellnodeidlist) {
            /// Only check if particles have overlap with neighbouring cells that are on another MPI domain and have not already been sent to
            if (sent_mpi_domain[cellnodeID] == 1) continue;
            //NNDataIn[nexport].Index=i;
            NNDataIn[nexport].ToTask=cellnodeID;
            NNDataIn[nexport].FromTask=ThisTask;
            NNDataIn[nexport].R2=rdist[i]*rdist[i];
            //NNDataIn[nexport].V2=vdist2[i];
            for (int k=0;k<3;k++) {
                NNDataIn[nexport].Pos[k]=pdata[i].gcm[k];
            }
            nexport++;
            nsend_local[cellnodeID]++;
            sent_mpi_domain[cellnodeID]++;
        }
    }

    //sort the export data such that all particles to be passed to thread j are together in ascending thread number
    if (nexport>0) std::sort(NNDataIn, NNDataIn + nexport, nn_export_cmp_vec);

    //then store the offset in the export data for the jth Task in order to send data.
    noffset[0] = 0; for(auto j = 1; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of items to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (auto j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            MPI_Sendrecv(&NNDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct nndata_in), MPI_BYTE,
                recvTask, TAG_NN_A+ichunk,
                &NNDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct nndata_in), MPI_BYTE, 
                recvTask, TAG_NN_A+ichunk, 
                MPI_COMM_WORLD, &status);
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }
}

/// @brief Mirror to @ref MPIGetHaloSearchExportNum, use exported positions, 
/// run ball search to find number of all local particles that need to be
/// imported back to exported positions's thread so that a proper search can be made.
void MPIGetHaloSearchImportNum(const Int_t nbodies, KDTree *tree, Particle *Part)
{
    Int_t i, j, nexport=0;
    Int_t nsend_local[NProcs],nbuffer[NProcs];
    Int_t oldnsend[NProcs*NProcs];
    Int_t *nn=new Int_t[nbodies];
    Double_t *nnr2=new Double_t[nbodies];
    for(j=0;j<NProcs;j++)
    {
        nbuffer[j]=0;
        for (int k=0;k<j;k++)nbuffer[j]+=mpi_nsend[ThisTask+k*NProcs];//offset on "receiver" end
    }
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (j=0;j<NProcs;j++) {
        for (i=0;i<nbodies;i++) nn[i]=-1;
        if (j==ThisTask) continue;
        if (mpi_nsend[ThisTask+j*NProcs]==0) continue;
            //search local list and tag all local particles that need to be exported back (or imported) to the exported particles thread
            for (i=nbuffer[j];i<nbuffer[j]+mpi_nsend[ThisTask+j*NProcs];i++) {
                tree->SearchBallPos(NNDataGet[i].Pos, NNDataGet[i].R2, j, nn, nnr2);
            }
            for (i=0;i<nbodies;i++) {
                if (nn[i]!=-1) {
                    nexport++;
                    nsend_local[j]++;
                }
            }

    }
    //must store old mpi nsend for accessing NNDataGet properly.
    for (j=0;j<NProcs;j++) for (int k=0;k<NProcs;k++) oldnsend[k+j*NProcs]=mpi_nsend[k+j*NProcs];
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    NImport=0;
    for (j=0;j<NProcs;j++)NImport+=mpi_nsend[ThisTask+j*NProcs];
    NExport=nexport;
    for (j=0;j<NProcs;j++) for (int k=0;k<NProcs;k++) mpi_nsend[k+j*NProcs]=oldnsend[k+j*NProcs];
    delete[] nn;
    delete[] nnr2;
}

/// @brief Mirror to \ref MPIBuildHaloSearchExportList, use exported particles, 
/// run ball search to find all local particles that need to be
/// imported back to exported particle's thread so that a proper NN search can be made.
Int_t MPIBuildHaloSearchImportList(Options &opt, const Int_t nbodies, KDTree *tree, Particle *Part){
    Int_t i, j, nexport=0,ncount;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    Int_t *nn=new Int_t[nbodies];
    Double_t *nnr2=new Double_t[nbodies];
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(Particle);
    MPI_Status status;
    for(j=0;j<NProcs;j++)
    {
        nbuffer[j]=0;
        for (int k=0;k<j;k++)nbuffer[j]+=mpi_nsend[ThisTask+k*NProcs];//offset on "receiver" end
    }

    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (j=0;j<NProcs;j++) {
        for (i=0;i<nbodies;i++) nn[i]=-1;
        if (j==ThisTask) continue;
        if (mpi_nsend[ThisTask+j*NProcs]==0) continue;
        //search local list and tag all local particles that need to be exported back (or imported) to the exported particles thread
        for (i=nbuffer[j];i<nbuffer[j]+mpi_nsend[ThisTask+j*NProcs];i++) {
            tree->SearchBallPos(NNDataGet[i].Pos, NNDataGet[i].R2, j, nn, nnr2);
        }
        for (i=0;i<nbodies;i++) {
            if (nn[i]==-1) continue;
            for (int k=0;k<3;k++) {
                PartDataIn[nexport].SetPosition(k,Part[i].GetPosition(k));
                PartDataIn[nexport].SetVelocity(k,Part[i].GetVelocity(k));
            }
            nexport++;
            nsend_local[j]++;
        }
    }
    //sort the export data such that all particles to be passed to thread j are together in ascending thread number
    //qsort(NNDataReturn, nexport, sizeof(struct nndata_in), nn_export_cmp);

    //then store the offset in the export particle data for the jth Task in order to send data.
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            MPI_Sendrecv(&PartDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(Particle), MPI_BYTE,
                recvTask, TAG_NN_B+ichunk,
                &PartDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(Particle), MPI_BYTE, 
                recvTask, TAG_NN_B+ichunk, 
                MPI_COMM_WORLD, &status);
            // MPISendReceiveHydroInfoBetweenThreads(opt, cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], recvTask, TAG_FOF_B_HYDRO, mpi_comm);
            // MPISendReceiveStarInfoBetweenThreads(opt, cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], recvTask, TAG_FOF_B_STAR, mpi_comm);
            // MPISendReceiveBHInfoBetweenThreads(opt, cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], recvTask, TAG_FOF_B_BH, mpi_comm);
            // MPISendReceiveExtraDMInfoBetweenThreads(opt, cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], recvTask, TAG_FOF_B_EXTRA_DM, mpi_comm);
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }

    ncount=0;for (int k=0;k<NProcs;k++)ncount+=mpi_nsend[ThisTask+k*NProcs];
    delete[] nn;
    delete[] nnr2;
    return ncount;
}


/// @brief Similar to \ref MPIBuildParticleExportList, however this is for associated baryon search 
/// where particles have been moved from original mpi domains and their group id accessed through 
/// the id array and their stored id and length in numingroup
void MPIBuildParticleExportBaryonSearchList(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t *ids, Int_t *numingroup, Double_t rdist){
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    Double_t xsearch[3][2];
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(fofdata_in);
    MPI_Status status;
    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    ///\todo would like to add openmp to this code. In particular, loop over nbodies but issue is nexport.
    ///This would either require making a FoFDataIn[nthreads][NExport] structure so that each omp thread
    ///can only access the appropriate memory and adjust nsend_local.\n
    ///\em Or outer loop is over threads, inner loop over nbodies and just have a idlist of size Nlocal that tags particles
    ///which must be exported. Then its a much quicker follow up loop (no if statement) that stores the data
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=0;i<nbodies;i++) {
        for (int k=0;k<3;k++) {xsearch[k][0]=Part[i].GetPosition(k)-rdist;xsearch[k][1]=Part[i].GetPosition(k)+rdist;}
        for (j=0;j<NProcs;j++) {
            if (j!=ThisTask) {
                //determine if search region is not outside of this processors domain
                if(MPIInDomain(xsearch,mpi_domain[j].bnd))
                {
                    //FoFDataIn[nexport].Part=Part[i];
                    FoFDataIn[nexport].Index = i;
                    FoFDataIn[nexport].Task = j;
                    FoFDataIn[nexport].iGroup = pfof[ids[Part[i].GetID()]];//set group id
                    FoFDataIn[nexport].iGroupTask = ThisTask;//and the task of the group
                    FoFDataIn[nexport].iLen = numingroup[pfof[ids[Part[i].GetID()]]];
                    nexport++;
                    nsend_local[j]++;
                }
            }
        }
    }
    if (nexport>0) {
        //sort the export data such that all particles to be passed to thread j are together in ascending thread number
        std::sort(FoFDataIn, FoFDataIn + nexport, fof_export_cmp_vec);
        for (i=0;i<nexport;i++) PartDataIn[i] = Part[FoFDataIn[i].Index];
    }
    //then store the offset in the export particle data for the jth Task in order to send data.
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    //and then gather the number of particles to be sent from mpi thread m to mpi thread n in the mpi_nsend[NProcs*NProcs] array via [n+m*NProcs]
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    NImport=0;for (j=0;j<NProcs;j++)NImport+=mpi_nsend[ThisTask+j*NProcs];
    //now send the data.
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];

    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            MPI_Sendrecv(&FoFDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct fofdata_in), MPI_BYTE,
                recvTask, TAG_FOF_A+ichunk,
                &FoFDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct fofdata_in), MPI_BYTE, 
                recvTask, TAG_FOF_A+ichunk, 
                MPI_COMM_WORLD, &status);
            MPI_Sendrecv(&PartDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(Particle), MPI_BYTE,
                recvTask, TAG_FOF_B+ichunk,
                &PartDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(Particle), MPI_BYTE, 
                recvTask, TAG_FOF_B+ichunk, 
                MPI_COMM_WORLD, &status);
            MPISendReceiveHydroInfoBetweenThreads(opt, cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], recvTask, TAG_FOF_B_HYDRO, mpi_comm);
            MPISendReceiveStarInfoBetweenThreads(opt, cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], recvTask, TAG_FOF_B_STAR, mpi_comm);
            MPISendReceiveBHInfoBetweenThreads(opt, cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], recvTask, TAG_FOF_B_BH, mpi_comm);
            MPISendReceiveExtraDMInfoBetweenThreads(opt, cursendchunksize,  &PartDataIn[noffset[recvTask]+sendoffset], currecvchunksize, &PartDataGet[nbuffer[recvTask]+recvoffset], recvTask, TAG_FOF_B_EXTRA_DM, mpi_comm);
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }
}

//@}

/// \name FOF related mpi routines
//@{
///Set fof task id of particle
short_mpi_t* MPISetTaskID(const Int_t nbodies){
    short_mpi_t *foftask=new short_mpi_t[nbodies];
    for (Int_t i=0;i<nbodies;i++) foftask[i]=ThisTask;
    return foftask;
}

///Offset pfof array based so that local group numbers do not overlap
///\todo alter to that this now ranks threads so that group ids are larger if thread has more particles. This ensures that mpi threads
///send particles to thread with the fewest particles when linking across mpi domains during a FOF search
void MPIAdjustLocalGroupIDs(const Int_t nbodies, Int_t *pfof){
    PriorityQueue *pq=new PriorityQueue(NProcs);
    for (int j=0;j<NProcs;j++) pq->Push(j,mpi_nlocal[j]);
    Int_t rankorder[NProcs];
    for (int j=0;j<NProcs;j++) {rankorder[NProcs-1-j]=pq->TopQueue();pq->Pop();}
    Int_t offset=0;
    for (int j=0;j<NProcs;j++){if(rankorder[j]==ThisTask) break; offset+=mpi_nlocal[rankorder[j]];}
    //Int_t offset=nbodies*ThisTask;
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) pfof[i]+=offset;
    mpi_maxgid=0;
    for (int j=0;j<NProcs;j++){mpi_maxgid+=mpi_nlocal[rankorder[j]];}
    mpi_gidoffset=0;
    for (int j=0;j<NProcs;j++){if(rankorder[j]==ThisTask) break; mpi_gidoffset+=mpi_ngroups[rankorder[j]];}
    delete pq;
}

//My idea is this for doing the stiching. First generate export list of particle data and another seperate data structure for the FOF data
//Next, when examining local search using export particles (since iGroup=0 is unlinked) if a export particle current iGroup is larger
//Then adjust the local particle and all members of its group so long as its group is NOT group zero. Calculate the number of new links
//and determine the total number of new links across all mpi threads.
//If that number is not zero, then groups have been found that are across processor domains.
//One has to iterate the check across the domains till no more new links have been found. That is must update the export particles Group ids
//then begin the check anew.

//Couple of key things to think about are, one I really shouldn't have to run the check again to find the particles that meet the conditions across
//threads since that has NOT changed. must figure out a way to store relevant particles. Otherwise, continuously checking, seems a waste of cpu cycles.
//second I must pass along head, tail, next, and length information, maybe by using a plist structure so that it is easy to alter the particles locally to new group id
//Also must determine optimal way of setting which processor the group should end up on. Best way might be to use the length of the group locally since
//that would minimize the broadcasts.

/*! Particles that have been marked for export may have had their fof information updated so need to update this info
*/
void MPIUpdateExportList(const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_tree_t *&Len){
    Int_t i, j, nexport;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(fofdata_in);
    int nsendchunks,nrecvchunks,numsendrecv;
    int sendoffset,recvoffset;
    int cursendchunksize,currecvchunksize;
    MPI_Status status;

    nexport=0;
    for (j=0;j<NProcs;j++) {nexport+=mpi_nsend[j+ThisTask*NProcs];nsend_local[j]=mpi_nsend[j+ThisTask*NProcs];}
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    for (i=0;i<nexport;i++) {
        FoFDataIn[i].iGroup = pfof[Part[FoFDataIn[i].Index].GetID()];
        FoFDataIn[i].iGroupTask=mpi_foftask[Part[FoFDataIn[i].Index].GetID()];
        FoFDataIn[i].iLen=Len[FoFDataIn[i].Index];
    }

    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        nbuffer[recvTask] = 0;
        for (int k=0;k<recvTask;k++) nbuffer[recvTask]+=mpi_nsend[sendTask+k*NProcs]; //offset on local receiving buffer
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk = 0; ichunk < numsendrecv; ichunk++)
        {
            MPI_Sendrecv(&FoFDataIn[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct fofdata_in), MPI_BYTE,
                recvTask, TAG_FOF_A+ichunk,
                &FoFDataGet[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct fofdata_in), MPI_BYTE, 
                recvTask, TAG_FOF_A+ichunk, 
                MPI_COMM_WORLD, &status);
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }

}

/*! This routine searches the local particle list using the positions of the exported particles to see if any local particles
    met the linking criterion and any other FOF criteria of said exported particle. If that is the case, then the group id of the local particle
    and all other particles that belong to the same group are adjusted if the group id of the exported particle is smaller. This routine returns
    the number of links found between the local particles and all other exported particles from all other mpi domains.
    \todo need to update lengths if strucden flag used to limit particles for which real velocity density calculated
*/
Int_t MPILinkAcross(const Int_t nbodies, KDTree *&tree, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Int_tree_t *&Head, Int_tree_t *&Next, Double_t rdist2){
    Int_t i,j,k;
    Int_t links=0;
    Int_t *nn=new Int_t[nbodies];
    Int_t nt,ss,oldlen;
    Coordinate x;
    for (i=0;i<NImport;i++) {
        for (j=0;j<3;j++) x[j]=PartDataGet[i].GetPosition(j);
        //find all particles within a search radius of the imported particle
        nt=tree->SearchBallPosTagged(x, rdist2, nn);
        for (Int_t ii=0;ii<nt;ii++) {
            k=nn[ii];
            //if the imported particle does not belong to a group
            if (FoFDataGet[i].iGroup==0)
            {
                //if the current local particle's group head is zero and the exported particle group is zero
                //update the local particle's group id and the task to which it belongs
                //then one should link it and to make global decision base whether this task handles
                //the change on the PID of the particle
                if (pfof[Part[Head[k]].GetID()]==0&&Part[Head[k]].GetPID() > PartDataGet[i].GetPID()) {
                    pfof[Part[k].GetID()]=mpi_maxgid+mpi_gidoffset;///some unique identifier based on this task
                    mpi_gidoffset++;//increase unique identifier
                    Len[k]=1;
                    mpi_foftask[Part[k].GetID()]=FoFDataGet[i].iGroupTask;
                    links++;
                }
                //if the local particle does belong to a group, let the task from which the imported particle came from
                //handle the change
            }
            //if imported particle has already been linked
            else {
                //check to see if local particle has already been linked
                if (pfof[Part[Head[k]].GetID()]>0)  {
                    //as iGroups and pfof have been rank ordered globally
                    //proceed to link local particle to imported particle if
                    //its group id is larger
                    if(pfof[Part[Head[k]].GetID()] > FoFDataGet[i].iGroup) {
                        ss = Head[k];
                        oldlen=Len[k];
                        do{
                            pfof[Part[ss].GetID()]=FoFDataGet[i].iGroup;
                            mpi_foftask[Part[ss].GetID()]=FoFDataGet[i].iGroupTask;
                            Len[ss]=FoFDataGet[i].iLen+oldlen;
                        }while((ss = Next[ss]) >= 0);
                        FoFDataGet[i].iLen+=oldlen;
                        ss = Head[k];
                        links++;
                    }
                    //otherwise, let the task from which this imported particle came from
                    //handle the change
                }
                //if not in local group, add the particle to the imported particles group
                else {
                    pfof[Part[k].GetID()]=FoFDataGet[i].iGroup;
                    Len[k]=FoFDataGet[i].iLen;
                    mpi_foftask[Part[k].GetID()]=FoFDataGet[i].iGroupTask;
                    FoFDataGet[i].iLen+=1;
                    links++;
                }
            }
        }
    }
    delete[] nn;
    return links;
}
///link particles belonging to the same group across mpi domains using comparison function
Int_t MPILinkAcross(const Int_t nbodies, KDTree *&tree, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Int_tree_t *&Head, Int_tree_t *&Next, Double_t rdist2, FOFcompfunc &cmp, Double_t *params){
    Int_t i,k;
    Int_t links=0;
    Int_t *nn=new Int_t[nbodies];
    Int_t nt;
    for (i=0;i<NImport;i++) {
        nt=tree->SearchCriterionTagged(PartDataGet[i], cmp, params, nn);
        for (Int_t ii=0;ii<nt;ii++) {
            k=nn[ii];
            if (FoFDataGet[i].iGroup==0)
            {
                if (pfof[Part[Head[k]].GetID()]==0&&Part[Head[k]].GetPID() > PartDataGet[i].GetPID()) {
                    pfof[Part[k].GetID()]=mpi_maxgid+mpi_gidoffset;///some unique identifier based on this task
                    mpi_gidoffset++;//increase unique identifier
                    Len[k]=1;
                    mpi_foftask[Part[k].GetID()]=FoFDataGet[i].iGroupTask;
                    links++;
                }
            }
            else {
            if (pfof[Part[Head[k]].GetID()]>0)  {
                if(pfof[Part[Head[k]].GetID()] > FoFDataGet[i].iGroup) {
                    Int_t ss = Head[k];
                    do{
                        pfof[Part[ss].GetID()]=FoFDataGet[i].iGroup;
                        mpi_foftask[Part[ss].GetID()]=FoFDataGet[i].iGroupTask;
                        Len[ss]=FoFDataGet[i].iLen;
                    }while((ss = Next[ss]) >= 0);
                    ss = Head[k];
                    links++;
                }
            }
            else {
                pfof[Part[k].GetID()]=FoFDataGet[i].iGroup;
                Len[k]=FoFDataGet[i].iLen;
                mpi_foftask[Part[k].GetID()]=FoFDataGet[i].iGroupTask;
                links++;
            }
            }
        }
    }
    delete[] nn;
    return links;
}

///link particles belonging to the same group across mpi domains given a type check function
Int_t MPILinkAcross(const Int_t nbodies, KDTree *&tree, Particle *Part, Int_t *&pfof, Int_tree_t *&Len, Int_tree_t *&Head, Int_tree_t *&Next, Double_t rdist2, FOFcheckfunc &check, Double_t *params){
    Int_t i,j,k;
    Int_t links=0;
    Int_t *nn=new Int_t[nbodies];
    Int_t nt;
    Coordinate x;
    for (i=0;i<NImport;i++) {
        //if exported particle not in a group, do nothing
        if (FoFDataGet[i].iGroup==0) continue;
        for (j=0;j<3;j++) x[j]=PartDataGet[i].GetPosition(j);
        nt=tree->SearchBallPosTagged(x, rdist2, nn);
        for (Int_t ii=0;ii<nt;ii++) {
            k=nn[ii];
            //check that at least on of the particles meets the type criterion
            if (check(Part[k],params)!=0 && check(PartDataGet[i],params)!=0) continue;
            //if local particle in a group
            if (pfof[Part[Head[k]].GetID()]>0)  {
                //only change if both particles are appropriate type and group ids indicate local needs to be exported
                if (!(check(Part[k],params)==0 && check(PartDataGet[i],params)==0)) continue;
                if(pfof[Part[Head[k]].GetID()] > FoFDataGet[i].iGroup) {
                    Int_t ss = Head[k];
                    do{
                        pfof[Part[ss].GetID()]=FoFDataGet[i].iGroup;
                        mpi_foftask[Part[ss].GetID()]=FoFDataGet[i].iGroupTask;
                        Len[ss]=FoFDataGet[i].iLen;
                    }while((ss = Next[ss]) >= 0);
                    ss = Head[k];
                    links++;
                }
            }
            //if local particle not in a group and export is appropriate type, link
            else {
                if (check(PartDataGet[i],params)!=0) continue;
                pfof[Part[k].GetID()]=FoFDataGet[i].iGroup;
                Len[k]=FoFDataGet[i].iLen;
                mpi_foftask[Part[k].GetID()]=FoFDataGet[i].iGroupTask;
                links++;
            }
        }
    }
    delete[] nn;
    return links;
}
/*!
    Group particles belong to a group to a particular mpi thread so that locally easy to determine
    the maximum group size and reoder the group ids according to descending group size.
    return the new local number of particles
*/
Int_t MPIGroupExchange(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof){
    Int_t i, j, nexport,nimport,nlocal;
    Int_t nsend_local[NProcs],noffset_import[NProcs],noffset_export[NProcs],nbuffer[NProcs];
    //int sendTask,recvTask;
    int maxchunksize = 2147483648 / NProcs / sizeof(fofid_in);
    MPI_Status status;
    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    vr::Timer local_timer;
    int task;
    FoFGroupDataExport=NULL;
    FoFGroupDataLocal=NULL;
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    //first determine how big a local array is needed to store linked particles and broadcast information to create new nsend array
    nlocal=0;
    for (i=0;i<nbodies;i++) {
        if (mpi_foftask[i]!=ThisTask)
            nsend_local[mpi_foftask[i]]++;
    }
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    nexport=nimport=0;
    for (j=0;j<NProcs;j++){
        nimport+=mpi_nsend[ThisTask+j*NProcs];
        nexport+=mpi_nsend[j+ThisTask*NProcs];
    }
    //declare array for local storage of the appropriate size
    nlocal=nbodies-nexport+nimport;
    NImport=nimport;
    if (nexport >0) FoFGroupDataExport=new fofid_in[nexport];
    LOG(trace) <<" Exchanging ... nimport = "<<nimport<<", nexport="<<nexport<<" old nlocal="<<nbodies<<" new nlocal="<<nlocal;

    Int_t *storeval=new Int_t[nbodies];
    Noldlocal=nbodies-nexport;
    //store type in temporary array, then use type to store what task particle belongs to and sort values
    for (i=0;i<nbodies;i++) storeval[i]=Part[i].GetType();
    for (i=0;i<nbodies;i++) Part[i].SetType((mpi_foftask[i]!=ThisTask));
    // qsort(Part,nbodies,sizeof(Particle),TypeCompare);
    std::sort(Part, Part + nbodies, TypeCompareVec);
    for (i=0;i<nbodies;i++) Part[i].SetType(storeval[Part[i].GetID()]);
    //now use array to rearrange data
    for (i=0;i<nbodies;i++) storeval[i]=mpi_foftask[Part[i].GetID()];
    for (i=0;i<nbodies;i++) mpi_foftask[i]=storeval[i];
    for (i=0;i<nbodies;i++) storeval[i]=pfof[Part[i].GetID()];
    for (i=0;i<nbodies;i++) pfof[i]=storeval[i];
    for (i=0;i<nbodies;i++) Part[i].SetID(i);
    //for sorting purposes to place untagged particles at the end. Was done by setting type
    //now via storeval and ids
    for (i=0;i<nbodies;i++) storeval[i]=-pfof[Part[i].GetID()];
    for (i=0;i<nbodies;i++) Part[i].SetID(storeval[i]);
    if (nimport>0) FoFGroupDataLocal=new fofid_in[nimport];
    delete[] storeval;

    //determine offsets in arrays so that data contiguous with regards to processors for broadcasting
    //offset on transmitter end
    noffset_export[0]=0;
    for (j=1;j<NProcs;j++) noffset_export[j]=noffset_export[j-1]+mpi_nsend[(j-1)+ThisTask*NProcs];
    //offset on receiver end
    for (j=0;j<NProcs;j++) {
        noffset_import[j]=0;
        if (j!=ThisTask) for (int k=0;k<j;k++)noffset_import[j]+=mpi_nsend[ThisTask+k*NProcs];
    }
    for (j=0;j<NProcs;j++) nbuffer[j]=0;
    for (i=nbodies-nexport;i<nbodies;i++) {
        //if particle belongs to group that should belong on a different mpi thread, store for broadcasting
        task=mpi_foftask[i];
        if (task!=ThisTask) {
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].p=Part[i];
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].Index = i;
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].Task = task;
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].iGroup = pfof[i];
            //now that have all the particles that need broacasting, if extra information stored
            //then must also fill up appropriate hydro/star/bh buffers for communication.
        }
        nbuffer[task]++;
    }

    // if using mesh, mpi tasks of cells need to be updated
    if (opt.impiusemesh) {
        LOG(debug)<<" Updating mpi mesh ... ";
        // collect mesh data for exported particles
        vector<unordered_set<int>> newcellinfo(opt.numcells);
        for (i=0;i<nexport;i++) 
        {
            Coordinate x(FoFGroupDataExport[i].p.GetPosition());
            vector<unsigned int> ix(3);
            unsigned long long index;
            for (j=0; j<3; j++) ix[j] = floor(x[j]*opt.icellwidth[j]);
            index = ix[0]*opt.numcellsperdim*opt.numcellsperdim+ix[1]*opt.numcellsperdim+ix[2];
            newcellinfo[index].insert(FoFGroupDataExport[i].Task);
        }
        vector<int> newcellindex, newcelltask;
        for (auto i=0;i<opt.numcells;i++) 
        {
            if (newcellinfo[i].size() == 0) continue;
            for (auto &t:newcellinfo[i]) 
            {
                newcellindex.push_back(i);
                newcelltask.push_back(t);
            }
        }

        // now have for each cell set of tasks that the cell will belong to
        // must aggregate this information for all tasks
        // first determine how much is going to be sent by each task
        // and construct offsets and total
        int num = newcellindex.size(), mpi_num = 0;
        vector<int> mpi_sizes(NProcs), mpi_offsets(NProcs,0);
        MPI_Allgather(&num, 1, MPI_INTEGER, mpi_sizes.data(), 1, MPI_INTEGER, MPI_COMM_WORLD);
        mpi_num = mpi_sizes[0];
        for (auto i=1;i<NProcs;i++) {
            mpi_num += mpi_sizes[i];
            mpi_offsets[i] = mpi_sizes[i-1] + mpi_offsets[i-1];
        }
        // now collect information 
        vector<int> mpi_newcellindex(mpi_num), mpi_newcelltask(mpi_num);
        MPI_Allgatherv(newcellindex.data(), num, MPI_INTEGER, mpi_newcellindex.data(), mpi_sizes.data(), mpi_offsets.data(), MPI_INTEGER, MPI_COMM_WORLD);
        MPI_Allgatherv(newcelltask.data(), num, MPI_INTEGER, mpi_newcelltask.data(), mpi_sizes.data(), mpi_offsets.data(), MPI_INTEGER, MPI_COMM_WORLD);
        // process information 
        opt.newcellnodeids.resize(opt.numcells);
        for (auto i = 0; i < NProcs; i++) {
            auto istart = mpi_offsets[i];
            auto iend = istart + mpi_sizes[i];
            for (auto j = istart; j < iend; j++) {
                auto icell = mpi_newcellindex[j];
                auto itask = mpi_newcelltask[j];
                opt.newcellnodeids[icell].push_back(itask);
            }
        }
        LOG(debug) <<" Finished updating mpi mesh in "<<local_timer;
    }

    //now if there is extra information, strip off all the data from the FoFGroupDataExport
    //and store it explicitly into a buffer. Here are the buffers
    vector<Int_t> indices_gas_send;
    vector<float> propbuff_gas_send;
    vector<Int_t> indices_star_send;
    vector<float> propbuff_star_send;
    vector<Int_t> indices_bh_send;
    vector<float> propbuff_bh_send;
    vector<Int_t> indices_extra_dm_send;
    vector<float> propbuff_extra_dm_send;

    //now send the data.
    vr::Timer send_timer;
    LOG(debug)<<" Sending FOF data ... ";
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        vr::Timer comm_timer;
        LOG(trace) << "Send receive pair (" <<sendTask<<" "<<recvTask<<") "<<numsendrecv<<" chunks";
        for (auto ichunk=0;ichunk<numsendrecv;ichunk++)
        {
            // sending hydro, star and bh info            
            MPIFillFOFBuffWithHydroInfo(opt, cursendchunksize, &FoFGroupDataExport[noffset_export[recvTask]+sendoffset], Part, indices_gas_send, propbuff_gas_send);
            MPIFillFOFBuffWithStarInfo(opt, cursendchunksize, &FoFGroupDataExport[noffset_export[recvTask]+sendoffset], Part, indices_star_send, propbuff_star_send);
            MPIFillFOFBuffWithBHInfo(opt, cursendchunksize, &FoFGroupDataExport[noffset_export[recvTask]+sendoffset], Part, indices_bh_send, propbuff_bh_send);
            MPIFillFOFBuffWithExtraDMInfo(opt, cursendchunksize, &FoFGroupDataExport[noffset_export[recvTask]+sendoffset], Part, indices_extra_dm_send, propbuff_extra_dm_send);
            MPI_Sendrecv(&FoFGroupDataExport[noffset_export[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct fofid_in), MPI_BYTE,
                recvTask, TAG_FOF_C+ichunk,
                &FoFGroupDataLocal[noffset_import[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct fofid_in),
                MPI_BYTE, recvTask, TAG_FOF_C+ichunk, MPI_COMM_WORLD, &status);
            MPISendReceiveFOFHydroInfoBetweenThreads(opt, &FoFGroupDataLocal[noffset_import[recvTask]+recvoffset], indices_gas_send, propbuff_gas_send, recvTask, TAG_FOF_C+ichunk, mpi_comm);
            MPISendReceiveFOFStarInfoBetweenThreads(opt, &FoFGroupDataLocal[noffset_import[recvTask]+recvoffset], indices_star_send, propbuff_star_send, recvTask, TAG_FOF_C+ichunk, mpi_comm);
            MPISendReceiveFOFBHInfoBetweenThreads(opt, &FoFGroupDataLocal[noffset_import[recvTask]+recvoffset], indices_bh_send, propbuff_bh_send, recvTask, TAG_FOF_C+ichunk, mpi_comm);
            MPISendReceiveFOFExtraDMInfoBetweenThreads(opt, &FoFGroupDataLocal[noffset_import[recvTask]+recvoffset], indices_extra_dm_send, propbuff_extra_dm_send, recvTask, TAG_FOF_C+ichunk, mpi_comm);

            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
            LOG(trace) << "Finished send receive pair (" <<sendTask<<" "<<recvTask<<") in "<<comm_timer;
        }
    }
    LOG(debug)<<"Finished sending FOF information in "<<send_timer;
    Nlocal=nlocal;
    return nlocal;
}

/*!
    The baryon equivalent of \ref MPIGroupExchange. Here assume baryons are searched afterwards
*/
Int_t MPIBaryonGroupExchange(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof){
    Int_t i, j, nexport,nimport,nlocal;
    Int_t nsend_local[NProcs],noffset_import[NProcs],noffset_export[NProcs],nbuffer[NProcs];
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(fofid_in);
    int nsendchunks,nrecvchunks,numsendrecv;
    int sendoffset,recvoffset;
    int cursendchunksize,currecvchunksize;
    MPI_Status status;
    int task;
    FoFGroupDataExport=NULL;
    FoFGroupDataLocal=NULL;
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    //first determine how big a local array is needed to store linked particles and broadcast information to create new nsend array
    nlocal=0;
    for (i=0;i<nbodies;i++) {
        if (mpi_foftask[i]!=ThisTask)
            nsend_local[mpi_foftask[i]]++;
    }
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    nexport=nimport=0;
    for (j=0;j<NProcs;j++){
        nimport+=mpi_nsend[ThisTask+j*NProcs];
        nexport+=mpi_nsend[j+ThisTask*NProcs];
    }
    //declare array for local storage of the appropriate size
    nlocal=nbodies-nexport+nimport;
    NImport=nimport;
    if (nexport >0) FoFGroupDataExport=new fofid_in[nexport];
    else FoFGroupDataExport=new fofid_in[1];

    Nmemlocalbaryon=Nlocalbaryon[0];
    for (i=0;i<nbodies;i++) Part[i].SetID(i);
    Int_t *storeval=new Int_t[nbodies];
    //if trying to reduce memory allocation,  if nlocal < than the memory allocated adjust local list so that all particles to be exported are near the end.
    //and allocate the appropriate memory for pfof and mpi_idlist. otherwise, need to actually copy the particle data into FoFGroupDataLocal and proceed
    //as normal, storing info, send info, delete particle array, allocate a new array large enough to store info and copy over info
    ///\todo eventually I should replace arrays with vectors so that the size can change, removing the need to free and allocate
    if (nlocal<=Nmemlocalbaryon) {
        Noldlocal=nbodies-nexport;
        for (i=0;i<nbodies;i++) storeval[i]=Part[i].GetType();
        for (i=0;i<nbodies;i++) Part[i].SetType((mpi_foftask[i]!=ThisTask));
        //qsort(Part,nbodies,sizeof(Particle),TypeCompare);
        std::sort(Part, Part + nbodies, TypeCompareVec);
        for (i=0;i<nbodies;i++) Part[i].SetType(storeval[Part[i].GetID()]);
        //now use array to rearrange data
        for (i=0;i<nbodies;i++) storeval[i]=mpi_foftask[Part[i].GetID()];
        for (i=0;i<nbodies;i++) mpi_foftask[i]=storeval[i];
        for (i=0;i<nbodies;i++) storeval[i]=pfof[Part[i].GetID()];
        for (i=0;i<nbodies;i++) pfof[i]=storeval[i];
        for (i=0;i<nbodies;i++) Part[i].SetID(i);
        //now via storeval and ids
        for (i=0;i<nbodies;i++) storeval[i]=-pfof[Part[i].GetID()];
        for (i=0;i<nbodies;i++) Part[i].SetID(storeval[i]);
        if (nimport>0) FoFGroupDataLocal=new fofid_in[nimport];
    }
    //otherwise use FoFGroupDataLocal to store all the necessary data
    else {
        FoFGroupDataLocal=new fofid_in[nlocal];
        for (i=0;i<nbodies;i++) storeval[i]=Part[i].GetType();
        for (i=0;i<nbodies;i++) Part[i].SetType((mpi_foftask[i]!=ThisTask));
        //qsort(Part,nbodies,sizeof(Particle),TypeCompare);
        std::sort(Part, Part + nbodies,  TypeCompareVec);
        for (i=0;i<nbodies;i++) Part[i].SetType(storeval[Part[i].GetID()]);
        Int_t nn=nbodies-nexport;
        for (i=0;i<nn;i++) {
            FoFGroupDataLocal[i].p=Part[i];
            FoFGroupDataLocal[i].Index = i;
            FoFGroupDataLocal[i].Task = ThisTask;
            FoFGroupDataLocal[i].iGroup = pfof[Part[i].GetID()];
        }
        for (i=nn;i<nbodies;i++) storeval[i]=mpi_foftask[Part[i].GetID()];
        for (i=nn;i<nbodies;i++) mpi_foftask[i]=storeval[i];
        for (i=nn;i<nbodies;i++) storeval[i]=pfof[Part[i].GetID()];
        for (i=nn;i<nbodies;i++) pfof[i]=storeval[i];
        for (i=nn;i<nbodies;i++) Part[i].SetID(i);
    }
    delete[] storeval;
    //determine offsets in arrays so that data contiguous with regards to processors for broadcasting
    //offset on transmitter end
    noffset_export[0]=0;
    for (j=1;j<NProcs;j++) noffset_export[j]=noffset_export[j-1]+mpi_nsend[(j-1)+ThisTask*NProcs];
    //offset on receiver end
    for (j=0;j<NProcs;j++) {
        if (nlocal<Nlocalbaryon[0]) noffset_import[j]=0;
        else noffset_import[j]=nbodies-nexport;
        if (j!=ThisTask) for (int k=0;k<j;k++)noffset_import[j]+=mpi_nsend[ThisTask+k*NProcs];
    }
    for (j=0;j<NProcs;j++) nbuffer[j]=0;
    for (i=nbodies-nexport;i<nbodies;i++) {
        //if particle belongs to group that should belong on a different mpi thread, store for broadcasting
        task=mpi_foftask[i];
        if (task!=ThisTask) {
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].p=Part[i];
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].Index = i;
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].Task = task;
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].iGroup = pfof[i];
        }
        nbuffer[task]++;
    }
    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk=0;ichunk<numsendrecv;ichunk++)
        {
            MPI_Sendrecv(&FoFGroupDataExport[noffset_export[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct fofid_in), MPI_BYTE,
                recvTask, TAG_FOF_C+ichunk,
                &FoFGroupDataLocal[noffset_import[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct fofid_in), MPI_BYTE, 
                recvTask, TAG_FOF_C+ichunk, 
                MPI_COMM_WORLD, &status);

            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }
    Nlocalbaryon[0]=nlocal;
    return nlocal;
}

///Determine the local number of groups and their sizes (groups must be local to an mpi thread)
Int_t MPICompileGroups(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t minsize){
    Int_t i,start,ngroups;
    Int_t *numingroup,**plist;
    ngroups=0;

    //if not using mpi mesh, need to update mpi boundaries based on these exported particles
    //note that must ensure that periodicity account for
    if (!opt.impiusemesh)
    {
        MPI_Domain localdomain;
        for (auto j=0; j<3; j++)
        {
            localdomain.bnd[j][0] = mpi_domain[ThisTask].bnd[j][0];
            localdomain.bnd[j][1] = mpi_domain[ThisTask].bnd[j][1];
        }
        for (i=Noldlocal;i<nbodies;i++) {
            Coordinate x(FoFGroupDataLocal[i-Noldlocal].p.GetPosition());
            // adjust for period based on local mpi boundary
            for (auto j=0; j<3; j++)
            {
                if (x[j] < mpi_domain[ThisTask].bnd[j][0])
                    localdomain.bnd[j][0] = x[j];
                else if (x[j] > mpi_domain[ThisTask].bnd[j][1])
                    localdomain.bnd[j][1] = x[j];
            }
        }
        //now update the mpi_domains again
        MPI_Allgather(&localdomain, sizeof(MPI_Domain), MPI_BYTE, mpi_domain, sizeof(MPI_Domain), MPI_BYTE, MPI_COMM_WORLD);
    }
    for (i=Noldlocal;i<nbodies;i++) {
        Part[i]=FoFGroupDataLocal[i-Noldlocal].p;
        //note that before used type to sort particles
        //now use id
        Part[i].SetID(-FoFGroupDataLocal[i-Noldlocal].iGroup);
    }
    //used to use ID store store group id info
    std::sort(Part, Part + nbodies, IDCompareVec);
    //determine the # of groups, their size and the current group ID
    for (i=0,start=0;i<nbodies;i++) {
        if (Part[i].GetID()!=Part[start].GetID()) {
            //if group is too small set type to zero, which currently is used to store the group id
            if ((i-start)<minsize) for (Int_t j=start;j<i;j++) Part[j].SetID(0);
            else ngroups++;
            start=i;
        }
        if (Part[i].GetID()>=0) break;
    }
    //again resort to move untagged particles to the end.
    std::sort(Part, Part + nbodies, IDCompareVec);
    for (i=nbodies-NExport; i< nbodies; i++) Part[i].SetID(0);
    //now adjust pfof and ids.
    for (i=0;i<nbodies;i++) {
        pfof[i]=-Part[i].GetID();
        Part[i].SetID(i);
    }    
    numingroup=new Int_t[ngroups+1];
    plist=new Int_t*[ngroups+1];
    ngroups=1;//offset as group zero is untagged
    for (i=0,start=0;i<nbodies;i++) {
        if (pfof[i]!=pfof[start]) {
            numingroup[ngroups]=i-start;
            plist[ngroups]=new Int_t[numingroup[ngroups]];
            for (Int_t j=start,count=0;j<i;j++) plist[ngroups][count++]=j;
            ngroups++;
            start=i;
        }
        if (pfof[i]==0) break;
    }
    ngroups--;

    //reorder groups ids according to size
    ReorderGroupIDs(ngroups,ngroups,numingroup,pfof,plist);
    for (i=1;i<=ngroups;i++) delete[] plist[i];
    delete[] plist;
    delete[] numingroup;
    //broadcast number of groups so that ids can be properly offset
    MPI_Allgather(&ngroups, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
    if(FoFGroupDataLocal!=NULL) delete[] FoFGroupDataLocal;
    if(FoFGroupDataExport!=NULL) delete[] FoFGroupDataExport;
    return ngroups;
}

///Similar to \ref MPICompileGroups but optimised for separate baryon search
///\todo need to update to reflect vector implementation
Int_t MPIBaryonCompileGroups(Options &opt, const Int_t nbodies, Particle *Part, Int_t *&pfof, Int_t minsize, int iorder){
    Int_t i,start,ngroups;
    Int_t *numingroup, **plist;
    ngroups=0;

    //if minimizing memory load when using mpi (by adding extra routines to determine memory required)
    //first check to see if local memory is enough to contained expected number of particles
    //if local mem is enough, copy data from the FoFGroupDataLocal
    if(Nmemlocalbaryon>nbodies) {
        for (i=Noldlocal;i<nbodies;i++) {
            Part[i]=FoFGroupDataLocal[i-Noldlocal].p;
            Part[i].SetID(-FoFGroupDataLocal[i-Noldlocal].iGroup);
        }
        //now use ID
        //qsort(Part,nbodies,sizeof(Particle),IDCompare);
        std::sort(Part, Part + nbodies, IDCompareVec);
        //determine the # of groups, their size and the current group ID
        for (i=0,start=0;i<nbodies;i++) {
            if (Part[i].GetID()!=Part[start].GetID()) {
                //if group is too small set type to zero, which currently is used to store the group id
                if ((i-start)<minsize) for (Int_t j=start;j<i;j++) Part[j].SetID(0);
                else ngroups++;
                start=i;
            }
            if (Part[i].GetID()==0) break;
        }

        //again resort to move untagged particles to the end.
        //qsort(Part,nbodies,sizeof(Particle),IDCompare);
        std::sort(Part, Part + nbodies, IDCompareVec);
        //now adjust pfof and ids.
        for (i=0;i<nbodies;i++) {pfof[i]=-Part[i].GetID();Part[i].SetID(i);}
        numingroup=new Int_t[ngroups+1];
        plist=new Int_t*[ngroups+1];
        ngroups=1;//offset as group zero is untagged
        for (i=0,start=0;i<nbodies;i++) {
            if (pfof[i]!=pfof[start]) {
                numingroup[ngroups]=i-start;
                plist[ngroups]=new Int_t[numingroup[ngroups]];
                for (Int_t j=start,count=0;j<i;j++) plist[ngroups][count++]=j;
                ngroups++;
                start=i;
            }
            if (pfof[i]==0) break;
        }
        ngroups--;
    }
    else {
        //sort local list
        //qsort(FoFGroupDataLocal, nbodies, sizeof(struct fofid_in), fof_id_cmp);
        std::sort(FoFGroupDataLocal, FoFGroupDataLocal + nbodies, fof_id_cmp_vec);
        //determine the # of groups, their size and the current group ID
        for (i=0,start=0;i<nbodies;i++) {
            if (FoFGroupDataLocal[i].iGroup!=FoFGroupDataLocal[start].iGroup) {
                if ((i-start)<minsize){
                    for (Int_t j=start;j<i;j++) FoFGroupDataLocal[j].iGroup=0;
                }
                else ngroups++;
                start=i;
            }
            if (FoFGroupDataLocal[i].iGroup==0) break;
        }
        //now sort again which will put particles group then id order, and determine size of groups and their current group id;
        //qsort(FoFGroupDataLocal, nbodies, sizeof(struct fofid_in), fof_id_cmp);
        std::sort(FoFGroupDataLocal, FoFGroupDataLocal + nbodies, fof_id_cmp_vec);
        numingroup=new Int_t[ngroups+1];
        plist=new Int_t*[ngroups+1];
        ngroups=1;//offset as group zero is untagged
        for (i=0,start=0;i<nbodies;i++) {
            if (FoFGroupDataLocal[i].iGroup!=FoFGroupDataLocal[start].iGroup) {
                numingroup[ngroups]=i-start;
                plist[ngroups]=new Int_t[numingroup[ngroups]];
                for (Int_t j=start,count=0;j<i;j++) plist[ngroups][count++]=j;
                ngroups++;
                start=i;
            }
            if (FoFGroupDataLocal[i].iGroup==0) break;
        }
        ngroups--;
        for (i=0;i<nbodies;i++) pfof[i]=FoFGroupDataLocal[i].iGroup;
        //and store the particles global ids
        for (i=0;i<nbodies;i++) {
            Part[i]=FoFGroupDataLocal[i].p;
            Part[i].SetID(i);
        }
    }
    //reorder groups ids according to size if required.
    if (iorder) ReorderGroupIDs(ngroups,ngroups,numingroup,pfof,plist);
    for (i=1;i<=ngroups;i++) if (numingroup[i]>0) delete[] plist[i];
    delete[] plist;
    delete[] numingroup;

    //broadcast number of groups so that ids can be properly offset
    MPI_Allgather(&ngroups, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
    if(FoFGroupDataLocal!=NULL) delete[] FoFGroupDataLocal;
    if(FoFGroupDataExport!=NULL) delete[] FoFGroupDataExport;
    return ngroups;
}

///Determine which exported dm particle is closest in phase-space to a local baryon particle and assign that particle to the group of that dark matter particle if is closest particle
Int_t MPISearchBaryons(const Int_t nbaryons, Particle *&Pbaryons, Int_t *&pfofbaryons, Int_t *numingroup, Double_t *localdist, Int_t nsearch, Double_t *param, Double_t *period)
{
    Double_t D2, dval, rval;
    Coordinate x1;
    Particle p1;
    Int_t  i, j, k, pindex,nexport=0;
    int tid;
    FOFcompfunc fofcmp=FOF6d;
    Int_t *nnID;
    Double_t *dist2;
    if (NImport>0) {
    //now dark matter particles associated with a group existing on another mpi domain are local and can be searched.
    KDTree *mpitree =  new KDTree(PartDataGet,NImport,nsearch/2,mpitree->TPHYS,mpitree->KEPAN,100,0,0,0,period);
    if (nsearch>NImport) nsearch=NImport;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,j,k,tid,p1,pindex,x1,D2,dval,rval,nnID,dist2)
{
    nnID=new Int_t[nsearch];
    dist2=new Double_t[nsearch];
#pragma omp for reduction(+:nexport)
#endif
    for (i=0;i<nbaryons;i++)
    {
#ifdef USEOPENMP
        tid=omp_get_thread_num();
#else
        tid=0;
#endif
        p1=Pbaryons[i];
        x1=Coordinate(p1.GetPosition());
        rval=MAXVALUE;
        dval=localdist[i];
        mpitree->FindNearestPos(x1, nnID, dist2,nsearch);
        if (dist2[0]<param[6]) {
        for (j=0;j<nsearch;j++) {
            D2=0;
            pindex=PartDataGet[nnID[j]].GetID();
            if (numingroup[pfofbaryons[i]]<FoFDataGet[pindex].iLen) {
                if (fofcmp(p1,PartDataGet[nnID[j]],param)) {
                    for (k=0;k<3;k++) {
                        D2+=(p1.GetPosition(k)-PartDataGet[nnID[j]].GetPosition(k))*(p1.GetPosition(k)-PartDataGet[nnID[j]].GetPosition(k))/param[6]+(p1.GetVelocity(k)-PartDataGet[nnID[j]].GetVelocity(k))*(p1.GetVelocity(k)-PartDataGet[nnID[j]].GetVelocity(k))/param[7];
                    }
#ifdef GASON
                    D2+=p1.GetU()/param[7];
#endif
                    if (dval>D2) {dval=D2;pfofbaryons[i]=FoFDataGet[pindex].iGroup;rval=dist2[j];mpi_foftask[i]=FoFDataGet[pindex].iGroupTask;}
                }
            }
        }
        }
        nexport+=(mpi_foftask[i]!=ThisTask);
    }
    delete[] nnID;
    delete[] dist2;
#ifdef USEOPENMP
}
#endif
    }
    return nexport;
}

Int_t MPIBaryonExchange(Options &opt, const Int_t nbaryons, Particle *Pbaryons, Int_t *pfofbaryons){
    Int_t i, j, nexport,nimport,nlocal;
    Int_t nsend_local[NProcs],noffset_import[NProcs],noffset_export[NProcs],nbuffer[NProcs];
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(fofid_in);
    int nsendchunks,nrecvchunks,numsendrecv;
    int sendoffset,recvoffset;
    int cursendchunksize,currecvchunksize;
    MPI_Status status;
    int task;
    //initial containers to send info across threads
    FoFGroupDataExport=NULL;
    FoFGroupDataLocal=NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    //first determine how big a local array is needed to store tagged baryonic particles
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    nlocal=0;
    for (i=0;i<nbaryons;i++) {
        if (mpi_foftask[i]!=ThisTask)
            nsend_local[mpi_foftask[i]]++;
    }
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    nexport=nimport=0;
    for (j=0;j<NProcs;j++){
        nimport+=mpi_nsend[ThisTask+j*NProcs];
        nexport+=mpi_nsend[j+ThisTask*NProcs];
    }
    //declare array for local storage of the appropriate size
    nlocal=nbaryons-nexport+nimport;
    //store import number
    NImport=nimport;
    FoFGroupDataExport=new fofid_in[nexport+1];//+1 just buffer to ensure that if nothing broadcast, easy to allocate and deallocate memory

    //if trying to reduce memory allocation, then need to check amount stored locally and how much that needs to be adjusted by
    //if nlocal < than the memory allocated adjust local list so that all particles to be exported are near the end.
    //and allocate the appropriate memory for pfofbaryons and mpi_idlist. otherwise, need to actually copy the particle data into FoFGroupDataLocal and proceed
    //as normal, storing info, send info, delete particle array, allocate a new array large enough to store info and copy over info
    ///\todo eventually I should replace arrays with vectors so that the size can change, removing the need to free and allocate
    Int_t *storeval=new Int_t[nbaryons];
    if (nlocal<Nmemlocal) {
        Noldlocal=nbaryons-nexport;
        for (i=0;i<nbaryons;i++) storeval[i]=Pbaryons[i].GetType();
        for (i=0;i<nbaryons;i++) Pbaryons[i].SetType((mpi_foftask[i]!=ThisTask));
        //qsort(Pbaryons,nbaryons,sizeof(Particle),TypeCompare);
        std::sort(Pbaryons, Pbaryons + nbaryons, TypeCompareVec);
        for (i=0;i<nbaryons;i++) Pbaryons[i].SetType(storeval[Pbaryons[i].GetID()]);
        //now use array to rearrange data
        for (i=0;i<nbaryons;i++) storeval[i]=mpi_foftask[Pbaryons[i].GetID()];
        for (i=0;i<nbaryons;i++) mpi_foftask[i]=storeval[i];
        for (i=0;i<nbaryons;i++) storeval[i]=pfofbaryons[Pbaryons[i].GetID()];
        for (i=0;i<nbaryons;i++) pfofbaryons[i]=storeval[i];
        for (i=0;i<nbaryons;i++) Pbaryons[i].SetID(i);
        //for sorting purposes to place untagged particles at the end.
        for (i=0;i<nbaryons;i++) storeval[i]=-pfofbaryons[Pbaryons[i].GetID()];
        for (i=0;i<nbaryons;i++) Pbaryons[i].SetID(storeval[i]);
        if (nimport>0) FoFGroupDataLocal=new fofid_in[nimport];
    }
    //otherwise use FoFGroupDataLocal to store all the necessary data
    else {
        FoFGroupDataLocal=new fofid_in[nlocal];
        for (i=0;i<nbaryons;i++) storeval[i]=Pbaryons[i].GetType();
        for (i=0;i<nbaryons;i++) Pbaryons[i].SetType((mpi_foftask[i]!=ThisTask));
        // qsort(Pbaryons,nbaryons,sizeof(Particle),TypeCompare);
        std::sort(Pbaryons, Pbaryons + nbaryons, TypeCompareVec);
        for (i=0;i<nbaryons;i++) Pbaryons[i].SetType(storeval[Pbaryons[i].GetID()]);
        Int_t nn=nbaryons-nexport;
        for (i=0;i<nn;i++) {
            FoFGroupDataLocal[i].p=Pbaryons[i];
            FoFGroupDataLocal[i].Index = i;
            FoFGroupDataLocal[i].Task = ThisTask;
            FoFGroupDataLocal[i].iGroup = pfofbaryons[Pbaryons[i].GetID()];
        }
        for (i=nn;i<nbaryons;i++) storeval[i]=mpi_foftask[Pbaryons[i].GetID()];
        for (i=nn;i<nbaryons;i++) mpi_foftask[i]=storeval[i];
        for (i=nn;i<nbaryons;i++) storeval[i]=pfofbaryons[Pbaryons[i].GetID()];
        for (i=nn;i<nbaryons;i++) pfofbaryons[i]=storeval[i];
        for (i=nn;i<nbaryons;i++) Pbaryons[i].SetID(i);
    }
    delete[] storeval;

    //determine offsets in arrays so that data contiguous with regards to processors for broadcasting
    //offset on transmitter end
    noffset_export[0]=0;
    for (j=1;j<NProcs;j++) noffset_export[j]=noffset_export[j-1]+mpi_nsend[(j-1)+ThisTask*NProcs];
    for (j=0;j<NProcs;j++) {
        if (nlocal<Nmemlocal) noffset_import[j]=0;
        else noffset_import[j]=nbaryons-nexport;
        if (j!=ThisTask) for (int k=0;k<j;k++)noffset_import[j]+=mpi_nsend[ThisTask+k*NProcs];
    }
    for (j=0;j<NProcs;j++) nbuffer[j]=0;
    for (i=nbaryons-nexport;i<nbaryons;i++) {
        //if particle belongs to group that should belong on a different mpi thread, store for broadcasting
        task=mpi_foftask[i];
        if (task!=ThisTask) {
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].p=Pbaryons[i];
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].Index = i;
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].Task = task;
            FoFGroupDataExport[noffset_export[task]+nbuffer[task]].iGroup = pfofbaryons[i];
        }
        nbuffer[task]++;
    }
    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk=0;ichunk<numsendrecv;ichunk++)
        {
            MPI_Sendrecv(&FoFGroupDataExport[noffset_export[recvTask]+sendoffset],
                cursendchunksize * sizeof(struct fofid_in), MPI_BYTE,
                recvTask, TAG_FOF_C+ichunk,
                &FoFGroupDataLocal[noffset_import[recvTask]+recvoffset],
                currecvchunksize * sizeof(struct fofid_in), MPI_BYTE, 
                recvTask, TAG_FOF_C+ichunk, 
                MPI_COMM_WORLD, &status);

            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }
    Nlocalbaryon[0]=nlocal;
    return nlocal;
}
//@}


/// \name FOF routines related to modifying group ids
//@{

///This alters the group ids by an offset determined by the number of groups on all previous mpi threads so that the group
///has a unique group id. Prior to this, group ids are determined locally.
inline void MPIAdjustGroupIDs(const Int_t nbodies,Int_t *pfof) {
    Int_t noffset=0;
    for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
    for (Int_t i=0;i<nbodies;i++) if (pfof[i]>0) pfof[i]+=noffset;
}

///Collect FOF from all
void MPICollectFOF(const Int_t nbodies, Int_t *&pfof){
    Int_t sendTask,recvTask;
    MPI_Status status;
    //if using mpi, offset the pfof so that group ids are no unique, before just local to thread
    MPIAdjustGroupIDs(Nlocal,pfof);
    //now send the data from all MPI threads to thread zero
    //first must send how much data is local to a processor
    Int_t nsend_local[NProcs];
    for (int i=0;i<NProcs;i++) nsend_local[i]=0;
    if (ThisTask!=0)nsend_local[0]=Nlocal;
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    recvTask=0;
    //next copy task zero pfof into global mpi_pfof to the appropriate indices
    if (ThisTask==0) {
        for (Int_t i=0;i<Nlocal;i++) mpi_pfof[mpi_indexlist[i]]=pfof[i];
    }
    //then for each processor copy their values into local pfof and mpi_indexlist
    //note mpi_idlist contains id values whereas indexlist contains the index order of how particles were loaded
    //this requires determining the largest size needed
    if (ThisTask==0) {
        Int_t maxnlocal=0;
        for(int j=1;j<NProcs;j++) if (maxnlocal<mpi_nsend[j*NProcs]) maxnlocal=mpi_nsend[j*NProcs];
        delete[] pfof;
        delete[] mpi_indexlist;
        pfof=new Int_t[maxnlocal];
        mpi_indexlist=new Int_t[maxnlocal];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //now for each mpi task, copy appropriate data to mpi thread 0 local buffers
    for(int j=1;j<NProcs;j++)
    {
        sendTask=j;
        recvTask=0;
        if (ThisTask==sendTask) {
            MPI_Ssend(pfof, Nlocal , MPI_Int_t,0, TAG_FOF_D, MPI_COMM_WORLD);
            MPI_Ssend(mpi_indexlist, Nlocal , MPI_Int_t,0, TAG_FOF_E, MPI_COMM_WORLD);
        }
        if(ThisTask==0) {
            MPI_Recv(pfof, mpi_nsend[sendTask*NProcs], MPI_Int_t, sendTask, TAG_FOF_D, MPI_COMM_WORLD, &status);
            MPI_Recv(mpi_indexlist, mpi_nsend[sendTask*NProcs], MPI_Int_t, sendTask, TAG_FOF_E, MPI_COMM_WORLD, &status);
            for (Int_t i=0;i<mpi_nsend[sendTask*NProcs];i++) mpi_pfof[mpi_indexlist[i]]=pfof[i];
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
//@}

/// \name Routines related to distributing the grid cells used to calculate the coarse-grained mean field
//@{
/*! Collects all the grid data
*/
void MPIBuildGridData(const Int_t ngrid, GridCell *grid, Coordinate *gvel, Matrix *gveldisp){
    Int_t i, j;
    Int_t nsend_local[NProcs],noffset[NProcs];
    Int_t sendTask,recvTask;
    MPI_Status status;

    for (j=0;j<NProcs;j++) nsend_local[j]=Ngridlocal;
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    noffset[0]=0;
    for (j=1;j<NProcs;j++) noffset[j]=noffset[j]+mpi_nsend[ThisTask+j*NProcs];
    for (i=0;i<Ngridlocal;i++) {
        for (j=0;j<3;j++) mpi_grid[noffset[ThisTask]+i].xm[j]=grid[i].xm[j];
        mpi_gvel[noffset[ThisTask]+i]=gvel[i];
        mpi_gveldisp[noffset[ThisTask]+i]=gveldisp[i];
    }

    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        MPI_Sendrecv(grid,
            Ngridlocal* sizeof(GridCell), MPI_BYTE,
            recvTask, TAG_GRID_A,
            &mpi_grid[noffset[recvTask]],
            mpi_nsend[sendTask+recvTask * NProcs] * sizeof(GridCell),
            MPI_BYTE, recvTask, TAG_GRID_A, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(gvel,
            Ngridlocal* sizeof(Coordinate), MPI_BYTE,
            recvTask, TAG_GRID_B,
            &mpi_gvel[noffset[recvTask]],
            mpi_nsend[sendTask+recvTask * NProcs] * sizeof(Coordinate),
            MPI_BYTE, recvTask, TAG_GRID_B, MPI_COMM_WORLD, &status);
        MPI_Sendrecv(gveldisp,
            Ngridlocal* sizeof(Matrix), MPI_BYTE,
            recvTask, TAG_GRID_C,
            &mpi_gveldisp[noffset[recvTask]],
            mpi_nsend[sendTask+recvTask * NProcs] * sizeof(Matrix),
            MPI_BYTE, recvTask, TAG_GRID_C, MPI_COMM_WORLD, &status);
    }
}
//@}

/// \name config updates for MPI
//@{
///Update config option for particle types present
void MPIUpdateUseParticleTypes(Options &opt)
{
    MPI_Bcast(&(opt.iusestarparticles),sizeof(opt.iusestarparticles),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.iusesinkparticles),sizeof(opt.iusesinkparticles),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.iusewindparticles),sizeof(opt.iusewindparticles),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.iusetracerparticles),sizeof(opt.iusetracerparticles),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.iuseextradarkparticles),sizeof(opt.iuseextradarkparticles),MPI_BYTE,0,MPI_COMM_WORLD);
}
//@}

/// \name comparison functions used to assign particles to a specific mpi thread
//@{
///comprasion function used to sort particles for export so that all particles being exported to the same processor are in a contiguous block and well ordered
///\todo, I should remove this and determine the ordering before hand.
int fof_export_cmp(const void *a, const void *b)
{
  if(((struct fofdata_in *) a)->Task < (((struct fofdata_in *) b)->Task))
    return -1;
  if(((struct fofdata_in *) a)->Task > (((struct fofdata_in *) b)->Task))
    return +1;
  return 0;
}

///comprasion function used to sort particles for export so that all particles being exported to the same processor are in a contiguous block and well ordered
int nn_export_cmp(const void *a, const void *b)
{
  if(((struct nndata_in *) a)->ToTask < (((struct nndata_in *) b)->ToTask))
    return -1;
  if(((struct nndata_in *) a)->ToTask > (((struct nndata_in *) b)->ToTask))
    return +1;
  return 0;
}

///comprasion function used to sort grouped particles so that easy to determine total number of groups locally, size of groups, etc.
int fof_id_cmp(const void *a, const void *b)
{
  if(((struct fofid_in *) a)->iGroup > (((struct fofid_in *) b)->iGroup))
    return -1;

  if(((struct fofid_in *) a)->iGroup < (((struct fofid_in *) b)->iGroup))
    return +1;

  if(((struct fofid_in *) a)->p.GetType() < (((struct fofid_in *) b)->p.GetType()))
    return -1;

  if(((struct fofid_in *) a)->p.GetType() > (((struct fofid_in *) b)->p.GetType()))
    return +1;


  if(((struct fofid_in *) a)->p.GetID() < (((struct fofid_in *) b)->p.GetID()))
    return -1;

  if(((struct fofid_in *) a)->p.GetID() > (((struct fofid_in *) b)->p.GetID()))
    return +1;

  return 0;
}

bool fof_export_cmp_vec(const fofdata_in &a, const fofdata_in &b)
{
    if (a.Task < b.Task) return true;
    else return false ;
}

bool nn_export_cmp_vec(const nndata_in &a, const nndata_in &b)
{
    if(a.ToTask < b.ToTask) return true;
    else return false;
}

bool fof_id_cmp_vec(const fofid_in &a, const fofid_in &b)
{
    if (a.iGroup < b.iGroup) return true;
    else if (a.iGroup > b.iGroup) return false ;
    else {
        if (a.p.GetType() < b.p.GetType()) return true;
        else if (a.p.GetType() > b.p.GetType()) return false;
        else {
            if (a.p.GetID() < b.p.GetID()) return true;
            else return false;
        }
    }
}

//@}

///\name mesh MPI decomposition related functions
//@{
vector<int> MPIGetCellListInSearchUsingMesh(Options &opt, Double_t xsearch[3][2], bool ignorelocalcells)
{
    int ixstart,iystart,izstart,ixend,iyend,izend,index;
    vector<int> celllist;
    ixstart=floor(xsearch[0][0]*opt.icellwidth[0]);
    ixend=floor(xsearch[0][1]*opt.icellwidth[0]);
    iystart=floor(xsearch[1][0]*opt.icellwidth[1]);
    iyend=floor(xsearch[1][1]*opt.icellwidth[1]);
    izstart=floor(xsearch[2][0]*opt.icellwidth[2]);
    izend=floor(xsearch[2][1]*opt.icellwidth[2]);

    for (auto ix=ixstart;ix<=ixend;ix++){
        for (auto iy=iystart;iy<=iyend;iy++){
            for (auto iz=izstart;iz<=izend;iz++){
                index=0;
                if (iz<0) index+=opt.numcellsperdim+iz;
                else if (iz>=opt.numcellsperdim) index+=iz-opt.numcellsperdim;
                else index+=iz;
                if (iy<0) index+=(opt.numcellsperdim+iy)*opt.numcellsperdim;
                else if (iy>=opt.numcellsperdim) index+=(iy-opt.numcellsperdim)*opt.numcellsperdim;
                else index+=iy*opt.numcellsperdim;
                if (ix<0) index+=(opt.numcellsperdim+ix)*opt.numcellsperdim*opt.numcellsperdim;
                else if (ix>=opt.numcellsperdim) index+=(ix-opt.numcellsperdim)*opt.numcellsperdim*opt.numcellsperdim;
                else index+=ix*opt.numcellsperdim*opt.numcellsperdim;
                //if ignoring local cells and cell not local, add to cell list
                //or add regardless if not ignoring local cells
                if (ignorelocalcells && opt.cellnodeids[index]==ThisTask) continue;
                celllist.push_back(index);
            }
        }
    }
    return celllist;
}

vector<int> MPIGetCellNodeIDListInSearchUsingMesh(Options &opt, Double_t xsearch[3][2])
{
    int ixstart,iystart,izstart,ixend,iyend,izend,index;
    vector<int> cellnodeidlist;
    ixstart=floor(xsearch[0][0]*opt.icellwidth[0]);
    ixend=floor(xsearch[0][1]*opt.icellwidth[0]);
    iystart=floor(xsearch[1][0]*opt.icellwidth[1]);
    iyend=floor(xsearch[1][1]*opt.icellwidth[1]);
    izstart=floor(xsearch[2][0]*opt.icellwidth[2]);
    izend=floor(xsearch[2][1]*opt.icellwidth[2]);

    for (auto ix=ixstart;ix<=ixend;ix++){
        for (auto iy=iystart;iy<=iyend;iy++){
            for (auto iz=izstart;iz<=izend;iz++){
                index=0;
                if (iz<0) index+=opt.numcellsperdim+iz;
                else if (iz>=opt.numcellsperdim) index+=iz-opt.numcellsperdim;
                else index+=iz;
                if (iy<0) index+=(opt.numcellsperdim+iy)*opt.numcellsperdim;
                else if (iy>=opt.numcellsperdim) index+=(iy-opt.numcellsperdim)*opt.numcellsperdim;
                else index+=iy*opt.numcellsperdim;
                if (ix<0) index+=(opt.numcellsperdim+ix)*opt.numcellsperdim*opt.numcellsperdim;
                else if (ix>=opt.numcellsperdim) index+=(ix-opt.numcellsperdim)*opt.numcellsperdim*opt.numcellsperdim;
                else index+=ix*opt.numcellsperdim*opt.numcellsperdim;
                //if ignoring local cells and cell not local, add to cell list
                //or add regardless if not ignoring local cells
                if (opt.cellnodeids[index] != ThisTask) {
                    cellnodeidlist.push_back(opt.cellnodeids[index]);
                }
                //and also check to see any cells that have been newly associated
                // mpi mpi domains. if newcellnodeids has zero size, no cells
                // have been newly associated to mpi domains
                if (opt.newcellnodeids.size() == 0) continue;
                if (opt.newcellnodeids[index].size() == 0) continue;
                for (auto &c:opt.newcellnodeids[index]) {
                    if (c!=ThisTask) cellnodeidlist.push_back(c);
                }
            }
        }
    }
    return cellnodeidlist;
}
//@}

//@{
///Find local particle that originated from foreign swift tasks
#ifdef SWIFTINTERFACE
void MPISwiftExchange(vector<Particle> &Part){
    Int_t nbodies = Part.size();
    Int_t i, j, nexport=0,nimport=0;
    Int_t nsend_local[NProcs],noffset[NProcs],nbuffer[NProcs];
    int sendTask,recvTask;
    int maxchunksize=2147483648/NProcs/sizeof(Particle);
    int nsendchunks,nrecvchunks,numsendrecv;
    int sendoffset,recvoffset;
    int cursendchunksize,currecvchunksize;
    MPI_Status status;
    Particle *PartBufSend, *PartBufRecv;
    for (j=0;j<NProcs;j++) nsend_local[j]=0;
    for (i=0;i<nbodies;i++) {
        if (Part[i].GetSwiftTask() != ThisTask) {
            nexport++;
            nsend_local[Part[i].GetSwiftTask()]++;
        }
    }
    for(j = 1, noffset[0] = 0; j < NProcs; j++) noffset[j]=noffset[j-1] + nsend_local[j-1];
    MPI_Allgather(nsend_local, NProcs, MPI_Int_t, mpi_nsend, NProcs, MPI_Int_t, MPI_COMM_WORLD);
    for (j=0;j<NProcs;j++)nimport+=mpi_nsend[ThisTask+j*NProcs];
    ///\todo need to copy information and see what is what

    if (nexport >0) {
        PartBufSend=new Particle[nexport];
        for (i=0;i<nexport;i++) {
            #ifdef GASON
            Part[i].SetHydroProperties();
            #endif
            #ifdef STARON
            Part[i].SetStarProperties();
            #endif
            #ifdef BHON
            Part[i].SetBHProperties();
            #endif
            #ifdef EXTRADMON
            Part[i].SetExtraDMProperties();
            #endif
            PartBufSend[i]=Part[i+nbodies-nexport];
            PartBufSend[i].SetID(Part[i+nbodies-nexport].GetSwiftTask());
        }
        std::sort(PartBufSend, PartBufSend + nexport, IDCompareVec);
    }
    if (nimport > 0) PartBufRecv = new Particle[nimport];

    //now send the data.
    auto commpair = MPIGenerateCommPairs(mpi_nsend);
    for(auto [task1, task2]:commpair)
    {
        if (ThisTask != task1 && ThisTask != task2) continue;
        auto [sendTask,recvTask] = MPISetSendRecvTask(task1, task2);
        auto [numsendrecv, cursendchunksize, currecvchunksize, sendoffset, recvoffset] = MPIInitialzeCommChunks(
            mpi_nsend[recvTask + sendTask * NProcs], 
            mpi_nsend[sendTask + recvTask * NProcs], 
            maxchunksize);
        for (auto ichunk=0;ichunk<numsendrecv;ichunk++)
        {
            MPI_Sendrecv(&PartBufSend[noffset[recvTask]+sendoffset],
                cursendchunksize * sizeof(Particle), MPI_BYTE,
                recvTask, TAG_SWIFT_A+ichunk,
                &PartBufRecv[nbuffer[recvTask]+recvoffset],
                currecvchunksize * sizeof(Particle), MPI_BYTE, 
                recvTask, TAG_SWIFT_A+ichunk, 
                MPI_COMM_WORLD, &status);
            MPIUpdateCommChunks(mpi_nsend[recvTask + sendTask * NProcs], mpi_nsend[sendTask + recvTask * NProcs], cursendchunksize, currecvchunksize, sendoffset, recvoffset);
        }
    }

    Part.resize(nbodies-nexport+nimport);
    if (nexport > 0) delete[] PartBufSend;
    if (nimport > 0) {
        for (i=0;i<nimport;i++) Part[i+nbodies-nexport]=PartBufRecv[i];
        delete[] PartBufRecv;
    }
}

#endif

//@}

#endif
