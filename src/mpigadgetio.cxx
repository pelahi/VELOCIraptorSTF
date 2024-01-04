/*! \file mpigadgetio.cxx
 *  \brief this file contains routines used with MPI compilation and gadget io and domain construction.
 */

#ifdef USEMPI

//-- For MPI

#include "stf.h"

#include "gadgetitems.h"
#include "endianutils.h"

/// \name Gadget Domain decomposition
//@{

/*!
    Determine the domain decomposition.\n
    Here the domains are constructured in data units
    only ThisTask==0 should call this routine. It is tricky to get appropriate load balancing and correct number of particles per processor.\n

    I could use recursive binary splitting like kd-tree along most spread axis till have appropriate number of volumes corresponding
    to number of processors.

    NOTE: assume that cannot store data so position information is read Nsplit times to determine boundaries of subvolumes
    could also randomly subsample system and produce tree from that
    should store for each processor the node structure generated by the domain decomposition
    what I could do is read file twice, one to get extent and other to calculate entropy then decompose
    along some primary axis, then choose orthogonal axis, keep iterating till have appropriate number of subvolumes
    then store the boundaries of the subvolume. this means I don't store data but get at least reasonable domain decomposition

    NOTE: pkdgrav uses orthoganl recursive bisection along with kd-tree, gadget-2 uses peno-hilbert curve to map particles and oct-trees
    the question with either method is guaranteeing load balancing. for ORB achieved by splitting (sub)volume along a dimension (say one with largest spread or max entropy)
    such that either side of the cut has approximately the same number of particles (ie: median splitting). But for both cases, load balancing requires particle information
    so I must load the system then move particles about to ensure load balancing.

    Main thing first is get the dimensional extent of the system.
    then I could get initial splitting just using mid point between boundaries along each dimension.
    once have that initial splitting just load data then start shifting data around.
*/
void MPIDomainExtentGadget(Options &opt){
    #define SKIP2 Fgad[i].read((char*)&dummy, sizeof(dummy));
    Int_t i,m;
    int dummy;
    char   buf[200];
#ifdef GADGET2FORMAT
    char DATA[5];
#endif
    fstream *Fgad;
    struct gadget_header *header;

    if (ThisTask==0) {
    Fgad=new fstream[opt.num_files];
    header=new gadget_header[opt.num_files];
    for(i=0; i<opt.num_files; i++)
    {
        if(opt.num_files>1) sprintf(buf,"%s.%lld",opt.fname,i);
        else sprintf(buf,"%s",opt.fname);
        Fgad[i].open(buf,ios::in);
        if(!Fgad[i]) {
            LOG(info)<<"can't open file "<<buf;
            exit(0);
        }
        else LOG(info)<<"reading "<<buf;
#ifdef GADGET2FORMAT
        SKIP2;
        Fgad[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
        SKIP2;
        SKIP2;
#endif
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        Fgad[i].read((char*)&header[i], sizeof(gadget_header));
        Fgad[i].read((char*)&dummy, sizeof(dummy));
        //endian indep call
        header[i].Endian();
    }
    for (m=0;m<3;m++) {mpi_xlim[m][0]=0;mpi_xlim[m][1]=header[0].BoxSize;}
#ifdef MPIEXPANDLIM
    for (int j=0;j<3;j++) {
        Double_t dx=0.001*(mpi_xlim[j][1]-mpi_xlim[j][0]);
        mpi_xlim[j][0]-=dx;mpi_xlim[j][1]+=dx;
    }
#endif
    }

    //make sure limits have been found
    MPI_Barrier(MPI_COMM_WORLD);
    if (NProcs==1) {
        for (i=0;i<3;i++) {
            mpi_domain[ThisTask].bnd[i][0]=mpi_xlim[i][0];
            mpi_domain[ThisTask].bnd[i][1]=mpi_xlim[i][1];
        }
    }
}

///to update the decomposition based on gadget information
void MPIDomainDecompositionGadget(Options &opt){
    if (ThisTask==0) {
    }
}

///reads a gadget file to determine number of particles in each MPIDomain
void MPINumInDomainGadget(Options &opt)
{
    #define SKIP2 Fgad[i].read((char*)&dummy, sizeof(dummy));
    InitEndian();
    if (NProcs == 1) return;
    if (opt.cellnodeids.size() == 0) {
        MPIDomainExtentGadget(opt);
        MPIInitialDomainDecomposition(opt);
        MPIDomainDecompositionGadget(opt);
    }
    Int_t i,k,n;
    int dummy;
    FLOAT ctemp[3];
    char   buf[200];
#ifdef GADGET2FORMAT
    char DATA[5];
#endif
    fstream *Fgad;
    struct gadget_header *header;
    Int_t Nlocalold=Nlocal;
    int *ireadtask,*readtaskID;
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    std::vector<int> ireadfile(opt.num_files);
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);

    Int_t ibuf=0,*Nbuf, *Nbaryonbuf;
    Nbuf=new Int_t[NProcs];
    Nbaryonbuf=new Int_t[NProcs];
    for (int j=0;j<NProcs;j++) Nbuf[j]=0;
    for (int j=0;j<NProcs;j++) Nbaryonbuf[j]=0;

    //opening file
    Fgad=new fstream[opt.num_files];
    header=new gadget_header[opt.num_files];
    if (ireadtask[ThisTask]>=0) {
        MPISetFilesRead(opt,ireadfile,ireadtask);
        for(i=0; i<opt.num_files; i++) if(ireadfile[i])
        {
            if(opt.num_files>1) sprintf(buf,"%s.%lld",opt.fname,i);
            else sprintf(buf,"%s",opt.fname);
            Fgad[i].open(buf,ios::in);
		    LOG(info) <<"Reading file "<<buf;
        #ifdef GADGET2FORMAT
            SKIP2;
            Fgad[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
            SKIP2;
            SKIP2;
            LOG(trace) "Reading... %s\n",DATA);
        #endif
            Fgad[i].read((char*)&dummy, sizeof(dummy));
            Fgad[i].read((char*)&header[i], sizeof(gadget_header));
            Fgad[i].read((char*)&dummy, sizeof(dummy));
            //endian indep call
            header[i].Endian();

        #ifdef GADGET2FORMAT
            Fgad[i].read((char*)&dummy, sizeof(dummy));
            Fgad[i].read((char*)&DATA[0],sizeof(char)*4);DATA[4] = '\0';
            Fgad[i].read((char*)&dummy, sizeof(dummy));
            Fgad[i].read((char*)&dummy, sizeof(dummy));
        #endif
            Fgad[i].read((char*)&dummy, sizeof(dummy));
            for(k=0;k<NGTYPE;k++)
            {
            	LOG(trace)<<k<<" "<<header[i].npart[k];
                for(n=0;n<header[i].npart[k];n++)
                {
                    Fgad[i].read((char*)&ctemp[0], sizeof(FLOAT)*3);
#ifdef PERIODWRAPINPUT
                    PeriodWrapInput<FLOAT>(header[i].BoxSize, ctemp);
#endif
                    ibuf=MPIGetParticlesProcessor(opt, ctemp[0],ctemp[1],ctemp[2]);
                    if (opt.partsearchtype==PSTALL) {
                        Nbuf[ibuf]++;
                    }
                    else if (opt.partsearchtype==PSTDARK) {
                        if (!(k==GGASTYPE||k==GSTARTYPE||k==GBHTYPE)) {
                            Nbuf[ibuf]++;
                        }
                        else {
                            if (opt.iBaryonSearch) {
                                Nbaryonbuf[ibuf]++;
                            }
                        }
                    }
                    else if (opt.partsearchtype==PSTSTAR) {
                        if (k==STARTYPE) {
                            Nbuf[ibuf]++;
                        }
                    }
                    else if (opt.partsearchtype==PSTGAS) {
                        if (k==GASTYPE) {
                            Nbuf[ibuf]++;
                        }
                    }
                }
            }
            Fgad[i].close();
        }
    }
    //now having read number of particles, run all gather
    Int_t mpi_nlocal[NProcs];
    MPI_Allreduce(Nbuf,mpi_nlocal,NProcs,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
    Nlocal=mpi_nlocal[ThisTask];
    if (opt.iBaryonSearch) {
        MPI_Allreduce(Nbaryonbuf,mpi_nlocal,NProcs,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
        Nlocalbaryon[0]=mpi_nlocal[ThisTask];
    }
}

//@}

#endif
