/*! \file ramsesio.cxx
 *  \brief this file contains routines for ramses snapshot file io
 *
 * \todo need to check if amr file quantity ngrid_current is actually the number of cells in the file as
 * an example fortran code I have been given seems to also use the ngridlevel array, which stores the number of cells
 * at a given resolution level.
 * \todo need to add in ability for multiple read threads and sends between read threads
 *
 *
 * Edited by:    Rodrigo Ca\~nas
 *               rodrigo.canas@icrar.org
 *
 * Last edited:  7 - Jun - 2017
 *
 *
 */

//-- RAMSES SPECIFIC IO

#include "stf.h"

#include "ramsesitems.h"
#include "endianutils.h"


int RAMSES_fortran_read(fstream &F, int &i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)&i,sizeof(int)); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, RAMSESFLOAT &f){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)&f,sizeof(RAMSESFLOAT)); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, int *i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)i,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, unsigned int *i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy)); byteoffset += sizeof(int);
    F.read((char*)i,dummy);               byteoffset += dummy;
    F.read((char*)&dummy, sizeof(dummy)); byteoffset += sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, long long *i){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)i,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, RAMSESFLOAT *f){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    F.read((char*)f,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    return byteoffset;
}
int RAMSES_fortran_read(fstream &F, char *f){
    int dummy,byteoffset=0;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(char);
    F.read((char*)f,dummy); byteoffset+=dummy;
    F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(char);
    return byteoffset;
}
int RAMSES_fortran_skip(fstream &F, int nskips){
    int dummy,byteoffset=0;
    for (int i=0;i<nskips;i++) {
        F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
        F.seekg(dummy,ios::cur); byteoffset+=dummy;
        F.read((char*)&dummy, sizeof(dummy));byteoffset+=sizeof(int);
    }
    return byteoffset;
}

Int_t RAMSES_get_nbodies(char *fname, int ptype, Options &opt)
{
    char buf[2000],buf1[2000],buf2[2000];
    double * dummy_age, * dummy_mass;
    char *dummy_family2;
    int *dummy_family;
    double dmp_mass;
    double OmegaM, OmegaB;
    int totalghost = 0;
    int totalstars = 0;
    int totaldm    = 0;
    int alltotal   = 0;
    int ghoststars = 0;
    string stringbuf;
    int ninputoffset = 0;
    sprintf(buf1,"%s/amr_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/amr_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find AMR data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }
    //if gas searched in some fashion then load amr/hydro data
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
    sprintf(buf1,"%s/hydro_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/hydro_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find Hydro data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }
    }
    sprintf(buf1,"%s/part_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/part_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find Particle data \nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }

    fstream Framses;
    RAMSES_Header ramses_header_info;
    //buffers to load data
    int intbuff[NRAMSESTYPE];
    long long longbuff[NRAMSESTYPE];
    int i,j,k,ireaderror=0;
    Int_t nbodies=0;
    //IntType inttype;
    int dummy;

    int nusetypes,usetypes[NRAMSESTYPE];

    if (ptype==PSTALL) {nusetypes=4;usetypes[0]=RAMSESGASTYPE;usetypes[1]=RAMSESDMTYPE;usetypes[2]=RAMSESSTARTYPE;usetypes[3]=RAMSESBHTYPE;}
    else if (ptype==PSTDARK) {nusetypes=1;usetypes[0]=RAMSESDMTYPE;}
    else if (ptype==PSTGAS) {nusetypes=1;usetypes[0]=RAMSESGASTYPE;}
    else if (ptype==PSTSTAR) {nusetypes=1;usetypes[0]=RAMSESSTARTYPE;}
    else if (ptype==PSTBH) {nusetypes=1;usetypes[0]=RAMSESBHTYPE;}


    for (j = 0; j < NRAMSESTYPE; j++) ramses_header_info.npartTotal[j] = 0;

    //Open the specified file and the specified dataset in the file.
    //first open amr data
    sprintf(buf1,"%s/amr_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/amr_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    Framses.open(buf, ios::binary|ios::in);
    //read header info
    //this is based on a ramsestotipsy fortran code that does not detail everything in the header block but at least its informative. The example has
    /*
    read(10)ncpu
    read(10)ndim
    read(10)nx,ny,nz
    read(10)nlevelmax
    read(10)ngridmax
    read(10)nboundary
    read(10)ngrid_current
    read(10)boxlen
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)
    read(10)msph
    close(10)
    */

    // Number of files
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.num_files, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));
    opt.num_files=ramses_header_info.num_files;

    // Number of dimensions
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.ndim, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    //
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nx, sizeof(int));
    Framses.read((char*)&ramses_header_info.ny, sizeof(int));
    Framses.read((char*)&ramses_header_info.nz, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    // Maximum refinement level
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nlevelmax, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    //
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.ngridmax, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    //
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nboundary, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));

    //this would be number of active grids but ignore for now
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.seekg(dummy,ios::cur);
    Framses.read((char*)&dummy, sizeof(dummy));

    // Boxsize
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.BoxSize, sizeof(RAMSESFLOAT));
    Framses.read((char*)&dummy, sizeof(dummy));

    //now skip 10 blocks
    for (i=0;i<10;i++) {
        Framses.read((char*)&dummy, sizeof(dummy));
        //skip block size given by dummy
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
    }
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.mass[RAMSESGASTYPE], sizeof(RAMSESFLOAT));
    Framses.read((char*)&dummy, sizeof(dummy));

    Framses.close();

    //reopen to get number of amr cells might need to alter to read grid information and what cells have no so-called son cells
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
    for (i=0;i<ramses_header_info.num_files;i++) {
        sprintf(buf1,"%s/amr_%s.out%05d",fname,opt.ramsessnapname,i+1);
        sprintf(buf2,"%s/amr_%s.out",fname,opt.ramsessnapname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Framses.open(buf, ios::binary|ios::in);
        for (j=0;j<6;j++) {
            Framses.read((char*)&dummy, sizeof(dummy));
            Framses.seekg(dummy,ios::cur);
            Framses.read((char*)&dummy, sizeof(dummy));
        }
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&ramses_header_info.npart[RAMSESGASTYPE], sizeof(int));
        Framses.read((char*)&dummy, sizeof(dummy));
        ramses_header_info.npartTotal[RAMSESGASTYPE]+=ramses_header_info.npart[RAMSESGASTYPE];
    }

    //now hydro header data
    sprintf(buf1,"%s/hydro_%s.out00001",fname,opt.ramsessnapname);
    sprintf(buf2,"%s/hydro_%s.out",fname,opt.ramsessnapname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    Framses.open(buf, ios::binary|ios::in);

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.seekg(dummy,ios::cur);
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.nvarh, sizeof(int));
    Framses.read((char*)&dummy, sizeof(dummy));
    for (i=0;i<3;i++) {
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));
    }

    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.read((char*)&ramses_header_info.gamma_index, sizeof(RAMSESFLOAT));
    Framses.read((char*)&dummy, sizeof(dummy));
    Framses.close();
    }


    //
    // Compute Mass of DM particles in RAMSES code units
    //
    fstream Finfo;
    sprintf(buf1,"%s/info_%s.txt", fname,opt.ramsessnapname);
    Finfo.open(buf1, ios::in);
    Finfo>>stringbuf>>stringbuf>>opt.num_files;
    getline(Finfo,stringbuf);//ndim
    getline(Finfo,stringbuf);//lmin
    getline(Finfo,stringbuf);//lmax
    getline(Finfo,stringbuf);//ngridmax
    getline(Finfo,stringbuf);//nstep
    getline(Finfo,stringbuf);//blank
    getline(Finfo,stringbuf);//box
    getline(Finfo,stringbuf);//time
    getline(Finfo,stringbuf);//a
    getline(Finfo,stringbuf);//hubble

    getline(Finfo,stringbuf);//Mysterious mis-reading one line

    Finfo>>stringbuf>>stringbuf>>OmegaM;
    getline(Finfo,stringbuf);
    getline(Finfo,stringbuf);
    getline(Finfo,stringbuf);
    Finfo>>stringbuf>>stringbuf>>OmegaB;
    Finfo.close();
    dmp_mass = 1.0 / (opt.Neff*opt.Neff*opt.Neff) * (OmegaM - OmegaB) / OmegaM;

    //now particle info
    for (i=0;i<ramses_header_info.num_files;i++)
    {
        sprintf(buf1,"%s/part_%s.out%05d",fname,opt.ramsessnapname,i+1);
        sprintf(buf2,"%s/part_%s.out",fname,opt.ramsessnapname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Framses.open(buf, ios::binary|ios::in);

        ramses_header_info.npart[RAMSESDMTYPE]   = 0;
        ramses_header_info.npart[RAMSESSTARTYPE] = 0;
        ramses_header_info.npart[RAMSESSINKTYPE] = 0;

        //number of cpus
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));

        //number of dimensions
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));

        // Total number of LOCAL particles
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&ramses_header_info.npartlocal, sizeof(int));
        Framses.read((char*)&dummy, sizeof(dummy));

        // Random seeds
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));

        // Total number of Stars over all processors
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&ramses_header_info.nstarTotal, sizeof(int));
        Framses.read((char*)&dummy, sizeof(dummy));

        // Total mass of stars
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));

        // Total lost mass of stars
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));

        // Number of sink particles over the whole simulation (all are included in
        // all processors)
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&ramses_header_info.npartTotal[RAMSESSINKTYPE], sizeof(int));
        Framses.read((char*)&dummy, sizeof(dummy));

        //to determine how many particles of each type, need to look at the mass
        // Skip pos, vel, mass
        for (j = 0; j < 6; j++)
        {
            Framses.read((char*)&dummy, sizeof(dummy));
            Framses.seekg(dummy,ios::cur);
            Framses.read((char*)&dummy, sizeof(dummy));
        }
        //allocate memory to store masses and ages
        dummy_mass = new double [ramses_header_info.npartlocal];
        dummy_age  = new double [ramses_header_info.npartlocal];
	dummy_family = new int [ramses_header_info.npartlocal];
	dummy_family2= new char[ramses_header_info.npartlocal];

        // Read Mass
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy_mass[0], dummy);
        Framses.read((char*)&dummy, sizeof(dummy));

        // Skip Id
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));

        // Skip level
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.seekg(dummy,ios::cur);
        Framses.read((char*)&dummy, sizeof(dummy));

	// Read Family
	Framses.read((char*)&dummy, sizeof(dummy));
	Framses.read((char*)&dummy_family2[0], dummy);
	Framses.read((char*)&dummy, sizeof(dummy));

	for(int ii2=0; ii2<ramses_header_info.npartlocal; ii2++) dummy_family[ii2] = dummy_family2[ii2] + 0;

	// Skip Tag
	Framses.read((char*)&dummy, sizeof(dummy));
	Framses.seekg(dummy,ios::cur);
	Framses.read((char*)&dummy, sizeof(dummy));

        // Read Birth epoch
        //necessary to separate ghost star particles with negative ages from real one
        Framses.read((char*)&dummy, sizeof(dummy));
        Framses.read((char*)&dummy_age[0], dummy);
        Framses.read((char*)&dummy, sizeof(dummy));

        ghoststars = 0;
        for (j = 0; j < ramses_header_info.npartlocal; j++)
        {
            if ((fabs((dummy_mass[j]-dmp_mass)/dmp_mass) < 1e-5) && (dummy_family[j]==1))
                ramses_header_info.npart[RAMSESDMTYPE]++;
	    if (dummy_family[j] == 2) ramses_header_info.npart[RAMSESSTARTYPE] ++;
	    if (dummy_family[j] == 3) ramses_header_info.npart[RAMSESSINKTYPE] ++;
            //else
            //    if (dummy_age[j] != 0.0)
            //        ramses_header_info.npart[RAMSESSTARTYPE]++;
            //    else
            //    ghoststars++;
        }
        delete [] dummy_age;
        delete [] dummy_mass;
	delete [] dummy_family;
	delete [] dummy_family2;
        Framses.close();

        totalghost += ghoststars;
        totalstars += ramses_header_info.npart[RAMSESSTARTYPE];
        totaldm    += ramses_header_info.npart[RAMSESDMTYPE];
        alltotal   += ramses_header_info.npartlocal;

        //now with information loaded, set totals
        ramses_header_info.npartTotal[RAMSESDMTYPE]+=ramses_header_info.npart[RAMSESDMTYPE];
        ramses_header_info.npartTotal[RAMSESSTARTYPE]+=ramses_header_info.npart[RAMSESSTARTYPE];
        ramses_header_info.npartTotal[RAMSESSINKTYPE]+=ramses_header_info.npart[RAMSESSINKTYPE];
    }
    for(j=0, nbodies=0; j<nusetypes; j++) {
        k=usetypes[j];
        nbodies+=ramses_header_info.npartTotal[k];
        //nbodies+=((long long)(ramses_header_info.npartTotalHW[k]) << 32);
    }

    for (j=0;j<NPARTTYPES;j++) opt.numpart[j]=0;
    if (ptype==PSTALL || ptype==PSTDARK) opt.numpart[DARKTYPE]=ramses_header_info.npartTotal[RAMSESDMTYPE];
    if (ptype==PSTALL || ptype==PSTGAS) opt.numpart[GASTYPE]=ramses_header_info.npartTotal[RAMSESGASTYPE];
    if (ptype==PSTALL || ptype==PSTSTAR) opt.numpart[STARTYPE]=ramses_header_info.npartTotal[RAMSESSTARTYPE];
    if (ptype==PSTALL || ptype==PSTBH) opt.numpart[BHTYPE]=ramses_header_info.npartTotal[RAMSESSINKTYPE];
    return nbodies;

}

/// Reads a ramses file. If cosmological simulation uses cosmology (generally
/// assuming LCDM or small deviations from this) to estimate the mean interparticle
/// spacing and scales physical linking length passed by this distance. Also reads
/// header and overrides passed cosmological parameters with ones stored in header.
void ReadRamses(Options &opt, vector<Particle> &Part, const Int_t nbodies, Particle *&Pbaryons, Int_t nbaryons)
{
    char buf[2000],buf1[2000],buf2[2000];
    string stringbuf,orderingstring;
    fstream Finfo;
    fstream *Famr;
    fstream *Fhydro;
    fstream *Fpart, *Fpartvel,*Fpartid,*Fpartmass, *Fpartlevel, *Fpartage, *Fpartmet, *Fpartfam;
    RAMSES_Header *header;
    int intbuff[NRAMSESTYPE];
    long long longbuff[NRAMSESTYPE];
    int i,j,k,n,idim,ivar,igrid,ireaderror=0;
    Int_t count2,bcount2;
    //IntType inttype;
    int dummy,byteoffset;
    Double_t MP_DM=MAXVALUE,LN,N_DM,MP_B=0;
    double z,aadjust,Hubble,Hubbleflow;
    Double_t mscale,lscale,lvscale,rhoscale;
    Double_t mtemp,utemp,rhotemp,Ztemp,Ttemp;
    Coordinate xpos,vpos;
    RAMSESFLOAT xtemp[3],vtemp[3];
    RAMSESIDTYPE idval;
    int typeval;
    RAMSESFLOAT ageval,metval;
    int *ngridlevel,*ngridbound,*ngridfile;
    int lmin=1000000,lmax=0;
    double dmp_mass;

    char *famtempchunk2;
    int *famtempchunk;

    int ninputoffset = 0;
    int ifirstfile=0,*ireadfile,ibuf=0;
    Int_t ibufindex;
    int *ireadtask,*readtaskID;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
    ireadfile=new int[opt.num_files];
    for (i=0;i<opt.num_files;i++) ireadfile[i]=1;
#else
    MPI_Bcast (&(opt.num_files), sizeof(opt.num_files), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);
#endif
    ///\todo because of the stupid fortran format, easier if chunksize is BIG so that
    ///number of particles local to a file are smaller
    Int_t chunksize=RAMSESCHUNKSIZE,nchunk;
    RAMSESFLOAT *xtempchunk, *vtempchunk, *mtempchunk, *sphtempchunk, *agetempchunk, *mettempchunk, *hydrotempchunk;
    RAMSESIDTYPE *idvalchunk, *levelchunk;
    int *icellchunk;

    Famr       = new fstream[opt.num_files];
    Fhydro     = new fstream[opt.num_files];
    Fpart      = new fstream[opt.num_files];
    Fpartvel   = new fstream[opt.num_files];
    Fpartmass  = new fstream[opt.num_files];
    Fpartid    = new fstream[opt.num_files];
    Fpartlevel = new fstream[opt.num_files];
    Fpartage   = new fstream[opt.num_files];
    Fpartmet   = new fstream[opt.num_files];
    Fpartfam   = new fstream[opt.num_files];
    header     = new RAMSES_Header[opt.num_files];

    Particle *Pbuf;
#ifdef USEMPI
    MPI_Status status;
    MPI_Comm mpi_comm_read;

    vector<Particle> *Preadbuf;
    //for parallel io
    Int_t BufSize=opt.mpiparticlebufsize;
    Int_t Nlocalbuf,*Nbuf, *Nreadbuf,*nreadoffset;
    Int_t *Nlocalthreadbuf,Nlocaltotalbuf;
    int *irecv, sendTask,recvTask,irecvflag, *mpi_irecvflag;
    MPI_Request *mpi_request;
    Int_t inreadsend,totreadsend;
    Int_t *mpi_nsend_baryon;
    Int_t *mpi_nsend_readthread;
    Int_t *mpi_nsend_readthread_baryon;
    if (opt.iBaryonSearch) mpi_nsend_baryon=new Int_t[NProcs*NProcs];
    if (opt.nsnapread>1) {
        mpi_nsend_readthread=new Int_t[opt.nsnapread*opt.nsnapread];
        if (opt.iBaryonSearch) mpi_nsend_readthread_baryon=new Int_t[opt.nsnapread*opt.nsnapread];
    }

    Nbuf=new Int_t[NProcs];
    nreadoffset=new Int_t[opt.nsnapread];
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);
    MPI_Comm_split(MPI_COMM_WORLD, (ireadtask[ThisTask]>=0), ThisTask, &mpi_comm_read);
    if (ireadtask[ThisTask]>=0)
    {
        //to temporarily store data
        Pbuf=new Particle[BufSize*NProcs];
        Nreadbuf=new Int_t[opt.nsnapread];
        for (int j=0;j<NProcs;j++) Nbuf[j]=0;
        for (int j=0;j<opt.nsnapread;j++) Nreadbuf[j]=0;
        if (opt.nsnapread>1){
            Preadbuf=new vector<Particle>[opt.nsnapread];
            for (int j=0;j<opt.nsnapread;j++) Preadbuf[j].reserve(BufSize);
        }
        //to determine which files the thread should read
        ireadfile=new int[opt.num_files];
        ifirstfile=MPISetFilesRead(opt,ireadfile,ireadtask);
        inreadsend=0;
        for (int j=0;j<opt.num_files;j++) inreadsend+=ireadfile[j];
        MPI_Allreduce(&inreadsend,&totreadsend,1,MPI_Int_t,MPI_MIN,mpi_comm_read);
    }
    else {
        Nlocalthreadbuf=new Int_t[opt.nsnapread];
        irecv=new int[opt.nsnapread];
        mpi_irecvflag=new int[opt.nsnapread];
        for (i=0;i<opt.nsnapread;i++) irecv[i]=1;
        mpi_request=new MPI_Request[opt.nsnapread];
    }
    Nlocal=0;
    if (opt.iBaryonSearch) Nlocalbaryon[0]=0;

    if (ireadtask[ThisTask]>=0) {
#endif

    //first read cosmological information
    sprintf(buf1,"%s/info_%s.txt",opt.fname,opt.ramsessnapname);
    Finfo.open(buf1, ios::in);
    getline(Finfo,stringbuf);//nfiles
    getline(Finfo,stringbuf);//ndim
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].levelmin;
    getline(Finfo,stringbuf);//lmax
    getline(Finfo,stringbuf);//ngridmax
    getline(Finfo,stringbuf);//nstep
    getline(Finfo,stringbuf);//blank

    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].BoxSize;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].time;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].aexp;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].HubbleParam;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegam;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].OmegaLambda;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegak;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].Omegab;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_l;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_d;
    Finfo>>stringbuf>>stringbuf>>header[ifirstfile].scale_t;

    //convert boxsize to comoving kpc/h
    header[ifirstfile].BoxSize*=header[ifirstfile].scale_l/3.086e21/header[ifirstfile].aexp*header[ifirstfile].HubbleParam/100.0;
    getline(Finfo,stringbuf);
    Finfo>>stringbuf>>orderingstring;
    getline(Finfo,stringbuf);
    Finfo.close();

    opt.p            = header[ifirstfile].BoxSize;
    opt.a            = header[ifirstfile].aexp;
    opt.Omega_m      = header[ifirstfile].Omegam;
    opt.Omega_Lambda = header[ifirstfile].OmegaLambda;
    opt.Omega_b      = header[ifirstfile].Omegab;
    opt.h            = header[ifirstfile].HubbleParam/100.0;
    opt.Omega_cdm    = opt.Omega_m-opt.Omega_b;
    //set hubble unit to km/s/kpc
    opt.H = 0.1;
    //set Gravity to value for kpc (km/s)^2 / solar mass
    opt.G = 4.30211349e-6;
    //and for now fix the units
    opt.lengthtokpc=opt.velocitytokms=opt.masstosolarmass=1.0;

    //Hubble flow
    if (opt.comove) aadjust=1.0;
    else aadjust=opt.a;
    CalcOmegak(opt);
    Hubble=GetHubble(opt, aadjust);
    CalcCriticalDensity(opt, aadjust);
    CalcBackgroundDensity(opt, aadjust);
    CalcVirBN98(opt,aadjust);
    //if opt.virlevel<0, then use virial overdensity based on Bryan and Norman 1997 virialization level is given by
    if (opt.virlevel<0) opt.virlevel=opt.virBN98;
    PrintCosmology(opt);

    //adjust length scale so that convert from 0 to 1 (box units) to kpc comoving
    //to scale mpi domains correctly need to store in opt.lengthinputconversion the box size in comoving little h value
    //opt.lengthinputconversion= opt.p*opt.h/opt.a;
    opt.lengthinputconversion = header[ifirstfile].BoxSize;
    //adjust velocity scale to that ramses is converted to km/s from which you can convert again;
    opt.velocityinputconversion = header[ifirstfile].scale_l/header[ifirstfile].scale_t*1e-5;

    //convert mass from code units to Solar Masses
    mscale   = header[ifirstfile].scale_d * (header[ifirstfile].scale_l * header[ifirstfile].scale_l * header[ifirstfile].scale_l) / 1.988e+33;

    //convert length from code units to kpc (this lscale includes the physical box size)
    lscale   = header[ifirstfile].scale_l/3.086e21;

    //convert velocity from code units to km/s
    lvscale  = header[ifirstfile].scale_l/header[ifirstfile].scale_t*1e-5;

    //convert density to Msun/kpc^3
    rhoscale = mscale/(lscale*lscale*lscale);

    //ignore hubble flow
    Hubbleflow=0.;

    //for (int j=0;j<NPARTTYPES;j++) nbodies+=opt.numpart[j];
    cout<<"Particle system contains "<<nbodies<<" particles (of interest) at is at time "<<opt.a<<" in a box of size "<<opt.p<<endl;

    //number of DM particles
    //NOTE: this assumes a uniform box resolution. However this is not used in the rest of this function
    N_DM = opt.Neff*opt.Neff*opt.Neff;

    //interparticle spacing (assuming a uniform resolution box)
    LN   = (lscale/(double)opt.Neff);
    opt.ellxscale = LN;

    //grab from the first particle file the dimensions of the arrays and also the number of cpus (should be number of files)
    sprintf(buf1,"%s/part_%s.out00001",opt.fname,opt.ramsessnapname);
    Fpart[ifirstfile].open(buf1, ios::binary|ios::in);
    RAMSES_fortran_read(Fpart[ifirstfile],header[ifirstfile].nfiles);
    RAMSES_fortran_read(Fpart[ifirstfile],header[ifirstfile].ndim);
    //adjust the number of files
    opt.num_files=header[ifirstfile].nfiles;
    Fpart[ifirstfile].close();
#ifdef USEMPI
    //now read tasks prepped and can read files to send information
    }
#endif

#ifdef USEMPI
    if (ireadtask[ThisTask]>=0)
    {
#endif
      dmp_mass = 1.0 / (opt.Neff*opt.Neff*opt.Neff) * opt.Omega_cdm / opt.Omega_m;
#ifdef USEMPI
    }
    MPI_Bcast (&dmp_mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    //if not only gas being searched open particle data
    count2=bcount2=0;
    if (opt.partsearchtype!=PSTGAS) {
#ifdef USEMPI
    if (ireadtask[ThisTask]>=0) {
        inreadsend=0;
#endif
    //read particle files consists of positions,velocities, mass, id, and level (along with ages and met if some flags set)
    //
    for (i=0;i<opt.num_files;i++) {
    if (ireadfile[i]) {
        sprintf(buf1,"%s/part_%s.out%05d",opt.fname,opt.ramsessnapname,i+1);
        sprintf(buf2,"%s/part_%s.out",opt.fname,opt.ramsessnapname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Fpart[i].open(buf, ios::binary|ios::in);
        Fpartvel[i].open(buf, ios::binary|ios::in);
        Fpartmass[i].open(buf, ios::binary|ios::in);
        Fpartid[i].open(buf, ios::binary|ios::in);
        Fpartlevel[i].open(buf, ios::binary|ios::in);
        Fpartage[i].open(buf, ios::binary|ios::in);
        Fpartmet[i].open(buf, ios::binary|ios::in);
        Fpartfam[i].open(buf, ios::binary|ios::in);

        //skip header information in each file save for number in the file
        //@{
        byteoffset=0;
        // ncpus
        byteoffset+=RAMSES_fortran_skip(Fpart[i]);
        // ndims
        byteoffset+=RAMSES_fortran_skip(Fpart[i]);
        // store number of particles locally in file
        byteoffset+=RAMSES_fortran_read(Fpart[i],header[i].npartlocal);
        // skip local seeds, nstartot, mstartot, mstarlost, nsink
        byteoffset+=RAMSES_fortran_skip(Fpart[i],5);

        //byteoffset now stores size of header offset for particles
        Fpartvel[i].seekg(byteoffset,ios::cur);
        Fpartmass[i].seekg(byteoffset,ios::cur);
        Fpartid[i].seekg(byteoffset,ios::cur);
        Fpartlevel[i].seekg(byteoffset,ios::cur);
        Fpartage[i].seekg(byteoffset,ios::cur);
        Fpartmet[i].seekg(byteoffset,ios::cur);
	Fpartfam[i].seekg(byteoffset,ios::cur);
        //skip positions
        for(idim=0;idim<header[ifirstfile].ndim;idim++)
        {
            RAMSES_fortran_skip(Fpartvel[i]);
            RAMSES_fortran_skip(Fpartmass[i]);
            RAMSES_fortran_skip(Fpartid[i]);
            RAMSES_fortran_skip(Fpartlevel[i]);
            RAMSES_fortran_skip(Fpartage[i]);
            RAMSES_fortran_skip(Fpartmet[i]);
            RAMSES_fortran_skip(Fpartfam[i]);
        }
        //skip velocities
        for(idim=0;idim<header[ifirstfile].ndim;idim++)
        {
            RAMSES_fortran_skip(Fpartmass[i]);
            RAMSES_fortran_skip(Fpartid[i]);
            RAMSES_fortran_skip(Fpartlevel[i]);
            RAMSES_fortran_skip(Fpartage[i]);
            RAMSES_fortran_skip(Fpartmet[i]);
            RAMSES_fortran_skip(Fpartfam[i]);
        }
        //skip mass
        RAMSES_fortran_skip(Fpartid[i]);
        RAMSES_fortran_skip(Fpartlevel[i]);
        RAMSES_fortran_skip(Fpartage[i]);
        RAMSES_fortran_skip(Fpartmet[i]);
        RAMSES_fortran_skip(Fpartfam[i]);
        //skip ids;
        RAMSES_fortran_skip(Fpartlevel[i]);
        RAMSES_fortran_skip(Fpartage[i]);
        RAMSES_fortran_skip(Fpartmet[i]);
        RAMSES_fortran_skip(Fpartfam[i]);
        //skip levels
        RAMSES_fortran_skip(Fpartage[i]);
        RAMSES_fortran_skip(Fpartmet[i]);
        RAMSES_fortran_skip(Fpartfam[i]);
	//skip family
	RAMSES_fortran_skip(Fpartage[i]);
	RAMSES_fortran_skip(Fpartmet[i]);
	//skip tag
	RAMSES_fortran_skip(Fpartage[i]);
	RAMSES_fortran_skip(Fpartmet[i]);
        //skip ages
        RAMSES_fortran_skip(Fpartmet[i]);
        //@}
        //
        //data loaded into memory in chunks
        chunksize    = nchunk = header[i].npartlocal;
        ninputoffset = 0;
        xtempchunk   = new RAMSESFLOAT  [3*chunksize];
        vtempchunk   = new RAMSESFLOAT  [3*chunksize];
        mtempchunk   = new RAMSESFLOAT  [chunksize];
        idvalchunk   = new RAMSESIDTYPE [chunksize];
        levelchunk   = new RAMSESIDTYPE [chunksize];
        agetempchunk = new RAMSESFLOAT  [chunksize];
        mettempchunk = new RAMSESFLOAT  [chunksize];
	famtempchunk = new int [chunksize];
	famtempchunk2= new char [chunksize];
        for(idim=0;idim<header[ifirstfile].ndim;idim++)
        {
            RAMSES_fortran_read(Fpart[i],&xtempchunk[idim*nchunk]);
            RAMSES_fortran_read(Fpartvel[i],&vtempchunk[idim*nchunk]);
        }
        RAMSES_fortran_read(Fpartmass[i],  mtempchunk);
        RAMSES_fortran_read(Fpartid[i],    idvalchunk);
        RAMSES_fortran_read(Fpartlevel[i], levelchunk);
        RAMSES_fortran_read(Fpartage[i],   agetempchunk);
        RAMSES_fortran_read(Fpartmet[i],   mettempchunk);
	RAMSES_fortran_read(Fpartfam[i],   famtempchunk2);
	for(int nn=0; nn<chunksize; nn++) famtempchunk[nn]=famtempchunk2[nn]+0;
        //RAMSES_fortran_read(Fpartid[i],idvalchunk);
        for (int nn=0;nn<nchunk;nn++)
        {
            if (fabs((mtempchunk[nn]-dmp_mass)/dmp_mass) > 1e-5 && (agetempchunk[nn] == 0.0))
            {
              //  GHOST PARTIRCLE!!!
            }
            else
            {
                xtemp[0] = xtempchunk[nn];
                xtemp[1] = xtempchunk[nn+nchunk];
                xtemp[2] = xtempchunk[nn+2*nchunk];

                vtemp[0] = vtempchunk[nn];
                vtemp[1] = vtempchunk[nn+nchunk];
                vtemp[2] = vtempchunk[nn+2*nchunk];

                idval = idvalchunk[nn];

                ///Need to check this for correct 'endianness'
//             for (int kk=0;kk<3;kk++) {xtemp[kk]=LittleRAMSESFLOAT(xtemp[kk]);vtemp[kk]=LittleRAMSESFLOAT(vtemp[kk]);}
#ifndef NOMASS
            mtemp=mtempchunk[nn];
#else
            mtemp=1.0;
#endif
            ageval = agetempchunk[nn];
	    typeval = 100;

            if (fabs((mtemp-dmp_mass)/dmp_mass) < 1e-5) typeval = DARKTYPE;
	    if (famtempchunk[nn]==2) typeval = STARTYPE;
	    if (famtempchunk[nn]==3) typeval = BHTYPE;
            //else typeval = STARTYPE;
/*
            if (ageval==0 && idval>0) typeval=DARKTYPE;
            else if (idval>0) typeval=STARTYPE;
            else typeval=BHTYPE;
*/
#ifdef USEMPI
            //determine processor this particle belongs on based on its spatial position
            if(opt.partsearchtype==PSTALL){
	    	ibuf=MPIGetParticlesProcessor(opt, xtemp[0],xtemp[1],xtemp[2]);
		ibufindex=ibuf*BufSize+Nbuf[ibuf];
	    }
	    else{
		if((opt.partsearchtype==PSTDARK && typeval==DARKTYPE) ||
		(opt.partsearchtype==PSTSTAR && typeval==STARTYPE))
		{
	    	  ibuf=MPIGetParticlesProcessor(opt, xtemp[0],xtemp[1],xtemp[2]);
		  ibufindex=ibuf*BufSize+Nbuf[ibuf];
		}
	    }
#endif
            //reset hydro quantities of buffer
#ifdef USEMPI
#ifdef GASON
            Pbuf[ibufindex].SetU(0);
#ifdef STARON
            Pbuf[ibufindex].SetSFR(0);
            Pbuf[ibufindex].SetZmet(0);
#endif
#endif
#ifdef STARON
            Pbuf[ibufindex].SetZmet(0);
            Pbuf[ibufindex].SetTage(0);
#endif
#ifdef BHON
#endif
#endif

            if (opt.partsearchtype==PSTALL) {
#ifdef USEMPI
                Pbuf[ibufindex]=Particle(mtemp*mscale,
                    xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                    vtemp[0]*opt.velocityinputconversion+Hubbleflow*xtemp[0],
                    vtemp[1]*opt.velocityinputconversion+Hubbleflow*xtemp[1],
                    vtemp[2]*opt.velocityinputconversion+Hubbleflow*xtemp[2],
                    count2,typeval);
                Pbuf[ibufindex].SetPID(idval);
#ifdef EXTRAINPUTINFO
                if (opt.iextendedoutput)
                {
                    Part[ibufindex].SetInputFileID(i);
                    Part[ibufindex].SetInputIndexInFile(nn+ninputoffset);
                }
#endif
                Nbuf[ibuf]++;
                MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
                Part[count2]=Particle(mtemp*mscale,
                    xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                    vtemp[0]*opt.velocityinputconversion+Hubbleflow*xtemp[0],
                    vtemp[1]*opt.velocityinputconversion+Hubbleflow*xtemp[1],
                    vtemp[2]*opt.velocityinputconversion+Hubbleflow*xtemp[2],
                    count2,typeval);
                Part[count2].SetPID(idval);
#ifdef EXTRAINPUTINFO
                if (opt.iextendedoutput)
                {
                  Part[count2].SetInputFileID(i);
                  Part[count2].SetInputIndexInFile(nn+ninputoffset);
                }
#endif
#endif
                count2++;
            }
            else if (opt.partsearchtype==PSTDARK) {
                if (!(typeval==STARTYPE||typeval==BHTYPE)) {
#ifdef USEMPI
                    Pbuf[ibufindex]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.velocityinputconversion+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.velocityinputconversion+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.velocityinputconversion+Hubbleflow*xtemp[2],
                        count2,DARKTYPE);
                    Pbuf[ibufindex].SetPID(idval);
#ifdef EXTRAINPUTINFO
                    if (opt.iextendedoutput)
                    {
                        Pbuf[ibufindex].SetInputFileID(i);
                        Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
                    }
#endif
                    //ensure that store number of particles to be sent to other reading threads
                    Nbuf[ibuf]++;
                    MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
                    Part[count2]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.velocityinputconversion+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.velocityinputconversion+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.velocityinputconversion+Hubbleflow*xtemp[2],
                        count2,typeval);
                    Part[count2].SetPID(idval);
#ifdef EXTRAINPUTINFO
                    if (opt.iextendedoutput)
                    {
                        Part[count2].SetInputFileID(i);
                        Part[count2].SetInputIndexInFile(nn+ninputoffset);
                    }
#endif
#endif
                    count2++;
                }
                else if (opt.iBaryonSearch) {
#ifdef USEMPI
                    Pbuf[ibufindex]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.velocityinputconversion+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.velocityinputconversion+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.velocityinputconversion+Hubbleflow*xtemp[2],
                        count2);
                    Pbuf[ibufindex].SetPID(idval);
#ifdef EXTRAINPUTINFO
                    if (opt.iextendedoutput)
                    {
                        Pbuf[ibufindex].SetInputFileID(i);
                        Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
                    }
#endif
                    if (typeval==STARTYPE) Pbuf[ibufindex].SetType(STARTYPE);
                    else if (typeval==BHTYPE) Pbuf[ibufindex].SetType(BHTYPE);
                    //ensure that store number of particles to be sent to the reading threads
                    Nbuf[ibuf]++;
                    if (ibuf==ThisTask) {
                        if (k==RAMSESSTARTYPE) Nlocalbaryon[2]++;
                        else if (k==RAMSESSINKTYPE) Nlocalbaryon[3]++;
                    }
                    MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocalbaryon[0], Pbaryons, Nreadbuf, Preadbuf);
#else
                    Pbaryons[bcount2]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.velocityinputconversion+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.velocityinputconversion+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.velocityinputconversion+Hubbleflow*xtemp[2],
                        count2,typeval);
                    Pbaryons[bcount2].SetPID(idval);
#ifdef EXTRAINPUTINFO
                    if (opt.iextendedoutput)
                    {
                        Part[bcount2].SetInputFileID(i);
                        Part[bcount2].SetInputIndexInFile(nn+ninputoffset);
                    }
#endif
#endif
                    bcount2++;
                }
            }
            else if (opt.partsearchtype==PSTSTAR) {
                if (typeval==STARTYPE) {
#ifdef USEMPI
                    //if using MPI, determine proccessor and place in ibuf, store particle in particle buffer and if buffer full, broadcast data
                    //unless ibuf is 0, then just store locally
                    Pbuf[ibufindex]=Particle(mtemp*mscale,
                        xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                        vtemp[0]*opt.velocityinputconversion+Hubbleflow*xtemp[0],
                        vtemp[1]*opt.velocityinputconversion+Hubbleflow*xtemp[1],
                        vtemp[2]*opt.velocityinputconversion+Hubbleflow*xtemp[2],
                        count2,STARTYPE);
                    //ensure that store number of particles to be sent to the reading threads
                    Pbuf[ibufindex].SetPID(idval);
#ifdef EXTRAINPUTINFO
                    if (opt.iextendedoutput)
                    {
                        Pbuf[ibufindex].SetInputFileID(i);
                        Pbuf[ibufindex].SetInputIndexInFile(nn+ninputoffset);
                    }
#endif
                    Nbuf[ibuf]++;
                    MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
                Part[count2]=Particle(mtemp*mscale,
                    xtemp[0]*lscale,xtemp[1]*lscale,xtemp[2]*lscale,
                    vtemp[0]*opt.velocityinputconversion+Hubbleflow*xtemp[0],
                    vtemp[1]*opt.velocityinputconversion+Hubbleflow*xtemp[1],
                    vtemp[2]*opt.velocityinputconversion+Hubbleflow*xtemp[2],
                    count2,typeval);
                Part[count2].SetPID(idval);
#ifdef EXTRAINPUTINFO
                  if (opt.iextendedoutput)
                  {
                      Part[count2].SetInputFileID(i);
                      Part[count2].SetInputIndexInFile(nn+ninputoffset);
	              }
#endif
#endif
                    count2++;
                }
            }
        }//end of ghost particle check
        }//end of loop over chunk
        delete[] xtempchunk;
        delete[] vtempchunk;
        delete[] mtempchunk;
        delete[] idvalchunk;
        delete[] agetempchunk;
        delete[] levelchunk;
        delete[] mettempchunk;
	delete[] famtempchunk;
	delete[] famtempchunk2;
        Fpart[i].close();
        Fpartvel[i].close();
        Fpartmass[i].close();
        Fpartid[i].close();
        Fpartlevel[i].close();
        Fpartage[i].close();
        Fpartmet[i].close();
	Fpartfam[i].close();

#ifdef USEMPI

        //send information between read threads
        if (opt.nsnapread>1&&inreadsend<totreadsend){
            MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
            MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
            inreadsend++;
            for(ibuf = 0; ibuf < opt.nsnapread; ibuf++) Nreadbuf[ibuf]=0;
        }
#endif
    }//end of whether reading a file
    }//end of loop over file
#ifdef USEMPI
    //once finished reading the file if there are any particles left in the buffer broadcast them
    for(ibuf = 0; ibuf < NProcs; ibuf++) if (ireadtask[ibuf]<0)
    {
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
        if (Nbuf[ibuf]>0) {
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            Nbuf[ibuf]=0;
            //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
        }
    }
    if (opt.nsnapread>1){
        MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
        MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
    }
    }//end of ireadtask[ibuf]>0
#endif
#ifdef USEMPI
    else {
        MPIReceiveParticlesFromReadThreads(opt,Pbuf,Part.data(),readtaskID, irecv, mpi_irecvflag, Nlocalthreadbuf, mpi_request,Pbaryons);
    }
#endif
    }

    //if gas searched in some fashion then load amr/hydro data
    if (opt.partsearchtype==PSTGAS||opt.partsearchtype==PSTALL||(opt.partsearchtype==PSTDARK&&opt.iBaryonSearch)) {
#ifdef USEMPI
    if (ireadtask[ThisTask]>=0) {
    inreadsend=0;
#endif
    for (i=0;i<opt.num_files;i++) if (ireadfile[i]) {
        sprintf(buf1,"%s/amr_%s.out%05d",opt.fname,opt.ramsessnapname,i+1);
        sprintf(buf2,"%s/amr_%s.out",opt.fname,opt.ramsessnapname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Famr[i].open(buf, ios::binary|ios::in);
        sprintf(buf1,"%s/hydro_%s.out%05d",opt.fname,opt.ramsessnapname,i+1);
        sprintf(buf2,"%s/hydro_%s.out",opt.fname,opt.ramsessnapname);
        if (FileExists(buf1)) sprintf(buf,"%s",buf1);
        else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
        Fhydro[i].open(buf, ios::binary|ios::in);
        //read some of the amr header till get to number of cells in current file
        //@{
        byteoffset=0;
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].ndim);
        header[i].twotondim=pow(2,header[i].ndim);
        Famr[i].read((char*)&dummy, sizeof(dummy));
        Famr[i].read((char*)&header[i].nx, sizeof(int));
        Famr[i].read((char*)&header[i].ny, sizeof(int));
        Famr[i].read((char*)&header[i].nz, sizeof(int));
        Famr[i].read((char*)&dummy, sizeof(dummy));
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].nlevelmax);
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].ngridmax);
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].nboundary);
        byteoffset+=RAMSES_fortran_read(Famr[i],header[i].npart[RAMSESGASTYPE]);

        //then skip the rest
        for (j=0;j<14;j++) RAMSES_fortran_skip(Famr[i]);
        if (lmin>header[i].nlevelmax) lmin=header[i].nlevelmax;
        if (lmax<header[i].nlevelmax) lmax=header[i].nlevelmax;
        //@}
        //read header info from hydro files
        //@{
        RAMSES_fortran_skip(Fhydro[i]);
        RAMSES_fortran_read(Fhydro[i],header[i].nvarh);
        RAMSES_fortran_skip(Fhydro[i]);
        RAMSES_fortran_skip(Fhydro[i]);
        RAMSES_fortran_skip(Fhydro[i]);
        RAMSES_fortran_read(Fhydro[i],header[i].gamma_index);
        //@}
    }


    for (i=0;i<opt.num_files;i++) if (ireadfile[i]) {
        //then apparently read ngridlevels, which appears to be an array storing the number of grids at a given level
        ngridlevel=new int[header[i].nlevelmax];
        ngridfile=new int[(1+header[i].nboundary)*header[i].nlevelmax];
        RAMSES_fortran_read(Famr[i],ngridlevel);
        for (j=0;j<header[i].nlevelmax;j++) ngridfile[j]=ngridlevel[j];
        //skip some more
        RAMSES_fortran_skip(Famr[i]);
        //if nboundary>0 then need two skip twice then read ngridbound
        if(header[i].nboundary>0) {
            ngridbound=new int[header[i].nboundary*header[i].nlevelmax];
            RAMSES_fortran_skip(Famr[i]);
            RAMSES_fortran_skip(Famr[i]);
            //ngridbound is an array of some sort but I don't see what it is used for
            RAMSES_fortran_read(Famr[i],ngridbound);
            for (j=0;j<header[i].nlevelmax;j++) ngridfile[header[i].nlevelmax+j]=ngridbound[j];
        }
        //skip some more
        RAMSES_fortran_skip(Famr[i],2);
        //if odering list in info is bisection need to skip more
        if (orderingstring==string("bisection")) RAMSES_fortran_skip(Famr[i],5);
        else RAMSES_fortran_skip(Famr[i],4);

        ninputoffset=0;

        for (k=0;k<header[i].nboundary+1;k++) {
            for (j=0;j<header[i].nlevelmax;j++) {
                //first read amr for positions
                chunksize=nchunk=ngridfile[k*header[i].nlevelmax+j];
                if (chunksize>0) {
                    xtempchunk=new RAMSESFLOAT[3*chunksize];
                    //store son value in icell
                    icellchunk=new int[header[i].twotondim*chunksize];
                    //skip grid index, next index and prev index.
                    RAMSES_fortran_skip(Famr[i],3);
                    //now read grid centre
                    for (idim=0;idim<header[i].ndim;idim++) {
                        RAMSES_fortran_read(Famr[i],&xtempchunk[idim*chunksize]);
                    }
                    //skip father index, then neighbours index
                    RAMSES_fortran_skip(Famr[i],1+2*header[i].ndim);
                    //read son index to determine if a cell in a specific grid is at the highest resolution and needs to be represented by a particle
                    for (idim=0;idim<header[i].twotondim;idim++) {
                        RAMSES_fortran_read(Famr[i],&icellchunk[idim*chunksize]);
                    }
                    //skip cpu map and refinement map (2^ndim*2)
                    RAMSES_fortran_skip(Famr[i],2*header[i].twotondim);
                }
                RAMSES_fortran_skip(Fhydro[i]);
                //then read hydro for other variables (first is density, then velocity, then pressure, then metallicity )
                if (chunksize>0) {
                    hydrotempchunk=new RAMSESFLOAT[chunksize*header[i].twotondim*header[i].nvarh];
                    //first read velocities (for 2 cells per number of dimensions (ie: cell corners?))
                    for (idim=0;idim<header[i].twotondim;idim++) {
                        for (ivar=0;ivar<header[i].nvarh;ivar++) {
                            RAMSES_fortran_read(Fhydro[i],&hydrotempchunk[idim*chunksize*header[i].nvarh+ivar*chunksize]);
                            for (igrid=0;igrid<chunksize;igrid++) {
                                //once we have looped over all the hydro data then can start actually storing it into the particle structures
                                if (ivar==header[i].nvarh-1) {
                                    //if cell has no internal cells or at maximum level produce a particle
                                    if (icellchunk[idim*chunksize+igrid]==0 || j==header[i].nlevelmax-1) {
                                        //first suggestion is to add some jitter to the particle positions
                                        double dx = pow(0.5, j);
                                        int ix, iy, iz;
                                        //below assumes three dimensions with 8 corners (? maybe cells) per grid
                                        iz = idim/4;
                                        iy = (idim - (4*iz))/2;
                                        ix = idim - (2*iy) - (4*iz);
                                        // Calculate absolute coordinates + jitter, and generate particle
                                        xpos[0] = ((((float)rand()/(float)RAND_MAX) * header[i].BoxSize * dx) +(header[i].BoxSize * (xtempchunk[igrid] + (double(ix)-0.5) * dx )) - (header[i].BoxSize*dx/2.0)) ;
                                        xpos[1] = ((((float)rand()/(float)RAND_MAX) * header[i].BoxSize * dx) +(header[i].BoxSize * (xtempchunk[igrid+1*chunksize] + (double(iy)-0.5) * dx )) - (header[i].BoxSize*dx/2.0)) ;
                                        xpos[2] = ((((float)rand()/(float)RAND_MAX) * header[i].BoxSize * dx) +(header[i].BoxSize * (xtempchunk[igrid+2*chunksize] + (double(iz)-0.5) * dx )) - (header[i].BoxSize*dx/2.0)) ;
#ifdef USEMPI
                                        //determine processor this particle belongs on based on its spatial position
                                        ibuf=MPIGetParticlesProcessor(opt, xpos[0],xpos[1],xpos[2]);
                                        ibufindex=ibuf*BufSize+Nbuf[ibuf];
#endif

                                        vpos[0]=hydrotempchunk[idim*chunksize*header[i].nvarh+1*chunksize+igrid];
                                        vpos[1]=hydrotempchunk[idim*chunksize*header[i].nvarh+2*chunksize+igrid];
                                        vpos[2]=hydrotempchunk[idim*chunksize*header[i].nvarh+3*chunksize+igrid];
                                        mtemp=dx*dx*dx*hydrotempchunk[idim*chunksize*header[i].nvarh+0*chunksize+igrid];
                                        //the self energy P/rho is given by
                                        utemp=hydrotempchunk[idim*chunksize*header[i].nvarh+4*chunksize+igrid]/hydrotempchunk[idim*chunksize*header[i].nvarh+0*chunksize+igrid]/(header[i].gamma_index-1.0);
                                        rhotemp=hydrotempchunk[idim*chunksize*header[i].nvarh+0*chunksize+igrid]*rhoscale;
                                        Ztemp=hydrotempchunk[idim*chunksize*header[i].nvarh+5*chunksize+igrid];
                                        if (opt.partsearchtype==PSTALL) {
#ifdef USEMPI
                                            Pbuf[ibufindex]=Particle(mtemp*mscale,
                                                xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
                                                vpos[0]*opt.velocityinputconversion+Hubbleflow*xpos[0],
                                                vpos[1]*opt.velocityinputconversion+Hubbleflow*xpos[1],
                                                vpos[2]*opt.velocityinputconversion+Hubbleflow*xpos[2],
                                                count2,GASTYPE);
                                            Pbuf[ibufindex].SetPID(idval);
#ifdef GASON
                                            Pbuf[ibufindex].SetU(utemp);
                                            Pbuf[ibufindex].SetSPHDen(rhotemp);
#ifdef STARON
                                            Pbuf[ibufindex].SetZmet(Ztemp);
#endif
#endif
#ifdef EXTRAINPUTINFO
                                            if (opt.iextendedoutput)
                                            {
                                                Pbuf[ibufindex].SetInputFileID(i);
                                                Pbuf[ibufindex].SetInputIndexInFile(idim*chunksize*header[i].nvarh+0*chunksize+igrid+ninputoffset);
                                            }
#endif
                                            //ensure that store number of particles to be sent to the threads involved with reading snapshot files
                                            Nbuf[ibuf]++;
                                            MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocal, Part.data(), Nreadbuf, Preadbuf);
#else
                                            Part[count2]=Particle(mtemp*mscale,
                                                xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
                                                vpos[0]*opt.velocityinputconversion+Hubbleflow*xpos[0],
                                                vpos[1]*opt.velocityinputconversion+Hubbleflow*xpos[1],
                                                vpos[2]*opt.velocityinputconversion+Hubbleflow*xpos[2],
                                                count2,GASTYPE);
                                            Part[count2].SetPID(idval);
#ifdef GASON
                                            Part[count2].SetU(utemp);
                                            Part[count2].SetSPHDen(rhotemp);
#ifdef STARON
                                            Part[count2].SetZmet(Ztemp);
#endif
#endif
#ifdef EXTRAINPUTINFO
                                            if (opt.iextendedoutput)
                                            {
                                                Part[count2].SetInputFileID(i);
                                                Part[count2].SetInputIndexInFile(idim*chunksize*header[i].nvarh+0*chunksize+igrid+ninputoffset);
                                            }
#endif

#endif
                                            count2++;
                                        }
                                        else if (opt.partsearchtype==PSTDARK&&opt.iBaryonSearch) {
#ifdef USEMPI
                                            Pbuf[ibufindex]=Particle(mtemp*mscale,
                                                xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
                                                vpos[0]*opt.velocityinputconversion+Hubbleflow*xpos[0],
                                                vpos[1]*opt.velocityinputconversion+Hubbleflow*xpos[1],
                                                vpos[2]*opt.velocityinputconversion+Hubbleflow*xpos[2],
                                                count2,GASTYPE);
                                            Pbuf[ibufindex].SetPID(idval);
#ifdef GASON
                                            Pbuf[ibufindex].SetU(utemp);
                                            Pbuf[ibufindex].SetSPHDen(rhotemp);
#ifdef STARON
                                            Pbuf[ibufindex].SetZmet(Ztemp);
#endif
#endif
#ifdef EXTRAINPUTINFO
                                            if (opt.iextendedoutput)
                                            {
                                                Pbuf[ibufindex].SetInputFileID(i);
                                                Pbuf[ibufindex].SetInputIndexInFile(idim*chunksize*header[i].nvarh+0*chunksize+igrid+ninputoffset);
                                            }
#endif
                                            //ensure that store number of particles to be sent to the reading threads
                                            if (ibuf==ThisTask) {
                                                Nlocalbaryon[1]++;
                                            }
                                            MPIAddParticletoAppropriateBuffer(opt, ibuf, ibufindex, ireadtask, BufSize, Nbuf, Pbuf, Nlocalbaryon[0], Pbaryons, Nreadbuf, Preadbuf);
#else
                                            Pbaryons[bcount2]=Particle(mtemp*mscale,
                                                xpos[0]*lscale,xpos[1]*lscale,xpos[2]*lscale,
                                                vpos[0]*opt.velocityinputconversion+Hubbleflow*xpos[0],
                                                vpos[1]*opt.velocityinputconversion+Hubbleflow*xpos[1],
                                                vpos[2]*opt.velocityinputconversion+Hubbleflow*xpos[2],
                                                count2,GASTYPE);
                                            Pbaryons[bcount2].SetPID(idval);
#ifdef GASON
                                            Pbaryons[bcount2].SetU(utemp);
                                            Pbaryons[bcount2].SetSPHDen(rhotemp);
#ifdef STARON
                                            Pbaryons[bcount2].SetZmet(Ztemp);
#endif
#endif
#ifdef EXTRAINPUTINFO
                                            if (opt.iextendedoutput)
                                            {
                                                Pbaryons[bcount2].SetInputFileID(i);
                                                Pbaryons[bcount2].SetInputIndexInFile(idim*chunksize*header[i].nvarh+0*chunksize+igrid+ninputoffset);
                                            }
#endif

#endif
                                        bcount2++;
                                    }
                                }
                                }
                            }
                        }
                    }
                }
                if (chunksize>0) {
                    delete[] xtempchunk;
                    delete[] hydrotempchunk;
                }
            }
        }
        Famr[i].close();
#ifdef USEMPI
        //send information between read threads
        if (opt.nsnapread>1&&inreadsend<totreadsend){
            MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
            MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
            inreadsend++;
            for(ibuf = 0; ibuf < opt.nsnapread; ibuf++) Nreadbuf[ibuf]=0;
        }
#endif
    }
#ifdef USEMPI
    //once finished reading the file if there are any particles left in the buffer broadcast them
    for(ibuf = 0; ibuf < NProcs; ibuf++) if (ireadtask[ibuf]<0)
    {
        MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t, ibuf, ibuf+NProcs, MPI_COMM_WORLD);
        if (Nbuf[ibuf]>0) {
            MPI_Ssend(&Pbuf[ibuf*BufSize], sizeof(Particle)*Nbuf[ibuf], MPI_BYTE, ibuf, ibuf, MPI_COMM_WORLD);
            MPISendHydroInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
            MPISendStarInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
            MPISendBHInfoFromReadThreads(opt, Nbuf[ibuf], &Pbuf[ibuf*BufSize], ibuf);
            Nbuf[ibuf]=0;
            //last broadcast with Nbuf[ibuf]=0 so that receiver knows no more particles are to be broadcast
            MPI_Ssend(&Nbuf[ibuf],1,MPI_Int_t,ibuf,ibuf+NProcs,MPI_COMM_WORLD);
        }
    }
    if (opt.nsnapread>1){
        MPI_Allgather(Nreadbuf, opt.nsnapread, MPI_Int_t, mpi_nsend_readthread, opt.nsnapread, MPI_Int_t, mpi_comm_read);
        MPISendParticlesBetweenReadThreads(opt, Preadbuf, Part.data(), ireadtask, readtaskID, Pbaryons, mpi_comm_read, mpi_nsend_readthread, mpi_nsend_readthread_baryon);
    }
    }//end of reading task
#endif
#ifdef USEMPI
    //if not reading information than waiting to receive information
    else {
        MPIReceiveParticlesFromReadThreads(opt,Pbuf,Part.data(),readtaskID, irecv, mpi_irecvflag, Nlocalthreadbuf, mpi_request,Pbaryons);
    }
#endif
    }//end of check if gas loaded
	
    //update info
    opt.p*=opt.a/opt.h;
#ifdef HIGHRES
    opt.zoomlowmassdm=MP_DM*mscale;
#endif

#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
    //update cosmological data and boundary in code units
    MPI_Bcast(&(opt.p),sizeof(opt.p),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.a),sizeof(opt.a),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_cdm),sizeof(opt.Omega_cdm),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_b),sizeof(opt.Omega_b),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_m),sizeof(opt.Omega_m),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.Omega_Lambda),sizeof(opt.Omega_Lambda),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.h),sizeof(opt.h),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhocrit),sizeof(opt.rhocrit),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.rhobg),sizeof(opt.rhobg),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virlevel),sizeof(opt.virlevel),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.virBN98),sizeof(opt.virBN98),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.ellxscale),sizeof(opt.ellxscale),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.lengthinputconversion),sizeof(opt.lengthinputconversion),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.velocityinputconversion),sizeof(opt.velocityinputconversion),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.massinputconversion),sizeof(opt.massinputconversion),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&(opt.G),sizeof(opt.G),MPI_BYTE,0,MPI_COMM_WORLD);
#ifdef NOMASS
    MPI_Bcast(&(opt.MassValue),sizeof(opt.MassValue),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    MPI_Bcast(&(Ntotal),sizeof(Ntotal),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&opt.zoomlowmassdm,sizeof(opt.zoomlowmassdm),MPI_BYTE,0,MPI_COMM_WORLD);
#endif
    //store how to convert input internal energies to physical output internal energies
    //as we already convert ramses units to sensible output units, nothing to do.
    opt.internalenergyinputconversion = 1.0;

    //a bit of clean up
#ifdef USEMPI
    MPI_Comm_free(&mpi_comm_read);
    if (opt.iBaryonSearch) delete[] mpi_nsend_baryon;
    if (opt.nsnapread>1) {
        delete[] mpi_nsend_readthread;
        if (opt.iBaryonSearch) delete[] mpi_nsend_readthread_baryon;
        if (ireadtask[ThisTask]>=0) delete[] Preadbuf;
    }
    delete[] Nbuf;
    if (ireadtask[ThisTask]>=0) {
        delete[] Nreadbuf;
        delete[] Pbuf;
        delete[] ireadfile;
    }
    delete[] ireadtask;
    delete[] readtaskID;
#endif
}
