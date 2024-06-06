/*! \file main.cxx
 *  \brief Main program

    This code initialize MPI (if necessary) and reads parameters from the user input (see \ref ui.cxx)
    loads particle data from a N-body output (see \ref io.cxx) into an array of \ref NBody::Particle,
    analyzes the system, searches for structures then substructures (see \ref search.cxx),
    outputs the group ids (in two different fashions see \ref io.cxx) and can also analyze the structures
    (see \ref substructureproperties.cxx and \ref io.cxx).

    \todo remove array mpi_idlist and mpi_indexlist as these arrays are unnecessary
    \todo alter unbinding/sortbybindingenergy calls since seems a waste of cpu cycles
*/

#include <iterator>
#include <sstream>
#include <unistd.h>

#include "compilation_info.h"
#include "ioutils.h"
#include "logging.h"
#include "stf.h"
#include "timer.h"

using namespace std;
using namespace Math;
using namespace NBody;
using namespace velociraptor;

void finish_vr(Options &opt)
{
    //get memory useage
    LOG(info) << "Finished running VR";
    MEMORY_USAGE_REPORT(info);
}

void show_version_info(int argc, char *argv[])
{
	LOG_RANK0(info) << "VELOCIraptor git version: " << git_sha1();
	LOG_RANK0(info) << "VELOCIraptor compiled with: " << vr::get_cxx_flags() << " " << vr::get_macros();
	LOG_RANK0(info) << "VELOCIraptor was built on: " << vr::get_compilation_date();

	char cwd[PATH_MAX];
	std::ostringstream os;
	LOG_RANK0(info) << "VELOCIraptor running at: " << getcwd(cwd, PATH_MAX);
	std::copy(argv, argv + argc, std::ostream_iterator<char *>(os, " "));
	LOG_RANK0(info) << "VELOCIraptor started with command line: " << os.str();

    // MPI/OpenMP information
    std::ostringstream mpi_info;
    mpi_info << "VELOCIratptor MPI support: ";
    #ifdef USEMPI
    mpi_info << "yes. ";
    mpi_info << "\n" << "MPI comm world contains "<<NProcs << " MPI ranks";
    char hostname[NAME_MAX + 1];
    ::gethostname(hostname, NAME_MAX);
    std::vector<char> all_hostnames;
    if (ThisTask == 0) {
        all_hostnames.resize((NAME_MAX + 1)* NProcs);
    }
    MPI_Gather(hostname, NAME_MAX + 1, MPI_CHAR, all_hostnames.data(), NAME_MAX + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (ThisTask == 0) {
        std::vector<std::string> proper_hostnames;
        for (int rank = 0; rank != NProcs; rank++) {
            proper_hostnames.emplace_back(all_hostnames.data() + (NAME_MAX + 1) * rank);
        }
        std::set<std::string> s(proper_hostnames.begin(), proper_hostnames.end());
        std::vector<std::string> unique_hostnames(s.size());
        std::copy(s.begin(), s.end(), unique_hostnames.begin());
        mpi_info << " running in " << unique_hostnames.size() << " nodes: ";
        mpi_info << vr::printable_range(unique_hostnames);
    }
    #else
    mpi_info << "no.";
    #endif
	LOG_RANK0(info) << "VELOCIratptor MPI support: " << mpi_info.str();

	LOG_RANK0(info) << "VELOCIratptor OpenMP support: "
#ifdef USEOPENMP
	<< "yes, " << omp_get_max_threads() << " OpenMP threads, OpenMP version " << _OPENMP;
#else
	<< "no";
#endif
    // report binding
    report_binding();
}

int run(int argc,char **argv)
{

#ifdef USEMPI
    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
#else
    int ThisTask=0,NProcs=1;
    Int_t Nlocal,Ntotal;
    Int_t Ntotalbaryon[NBARYONTYPES], Nlocalbaryon[NBARYONTYPES];
#endif
    int nthreads;
#ifdef USEOPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads=1;
#endif

    show_version_info(argc, argv);

#ifdef SWIFTINTERFACE
    cout<<"Built with SWIFT interface enabled when running standalone VELOCIraptor. Should only be enabled when running VELOCIraptor as a library from SWIFT. Exiting..."<<endl;
    exit(0);
#endif

    gsl_set_error_handler_off();

    Options opt;
    //get arguments
    GetArgs(argc, argv, opt);
    cout.precision(10);

#ifdef USEMPI
    MPIInitWriteComm();
#ifdef USEADIOS
    //init adios
    adios_init_noxml(MPI_COMM_WORLD);
    //specify the buffer size (use the tot particle buffer size in bytes and convert to MB)
    adios_set_max_buffer_size(opt.mpiparticletotbufsize/1024/1024);
#endif
#endif

    //variables
    //number of particles, (also number of baryons if use dm+baryon search)
    //to store (point to) particle data
    Int_t nbodies,nbaryons,ndark;
    vector<Particle> Part;
    Particle *Pbaryons;
    KDTree *tree;

    //number in subset, number of grids used if iSingleHalo==0;
    Int_t nsubset, ngrid;
    //number of groups, temp num groups, and number of halos
    Int_t ngroup, ng, nhalos;

    //to store group value (pfof), and also arrays to parse particles
    //vector<Int_t> pfof, pfofbaryons;
    Int_t *pfof, *pfofbaryons;;
    Int_t *numingroup,**pglist;
    Int_t *pfofall;
    //to store information about the group
    PropData *pdata=NULL,*pdatahalos=NULL;

    //to store time and output time taken
    vr::Timer total_timer;

    Coordinate cm,cmvel;
    Double_t Mtot;
    char fname1[1000],fname4[1000];

#ifdef USEMPI
    mpi_nlocal=new Int_t[NProcs];
    mpi_domain=new MPI_Domain[NProcs];
    mpi_nsend=new Int_t[NProcs*NProcs];
    mpi_ngroups=new Int_t[NProcs];
    mpi_nhalos=new Int_t[NProcs];
    //store MinSize as when using mpi prior to stitching use min of 2;
    MinNumMPI=2;
    //if single halo, use minsize to initialize the old minimum number
    //else use the halominsize since if mpi and not single halo, halos localized to mpi domain for substructure search
    if (opt.iSingleHalo) MinNumOld=opt.MinSize;
    else MinNumOld=opt.HaloMinSize;
#endif

    //read particle information and allocate memory
    vr::Timer loading_timer;
    //for MPI determine total number of particles AND the number of particles assigned to each processor
    if (ThisTask==0) {
        LOG(info) << "Reading header";
        nbodies=ReadHeader(opt);
        if (opt.iBaryonSearch>0) {
            for (int i=0;i<NBARYONTYPES;i++) Ntotalbaryon[i]=Nlocalbaryon[i]=0;
            nbaryons=0;
            int pstemp=opt.partsearchtype;
            opt.partsearchtype=PSTGAS;
            nbaryons+=ReadHeader(opt);
            if (opt.iusestarparticles) {
                opt.partsearchtype=PSTSTAR;
                nbaryons+=ReadHeader(opt);
            }
            if (opt.iusesinkparticles) {
                opt.partsearchtype=PSTBH;
                nbaryons+=ReadHeader(opt);
            }
            opt.partsearchtype=pstemp;
        }
        else nbaryons=0;
    }

#ifdef USEMPI
    MPI_Bcast(&nbodies,1, MPI_Int_t,0,MPI_COMM_WORLD);
    //if MPI, possible desired particle types not present in file so update used particle types
    MPIUpdateUseParticleTypes(opt);
    if (opt.iBaryonSearch>0) MPI_Bcast(&nbaryons,1, MPI_Int_t,0,MPI_COMM_WORLD);
    //initial estimate need for memory allocation assuming that work balance is not greatly off
#endif
    LOG_RANK0(info) << "There are " << nbodies << " particles in total that require " << vr::memory_amount(nbodies * sizeof(Particle));
    if (opt.iBaryonSearch > 0) {
        LOG_RANK0(info) << "There are " << nbaryons << " baryon particles in total that require " << vr::memory_amount(nbaryons * sizeof(Particle));
    }


    //note that for nonmpi particle array is a contiguous block of memory regardless of whether a separate baryon search is required
#ifndef USEMPI
    Nlocal=nbodies;
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) {
        Part.resize(nbodies+nbaryons);
        Pbaryons=&(Part.data()[nbodies]);
        Nlocalbaryon[0]=nbaryons;
    }
    else {
        Part.resize(nbodies);
        Pbaryons=NULL;
        nbaryons=0;
    }
#else
    //for mpi however, it is not possible to have a simple contiguous block of memory IFF a separate baryon search is required.
    //for the simple reason that the local number of particles changes to ensure large fof groups are local to an mpi domain
    //however, when reading data, it is much simplier to have a contiguous block of memory, sort that memory (if necessary)
    //and then split afterwards the dm particles and the baryons
    if (NProcs==1) {Nlocal=Nmemlocal=nbodies;NExport=NImport=1;}
    else {
#ifdef MPIREDUCEMEM
		LOG_RANK0(trace) << "Running MPI with reduced mem";
        //if allocating reasonable amounts of memory, use MPIREDUCEMEM
        //this determines number of particles in the mpi domains
        MPINumInDomain(opt);
        LOG(info) << "There are " << Nlocal << " particles and have allocated enough memory for "
                  << Nmemlocal << " requiring " << vr::memory_amount(Nmemlocal * sizeof(Particle));
        if (opt.iBaryonSearch > 0) {
            LOG(info) << "There are " << Nlocalbaryon[0] << " baryon particles and have allocated enough memory for "
                      << Nmemlocalbaryon << " requiring " << vr::memory_amount(Nmemlocalbaryon * sizeof(Particle));
        }
#else
        //otherwise just base on total number of particles * some factor and initialise the domains
		LOG_RANK0(trace) << "Running MPI with some extra memory as buffer";
        MPIDomainExtent(opt);
        MPIDomainDecomposition(opt);
        Nlocal=nbodies/NProcs*MPIProcFac;
        Nmemlocal=Nlocal;
        Nlocalbaryon[0]=nbaryons/NProcs*MPIProcFac;
        Nmemlocalbaryon=Nlocalbaryon[0];
        NExport=NImport=Nlocal*MPIExportFac;
        LOG(info) << "Have allocated enough memory for " << Nmemlocal << " requiring "
                  << vr::memory_amount(Nmemlocal * sizeof(Particle));
        if (opt.iBaryonSearch > 0) {
            LOG(info) << "Have allocated enough memory for " << Nmemlocalbaryon << " baryons particles requiring "
                      << vr::memory_amount(Nmemlocalbaryon * sizeof(Particle));
        }
#endif
    }
    LOG(info) << "Will also require additional memory for FOF algorithms and substructure search. "
              << "Largest mem needed for preliminary FOF search. Rough estimate is " << vr::memory_amount(Nlocal * sizeof(Int_tree_t) * 8);
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) {
        Part.resize(Nmemlocal+Nmemlocalbaryon);
        Pbaryons=&(Part.data()[Nlocal]);
        nbaryons=Nlocalbaryon[0];
    }
    else {
        Part.resize(Nmemlocal);
        Pbaryons=NULL;
        nbaryons=0;
    }
#endif

    //now read particle data
    ReadData(opt, Part, nbodies, Pbaryons, nbaryons);
#ifdef USEMPI
    //if mpi and want separate baryon search then once particles are loaded into contigous block of memory and sorted according to type order,
    //allocate memory for baryons
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) {
        Pbaryons=new Particle[Nmemlocalbaryon];
        nbaryons=Nlocalbaryon[0];

        for (Int_t i=0;i<Nlocalbaryon[0];i++) Pbaryons[i]=Part[i+Nlocal];
        Part.resize(Nlocal);
    }
#endif

#ifdef USEMPI
    Ntotal=nbodies;
    nbodies=Nlocal;
    NExport=NImport=Nlocal*MPIExportFac;
    mpi_period=opt.p;
    MPI_Allgather(&nbodies, 1, MPI_Int_t, mpi_nlocal, 1, MPI_Int_t, MPI_COMM_WORLD);
    MPI_Allreduce(&nbodies, &Ntotal, 1, MPI_Int_t, MPI_SUM, MPI_COMM_WORLD);
    LOG(info) << Nlocal << '/' << Ntotal << " particles loaded in " << loading_timer;
#else
    LOG(info) << nbodies << " particles loaded in " << loading_timer;
#endif

    //write out the configuration used by velociraptor having read in the data (as input data can contain cosmological information)
    WriteVELOCIraptorConfig(opt);
    WriteSimulationInfo(opt);
    WriteUnitInfo(opt);

    //set filenames if they have been passed
#ifdef USEMPI
    if (opt.smname!=NULL) sprintf(fname4,"%s.%d",opt.smname,ThisTask);
#else
    if (opt.smname!=NULL) sprintf(fname4,"%s",opt.smname);
#endif

    //read local velocity data or calculate it
    //(and if STRUCDEN flag or HALOONLYDEN is set then only calculate the velocity density function for objects within a structure
    //as found by SearchFullSet)
#if defined (STRUCDEN) || defined (HALOONLYDEN)
#else
    if (opt.iSubSearch==1) {
        vr::Timer timer;
        if(FileExists(fname4)) ReadLocalVelocityDensity(opt, nbodies,Part);
        else  {
            GetVelocityDensity(opt, nbodies, Part.data());
            WriteLocalVelocityDensity(opt, nbodies,Part);
        }
        LOG(info) << "Local velocity density read/analised for " << Nlocal << " particles with "
                  << nthreads << " threads in " << timer;
    }
#endif

    //here adjust Efrac to Omega_cdm/Omega_m from what it was before if baryonic search is separate
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) opt.uinfo.Eratio*=opt.Omega_cdm/opt.Omega_m;

    //From here can either search entire particle array for "Halos" or if a single halo is loaded, then can just search for substructure
    if (!opt.iSingleHalo) {
        vr::Timer timer;
#ifndef USEMPI
        pfof=SearchFullSet(opt,nbodies,Part,ngroup);
        nhalos=ngroup;
#else
        //nbodies=Ntotal;
        ///\todo Communication Buffer size determination and allocation. For example, eventually need something like FoFDataIn = (struct fofdata_in *) CommBuffer;
        ///At the moment just using NExport
        NExport=Nlocal*MPIExportFac;
        //Now when MPI invoked this returns pfof after local linking and linking across and also reorders groups
        //according to size and localizes the particles belong to the same group to the same mpi thread.
        //after this is called Nlocal is adjusted to the local subset where groups are localized to a given mpi thread.
        pfof=SearchFullSet(opt,Nlocal,Part,ngroup);
        nbodies=Nlocal;
        nhalos=ngroup;
#endif
        LOG(info) << "Search over " << nbodies << " with " << nthreads << " took " << timer;
        //if compiled to determine inclusive halo masses, then for simplicity, I assume halo id order NOT rearranged!
        //this is not necessarily true if baryons are searched for separately.
        if (opt.iInclusiveHalo > 0 && opt.iInclusiveHalo < 3) {
            pdatahalos=new PropData[nhalos+1];
            Int_t *numinhalos=BuildNumInGroup(nbodies, nhalos, pfof);
            Int_t *sortvalhalos=new Int_t[nbodies];
            Int_t *originalID=new Int_t[nbodies];
            for (Int_t i=0;i<nbodies;i++) {sortvalhalos[i]=pfof[i]*(pfof[i]>0)+nbodies*(pfof[i]==0);originalID[i]=Part[i].GetID();Part[i].SetID(i);}
            Int_t *noffsethalos=BuildNoffset(nbodies, Part.data(), nhalos, numinhalos, sortvalhalos);
            ///here if inclusive halo flag is 3, then S0 masses are calculated after substructures are found for field objects
            ///and only calculate FOF masses. Otherwise calculate inclusive masses at this moment.
            GetInclusiveMasses(opt, nbodies, Part.data(), nhalos, pfof, numinhalos, pdatahalos, noffsethalos);
            sort(Part.begin(), Part.end(), IDCompareVec);
            delete[] numinhalos;
            delete[] sortvalhalos;
            delete[] noffsethalos;
            for (Int_t i=0;i<nbodies;i++) Part[i].SetID(originalID[i]);
            delete[] originalID;
        }
    }
    else {
        Coordinate *gvel;
        Matrix *gveldisp;
        GridCell *grid;
        ///\todo Scaling is still not MPI compatible
        if (opt.iScaleLengths) ScaleLinkingLengths(opt,nbodies,Part.data(),cm,cmvel,Mtot);
        opt.Ncell=opt.Ncellfac*nbodies;
        //build grid using leaf nodes of tree (which is guaranteed to be adaptive and have maximum number of particles in cell of tree bucket size)
        tree=InitializeTreeGrid(opt,nbodies,Part.data());
        ngrid=tree->GetNumLeafNodes();
        LOG(info) << "Given " << nbodies << " particles, and max cell size of " << opt.Ncell
                  << " there are " << ngrid << " leaf nodes or grid cells, with each node containing ~"
                  << nbodies / ngrid << " particles";
        grid=new GridCell[ngrid];
        //note that after this system is back in original order as tree has been deleted.
        FillTreeGrid(opt, nbodies, ngrid, tree, Part.data(), grid);
        //calculate cell quantities to get mean field
        gvel=GetCellVel(opt,nbodies,Part.data(),ngrid,grid);
        gveldisp=GetCellVelDisp(opt,nbodies,Part.data(),ngrid,grid,gvel);
        opt.HaloSigmaV=0;for (int j=0;j<ngrid;j++) opt.HaloSigmaV+=pow(gveldisp[j].Det(),1./3.);opt.HaloSigmaV/=(double)ngrid;

        //now that have the grid cell volume quantities and local volume density
        //can determine the logarithmic ratio between the particle velocity density and that predicted by the background velocity distribution
        GetDenVRatio(opt,nbodies,Part.data(),ngrid,grid,gvel,gveldisp);
        //WriteDenVRatio(opt,nbodies,Part);
        //and then determine how much of an outlier it is
        nsubset=GetOutliersValues(opt,nbodies,Part.data());
        //save the normalized denvratio and also determine how many particles lie above the threshold.
        //nsubset=WriteOutlierValues(opt, nbodies,Part);
        //Now check if any particles are above the threshold
        if (nsubset==0) {
            LOG(info) << "No particles found above threshold of " << opt.ellthreshold;
            LOG(info) << "Exiting";
            return 0;
        }
        else {
            LOG(info) << nsubset << " particles above threshold of " << opt.ellthreshold << " to be searched";
        }
#ifndef USEMPI
        pfof=SearchSubset(opt,nbodies,nbodies,Part.data(),ngroup);
#else
        //nbodies=Ntotal;
        ///\todo Communication Buffer size determination and allocation. For example, eventually need something like FoFDataIn = (struct fofdata_in *) CommBuffer;
        ///At the moment just using NExport
        NExport=Nlocal*MPIExportFac;
        mpi_foftask=MPISetTaskID(nbodies);

        //Now when MPI invoked this returns pfof after local linking and linking across and also reorders groups
        //according to size and localizes the particles belong to the same group to the same mpi thread.
        //after this is called Nlocal is adjusted to the local subset where groups are localized to a given mpi thread.
        pfof=SearchSubset(opt,Nlocal,Nlocal,Part.data(),ngroup);
        nbodies=Nlocal;
#endif
    }
    if (opt.iSubSearch) {
        LOG(info) << "Searching subset";
        vr::Timer timer;
        //if groups have been found (and localized to single MPI thread) then proceed to search for subsubstructures
        SearchSubSub(opt, nbodies, Part, pfof,ngroup,nhalos, pdatahalos);
        LOG(info) << "Search for substructures " << Nlocal << " with " << nthreads
                  << " threads finished in " << timer;
    }
    pdata=new PropData[ngroup+1];
    //if inclusive halo mass required
    if (opt.iInclusiveHalo > 0 && opt.iInclusiveHalo < 3 && ngroup>0) {
        CopyMasses(opt,nhalos,pdatahalos,pdata);
        delete[] pdatahalos;
    }

    //if only searching initially for dark matter groups, once found, search for associated baryonic structures if requried
    if (opt.iBaryonSearch>0) {
        vr::Timer timer;
        if (opt.partsearchtype==PSTDARK) {
            pfofall=SearchBaryons(opt, nbaryons, Pbaryons, nbodies, Part, pfof, ngroup,nhalos,opt.iseparatefiles,opt.iInclusiveHalo,pdata);
            pfofbaryons=&pfofall[nbodies];
        }
        //if FOF search overall particle types then running sub search over just dm and need to associate baryons to just dm particles must determine number of baryons, sort list, run search, etc
        //but only need to run search if substructure has been searched
        else if (opt.iSubSearch==1){
            nbaryons=0;
            ndark=0;
            for (Int_t i=0;i<nbodies;i++) {
                if (Part[i].GetType()==DARKTYPE)ndark++;
                else nbaryons++;
            }
            Pbaryons=NULL;
            SearchBaryons(opt, nbaryons, Pbaryons, ndark, Part, pfof, ngroup,nhalos,opt.iseparatefiles,opt.iInclusiveHalo,pdata);
        }
        LOG(info) << "Baryon search with " << nthreads << " threads finished in " << timer;
    }

    //get mpi local hierarchy
    Int_t *nsub,*parentgid, *uparentgid,*stype;
    nsub=new Int_t[ngroup+1];
    parentgid=new Int_t[ngroup+1];
    uparentgid=new Int_t[ngroup+1];
    stype=new Int_t[ngroup+1];
    Int_t nhierarchy=GetHierarchy(opt,ngroup,nsub,parentgid,uparentgid,stype);
    CopyHierarchy(opt,pdata,ngroup,nsub,parentgid,uparentgid,stype);

    //if a separate baryon search has been run, now just place all particles together
    if (opt.iBaryonSearch>0 && opt.partsearchtype!=PSTALL) {
        delete[] pfof;
        pfof=&pfofall[0];
        nbodies+=nbaryons;
        Nlocal=nbodies;
    }

    //output results
    //if want to ignore any information regard particles themselves as particle PIDS are meaningless
    //which might be useful for runs where not interested in tracking just halo catalogues (save for
    //approximate methods like PICOLA. Here it writes desired output and exits
    if(opt.ireducedhaloproperties){
        numingroup = BuildNumInGroup(Nlocal, ngroup, pfof);
        CalculateHaloProperties(opt,Nlocal,Part.data(),ngroup,pfof,numingroup,pdata);
        WriteProperties(opt,ngroup,pdata);
        delete[] numingroup;
        delete[] pdata;
        finish_vr(opt);
    }
    if(opt.inoidoutput){
        numingroup=BuildNumInGroup(Nlocal, ngroup, pfof);
        CalculateHaloProperties(opt,Nlocal,Part.data(),ngroup,pfof,numingroup,pdata);
        WriteProperties(opt,ngroup,pdata);
        if (opt.iprofilecalc) WriteProfiles(opt, ngroup, pdata);
        delete[] numingroup;
        delete[] pdata;
        finish_vr(opt);
    }

    //if want a simple tipsy still array listing particles group ids in input order
    if(opt.iwritefof) {
#ifdef USEMPI
        if (ThisTask==0) {
            mpi_pfof=new Int_t[Ntotal];
            //since pfof is a local subset, not all pfof values have been set, thus initialize them to zero.
            for (Int_t i=0;i<Ntotal;i++) mpi_pfof[i]=0;
        }
        MPICollectFOF(Ntotal, pfof);
        if (ThisTask==0) WriteFOF(opt,Ntotal,mpi_pfof);
#else
        WriteFOF(opt,nbodies,pfof);
#endif
    }
    numingroup=BuildNumInGroup(Nlocal, ngroup, pfof);

    //if separate files explicitly save halos, associated baryons, and subhalos separately
    if (opt.iseparatefiles) {
        if (nhalos>0) {
            pglist=SortAccordingtoBindingEnergy(opt,Nlocal,Part.data(),nhalos,pfof,numingroup,pdata);//alters pglist so most bound particles first
            WriteProperties(opt,nhalos,pdata);
            WriteGroupCatalog(opt, nhalos, numingroup, pglist, Part,ngroup-nhalos);
            //if baryons have been searched output related gas baryon catalogue
            if (opt.partsearchtype==PSTALL){
                WriteGroupPartType(opt, nhalos, numingroup, pglist, Part);
            }
            WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype);
            for (Int_t i=1;i<=nhalos;i++) delete[] pglist[i];
            delete[] pglist;
        }
        else {
#ifdef USEMPI
            //if calculating inclusive masses at end, must call SortAccordingtoBindingEnergy if
            //MPI as domain, despite having no groups might need to exchange particles
            if (opt.iInclusiveHalo==3) SortAccordingtoBindingEnergy(opt,Nlocal,Part.data(),nhalos,pfof,numingroup,pdata);
#endif
            WriteGroupCatalog(opt,nhalos,numingroup,NULL,Part);
            WriteHierarchy(opt,nhalos,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype);
            if (opt.partsearchtype==PSTALL){
                WriteGroupPartType(opt, nhalos, numingroup, NULL, Part);
            }
        }
    }
    Int_t indexii=0;
    ng=ngroup;
    //if separate files, alter offsets
    if (opt.iseparatefiles) {
        sprintf(fname1,"%s.sublevels",opt.outname);
        sprintf(opt.outname,"%s",fname1);
        //alter index point to just output sublevels (assumes no reordering and assumes no change in nhalos as a result of unbinding in SubSubSearch)
        indexii=nhalos;
        ng=ngroup-nhalos;
    }

    if (ng>0) {
        pglist=SortAccordingtoBindingEnergy(opt,Nlocal,Part.data(),ng,pfof,&numingroup[indexii],&pdata[indexii],indexii);//alters pglist so most bound particles first
        WriteProperties(opt,ng,&pdata[indexii]);
        WriteGroupCatalog(opt, ng, &numingroup[indexii], pglist, Part);
        if (opt.iseparatefiles) WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype,1);
        else WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype,-1);
        if (opt.partsearchtype==PSTALL){
            WriteGroupPartType(opt, ng, &numingroup[indexii], pglist, Part);
        }
        for (Int_t i=1;i<=ng;i++) delete[] pglist[i];
        delete[] pglist;
    }
    else {
#ifdef USEMPI
        //if calculating inclusive masses at end, must call SortAccordingtoBindingEnergy if
        //MPI as domain, despite having no groups might need to exchange particles
        if (opt.iInclusiveHalo==3) SortAccordingtoBindingEnergy(opt,Nlocal,Part.data(),ng,pfof,numingroup,pdata);
#endif
        WriteProperties(opt,ng,NULL);
        WriteGroupCatalog(opt,ng,&numingroup[indexii],NULL,Part);
        if (opt.iseparatefiles) WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype,1);
        else WriteHierarchy(opt,ngroup,nhierarchy,psldata->nsinlevel,nsub,parentgid,stype,-1);
        if (opt.partsearchtype==PSTALL){
            WriteGroupPartType(opt, ng, &numingroup[indexii], NULL, Part);
        }
    }

    if (opt.iprofilecalc) WriteProfiles(opt, ngroup, pdata);

#ifdef EXTENDEDHALOOUTPUT
    if (opt.iExtendedOutput) WriteExtendedOutput (opt, ngroup, Nlocal, pdata, Part, pfof);
#endif

    delete[] pfof;
    delete[] numingroup;
    delete[] pdata;
    delete psldata;


    delete[] nsub;
    delete[] parentgid;
    delete[] uparentgid;
    delete[] stype;

    LOG(info) << "VELOCIraptor finished in " << total_timer;

    finish_vr(opt);
    return 0;
}


int main(int argc, char *argv[])
{
#ifdef USEMPI
#ifndef USEOPENMP
    MPI_Init(&argc, &argv);
#else
    // if using hybrid then need to check that threads are available and use the correct initialization
    // Each thread will call MPI routines, but these calls will be coordinated to occur only one at a time within a process.
    int required = MPI_THREAD_FUNNELED;
    int provided;
    MPI_Init_thread(&argc, &argv, required, &provided);
    if (provided < required) {
        LOG_RANK0(error) << "This MPI implementation provides insufficient threading support. "
                         << "Required was " << required << " but provided was " << provided;
        MPI_Finalize();
        exit(9);
    }
#endif // USEOPENMP
#endif // USEMPI

    auto do_abort = [](const char *error_message) {
        LOG(error) << error_message;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#else
        exit(1);
#endif
    };

    try {
        run(argc, argv);
    } catch (const std::exception &e) {
        std::ostringstream os;
        os << "Exception while running VR, aborting now: " << e.what();
        do_abort(os.str().c_str());
    } catch (...) {
        do_abort("Unexpected exception while running VR, aborting now");
    }

#ifdef USEMPI
#ifdef USEADIOS
    adios_finalize(ThisTask);
#endif
    MPI_Finalize();
#endif

}
