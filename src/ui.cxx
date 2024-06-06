/*! \file ui.cxx
 *  \brief user interface
 */

#include "io.h"
#include "ioutils.h"
#include "logging.h"
#include "stf.h"

static void configure_logging(const Options &opt)
{
	// VR's verbosity goes from 0 (less verbose) to at least 2 (very verbose)
	// We thus consider 2 as "trace", 1 as "debug" and so forth
	constexpr int NUM_LOG_LEVELS = 5;
	constexpr int VR_MIN_LOG_LEVEL = 2;
	int as_log_level = -opt.iverbose + VR_MIN_LOG_LEVEL;
	vr::LogLevel level = static_cast<vr::LogLevel>(std::min(std::max(as_log_level, 0), NUM_LOG_LEVELS - 1));
	vr::init_logging(level);
}

static void log_options_summary(const Options &opt)
{
    LOG(info) << "CONFIG INFO SUMMARY --------------------------";
    const char *search_type;
    switch(opt.partsearchtype)
    {
            case PSTALL:
                search_type = "Searching all particles regardless of type";
                break;
            case PSTDARK:
                search_type = "Searching dark matter particles";
                break;
            case PSTGAS:
                search_type = "Searching gas particles";
                break;
            case PSTSTAR:
                search_type = "Searching star particles";
                break;
            case PSTBH:
                search_type = "Searching BH particles?! Really? Are there enough?";
                break;
            default:
                search_type = "Unknown search type";
    }
    LOG(info) << search_type;

    if (opt.fofbgtype == FOF6D) {
        LOG(info) << "Field objects found with 3d FOF";
    }
    else if (opt.fofbgtype == FOF6D) {
        LOG(info) <<"Field objects found with initial 3d FOF then use single dispersion to find 6d FOFs";
    }
    else if (opt.fofbgtype == FOF6DADAPTIVE) {
        LOG(info) << "Field objects found with initial 3d FOF then use adaptive dispersion to find 6d FOFs";
    }
    else if (opt.fofbgtype == FOFSTNOSUBSET) {
        LOG(info) << "Field objects found with initial 3d FOF then use phase-space stream finding algorithm";
    }
    if (opt.fofbgtype <= FOF6D && opt.iKeepFOF) {
        LOG(info) << "Field objects found with initial 3d FOF then use dispersion to find 6d FOFs and 3dFOF objects kept as base objects (consider them inter halo stellar mass or ICL for stellar only searches) ";
    }
    if (opt.iBaryonSearch == 1) {
        LOG(info) << "Baryons treated separately for substructure search";
    }
    else if (opt.iBaryonSearch == 2) {
        LOG(info) << "Baryons treated separately for substructure search and also FOF field search";
    }
    if (opt.iSingleHalo) {
        LOG(info) << "Field objects NOT searched for, assuming single Halo and subsearch using mean field first step";
    }
    LOG(info) << "Allowed potential to kinetic ratio when unbinding particles: " << opt.uinfo.Eratio;
    if (opt.HaloMinSize != opt.MinSize) {
        LOG(info) << "Field objects (aka Halos) have different minimum required size than substructures: " << opt.HaloMinSize << " vs " << opt.MinSize;
    }
    LOG(info) << "Units: L=" << opt.lengthinputconversion << ", M=" << opt.massinputconversion << ", V=" << opt.velocityinputconversion << ", G=" << opt.G;
    if (opt.ibinaryout) {
        LOG(info) << "Binary output";
    }
    if (opt.iseparatefiles) {
        LOG(info) << "Separate files output";
    }
    if (opt.icomoveunit) {
        LOG(info) << "Converting properties into comoving, little h units";
    }
    if (opt.iextrahalooutput) {
        LOG(info) << "Calculating extra halo properties";
    }
    if (opt.iextragasoutput) {
        LOG(info) << "Calculating extra gas properties";
    }
    if (opt.iextrastaroutput) {
        LOG(info) << "Calculating extra star properties";
    }
    if (opt.iaperturecalc) {
        std::ostringstream os;
        os << "Calculating quantities in apertures: ";
        for (auto apname: opt.aperture_names_kpc) {
            os << apname << ",";
        }
        LOG(info) << os.str();
    }
    if (opt.iprofilecalc) {
        std::ostringstream os;
        os << "Calculating radial profiles: ";
        for (auto binedge: opt.profile_bin_edges) {
            os << binedge << ",";
        }
        LOG(info) << os.str();
    }

    if (opt.iextendedoutput) {
        LOG(info) << "Writing extended output for particle extraction from input files";
    }
    if (opt.iSphericalOverdensityPartList) {
        LOG(info) << "Writing list of particles in SO regions of halos. ";
    }

    if (opt.iHaloCoreSearch) {
        LOG(info) << "Searching for 6dfof cores so as to disentangle mergers";
    }
    if (opt.iHaloCoreSearch && opt.iAdaptiveCoreLinking) {
        LOG(info) << "With adaptive linking lengths";
    }
    if (opt.iHaloCoreSearch && opt.iPhaseCoreGrowth) {
        LOG(info) << "With with phase-space tensor core assignment";
    }

    const char *input_type;
    switch (opt.inputtype)
    {
        case IOGADGET:
            input_type = "Gadget";
            break;
        case IOTIPSY:
            input_type = "Tipsy";
            break;
        case IORAMSES:
            input_type = "RAMSES";
            break;
#ifdef USEHDF
        case IOHDF:
            input_type = "HDF";
            break;
#endif
#ifdef USEXDR
        case IONCHILADA:
            input_type = "NCHILADA";
            break;
#endif
        default:
            input_type = "Unknown";
            break;
    }
    LOG(info) << input_type << " file particle input";

#ifdef USEMPI
    LOG(info) << "MPI input communication requires allocation of temporary particle buffer of size "
              << vr::memory_amount(opt.mpiparticletotbufsize);
#endif
    LOG(info) << "----------------------------------------------";
}

///routine to get arguments from command line
/// \todo alter interface as now need to be able to specify only smdata file (no grid, res, normalized res) and functionality to specify eps, background fof search, etc
void GetArgs(int argc, char *argv[], Options &opt)
{
#if defined(USEMPI) && defined(USEPARALLELHDF)
    opt.mpinprocswritesize=NProcs;
#endif
    int option;
    int NumArgs = 0;
    int configflag=0;
    while ((option = getopt(argc, argv, ":C:I:i:s:Z:o:G:S:B:t:")) != EOF)
    {
        switch(option)
        {
            case 'C':
                opt.pname = optarg;
                configflag=1;
                NumArgs += 2;
                break;
            case 'I':
                opt.inputtype = atoi(optarg);
                NumArgs += 2;
                break;
            case 'i':
                opt.fname = optarg;
                NumArgs += 2;
                break;
            case 's':
                opt.num_files = atoi(optarg);
                NumArgs += 2;
                break;
            case 'Z':
                opt.nsnapread = atoi(optarg);
                NumArgs += 2;
                break;
            case 'o':
                opt.outname = optarg;
                NumArgs += 2;
                break;
            case 'G':
                opt.gnsphblocks = atoi(optarg);
                NumArgs += 2;
                break;
            case 'S':
                opt.gnstarblocks = atoi(optarg);
                NumArgs += 2;
                break;
            case 'B':
                opt.gnbhblocks = atoi(optarg);
                NumArgs += 2;
                break;
            case 't':
                opt.ramsessnapname = optarg;
                NumArgs += 2;
                break;
            case '?':
                usage();
        }
    }
    if(configflag){
        LOG_RANK0(info) << "Reading config file " << opt.pname;
        GetParamFile(opt);
    }
    else {
        LOG_RANK0(warning) << "No configuration file passed, using default values";
    }
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    ConfigCheck(opt);
}

///Outputs the usage to stdout
void usage(void)
{
#ifndef USEMPI
    int ThisTask =0, NProcs =1;
#endif
    Options opt;
    if (ThisTask==0) {
    cerr<<"USAGE:\n";
    cerr<<"\n";
    cerr<<"-C <configuration file (overrides other options)> "<<endl;
    cerr<<"-I <input format [Gadget (Default) "<<IOGADGET<<", HDF (if implemented) "<<IOHDF<<", TIPSY "<<IOTIPSY<<", RAMSES "<<IORAMSES<<", NCHILADA "<<IONCHILADA<<">"<<endl;
    cerr<<"-i <input file> "<<endl;
    cerr<<"-s <number of files per output for gadget input 1 [default]>"<<endl;
    cerr<<"-Z <number of threads used in parallel read ("<<opt.nsnapread<<")>"<<endl;
    cerr<<"-o <output filename>"<<endl;
    cerr<<" ===== EXTRA OPTIONS FOR GADGET INPUT ====== "<<endl;
    cerr<<"-g <number of extra sph/gas blocks for gadget>"<<endl;
    cerr<<"-s <number of extra star blocks for gadget>"<<endl;
    cerr<<"-b <number of extra bh blocks for gadget>"<<endl;
    cerr<<" ===== EXTRA OPTIONS REQUIRED FOR RAMSES INPUT ====== "<<endl;
    cerr<<"-t <ramses snapnumber>"<<endl;

    }
#ifdef USEMPI
    MPI_Finalize();
#endif
    exit(1);
}

/*!
    \page configopt  Configuration options

    An example configuration file can be found the examples directory (see for instance <a href="../../examples/sample.cfg">sample</a>). within the repository.
    Only the keywords listed here will be used, all other words/characters are ignored.
    Note that if misspell a keyword it will not be used. One can check the options used by examining <b><em>foo</em>.configuration</b>,
    where <em>foo</em> is your base output filename \ref Options.outname. This file lists all the options. Since this file is
    always written <b><em>DO NOT</em></b> name your input configuration file <b><em>foo</em>.configuration</b>.

    There are numerous key words that can be passed to change the values stored in \ref Options. Here we list them
    \section outconfig Output related.
    See \ref outputs for types of output (and \ref io.cxx for implementation of the code)

    \arg <b> \e Output </b> Output base name. Overrides the name passed with the command line argument <b> \e -o </b>. Only implemented for completeness. \ref Options.outname \n
    \arg <b> \e Write_group_array_file </b> When producing output also produce a file which lists for every particle the group they belong to. Can be used with \b tipsy format or to tag every particle. \ref Options.iwritefof
    \arg <b> \e Output_den </b> A filename for storing the intermediate step of calculating local densities. This is particularly useful if the code is not compiled with \b STRUCDEN & \b HALOONLYDEN (see \ref STF-makeflags). \ref Options.smname \n

    \section searchconfig Parameters related to search type.
    See \ref io.cxx (and related ios like \ref gadgetio.cxx), \ref search.cxx, \ref fofalgo.h for extra details
    \arg <b> \e Particle_search_type </b> A integer describing what types of particles are searched \ref Options.partsearchtype. (See \ref SEARCHTYPES) > \n
        - \b 1 \e all particles are searched
        - \b 2 \e DarkMatter particles (which are typically defined as type 1,2,3 for gadget) are searched
        - \b 3 \e Star particles (which are typically defined as type 4 for gadget) are searched
        - \b 4 \e Gas particles (which are typically defined as type 0 for gadget) are searched
    \n
    \arg <b> \e Baryon_searchflag </b> 0/1/2 indicating gas/stellar search done separately from DM search \ref Options.iBaryonSearch. \n
        - \b 1 field search run as normal and then substructure search for baryons run using baryons identified in field search.
        - \b 2 field search also altered to treat baryons differently, allowing only DM particles to be used as head links (ie link dm-dm, dm-baryon, but not baryon-baryon nor baryon-dm).
        - \b 0 is do nothing special for baryon particles.
    \n
    \arg <b> \e Singlehalo_search </b> 0/1 flag indicates that no field search is going to be run and the entire volume will be treated as a background region. Useful if searching for substructures in non-cosmological simulations. But can also be co-opted for other searches using different outlier criteria and FOF algorithms. \ref Options.iSingleHalo \n

    \section localdensityconfig Parameters related to local density estimator.
    See \ref localfield.cxx, \ref bgfield.cxx & \ref localbgcomp.cxx for more details

    \arg <b> \e Nsearch_velocity </b> number of velocity neighbours used to calculate velocity density, adjust \ref Options.Nvel (suggested value is 32) \n
    \arg <b> \e Nsearch_physical </b> number of physical neighbours searched for Nv to calculate velocity density  \ref Options.Nsearch (suggested value is 256) \n
    \arg <b> \e Cell_fraction </b> fraction of a halo contained in a subvolume used to characterize the background  \ref Options.Ncellfac \n
    \arg <b> \e Grid_type </b> integer describing type of grid used to decompose volume for substructure search  \ref Options.gridtype (see \ref GRIDTYPES) \n
        - \b 1 \e standard physical shannon entropy, balanced KD tree volume decomposition into cells
        - \b 2 \e phase phase-space shannon entropy, balanced KD tree volume decomposition into cells
        - \b 3 \e simple simple physical balanced KD tree decomposition of volume into cells

    \section fofconfig Parameters related to FOF search
    See \ref search.cxx \ref fofalgo.h for more details

    \subsection fofsubconfig Configuration for substructure search
    \arg <b> \e Search_for_substructure </b> By default field objects are searched for internal substructures but can disable this by setting this to 0 \n
    \arg <b> \e Keep_FOF </b> if field 6DFOF search is done, allows to keep structures found in 3DFOF (can be interpreted as the inter halo stellar mass when only stellar search is used).\n
    \arg <b> \e FoF_search_type </b> There are several substructure FOF criteria implemented (see \ref FOFTYPES for more types and \ref fofalgo.h for implementation) \n
        - \b 1 \e standard phase-space based, well tested VELOCIraptor criterion.

    \arg <b> \e Outlier_threshold </b> threshold of sigma level of outliers to be searched (default is 2.5) \ref Options.ellthreshold \n
    \arg <b> \e Significance_level </b> minimum significance level of group (default is 1) \ref Options.siglevel \n
    \arg <b> \e Velocity_ratio </b> speed ratio used in linking particles \ref Options.Vratio \n
    \arg <b> \e Velocity_opening_angle </b> angle between velocities when linking (in units of \f$ \pi \f$) \ref Options.thetaopen \n
    \arg <b> \e Physical_linking_length </b> physical linking length used in fof (if cosmological gadget file then assumed to be in units of inter particle spacing \ref gadgetio.cxx, if loading in a single halo then based on average interparticle spacing calculated in \ref haloproperties.cxx). To be deprecated and replaced by Substructure_physical_linking_length \ref Options.ellphys \n
    \arg <b> \e Substructure_physical_linking_length </b> physical linking length used in fof (if cosmological gadget file then assumed to be in units of inter particle spacing \ref gadgetio.cxx, if loading in a single halo then based on average interparticle spacing calculated in \ref haloproperties.cxx) \ref Options.ellphys \n
    \arg <b> \e Minimum_size </b> Minimum group (substructure) size \ref Options.MinSize \n

    \arg <b> \e Iterative_searchflag </b> 1/0 use interactive substructure search which is designed to first identify spatially compact candidate outlier regions and then relaxes the criteria to find the more diffuse (in phase-space) regions associate with these candidate structures \ref Options.iiterflag \n
    \arg <b> \e Iterative_threshold_factor </b> factor multiplied with \ref Options.ellthreshold when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are \f$ \sim 1 \f$ \ref Options.ellfac \n
    \arg <b> \e Iterative_linking_length_factor </b> factor multiplied with \ref Options.ellphys when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are \f$ \sim 2 \f$ \ref Options.ellxfac \n
    \arg <b> \e Iterative_Vratio_length_factor </b> factor multiplied with \ref Options.Vratio when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are \f$ \sim 1 \f$ \ref Options.vfac \n
    \arg <b> \e Iterative_ThetaOp_length_factor </b> factor multiplied with \ref Options.thetaopen when using iterative method and identifying outlier regions associated with the initial candidate list of spatially compact outlier groups. Typical values are \f$ \sim 1 \f$ \ref Options.thetafac \n
    \arg <b> \e CMrefadjustsubsearch_flag </b> 1/0 flag indicating whether particles are moved to the rough CM velocity frame of the background before substructures are searched for. \ref Options.icmrefadjust \n

    \subsection foffieldconfig Configuration for field search
    \arg <b> \e FoF_Field_search_type </b> There are several FOF criteria implemented to search for so-called field objects (see \ref FOFTYPES for more types and \ref fofalgo.h for implementation) \ref Options.fofbgtype \n
        - \b 5 \e standard 3D FOF based algorithm
        - \b 4 \e standard 3D FOF based algorithm <b> FOLLOWED </b> by 6D FOF search using the velocity scale defined by the largest halo on <b> ONLY </b> particles in 3DFOF groups
        - \b 3 \e standard 3D FOF based algorithm <b> FOLLOWED </b> by 6D FOF search using the velocity scale for each 3DFOF group
    \arg <b> \e Minimum_halo_size </b> Allows field objects (or so-called halos) to require a different minimum size (typically would be <= \ref Options.MinSize. Default is -1 which sets it to \ref Options.MinSize) \ref Options.HaloMinSize \n
    \arg <b> \e Halo_linking_length_factor </b> allows one to use different physical linking lengths between field objects and substructures.  (Typically for 3DFOF searches of dark matter haloes, set to value such that this times \ref Options.ellphys = 0.2 the interparticle spacing when examining cosmological simulations ) \ref Options.ellhalophysfac \n
    \arg <b> \e Halo_velocity_linking_length_factor </b> allows one to use different velocity linking lengths between field objects and substructures when using 6D FOF searches.  (Since in such cases the general idea is to use the local velocity dispersion to define a scale, \f$ \geq5 \f$ times this value seems to correctly scale searches) \ref Options.ellhalovelfac \n
    \arg <b> \e Halo_6D_linking_length_factor </b> allows one to use different linking lengths between 3DFOF and 6DFOF field search. Typically values are \f$ \sim 1 \f$ \n
    \arg <b> \e Halo_6D_vel_linking_length_factor </b> scaling applied to dispersions used in 6DFOF field search. Typical values are \f$ \geq 1.25 \f$ \n

    \arg <b> \e Halo_core_search </b> 0/1/2 flag allows one to explicitly search for large 6D FOF cores that are indicative of a recent major merger. Since substructure is defined on the scale of the maximum cell size and major mergers typically result two or more phase-space dense regions that are \e larger than the cell size used in reasonable substructure searches, one can identify them using this search. The overall goal is to treat these objects differently than a substructure. However, if 2 is set, then smaller core is treated as substruture and all particles within the FOF envelop are assigned to the cores based on their phase-space distance to core particles \ref Options.iHaloCoreSearch \n
    \arg <b> \e Use_adaptive_core_search </b> 0/1 flag allows one to run complex adaptive phase-space search for large 6D FOF cores and then use these linking lengths to separate mergers. 0 is simple high density dispersively cold cores with v scale adaptive, 1 is adaptive with simple dispersions in both x and v.
    \arg <b> \e Use_phase_tensor_core_growth </b> 0/1 flag allows one to run complex phase-space growth of merger remnants (6D FOF cores found). 0 is assignment with simple x and v dispersion to nearest core particle, 1 is phase-space tensor distance assignemnt to CM of core.
    \arg <b> \e Halo_core_ellx_fac </b> scaling applied to linking length when identifying merger remnants. Typically values are \f$ \sim0.5 \f$  \ref Options.halocorexfac
    \arg <b> \e Halo_core_ellv_fac </b> scaling applied to local dispersion to define the velocity scale used to identify merger remnants. Typically values are \f$ \sim1 \f$  \ref Options.halocorevfac
    \arg <b> \e Halo_core_ncellfac </b> used to determine the minimum number of particles a merger remnants is composed of, specifically \f$ N_{\rm min}= f_{\rm ncell}* N_{\rm S} \f$. Typically values are 0.005
    \arg <b> \e Halo_core_adaptive_sigma_fac </b> used when running fully adaptive core search, specifies the width of the physical linking length in configuration space dispersion (think of this as how many sigma to include). Typically values are \f$ \sim2 \f$. This has been tested on hydrodynamnical simulations to separate galaxy mergers. \ref Options.halocoresigmafac
    \arg <b> \e Halo_core_num_loops </b> allows the core search to iterate, shrinking the velocity linking length to used till the number of cores identified decreases or this limit is reached. Allows apative search with larger linking length to be robust.  Typically values are \f$ \sim5 \f$ with loops only running often twice. \ref Options.halocorenumloops.
    \arg <b> \e Halo_core_loop_ellx_fac </b> Factor by which configuration linking length is decreased when running loops for core search.  Typically values are \f$ \sim0.75 \f$. \ref Options.halocorexfaciter
    \arg <b> \e Halo_core_loop_ellv_fac </b> Factor by which velocity linking length is decreased when running loops for core search.  Typically values are \f$ \sim0.75 \f$. \ref Options.halocorevfaciter
    \arg <b> \e Halo_core_loop_elln_fac </b> Factor by which min group size is changed when running loops for core search.  Typically values are \f$ \sim1.25 \f$. \ref Options.halocorenumfaciter
    \arg <b> \e Halo_core_phase_significance </b> Significance a core must be in terms of phase-space distance scaled by dispersions (sigma).

    \section unbindconfig Unbinding Parameters
    See \ref unbinding and \ref unbind.cxx for more details

    \arg <b> \e Unbind_flag  </b> whether or not substructures should be passed through an unbinding routine. \ref Options.uinfo & \ref UnbindInfo.unbindflag \n
    \arg <b> \e Unbind_Only_With_Baryons_flag  </b> whether or not substructures should be passed through an unbinding routine only after having baryons added. This applies to running velociraptor on a simulation containing all different types of particles where substructure is searched using initially dark matter particles before adding baryons in a phase-space search. \ref Options.uinfo & \ref UnbindInfo.iunbindwithbaryons \n
    \arg <b> \e Allowed_kinetic_potential_ratio  </b> ratio of kinetic to potential energy at which a particle is still considered bound, ie: particle is still bound
    if \f$ \alpha T+W<0 \f$, so \f$ \alpha=1 \f$ would be standard unbinding and \f$ \alpha<1 \f$ allows one to identify unbound tidal debris.
    Given that <b> VELOCIraptor </b> was designed to identify tidal streams, it makes little sense to have this set to 1 unless explicitly required.
    Note that the code still separates particles into bound and unbound. Typical values of \f$ \alpha\geq 0.2 \f$ seems to minimize the number of false positives
    in tidal debris while still identifying completely unbound tidal debris. \ref Options.uinfo & \ref UnbindInfo.Eratio \n
    \arg <b> \e Min_bound_mass_frac </b> Designed to demand a substructure still have a minimum amount of self-bound mass. \ref Options.uinfo & \ref UnbindInfo.minEfrac \n
    \arg <b> \e Bound_halos </b> 0/1/2 flag to make sure field objects such as haloes are self bound before (use 1) and also after (use 2) substructures have been identified
    and extracted from the halo. Demanding boundness after substructure search can have interesting consequences as it is possible that a multiple merger will appear as
    a single FOF halo, however all with all the cores removed, the FOF halo is actually an unbound structure. \ref Options.iBoundHalos \n
    \arg <b> \e Keep_background_potential </b> 1/0 flag When determining whether a structure is self-bound, the approach taken is to treat the candidate structure in isolation. Then determine the velocity reference frame to determine the kinetic energy of each particle and remove them. However, it is possible one wishes to keep the background particles when determining the potential, that is once one starts unbinding, don't treat the candidate structure in isolation but in a background sea. When finding tidal debris, it is useful to keep the background. \ref Options.uinfo & \ref UnbindInfo.bgpot \n
    \arg <b> \e Kinetic_reference_frame_type </b> specify kinetic frame when determining whether particle is bound.
    Default is to use the centre-of-mass velocity frame (0) but can also use region around minimum of the potential (1). \ref Options.uinfo & \ref UnbindInfo.cmvelreftype \n
    \arg <b> \e Min_npot_ref </b> Set the minimum number of particles used to calculate the velocity of the minimum of the potential (10). \ref Options.uinfo & \ref UnbindInfo.Npotref \n
    \arg <b> \e Frac_pot_ref </b> Set the fraction of particles used to calculate the velocity of the minimum of the potential (0.1). \ref Options.uinfo & \ref UnbindInfo.fracpotref \n
    \arg <b> \e Unbinding_type </b> Set the unbinding criteria, either just remove particles deemeed "unbound", that is those with \f$ \alpha T+W>0\f$, choosing \ref UPART. Or with \ref USYSANDPART
    removes "unbound" particles till system also has a true bound fraction > \ref UnbindInfo.minEfrac.
    \arg <b> \e Softening_length </b> Set the (simple plummer) gravitational softening length. \ref UnbindInfo.eps

    \section cosmoconfig Units & Cosmology
    \subsection unitconfig Units
    \arg <b> \e Length_unit </b> change the input length unit \ref Options.L \n
    \arg <b> \e Velocity_unit </b> change the input velocity unit \ref Options.V \n
    \arg <b> \e Mass_unit </b> change the input mass unit \ref Options.M \n
    \arg <b> \e Gravity </b> change the gravity unit. should be set such that \f$ v^2=GM/r \f$.   \ref Options.G \n
    \arg <b> \e Hubble_unit </b> change the value of Hubble expansion (from normal 100 km/s/Mpc to units used, that is velocity unit/length unit) \ref Options.G \n
    \arg <b> \e Mass_value </b> if not mass is stored using the option \b NOMASS (see \ref STF-makeflags) then this is the mass of the particles \ref Options.MassValue \n
    \arg <b> \e Length_unit_to_kpc </b> Specify the conversion factor from the output unit to kpc \ref Options.lengthtokpc \n
    \arg <b> \e Velocity_unit_to_kms </b> Specify the conversion factor from the output unit to kms \ref Options.velocitytokms \n
    \arg <b> \e Mass_unit_to_solarmass </b> Specify the conversion factor from the output unit to solar masses \ref Options.masstosolarmass \n

    \subsection cosmologyconfig Cosmology
    \arg <b> \e Period </b> if not defined in the input data one can pass the period in the input units of the data. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.p \n
    \arg <b> \e Scale_factor </b> scale factor time if not defined in the input data. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.a \n
    \arg <b> \e h_val </b> the "little h" value often used in cosmological simulation when defining expansion rates, distance scales, etc. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.h \n
    \arg <b> \e Omega_m </b> matter density in units of the critical density used in cosmological simulations. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.Omega_m \n
    \arg <b> \e Omega_Lambda </b> energy density of the cosmological constant (or can technically be used for a dark energy fluid) in units of the critical density used in cosmological simulations. This is not always necessary, for instance the gadget snapshot format has this in the header.  \ref Options.Omega_Lambda \n
    \arg <b> \e Omega_cdm </b> dark matter density in units of the critical density used in cosmological simulations. This is not always necessary since for simple pure DM cosmological simulations this is easily determined by the total mass in the simulation. <em> However, for multiple resolution or non-standard DM models, this should be provided. </em> \ref Options.Omega_cdm \n
    \arg <b> \e Omega_b </b> baryon density in units of the critical density used in cosmological simulations. This is not always necessary since for simple pure adiabatic cosmological simulations this is easily determined by the total gas mass in the simulation.  \ref Options.Omega_b \n
    \arg <b> \e Omega_r </b> radiation density in units of the critical density used in cosmological simulations. Typically 0 as negligible in most simulations.  \ref Options.Omega_r \n
    \arg <b> \e Omega_nu </b> neutrino density in units of the critical density used in cosmological simulations. Typically 0 as negligible in most simulations and may need alteration for consistent evolution across cosmic time.  \ref Options.Omega_nu \n
    \arg <b> \e Omega_DE </b> dark energy density in units of the critical density used in cosmological simulations. Typically 0 unless what evolution based on w_of_DE.  \ref Options.Omega_de \n
    \arg <b> \e w_of_DE </b> equation of state of the dark energy fluid, \f$ w=\frac{p}{\rho} \f$. This is not necessary unless one is using a cosmological simulation with a \f$ w\neq -1 \f$ \ref Options.w_de \n
    \arg <b> \e Virial_density </b> virial overdensity in units of the background matter density used in cosmological simulations. This is not absolutely necessary for gadget format as for this input, the Bryan & Norman 1998 virial density is calculated based on a LCDM cosmology as \ref Options.virlevel is set to -1 and if left, a virial density is calculated. If set and a gadget input is used, this overrides the Bryan & Norman calculation. \n
    \arg <b> \e Critical_density </b> critical density in input units used in cosmological simulations. This is not absolutely necessary for gadget format as it can be calculated from the quantities in the header for typical GR cosmologies. \ref Options.rhobg \n

    \section otherconfigs Other configuration options
    \arg <b> \e Effective_Resolution </b> If running a multiple resolution cosmological zoom simulation, simple method of scaling the linking length by using the period, ie: \f$ p/N_{\rm eff} \f$ \ref Options.Neff \n
    \arg <b> \e Snapshot_value </b> If halo ids need to be offset to some starting value based on the snapshot of the output, which is useful for some halo merger tree codes, one can specific a snapshot number, and all halo ids will be listed as internal haloid + \f$ sn\times10^{12}\f$. \ref Options.snapshotvalue \n
    \arg <b> \e Verbose </b> 2/1/0 flag indicating how talkative the code is (2 very verbose, 1 verbose, 0 quiet). \ref Options.iverbose \n

    \section propconfigs Property calculation options
    \arg <b> \e Inclusive_halo_mass </b> 1/0 flag indicating whether inclusive masses are calculated for field objects. \ref Options.iInclusiveHalo \n
    \arg <b> \e Extensive_halo_properties_output </b> 1/0 flag indicating whether to calculate/output even more halo properties. \ref Options.iextrahalooutput \n
    \arg <b> \e Extended_output </b> 1/0 flag indicating whether produce extended output for quick particle extraction from input catalog of particles in structures \ref Options.iextendedoutput \n
    \arg <b> \e Iterate_cm_flag </b> 1/0 flag indicating whether to use shrinking spheres to calculate the center of mass and velocity. \ref Options.iIterateCM \n
    \arg <b> \e Sort_by_binding_energy </b> 1/0 flag indicating whether to sort by particle binding energy or by potential (if 0). \ref Options.iSortByBindingEnergy \n

    \section ioconfigs I/O options
    \arg <b> \e Cosmological_input </b> 1/0 indicating that input simulation is cosmological or not. With cosmological input, a variety of length/velocity scales are set to determine such things as the virial overdensity, linking length. \ref Options.icosmologicalin \n
    \arg <b> \e Input_chunk_size </b> Amount of information to read from input file in one go (100000). \ref Options.inputbufsize \n
    \arg <b> \e Write_group_array_file </b> 0/1 flag indicating whether write a single large tipsy style group assignment file is written. \ref Options.iwritefof \n
    \arg <b> \e Separate_output_files </b> 1/0 flag indicating whether separate files are written for field and subhalo groups. \ref Options.iseparatefiles \n
    \arg <b> \e Binary_output </b> 3/2/1/0 flag indicating whether output is hdf, binary or ascii. \ref Options.ibinaryout, \ref OUTADIOS, \ref OUTHDF, \ref OUTBINARY, \ref OUTASCII \n
    \arg <b> \e Comoving_units </b> 1/0 flag indicating whether the properties output is in physical or comoving little h units. \ref Options.icomoveunit \n

    \section inputflags input flags related to varies input formats
    \arg <b> \e NSPH_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> SPH </b> blocks are read/in the file. \ref Options.gnsphblocks \n
    \arg <b> \e NStar_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> Star </b> blocks are read/in the file. \ref Options.gnstarblocks \n
    \arg <b> \e NBH_extra_blocks </b> If gadget snapshot is loaded one can specific the number of extra <b> Black hole </b> blocks are read/in the file. \ref Options.gnbhblocks \n

    \arg <b> \e HDF_name_convention </b> HDF dataset naming convection. See \ref hdfitems.h for what naming conventions are available and what names exist. Currently have \ref HDFNUMNAMETYPES. \ref Options.ihdfnameconvention \n
    \arg <b> \e Input_includes_dm_particle </b> If dm particle specific information is in the input file. \ref Options.iusedmparticles \n
    \arg <b> \e Input_includes_gas_particle </b> If gas particle specific information is in the input file. \ref Options.iusegasparticles \n
    \arg <b> \e Input_includes_star_particle </b> If star particle specific information is in the input file. \ref Options.iusestarparticles \n
    \arg <b> \e Input_includes_bh_particle </b> If bh/sink particle specific information is in the input file. \ref Options.iusesinkparticles \n
    \arg <b> \e Input_includes_wind_particle </b> If wind particle specific information is in the input file. \ref Options.iusewindparticles \n
    \arg <b> \e Input_includes_tracer_particle </b> If tracer particle specific information is in the input file. \ref Options.iusetracerparticles \n
    \arg <b> \e Input_includes_star_particle </b> If star particle specific information is in the input file. \ref Options.iusestarparticles \n

    \section mpiconfigs MPI specific options
    \arg <b> \e MPI_part_allocation_fac </b> Memory allocated in mpi mode to store particles is (1+factor)* the memory need for the initial mpi decomposition.
    This factor should be >0 and is mean to allow a little room for particles to be exchanged between mpi threads withouth having to require new memory allocations and copying
    of data. \ref Options.mpipartfac \n
    \arg <b> \e MPI_particle_total_buf_size </b> Total memory size in bytes used to store particles in temporary buffer such that
    particles are sent to non-reading mpi processes in one communication round in chunks of size buffer_size/NProcs/sizeof(Particle). \ref Options.mpiparticlebufsize \n

    */

inline void ConfigMessage(char *c) {
#ifndef USEMPI
    int ThisTask = 0;
#endif
    if (ThisTask==0) cerr<<c<<endl;
}

inline void ConfigExit(const char *error=nullptr) {
    if (error) {
        LOG_RANK0(error) << error;
    }
#ifdef USEMPI
    MPI_Abort(MPI_COMM_WORLD,8);
#else
    exit(8);
#endif
}

inline string ExtraFieldIndexName(unsigned int i){
    string s = "";
    s ="_index_"+to_string(i);
    return s;
}

inline vector<string> ExtraFieldCalculationsAllowedOptions()
{
    vector<string> list;
    list.push_back(("averagemassweighted"));
    list.push_back(("average"));
    list.push_back(("logaveragemassweighted"));
    list.push_back(("logaverage"));
    list.push_back(("totalmassweighted"));
    list.push_back(("total"));
    list.push_back(("stdmassweighted"));
    list.push_back(("std"));
    list.push_back(("logstdmassweighted"));
    list.push_back(("logstd"));
    list.push_back(("minmassweighted"));
    list.push_back(("min"));
    list.push_back(("maxmassweighted"));
    list.push_back(("max"));
    list.push_back(("aperture_total"));
    list.push_back(("aperture_average"));
    return list;
}

inline map<string,int> ExtraFieldCalculationsStringtoIntFlag()
{
    map<string, int> list;
    list[string("averagemassweighted")] = CALCAVERAGEMASSWEIGHT;
    list[string("average")] = CALCAVERAGE;
    list[string("logaveragemassweighted")] = CALCLOGAVERAGEMASSWEIGHT;
    list[string("logaverage")] = CALCLOGAVERAGE;
    list[string("totalmassweighted")] = CALCTOTALMASSWEIGHT;
    list[string("total")] = CALCTOTAL;
    list[string("stdmassweighted")] = CALCSTDMASSWEIGHT;
    list[string("std")] = CALCSTD;
    list[string("logstdmassweighted")] = CALCLOGSTDMASSWEIGHT;
    list[string("logstd")] = CALCLOGSTD;
    list[string("minmassweighted")] = CALCMINMASSWEIGHT;
    list[string("min")] = CALCMIN;
    list[string("maxmassweighted")] = CALCMAXMASSWEIGHT;
    list[string("max")] = CALCMAX;
    list[string("aperture_total")] = CALCQUANTITYAPERTURETOTAL;
    list[string("aperture_average")] = CALCQUANTITYAPERTUREAVERAGE;
    list[string("unkown")] = 0;
    return list;
}

inline map<int,string> ExtraFieldCalculationsIntFlagToString()
{
    map<int, string> list;
    list[CALCAVERAGEMASSWEIGHT] = string("averagemassweighted");;
    list[CALCAVERAGE] = string("average");
    list[CALCLOGAVERAGEMASSWEIGHT] = string("logaveragemassweighted");
    list[CALCLOGAVERAGE] = string("logaverage");
    list[CALCTOTALMASSWEIGHT] = string("totalmassweighted");
    list[CALCTOTAL] = string("total");
    list[CALCSTDMASSWEIGHT] = string("stdmassweighted");
    list[CALCSTD] = string("std");
    list[CALCLOGSTDMASSWEIGHT] = string("logstdmassweighted");
    list[CALCLOGSTD] = string("logstd");
    list[CALCMINMASSWEIGHT] = string("minmassweighted");
    list[CALCMIN] = string("min");
    list[CALCMAXMASSWEIGHT] = string("maxmassweighted");
    list[CALCMAX] = string("max");
    list[CALCQUANTITYAPERTURETOTAL] = string("aperture_total");
    list[CALCQUANTITYAPERTUREAVERAGE] = string("aperture_average");
    list[0] = string("unkown");
    return list;
}


inline void ExtraFieldCheck(string configentryname,
    vector<string> &names, vector<string>&output_names,
    vector<int> &calctypes, vector<unsigned int> &indices,
    vector<float> &conversions, vector<string> &units,
    vector<int> &pairindices,
    vector<string> &names_aperture, vector<string> &output_names_aperture,
    vector<int> &calctypes_aperture, vector<unsigned int> &indices_aperture,
    vector<float> &conversions_aperture, vector<string> &units_aperture
)
{
    if (names.size()==0) return;
    set<unsigned int> unique_int;
    set<string> unique_name, outputset;
    string outputfieldname;
    unsigned int entryindex, calctype;
    string scalctype;
    vector<int> stdfuncs, avefuncs;
    stdfuncs.push_back(CALCSTD);avefuncs.push_back(CALCAVERAGE);
    stdfuncs.push_back(CALCSTDMASSWEIGHT);avefuncs.push_back(CALCAVERAGEMASSWEIGHT);
    stdfuncs.push_back(CALCLOGSTD);avefuncs.push_back(CALCLOGAVERAGE);
    stdfuncs.push_back(CALCLOGSTDMASSWEIGHT);avefuncs.push_back(CALCLOGAVERAGEMASSWEIGHT);

    vector<string> newnames(names), newunits(units);
    vector<int> newcalctypes(calctypes);
    vector<unsigned int> newindices(indices);
    vector<float> newconversions(conversions);
    vector<string> allowedlist = ExtraFieldCalculationsAllowedOptions();
    string list;
    for (auto &x:allowedlist) list+= x+string(", ");
    map<string, int> calcstringtoint = ExtraFieldCalculationsStringtoIntFlag();
    map<int, string> calcinttostring = ExtraFieldCalculationsIntFlagToString();


    //check for unacceptable entries
    for (auto i=0;i<names.size();i++)
    {
        entryindex = indices[i];
        calctype = calctypes[i];
        scalctype = calcinttostring[calctype];
        outputfieldname = names[i]+ExtraFieldIndexName(entryindex)+string("_")+scalctype;
        if (calctype == 0) {
            LOG_RANK0(error) << "Unknown calculation requested for " << configentryname;
            LOG_RANK0(error) << "This is entry "  << i << " config options " << names[i] << "," << entryindex << "," << scalctype;
            LOG_RANK0(error) << "Please refer to documentation for list of calculations. Currently allowed:";
            LOG_RANK0(error) << list;
            ConfigExit();
        }
        if (outputset.count(outputfieldname) == 0) {
            outputset.insert(outputfieldname);
        }
        else{
            LOG_RANK0(warning) << "Duplicate entry found for " << configentryname;
            LOG_RANK0(warning) << "This is entry " << i << " config options " << names[i] << "," << entryindex << "," << calctype;
            LOG_RANK0(warning) << "This will be removed";
        }
    }

    //if standard deviation in list and average not, add it
    for (auto i=0;i<calctypes.size();i++) {
        for (auto j=0;j<stdfuncs.size();j++) {
            if (calctypes[i] == stdfuncs[j])
            {
                newnames.push_back(names[i]);
                newindices.push_back(indices[i]);
                newconversions.push_back(conversions[i]);
                newunits.push_back(units[i]);
                newcalctypes.push_back(avefuncs[j]);
            }
        }
    }
    if (names.size() != newnames.size()) {
        names = newnames;
        indices = newindices;
        units = newunits;
        conversions = newconversions;
        calctypes = newcalctypes;
    }
    newnames.clear();
    newindices.clear();
    newunits.clear();
    newconversions.clear();
    newcalctypes.clear();

    //remove duplicate entries
    outputset.clear();
    for (auto i=0;i<names.size();i++)
    {
        entryindex = indices[i];
        calctype = calctypes[i];
        scalctype = calcinttostring[calctype];
        outputfieldname = names[i]+ExtraFieldIndexName(entryindex)+scalctype;
        if (outputset.count(outputfieldname) == 0) {
            outputset.insert(outputfieldname);
            newnames.push_back(names[i]);
            newindices.push_back(indices[i]);
            newconversions.push_back(conversions[i]);
            newunits.push_back(units[i]);
            newcalctypes.push_back(calctypes[i]);
        }
    }
    //strip out aperture calculations and store them separately
    auto nentries = newnames.size();
    names.clear();
    calctypes.clear();
    indices.clear();
    units.clear();
    conversions.clear();
    for (auto i=0;i<nentries;i++) {
        if (newcalctypes[i] >0) {
            names.push_back(newnames[i]);
            indices.push_back(newindices[i]);
            units.push_back(newunits[i]);
            calctypes.push_back(newcalctypes[i]);
            conversions.push_back(newconversions[i]);
        }
        else {
            names_aperture.push_back(newnames[i]);
            indices_aperture.push_back(newindices[i]);
            units_aperture.push_back(newunits[i]);
            calctypes_aperture.push_back(newcalctypes[i]);
            conversions_aperture.push_back(newconversions[i]);
        }
    }
    pairindices.resize(names.size());
    for (auto i=0;i<names.size();i++)
    {
        entryindex = indices[i];
        calctype = calctypes[i];
        pairindices[i] = i;
        outputfieldname = names[i]+ExtraFieldIndexName(entryindex)
            +string("_")+calcinttostring[calctype];//+string("_")+units[i];
        output_names.push_back(outputfieldname);
    }
    for (auto i=0;i<calctypes.size();i++) {
        for (auto j=0;j<stdfuncs.size();j++) if (calctypes[i] == stdfuncs[j]) {
            for (auto k=0;k<calctypes.size();k++) {
                if (names[k] != names[i] || indices[k] != indices[i]) continue;
                if (calctypes[k] == avefuncs[j]) {
                    pairindices[i] = k;
                    break;
                }
            }
        }
    }
    for (auto i=0;i<names_aperture.size();i++)
    {
        entryindex = indices_aperture[i];
        calctype = calctypes_aperture[i];
        outputfieldname = names_aperture[i]+ExtraFieldIndexName(entryindex)
            +string("_")+calcinttostring[calctype];//+string("_")+units_aperture[i];
        output_names_aperture.push_back(outputfieldname);
    }
}

inline void ListDuplicateEntryCheck(string configentryname, int &num,
    vector<string> &names, vector<Double_t> &values)
{
    set<string> unique;
    if (num ==0) return;
    int iduplicates = 0;
    for (auto &val:names) {
        if(unique.count(val)) {
            LOG_RANK0(warning) << "Duplicate entry(ies) found in " << configentryname;
            LOG_RANK0(warning) << "Removing duplicate entries";
            iduplicates += 1;
        }
        else unique.insert(val);
    }
    if (iduplicates) {
        names = vector<string>(unique.begin(),unique.end());
        values.clear();
        for (auto &val:names) values.push_back(stof(val));
        num = values.size();

    }
}

inline vector<string> StripIndexOffName(
    vector<string> names1, vector<string> names2,
    vector<unsigned int> indices1, vector<unsigned int> indices2
)
{
    vector<string> result;
    for (auto iextra=0;iextra<names1.size();iextra++) {
        result.push_back(names1[iextra].substr(0,
            names1[iextra].size()-to_string(indices1[iextra]).size()));
    }
    for (auto iextra=0;iextra<names2.size();iextra++) {
        result.push_back(names2[iextra].substr(0,
            names2[iextra].size()-to_string(indices2[iextra]).size()));
    }
    return result;
}

///Read parameters from a parameter file. For list of currently implemented options see \ref configopt
///\todo still more parameters that can be adjusted
void GetParamFile(Options &opt)
{
#ifndef USEMPI
    int ThisTask =0, NProcs =1;
#endif
    string line,sep="=";
    string tag,val;
    char buff[1024],*pbuff,tbuff[1024],vbuff[1024];
    fstream paramfile,cfgfile;
    size_t pos;
    string dataline, token, delimiter = ",";
    map<string, int> calcstringtoint = ExtraFieldCalculationsStringtoIntFlag();
    map<int, string> calcinttostring = ExtraFieldCalculationsIntFlagToString();
    if (!FileExists(opt.pname)){
        if (ThisTask==0)
            LOG(error) << "Config file: " << opt.pname << " does not exist or can't be read, terminating";
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,9);
#else
            exit(9);
#endif
    }
    paramfile.open(opt.pname, ios::in);
    std::string::size_type j;
    //first find output name, determine the number of valid entries in
    if (paramfile.is_open())
    {
        while (paramfile.good()){
            getline(paramfile,line);
            //if line is not commented out or empty
            if (line[0]!='#'&&line.length()!=0) {
                j = line.find(sep);
                if (j != std::string::npos) {
                    //clean up string
                    tag=line.substr(0,j);
                    strcpy(buff, tag.c_str());
                    pbuff=strtok(buff," ");
                    strcpy(tbuff, pbuff);
                    val=line.substr(j+1,line.length()-(j+1));
                    strcpy(buff, val.c_str());
                    pbuff=strtok(buff," ");
                    if (pbuff==NULL) continue;
                    strcpy(vbuff, pbuff);
                    if (strcmp(tbuff, "Output")==0){
                        opt.outname=new char[1024];
                        strcpy(opt.outname,vbuff);
                        break;
                    }
                }
            }
        }
#ifndef SWIFTINTERFACE
        if (opt.outname==NULL) {
            if (ThisTask==0)
                cerr<<"No output name given, terminating"<<endl;
#ifdef USEMPI
            MPI_Finalize();
#endif
            exit(9);
        }
#endif
    }
    //sprintf(fname,"%s.cfg",opt.outname);
    //cfgfile.open(fname, ios::out);
    paramfile.clear();
    paramfile.seekg(0, ios::beg);
    if (paramfile.is_open())
    {
        while (paramfile.good()){
            getline(paramfile,line);
            //if line is not commented out or empty
            if (line[0]!='#'&&line.length()!=0) {
                j = line.find(sep);
                if (j != std::string::npos) {
                    //clean up string
                    tag=line.substr(0,j);
                    strcpy(buff, tag.c_str());
                    pbuff=strtok(buff," ");
                    strcpy(tbuff, pbuff);
                    val=line.substr(j+1,line.length()-(j+1));
                    strcpy(buff, val.c_str());
                    pbuff=strtok(buff," ");
                    if (pbuff==NULL) continue;
                    strcpy(vbuff, pbuff);
                    //store/read local density distribution funciton values
                    if (strcmp(tbuff, "Output_den")==0){
                        opt.smname=new char[1024];
                        sprintf(opt.smname,"%s.localden",opt.outname);
                    }
                    //config search type
                    else if (strcmp(tbuff, "Particle_search_type")==0)
                        opt.partsearchtype = atoi(vbuff);
                    else if (strcmp(tbuff, "FoF_search_type")==0)
                        opt.foftype = atoi(vbuff);
                    else if (strcmp(tbuff, "FoF_Field_search_type")==0)
                        opt.fofbgtype = atoi(vbuff);
                    else if (strcmp(tbuff, "Search_for_substructure")==0)
                        opt.iSubSearch = atoi(vbuff);
                    else if (strcmp(tbuff, "Keep_FOF")==0)
                        opt.iKeepFOF = atoi(vbuff);
                    else if (strcmp(tbuff, "Iterative_searchflag")==0)
                        opt.iiterflag = atoi(vbuff);
                    else if (strcmp(tbuff, "Baryon_searchflag")==0)
                        opt.iBaryonSearch = atoi(vbuff);
                    else if (strcmp(tbuff, "CMrefadjustsubsearch_flag")==0)
                        opt.icmrefadjust = atoi(vbuff);
                    else if (strcmp(tbuff, "Halo_core_search")==0)
                        opt.iHaloCoreSearch = atoi(vbuff);
                    else if (strcmp(tbuff, "Use_adaptive_core_search")==0)
                        opt.iAdaptiveCoreLinking = atof(vbuff);
                    else if (strcmp(tbuff, "Use_phase_tensor_core_growth")==0)
                        opt.iPhaseCoreGrowth = atof(vbuff);

                    //bg and fof parameters
                    else if (strcmp(tbuff, "Local_velocity_density_approximate_calculation")==0)
                        opt.iLocalVelDenApproxCalcFlag = atoi(vbuff);
                    else if (strcmp(tbuff, "Cell_fraction")==0)
                        opt.Ncellfac = atof(vbuff);
                    else if (strcmp(tbuff, "Grid_type")==0)
                        opt.gridtype = atoi(vbuff);
                    else if (strcmp(tbuff, "Nsearch_velocity")==0)
                        opt.Nvel = atoi(vbuff);
                    else if (strcmp(tbuff, "Nsearch_physical")==0)
                        opt.Nsearch = atoi(vbuff);
                    else if (strcmp(tbuff, "Outlier_threshold")==0)
                        opt.ellthreshold = atof(vbuff);
                    else if (strcmp(tbuff, "Significance_level")==0)
                        opt.siglevel = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_ratio")==0)
                        opt.Vratio = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_opening_angle")==0)
                        opt.thetaopen = atof(vbuff);
                    else if (strcmp(tbuff, "Substructure_physical_linking_length")==0)
                        opt.ellphys = atof(vbuff);
                    else if (strcmp(tbuff, "Physical_linking_length")==0)
                        opt.ellphys = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_linking_length")==0)
                        opt.ellvel = atof(vbuff);
                    else if (strcmp(tbuff, "Minimum_size")==0)
                        opt.MinSize = atoi(vbuff);
                    else if (strcmp(tbuff, "Minimum_halo_size")==0)
                        opt.HaloMinSize = atoi(vbuff);
                    else if (strcmp(tbuff, "Halo_linking_length_factor")==0)
                        opt.ellhalophysfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_3D_linking_length")==0)
                        opt.ellhalo3dxfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_velocity_linking_length_factor")==0)
                        opt.ellhalovelfac = atof(vbuff);
                    //specific to 6DFOF field search
                    else if (strcmp(tbuff, "Halo_6D_linking_length_factor")==0)
                        opt.ellhalo6dxfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_6D_vel_linking_length_factor")==0)
                        opt.ellhalo6dvfac = atof(vbuff);
                    //specific search for 6d fof core searches
                    else if (strcmp(tbuff, "Halo_core_ellx_fac")==0)
                        opt.halocorexfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_ellv_fac")==0)
                        opt.halocorevfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_ncellfac")==0)
                        opt.halocorenfac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_adaptive_sigma_fac")==0)
                        opt.halocoresigmafac = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_num_loops")==0)
                        opt.halocorenumloops = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_loop_ellx_fac")==0)
                        opt.halocorexfaciter = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_loop_ellv_fac")==0)
                        opt.halocorevfaciter = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_loop_elln_fac")==0)
                        opt.halocorenumfaciter = atof(vbuff);
                    else if (strcmp(tbuff, "Halo_core_phase_significance")==0)
                        opt.halocorephasedistsig = atof(vbuff);
                    //cleaning up substructures by merging if phase distance is small
                    else if (strcmp(tbuff, "Structure_phase_merge_dist")==0)
                        opt.coresubmergemindist = atof(vbuff);
                    else if (strcmp(tbuff, "Apply_phase_merge_to_host")==0)
                        opt.icoresubmergewithbg = atoi(vbuff);

                    //for changing factors used in iterative search
                    else if (strcmp(tbuff, "Iterative_threshold_factor")==0)
                        opt.ellfac = atof(vbuff);
                    else if (strcmp(tbuff, "Iterative_linking_length_factor")==0)
                        opt.ellxfac = atof(vbuff);
                    else if (strcmp(tbuff, "Iterative_Vratio_factor")==0)
                        opt.vfac = atof(vbuff);
                    else if (strcmp(tbuff, "Iterative_ThetaOp_factor")==0)
                        opt.thetafac = atof(vbuff);
                    //for changing effective resolution when rescaling linking lengh
                    else if (strcmp(tbuff, "Effective_resolution")==0)
                        opt.Neff = atof(vbuff);
                    //for changing effective resolution when rescaling linking lengh
                    else if (strcmp(tbuff, "Singlehalo_search")==0)
                        opt.iSingleHalo = atoi(vbuff);

                    //units, cosmology
                    else if (strcmp(tbuff, "Length_unit")==0)
                        opt.lengthinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_unit")==0)
                        opt.velocityinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Mass_unit")==0)
                        opt.massinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Hubble_unit")==0)
                        opt.H = atof(vbuff);
                    else if (strcmp(tbuff, "Gravity")==0)
                        opt.G = atof(vbuff);
                    else if (strcmp(tbuff, "Mass_value")==0)
                        opt.MassValue = atof(vbuff);
                    else if (strcmp(tbuff, "Period")==0)
                        opt.p = atof(vbuff);
                    else if (strcmp(tbuff, "Scale_factor")==0)
                        opt.a = atof(vbuff);
                    else if (strcmp(tbuff, "h_val")==0)
                        opt.h = atof(vbuff);
                    else if (strcmp(tbuff, "Critical_density")==0)
                        opt.rhobg = atof(vbuff);
                    else if (strcmp(tbuff, "Virial_density")==0)
                        opt.virlevel = atof(vbuff);
                    else if (strcmp(tbuff, "Omega_m")==0)
                        opt.Omega_m = atof(vbuff);
                    else if (strcmp(tbuff, "Omega_Lambda")==0)
                        opt.Omega_Lambda = atof(vbuff);
                    else if (strcmp(tbuff, "Omega_DE")==0)
                        opt.Omega_de = atof(vbuff);
                    else if (strcmp(tbuff, "Omega_k")==0)
                        opt.Omega_k = atof(vbuff);
                    else if (strcmp(tbuff, "Omega_cdm")==0)
                        opt.Omega_cdm= atof(vbuff);
                    else if (strcmp(tbuff, "Omega_b")==0)
                        opt.Omega_b= atof(vbuff);
                    else if (strcmp(tbuff, "Omega_r")==0)
                        opt.Omega_r= atof(vbuff);
                    else if (strcmp(tbuff, "Omega_nu")==0)
                        opt.Omega_nu= atof(vbuff);
                    else if (strcmp(tbuff, "w_of_DE")==0)
                        opt.w_de= atof(vbuff);
                    //so units can be specified to convert to kpc, km/s, solar mass
                    //not necessarily the code units data is reported but the
                    //translation of these code units to these units
                    else if (strcmp(tbuff, "Length_input_unit_conversion_to_output_unit")==0)
                        opt.lengthinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_input_unit_conversion_to_output_unit")==0)
                        opt.velocityinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Mass_input_unit_conversion_to_output_unit")==0)
                        opt.massinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Metallicity_input_unit_conversion_to_output_unit")==0)
                        opt.metallicityinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Star_formation_rate_input_unit_conversion_to_output_unit")==0)
                        opt.SFRinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Stellar_age_input_unit_conversion_to_output_unit")==0)
                        opt.stellarageinputconversion = atof(vbuff);
                    else if (strcmp(tbuff, "Stellar_age_input_is_cosmological_scalefactor")==0)
                        opt.istellaragescalefactor = atoi(vbuff);
                    else if (strcmp(tbuff, "Star_formation_rate_input_is_specific_star_formation_rate")==0)
                        opt.isfrisssfr = atoi(vbuff);
                    else if (strcmp(tbuff, "Length_unit_to_kpc")==0)
                        opt.lengthtokpc = atof(vbuff);
                    else if (strcmp(tbuff, "Velocity_to_kms")==0)
                        opt.velocitytokms = atof(vbuff);
                    else if (strcmp(tbuff, "Mass_to_solarmass")==0)
                        opt.masstosolarmass = atof(vbuff);
                    else if (strcmp(tbuff, "Metallicity_to_solarmetallicity")==0)
                        opt.metallicitytosolar = atof(vbuff);
                    else if (strcmp(tbuff, "Star_formation_rate_to_solarmassperyear")==0)
                        opt.SFRtosolarmassperyear = atof(vbuff);
                    else if (strcmp(tbuff, "Stellar_age_to_yr")==0)
                        opt.stellaragetoyrs = atof(vbuff);
                    else if (strcmp(tbuff, "Gas_star_forming_rate_threshold")==0)
                        opt.gas_sfr_threshold = atof(vbuff);
                    //unbinding
                    else if (strcmp(tbuff, "Unbind_flag")==0)
                        opt.uinfo.unbindflag = atoi(vbuff);
                    else if (strcmp(tbuff, "Unbind_Only_With_Baryons_flag")==0)
                        opt.uinfo.iunbindwithbaryons = atoi(vbuff);
                    else if (strcmp(tbuff, "Unbinding_type")==0)
                        opt.uinfo.unbindtype = atoi(vbuff);
                    else if (strcmp(tbuff, "Unbinding_use_thermal_energy")==0)
                        opt.uinfo.iuseinternalenergy = atoi(vbuff);
                    else if (strcmp(tbuff, "Bound_halos")==0)
                        opt.iBoundHalos = atoi(vbuff);
                    else if (strcmp(tbuff, "Allowed_kinetic_potential_ratio")==0)
                        opt.uinfo.Eratio = atof(vbuff);
                    else if (strcmp(tbuff, "Min_bound_mass_frac")==0)
                        opt.uinfo.minEfrac = atof(vbuff);
                    else if (strcmp(tbuff, "Keep_background_potential")==0)
                        opt.uinfo.bgpot = atoi(vbuff);
                    else if (strcmp(tbuff, "Kinetic_reference_frame_type")==0)
                        opt.uinfo.cmvelreftype = atoi(vbuff);
                    else if (strcmp(tbuff, "Min_npot_ref")==0)
                        opt.uinfo.Npotref = atoi(vbuff);
                    else if (strcmp(tbuff, "Frac_pot_ref")==0)
                        opt.uinfo.fracpotref = atof(vbuff);
                    else if (strcmp(tbuff, "Unbinding_max_unbound_removal_fraction_per_iteration")==0)
                        opt.uinfo.maxunbindfrac = atof(vbuff);
                    else if (strcmp(tbuff, "Unbinding_max_unbound_fraction")==0)
                        opt.uinfo.maxunboundfracforiterativeunbind = atof(vbuff);
                    else if (strcmp(tbuff, "Unbinding_max_unbound_fraction_allowed")==0)
                        opt.uinfo.maxallowedunboundfrac = atof(vbuff);
                    else if (strcmp(tbuff, "Softening_length")==0)
                        opt.uinfo.eps = atof(vbuff);
                    else if (strcmp(tbuff, "Approximate_potential_calculation")==0)
                        opt.uinfo.iapproxpot = atoi(vbuff);
                    else if (strcmp(tbuff, "Approximate_potential_calculation_particle_number_fraction")==0)
                        opt.uinfo.approxpotnumfrac = atof(vbuff);
                    else if (strcmp(tbuff, "Approximate_potential_calculation_min_particle_number")==0)
                        opt.uinfo.approxpotminnum = atoi(vbuff);
                    else if (strcmp(tbuff, "Approximate_potential_calculation_method")==0)
                        opt.uinfo.approxpotmethod = atoi(vbuff);

                    //property related
                    else if (strcmp(tbuff, "Reference_frame_for_properties")==0)
                        opt.iPropertyReferencePosition = atoi(vbuff);
                    else if (strcmp(tbuff, "Particle_type_for_reference_frames")==0)
                        opt.ParticleTypeForRefenceFrame = atoi(vbuff);
                    else if (strcmp(tbuff, "Iterate_cm_flag")==0)
                        opt.iIterateCM = atoi(vbuff);
                    else if (strcmp(tbuff, "Inclusive_halo_masses")==0)
                        opt.iInclusiveHalo = atoi(vbuff);
                    else if (strcmp(tbuff, "Spherical_overdenisty_calculation_limited_to_structure_types")==0) {
                        int stype = atoi(vbuff);
                        opt.SphericalOverdensitySeachMaxStructLevel = HALOSTYPE;
                        opt.SphericalOverdensitySeachMaxStructLevel += HALOCORESTYPE*stype;
                    }
                    else if (strcmp(tbuff, "Extensive_halo_properties_output")==0)
                        opt.iextrahalooutput = atoi(vbuff);
                    else if (strcmp(tbuff, "Extensive_gas_properties_output")==0)
                        opt.iextragasoutput = atoi(vbuff);
                    else if (strcmp(tbuff, "Extensive_star_properties_output")==0)
                        opt.iextrastaroutput = atoi(vbuff);
                    else if (strcmp(tbuff, "Extensive_interloper_properties_output")==0)
                        opt.iextrainterloperoutput = atoi(vbuff);
                    else if (strcmp(tbuff, "Calculate_aperture_quantities")==0)
                        opt.iaperturecalc = atoi(vbuff);
                    else if (strcmp(tbuff, "Number_of_apertures")==0)
                        opt.aperturenum = atoi(vbuff);
                    else if (strcmp(tbuff, "Aperture_values_in_kpc")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.aperture_names_kpc.push_back(token);
                            opt.aperture_values_kpc.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                    }
                    else if (strcmp(tbuff, "Number_of_projected_apertures")==0)
                        opt.apertureprojnum = atoi(vbuff);
                    else if (strcmp(tbuff, "Projected_aperture_values_in_kpc")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.aperture_proj_names_kpc.push_back(token);
                            opt.aperture_proj_values_kpc.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                    }
                    else if (strcmp(tbuff, "Calculate_radial_profiles")==0)
                        opt.iprofilecalc = atoi(vbuff);
                    else if (strcmp(tbuff, "Radial_profile_min_FOF_size")==0)
                        opt.profileminFOFsize = atoi(vbuff);
                    else if (strcmp(tbuff, "Radial_profile_min_size")==0)
                        opt.profileminsize = atoi(vbuff);
                    else if (strcmp(tbuff, "Number_of_radial_profile_bin_edges")==0)
                        opt.profilenbins = atoi(vbuff);
                    else if (strcmp(tbuff, "Radial_profile_norm")==0)
                        opt.iprofilenorm = atoi(vbuff);
                    else if (strcmp(tbuff, "Radial_profile_bin_edges")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.profile_bin_edges.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                    }
                    else if (strcmp(tbuff, "Number_of_overdensities")==0)
                        opt.SOnum = atoi(vbuff);
                    else if (strcmp(tbuff, "Overdensity_values_in_critical_density")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.SOthresholds_names_crit.push_back(token);
                            opt.SOthresholds_values_crit.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                    }
                    //other options
                    else if (strcmp(tbuff, "Verbose")==0)
                        opt.iverbose = atoi(vbuff);
                    else if (strcmp(tbuff, "Write_group_array_file")==0)
                        opt.iwritefof = atoi(vbuff);
                    else if (strcmp(tbuff, "Snapshot_value")==0)
                        opt.snapshotvalue = HALOIDSNVAL*atoi(vbuff);

                    //input related
                    else if (strcmp(tbuff, "Cosmological_input")==0)
                        opt.icosmologicalin = atoi(vbuff);
                    //input read related
                    else if (strcmp(tbuff, "Input_chunk_size")==0)
                        opt.inputbufsize = atol(vbuff);
                    else if (strcmp(tbuff, "MPI_particle_total_buf_size")==0)
                        opt.mpiparticletotbufsize = atol(vbuff);
                    //mpi memory related
                    else if (strcmp(tbuff, "MPI_part_allocation_fac")==0)
                        opt.mpipartfac = atof(vbuff);
                    else if (strcmp(tbuff, "MPI_number_of_tasks_per_write")==0)
                        opt.mpinprocswritesize = atoi(vbuff);
                    else if (strcmp(tbuff, "MPI_use_zcurve_mesh_decomposition")==0)
                        opt.impiusemesh = (atoi(vbuff)>0);
                    else if (strcmp(tbuff, "MPI_zcurve_mesh_decomposition_min_num_cells_per_dim")==0)
                        opt.minnumcellperdim = atoi(vbuff);
                    else if (strcmp(tbuff, "MPI_zcurve_mesh_decomposition_num_cells_per_dim")==0) {
                        opt.numcellsperdim = atoi(vbuff);
                        opt.numcells = opt.numcellsperdim*opt.numcellsperdim*opt.numcellsperdim;
                    }
                    // reduced halo properties 
                    else if (strcmp(tbuff, "Calculate_reduced_halo_properties")==0)
                        opt.ireducedhaloproperties = (atoi(vbuff)>0);

                    ///OpenMP related
                    else if (strcmp(tbuff, "OMP_run_fof")==0)
                        opt.iopenmpfof = atoi(vbuff);
                    else if (strcmp(tbuff, "OMP_fof_region_size")==0)
                        opt.openmpfofsize = atoi(vbuff);
                    else if (strcmp(tbuff, "Gas_internal_property_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.gas_internalprop_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.gas_internalprop_index.size() == 0) {
                            opt.gas_internalprop_index.resize(opt.gas_internalprop_names.size(),0);
                        }
                        if (opt.gas_internalprop_function.size() == 0) {
                            opt.gas_internalprop_function.resize(opt.gas_internalprop_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.gas_internalprop_input_output_unit_conversion_factors.size() == 0) {
                            opt.gas_internalprop_input_output_unit_conversion_factors.resize(
                                opt.gas_internalprop_names.size(),1.0);
                        }
                        if (opt.gas_internalprop_output_units.size() == 0) {
                            opt.gas_internalprop_output_units.resize(opt.gas_internalprop_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "Gas_internal_property_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_internalprop_index = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_internal_property_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_internalprop_function = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_internal_property_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_internalprop_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_internal_property_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_internalprop_input_output_unit_conversion_factors = tempvec;
                    }
		    else if (strcmp(tbuff, "Tcut_halogas")==0) {
                         opt.temp_max_cut = atof(vbuff);
		    }
		    else if (strcmp(tbuff, "Gas_temperature_input_unit_conversion_to_output_unit")==0) {
                         opt.temp_input_output_unit_conversion_factor = atof(vbuff);
		    }
		    else if (strcmp(tbuff, "Gas_hot_overdensity_normalisation")==0) {
                opt.hot_gas_overdensity_normalisation = atof(vbuff);
                // if this is define, then make sure the SO calculations are check to contain this value, otherwise modify it to contain them.
                auto iter = std::find(opt.SOthresholds_values_crit.begin(), opt.SOthresholds_values_crit.end(), opt.hot_gas_overdensity_normalisation);
                if(iter == opt.SOthresholds_values_crit.end()){
                    opt.SOthresholds_values_crit.push_back(opt.hot_gas_overdensity_normalisation);
                    opt.SOnum = opt.SOnum + 1;
                }
                else if (strcmp(tbuff, "Overdensity_values_in_critical_density")==0) {
                    pos=0;
                    dataline=string(vbuff);
                    while ((pos = dataline.find(delimiter)) != string::npos) {
                        token = dataline.substr(0, pos);
                        opt.SOthresholds_names_crit.push_back(token);
                        opt.SOthresholds_values_crit.push_back(stof(token));
                        dataline.erase(0, pos + delimiter.length());
                    }
                }
            }
                    else if (strcmp(tbuff, "Aperture_values_normalised_to_overdensity")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.aperture_hotgas_normalised_to_overdensity = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_internal_property_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.star_internalprop_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.star_internalprop_index.size() == 0) {
                            opt.star_internalprop_index.resize(opt.star_internalprop_names.size(),0);
                        }
                        if (opt.gas_internalprop_function.size() == 0) {
                            opt.star_internalprop_function.resize(opt.star_internalprop_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.star_internalprop_input_output_unit_conversion_factors.size() == 0) {
                            opt.star_internalprop_input_output_unit_conversion_factors.resize(
                                opt.star_internalprop_names.size(),1.0);
                        }
                        if (opt.star_internalprop_output_units.size() == 0) {
                            opt.star_internalprop_output_units.resize(opt.star_internalprop_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "Star_internal_property_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_internalprop_index = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_internal_property_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_internalprop_function = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_internal_property_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_internalprop_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_internal_property_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_internalprop_input_output_unit_conversion_factors = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_internal_property_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.bh_internalprop_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.bh_internalprop_index.size() == 0) {
                            opt.bh_internalprop_index.resize(opt.bh_internalprop_names.size(),0);
                        }
                        if (opt.bh_internalprop_function.size() == 0) {
                            opt.bh_internalprop_function.resize(opt.bh_internalprop_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.bh_internalprop_input_output_unit_conversion_factors.size() == 0) {
                            opt.bh_internalprop_input_output_unit_conversion_factors.resize(
                                opt.bh_internalprop_names.size(),1.0);
                        }
                        if (opt.bh_internalprop_output_units.size() == 0) {
                            opt.bh_internalprop_output_units.resize(opt.bh_internalprop_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "BH_internal_property_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_internalprop_index = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_internal_property_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_internalprop_function = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_internal_property_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_internalprop_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_internal_property_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_internalprop_input_output_unit_conversion_factors = tempvec;
                    }
                    else if (strcmp(tbuff, "Extra_DM_internal_property_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.extra_dm_internalprop_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.extra_dm_internalprop_index.size() == 0) {
                            opt.extra_dm_internalprop_index.resize(opt.extra_dm_internalprop_names.size(),0);
                        }
                        if (opt.extra_dm_internalprop_function.size() == 0) {
                            opt.extra_dm_internalprop_function.resize(opt.extra_dm_internalprop_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.extra_dm_internalprop_input_output_unit_conversion_factors.size() == 0) {
                            opt.extra_dm_internalprop_input_output_unit_conversion_factors.resize(
                                opt.extra_dm_internalprop_names.size(),1.0);
                        }
                        if (opt.extra_dm_internalprop_output_units.size() == 0) {
                            opt.extra_dm_internalprop_output_units.resize(opt.extra_dm_internalprop_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "Extra_DM_internal_property_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.extra_dm_internalprop_index = tempvec;
                    }
                    else if (strcmp(tbuff, "Extra_DM_internal_property_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.extra_dm_internalprop_function = tempvec;
                    }
                    else if (strcmp(tbuff, "Extra_DM_internal_property_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.extra_dm_internalprop_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "Extra_DM_internal_property_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.extra_dm_internalprop_input_output_unit_conversion_factors = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.gas_chem_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.gas_chem_index.size() == 0) {
                            opt.gas_chem_index.resize(opt.gas_chem_names.size(),0);
                        }
                        if (opt.gas_chem_function.size() == 0) {
                            opt.gas_chem_function.resize(opt.gas_chem_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.gas_chem_input_output_unit_conversion_factors.size() == 0) {
                            opt.gas_chem_input_output_unit_conversion_factors.resize(
                                opt.gas_chem_names.size(),1.0);
                        }
                        if (opt.gas_chem_output_units.size() == 0) {
                            opt.gas_chem_output_units.resize(opt.gas_chem_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_chem_index = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_chem_function = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_chem_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_chem_input_output_unit_conversion_factors = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_chemistry_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.star_chem_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.star_chem_index.size() == 0) {
                            opt.star_chem_index.resize(opt.star_chem_names.size(),0);
                        }
                        if (opt.star_chem_function.size() == 0) {
                            opt.star_chem_function.resize(opt.star_chem_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.star_chem_input_output_unit_conversion_factors.size() == 0) {
                            opt.star_chem_input_output_unit_conversion_factors.resize(
                                opt.star_chem_names.size(),1.0);
                        }
                        if (opt.star_chem_output_units.size() == 0) {
                            opt.star_chem_output_units.resize(opt.star_chem_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "Star_chemistry_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_chem_index = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_chemistry_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_chem_function = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_chemistry_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_chem_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_chemistry_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_chem_input_output_unit_conversion_factors = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_chemistry_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.bh_chem_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.bh_chem_index.size() == 0) {
                            opt.bh_chem_index.resize(opt.bh_chem_names.size(),0);
                        }
                        if (opt.bh_chem_function.size() == 0) {
                            opt.bh_chem_function.resize(opt.bh_chem_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.bh_chem_input_output_unit_conversion_factors.size() == 0) {
                            opt.bh_chem_input_output_unit_conversion_factors.resize(
                                opt.bh_chem_names.size(),1.0);
                        }
                        if (opt.bh_chem_output_units.size() == 0) {
                            opt.bh_chem_output_units.resize(opt.bh_chem_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "BH_chemistry_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_chem_index = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_chemistry_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_chem_function = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_chemistry_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_chem_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_chemistry_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_chem_input_output_unit_conversion_factors = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_production_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.gas_chemproduction_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.gas_chemproduction_index.size() == 0) {
                            opt.gas_chemproduction_index.resize(opt.gas_chemproduction_names.size(),0);
                        }
                        if (opt.gas_chemproduction_function.size() == 0) {
                            opt.gas_chemproduction_function.resize(opt.gas_chemproduction_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.gas_chemproduction_input_output_unit_conversion_factors.size() == 0) {
                            opt.gas_chemproduction_input_output_unit_conversion_factors.resize(
                                opt.gas_chemproduction_names.size(),1.0);
                        }
                        if (opt.gas_chemproduction_output_units.size() == 0) {
                            opt.gas_chemproduction_output_units.resize(opt.gas_chemproduction_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_production_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_chemproduction_index = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_prodution_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_chemproduction_function = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_production_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_chemproduction_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "Gas_chemistry_production_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.gas_chemproduction_input_output_unit_conversion_factors = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_chemistry_production_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.star_chemproduction_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.star_chemproduction_index.size() == 0) {
                            opt.star_chemproduction_index.resize(opt.star_chemproduction_names.size(),0);
                        }
                        if (opt.star_chemproduction_function.size() == 0) {
                            opt.star_chemproduction_function.resize(opt.star_chemproduction_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.star_chemproduction_input_output_unit_conversion_factors.size() == 0) {
                            opt.star_chemproduction_input_output_unit_conversion_factors.resize(
                                opt.star_chemproduction_names.size(),1.0);
                        }
                        if (opt.star_chemproduction_output_units.size() == 0) {
                            opt.star_chemproduction_output_units.resize(opt.star_chemproduction_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "Star_chemistry_production_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_chemproduction_index = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_chemistry_prodution_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_chemproduction_function = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_chemistry_production_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_chemproduction_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "Star_chemistry_production_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.star_chemproduction_input_output_unit_conversion_factors = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_chemistry_production_names")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            opt.bh_chemproduction_names.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        if (opt.bh_chemproduction_index.size() == 0) {
                            opt.bh_chemproduction_index.resize(opt.bh_chemproduction_names.size(),0);
                        }
                        if (opt.bh_chemproduction_function.size() == 0) {
                            opt.bh_chemproduction_function.resize(opt.bh_chemproduction_names.size(),CALCAVERAGEMASSWEIGHT);
                        }
                        if (opt.bh_chemproduction_input_output_unit_conversion_factors.size() == 0) {
                            opt.bh_chemproduction_input_output_unit_conversion_factors.resize(
                                opt.bh_chemproduction_names.size(),1.0);
                        }
                        if (opt.bh_chemproduction_output_units.size() == 0) {
                            opt.bh_chemproduction_output_units.resize(opt.bh_chemproduction_names.size(),string("unknownunits"));
                        }
                    }
                    else if (strcmp(tbuff, "BH_chemistry_production_index_in_file")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<unsigned int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stoi(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_chemproduction_index = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_chemistry_prodution_calculation_type")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<int> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(calcstringtoint[token]);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_chemproduction_function = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_chemistry_production_output_units")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<string> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(token);
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_chemproduction_output_units = tempvec;
                    }
                    else if (strcmp(tbuff, "BH_chemistry_production_input_output_unit_conversion_factors")==0) {
                        pos=0;
                        dataline=string(vbuff);
                        vector<float> tempvec;
                        while ((pos = dataline.find(delimiter)) != string::npos) {
                            token = dataline.substr(0, pos);
                            tempvec.push_back(stof(token));
                            dataline.erase(0, pos + delimiter.length());
                        }
                        opt.bh_chemproduction_input_output_unit_conversion_factors = tempvec;
                    }

                    //output related
                    else if (strcmp(tbuff, "Separate_output_files")==0)
                        opt.iseparatefiles = atoi(vbuff);
                    else if (strcmp(tbuff, "Binary_output")==0)
                        opt.ibinaryout = atoi(vbuff);
                    else if (strcmp(tbuff, "Comoving_units")==0)
                        opt.icomoveunit = atoi(vbuff);
                    else if (strcmp(tbuff, "Extended_output")==0)
                        opt.iextendedoutput = atoi(vbuff);
                    else if (strcmp(tbuff, "Spherical_overdensity_halo_particle_list_output")==0)
                        opt.iSphericalOverdensityPartList = atoi(vbuff);
                    else if (strcmp(tbuff, "Sort_by_binding_energy")==0)
                        opt.iSortByBindingEnergy = atoi(vbuff);
                    else if (strcmp(tbuff, "SUBFIND_like_output")==0)
                        opt.isubfindoutput = atoi(vbuff);
                    else if (strcmp(tbuff, "No_particle_ID_list_output")==0)
                        opt.inoidoutput = atoi(vbuff);

                    //gadget io related to extra info for sph, stars, bhs,
                    else if (strcmp(tbuff, "NSPH_extra_blocks")==0)
                        opt.gnsphblocks = atoi(vbuff);
                    else if (strcmp(tbuff, "NStar_extra_blocks")==0)
                        opt.gnstarblocks = atoi(vbuff);
                    else if (strcmp(tbuff, "NBH_extra_blocks")==0)
                        opt.gnbhblocks = atoi(vbuff);


                    //input related to info for stars, bhs, winds/tracers, etc
                    else if (strcmp(tbuff, "HDF_name_convention")==0)
                        opt.ihdfnameconvention = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_dm_particle")==0)
                        opt.iusedmparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_gas_particle")==0)
                        opt.iusegasparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_star_particle")==0)
                        opt.iusestarparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_bh_particle")==0)
                        opt.iusesinkparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_wind_particle")==0)
                        opt.iusewindparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_tracer_particle")==0)
                        opt.iusetracerparticles = atoi(vbuff);
                    else if (strcmp(tbuff, "Input_includes_extradm_particle")==0)
                        opt.iuseextradarkparticles = atoi(vbuff);
                }
            }
        }
        paramfile.close();
        //cfgfile.close();
    }
    configure_logging(opt);
}

void NOMASSCheck(Options &opt)
{
#ifdef NOMASS
    if (opt.MassValue<=0) {
        LOG_RANK0(error) << "Code compiled to not store mass per particle";
        LOG_RANK0(error) << "Valid value of particle mass must be extracted from input (HDF5, Gadget) or set in the config file through  Mass_value";
        LOG_RANK0(error) << "Currently value <=0. Either input did not contain a valid value and update to config file needed";
        ConfigExit();
    }
#endif
}

void ConfigCheck(Options &opt)
{
#ifndef USEMPI
    int ThisTask =0;
#endif
    //if code is on the fly, no point in checking fname
    // input type, number of input files, etc
    if (opt.iontheflyfinding == false) {
    if (opt.fname==NULL||opt.outname==NULL){
        ConfigExit("Must provide input and output file names");
    }
    if (opt.inputtype == IOHDF && opt.ihdfnameconvention == -1) {
        ConfigExit("HDF input passed but the naming convention is not set in the config file. Please set HDF_name_convetion");
    }
    if (opt.num_files<1){
        ConfigExit("Invalid number of input files (<1)");
    }
    if (opt.inputbufsize<1){
        ConfigExit("Invalid read buf size (<1)");
    }
    }

    if (opt.iBaryonSearch && !(opt.partsearchtype==PSTALL || opt.partsearchtype==PSTDARK))
    {
        ConfigExit("Conflict in config file: both gas/star/etc particle type search AND the separate baryonic (gas,star,etc) search flag are on. Check config");
    }
    if (opt.iBoundHalos && opt.iKeepFOF)
    {
        ConfigExit("Conflict in config file: Asking for Bound Field objects but also asking to keep the 3DFOF/then run 6DFOF. This is incompatible. Check config");
    }
    if (opt.HaloMinSize==-1) opt.HaloMinSize=opt.MinSize;

    if (opt.lengthtokpc<=0){
        ConfigExit("Invalid unit conversion, length unit to kpc is <=0 or was not set. Update config file");
    }
    //convert reference apertures
    opt.lengthtokpc30pow2 /= opt.lengthtokpc*opt.lengthtokpc;
    opt.lengthtokpc50pow2 /= opt.lengthtokpc*opt.lengthtokpc;
    if (opt.velocitytokms<=0){
        ConfigExit("Invalid unit conversion, velocity unit to km/s is <=0 or was not set. Update config file");
    }
    if (opt.masstosolarmass<=0){
        ConfigExit("Invalid unit conversion, mass unit to solar mass is <=0 or was not set. Update config file");
    }
#if defined(GASON) || defined(STARON) || defined(BHON)
    if (opt.SFRtosolarmassperyear<=0){
        ConfigExit("Invalid unit conversion, SFR unit to solar mass per year is <=0 or was not set. Update config file");
    }
    if (opt.metallicitytosolar<=0){
        ConfigExit("Invalid unit conversion, metallicity unit to solar is <=0 or was not set. Update config file");
    }
    if (opt.stellaragetoyrs<=0){
        ConfigExit("Invalid unit conversion, stellar age unit to year is <=0 or was not set. Update config file");
    }
#endif
#ifndef NOMASS
    opt.MassValue = 1.0;
#endif

    //check gravity and hubble unit
    double gravity, hubunit;
    gravity = CalcGravitationalConstant(opt);
    hubunit = CalcHubbleUnit(opt);
    if (opt.G<=0) opt.G = gravity;
    else {
        auto diff = fabs(gravity-opt.G)/opt.G;
        if (diff>1e-2) {
            LOG_RANK0(warning) << "Configuration provides gravitational constant that differs from the default by more than 1%";
            LOG_RANK0(warning) << "Expecation: " << gravity;
            LOG_RANK0(warning) << "Value passed: " << opt.G;
        }
    }
    if (opt.H<=0) opt.H = hubunit;
    else {
        auto diff = fabs(hubunit-opt.H)/opt.H;
        if (diff>1e-2) {
            LOG_RANK0(warning) << "Configuration provides hubble units that differs from the default by more than 1%";
            LOG_RANK0(warning) << "Expecation: " << hubunit;
            LOG_RANK0(warning) << "Value passed: " << opt.H;
        }
    }
    if (opt.ParticleTypeForRefenceFrame!=-1)
    {
        // map<int, string> typetostring = {{GASTYPE, "Gas/SPH"}, {DARKTYPE, "DM"}, {STARTYPE, "Star"}, {BHTYPE, "Black holes"}};
        map<int, string> typetostring;
        typetostring[GASTYPE]="gas"; typetostring[DARKTYPE]="dm"; typetostring[STARTYPE]="star"; typetostring[BHTYPE]="bh";
        if (
           (opt.ParticleTypeForRefenceFrame == GASTYPE && opt.iusegasparticles == 0) ||
           (opt.ParticleTypeForRefenceFrame == DARKTYPE && opt.iusedmparticles == 0) ||
           (opt.ParticleTypeForRefenceFrame == STARTYPE && opt.iusestarparticles == 0) ||
           (opt.ParticleTypeForRefenceFrame == BHTYPE && opt.iusesinkparticles == 0)
        )
        {
            LOG_RANK0(error) << "Using particle type for reference frame that is not loaded."; 
            LOG_RANK0(error) <<"Asking for "+to_string(opt.ParticleTypeForRefenceFrame)+"which corresponds to "+typetostring[opt.ParticleTypeForRefenceFrame];
            LOG_RANK0(error)<<"Update config to either load the appropriate paricle with Input_includes_"+typetostring[opt.ParticleTypeForRefenceFrame]+"_particle";
            LOG_RANK0(error)<<"Or change the reference particle type";
            ConfigExit();
        }
    }
#ifdef USEMPI
    if (opt.minnumcellperdim<MESH_MINCELLSPERDIM){
        LOG_RANK0(warning) << "MPI mesh too coarse, minimum number of cells per dimension from which to produce z-curve decomposition is "<<MESH_MINCELLSPERDIM<<". Resetting";
        opt.minnumcellperdim = MESH_MINCELLSPERDIM;
    }
    if (opt.numcellsperdim<MESH_MINCELLSPERDIM && opt.numcellsperdim > 0){
        LOG_RANK0(warning) << "MPI cells per dim set for mesh but too coarse, minimum number of cells per dimension from which to produce z-curve decomposition is "<<MESH_MINCELLSPERDIM<<". Resetting";
        opt.numcellsperdim = MESH_MINCELLSPERDIM;
    }
    if (opt.mpiparticletotbufsize<(long int)(sizeof(Particle)*NProcs) && opt.mpiparticletotbufsize!=-1){
        LOG_RANK0(error) << "Invalid input particle buffer send size, minimum input buffer size given particle byte size "
                         << vr::memory_amount(sizeof(Particle)) << " and have " << NProcs << " MPI processes is "
                         << vr::memory_amount(sizeof(Particle) * NProcs);
        ConfigExit();
    }
    //if total buffer size is -1 then calculate individual buffer size based on default mpi size
    else if (opt.mpiparticletotbufsize==-1) {
        opt.mpiparticletotbufsize=MPIPartBufSize*NProcs*sizeof(Particle);
        opt.mpiparticlebufsize=MPIPartBufSize;
    }
    else {
        opt.mpiparticlebufsize=opt.mpiparticletotbufsize/NProcs/sizeof(Particle);
    }
    if (opt.mpipartfac<0){
        ConfigExit("Invalid MPI particle allocation factor, must be >0");
    }
    else if (opt.mpipartfac>1){
        LOG_RANK0(warning) << "MPI Particle allocation factor is high (>1)";
    }
    if (opt.mpinprocswritesize<1){
#ifdef USEPARALLELHDF
        LOG_RANK0(warning) << "Number of MPI task writing collectively < 1. Setting to 1";
        opt.mpinprocswritesize = 1;
#endif
    }
    if (opt.mpinprocswritesize>NProcs){
#ifdef USEPARALLELHDF
        LOG_RANK0(warning) << "Number of MPI task writing collectively > NProcs. Setting to NProcs";
        opt.mpinprocswritesize = NProcs;
#endif
    }
#endif

#ifdef USEOPENMP
    if (opt.iopenmpfof == 1 && opt.openmpfofsize < ompfofsearchnum){
        LOG_RANK0(warning) << "OpenMP FOF search region is small, resetting to minimum of " << ompfofsearchnum;
        opt.openmpfofsize = ompfofsearchnum;
    }
#endif

#ifndef USEHDF
    if (opt.ibinaryout==OUTHDF){
        ConfigExit("Code not compiled with HDF output enabled. Recompile with this enabled or change Binary_output");
    }
    if (opt.isubfindoutput) {
        ConfigExit("Code not compiled with HDF output enabled. Recompile with this enabled to produce subfind like output or turn off subfind like output");
    }
#endif

#ifndef USEADIOS
    if (opt.ibinaryout==OUTADIOS){
        ConfigExit("Code not compiled with ADIOS output enabled. Recompile with this enabled or change Binary_output");
    }
#endif
    if (opt.iaperturecalc>0) {
        if (opt.aperturenum != opt.aperture_values_kpc.size()) {
            ConfigExit("Aperture calculations requested but mismatch between number stated and values provided. Check config");
        }
        if (opt.aperturenum == 0) {
            LOG_RANK0(warning) << "Aperture calculations requested but number of apertures is zero";
        }
        if (opt.apertureprojnum != opt.aperture_proj_values_kpc.size()) {
            ConfigExit("Projected aperture calculations requested but mismatch between number stated and values provided. Check config");
        }
        if (opt.apertureprojnum == 0) {
            LOG_RANK0(warning) << "Aperture calculations requested but number of projected apertures is zero";
        }
        if (opt.aperturenum>0) for (auto &x:opt.aperture_values_kpc) x/=opt.lengthtokpc;
        if (opt.apertureprojnum>0) for (auto &x:opt.aperture_proj_values_kpc) x/=opt.lengthtokpc;
    }
    if (opt.SOnum>0) {
        if (opt.SOnum != opt.SOthresholds_values_crit.size()) {
            ConfigExit("SO calculations requested but mismatch between number stated and values provided. Check config");
        }
    }
    if (opt.iprofilecalc>0) {
        if (opt.profileminsize < 0) {
            LOG_RANK0(warning) << "Radial profile calculations limited to objects of < 0 size! Ignoring and setting to 0";
            opt.profileminsize = 0;
        }
        if (opt.profileminFOFsize < 0) {
            LOG_RANK0(warning) << "Radial profile calculations limited to FOF objects of < 0 size! Ignoring and setting to 0";
            opt.profileminFOFsize = 0;
        }
        if (opt.profilenbins != opt.profile_bin_edges.size()) {
            ConfigExit("Radial profile calculations requested but mismatch between number of edges stated and number provided. Check config");
        }
        if (opt.profilenbins == 0) {
            ConfigExit("Radial profile calculations requested but number of bin edges is zero. Check config");
        }
        if (opt.iprofilebintype == PROFILERBINTYPELOG) {
            for (auto i=0;i<opt.profilenbins;i++) opt.profile_bin_edges[i]=pow(10.0,opt.profile_bin_edges[i]);
        }
        if (opt.iprofilenorm == PROFILERNORMR200CRIT) opt.profileradnormstring = "R_200crit";
        else if (opt.iprofilenorm == PROFILERNORMPHYS) opt.profileradnormstring = "Physical";
    }

    if (opt.gas_internalprop_names.size() != opt.gas_internalprop_function.size()){
        ConfigExit("GAS: # of Internal Property does not the # of function type entries. Check config");
    }
    if (opt.gas_internalprop_names.size() != opt.gas_internalprop_index.size()){
        ConfigExit("GAS: # of Internal Property names does not the # of index in file entries. Check config");
    }
    if (opt.gas_chem_names.size() != opt.gas_chem_function.size()){
        ConfigExit("GAS: # of Chemistry names does not the # of function type entries. Check config");
    }
    if (opt.gas_chem_names.size() != opt.gas_chem_index.size()){
        ConfigExit("GAS: # of Chemistry Production names does not the # of index in file entries. Check config");
    }
    if (opt.gas_chemproduction_names.size() != opt.gas_chemproduction_function.size()){
        ConfigExit("GAS: # of Chemistry Production names does not the # of function type entries. Check config");
    }
    if (opt.gas_chemproduction_names.size() != opt.gas_chemproduction_index.size()){
        ConfigExit("GAS: # of Chemistry Production names does not the # of index in file entries. Check config");
    }

    if (opt.star_internalprop_names.size() != opt.star_internalprop_function.size()){
        ConfigExit("STAR: # of Internal Property does not the # of function type entries. Check config.");
    }
    if (opt.star_internalprop_names.size() != opt.star_internalprop_index.size()){
        ConfigExit("STAR: # of Internal Property names does not the # of index in file entries. Check config.");
    }
    if (opt.star_chem_names.size() != opt.star_chem_function.size()){
        ConfigExit("STAR: # of Chemistry names does not the # of function type entries. Check config.");
    }
    if (opt.star_chem_names.size() != opt.star_chem_index.size()){
        ConfigExit("STAR: # of Chemistry Production names does not the # of index in file entries. Check config.");
    }
    if (opt.star_chemproduction_names.size() != opt.star_chemproduction_function.size()){
        ConfigExit("STAR: # of Chemistry Production names does not the # of function type entries. Check config.");
    }
    if (opt.star_chemproduction_names.size() != opt.star_chemproduction_index.size()){
        ConfigExit("STAR: # of Chemistry Production names does not the # of index in file entries. Check config.");
    }

    if (opt.bh_internalprop_names.size() != opt.bh_internalprop_function.size()){
        ConfigExit("BH: # of Internal Property does not the # of function type entries. Check config.");
    }
    if (opt.bh_internalprop_names.size() != opt.bh_internalprop_index.size()){
        ConfigExit("BH: # of Internal Property names does not the # of index in file entries. Check config.");
    }
    if (opt.bh_chem_names.size() != opt.bh_chem_function.size()){
        ConfigExit("BH: # of Chemistry names does not the # of function type entries. Check config.");
    }
    if (opt.bh_chem_names.size() != opt.bh_chem_index.size()){
        ConfigExit("BH: # of Chemistry Production names does not the # of index in file entries. Check config.");
    }
    if (opt.bh_chemproduction_names.size() != opt.bh_chemproduction_function.size()){
        ConfigExit("BH: # of Chemistry Production names does not the # of function type entries. Check config.");
    }
    if (opt.bh_chemproduction_names.size() != opt.bh_chemproduction_index.size()){
        ConfigExit("BH: # of Chemistry Production names does not the # of index in file entries. Check config.");
    }

    if (opt.extra_dm_internalprop_names.size() != opt.extra_dm_internalprop_function.size()){
        ConfigExit("Extra_DM: # of Internal Property does not the # of function type entries. Check config.");
    }
    if (opt.extra_dm_internalprop_names.size() != opt.extra_dm_internalprop_index.size()){
        ConfigExit("Extra_DM: # of Internal Property names does not the # of index in file entries. Check config.");
    }
    if (opt.uinfo.iapproxpot) {
        if (opt.uinfo.approxpotnumfrac <=0) {
            ConfigExit("Calculating approximate potential but fraction of particles <=0. Check config.");
        }
        if (opt.uinfo.approxpotminnum <=0) {
            ConfigExit("Calculating approximate potential but min number of particles to use <=0. Check config.");
        }
        if (opt.uinfo.approxpotmethod < POTAPPROXMETHODTREE || opt.uinfo.approxpotmethod > POTAPPROXMETHODRAND) {
            ConfigExit("In approximate potential but using invalid method for sampling particles. Use 0 for Tree and 1 for Rand. Check config.");
        }
    }

    set<string> uniqueval;
    set<string> outputset;
    string configentryname, outputfieldname, mainname;

    //clean up aperture list and spherical overdensity list to remove duplicates
    configentryname = "Overdensity_values_in_critical_density";
    ListDuplicateEntryCheck(configentryname, opt.SOnum,
        opt.SOthresholds_names_crit, opt.SOthresholds_values_crit);
    configentryname = "Aperture_values_in_kpc";
    ListDuplicateEntryCheck(configentryname, opt.aperturenum,
        opt.aperture_names_kpc, opt.aperture_values_kpc);
    configentryname = "Projected_aperture_values_in_kpc";
    ListDuplicateEntryCheck(configentryname, opt.apertureprojnum,
        opt.aperture_proj_names_kpc, opt.aperture_proj_values_kpc);

    //set output field names for extra properties
    configentryname = "Gas_internal_property_names";
    ExtraFieldCheck(configentryname,
        opt.gas_internalprop_names, opt.gas_internalprop_output_names,
        opt.gas_internalprop_function, opt.gas_internalprop_index,
        opt.gas_internalprop_input_output_unit_conversion_factors, opt.gas_internalprop_output_units,
        opt.gas_internalprop_index_paired_calc,
        opt.gas_internalprop_names_aperture, opt.gas_internalprop_output_names_aperture,
        opt.gas_internalprop_function_aperture, opt.gas_internalprop_index_aperture,
        opt.gas_internalprop_input_output_unit_conversion_factors_aperture, opt.gas_internalprop_output_units_aperture
    );
    configentryname = "Gas_chemistry_names";
    ExtraFieldCheck(configentryname,
        opt.gas_chem_names, opt.gas_chem_output_names,
        opt.gas_chem_function, opt.gas_chem_index,
        opt.gas_chem_input_output_unit_conversion_factors, opt.gas_chem_output_units,
        opt.gas_chem_index_paired_calc,
        opt.gas_chem_names_aperture, opt.gas_chem_output_names_aperture,
        opt.gas_chem_function_aperture, opt.gas_chem_index_aperture,
        opt.gas_chem_input_output_unit_conversion_factors_aperture, opt.gas_chem_output_units_aperture
    );
    configentryname = "Gas_chemistry_production_names";
    ExtraFieldCheck(configentryname,
        opt.gas_chemproduction_names, opt.gas_chemproduction_output_names,
        opt.gas_chemproduction_function, opt.gas_chemproduction_index,
        opt.gas_chemproduction_input_output_unit_conversion_factors, opt.gas_chemproduction_output_units,
        opt.gas_chemproduction_index_paired_calc,
        opt.gas_chemproduction_names_aperture, opt.gas_chemproduction_output_names_aperture,
        opt.gas_chemproduction_function_aperture, opt.gas_chemproduction_index_aperture,
        opt.gas_chemproduction_input_output_unit_conversion_factors_aperture, opt.gas_chemproduction_output_units_aperture
    );
    opt.gas_extraprop_aperture_calc = (opt.gas_internalprop_names_aperture.size() +
        opt.gas_chem_names_aperture.size()+opt.gas_chemproduction_names_aperture.size()>0);

    configentryname = "Star_internal_property_names";
    ExtraFieldCheck(configentryname,
        opt.star_internalprop_names, opt.star_internalprop_output_names,
        opt.star_internalprop_function, opt.star_internalprop_index,
        opt.star_internalprop_input_output_unit_conversion_factors, opt.star_internalprop_output_units,
        opt.star_internalprop_index_paired_calc,
        opt.star_internalprop_names_aperture, opt.star_internalprop_output_names_aperture,
        opt.star_internalprop_function_aperture, opt.star_internalprop_index_aperture,
        opt.star_internalprop_input_output_unit_conversion_factors_aperture, opt.star_internalprop_output_units_aperture
    );
    configentryname = "Star_chemistry_names";
    ExtraFieldCheck(configentryname,
        opt.star_chem_names, opt.star_chem_output_names,
        opt.star_chem_function, opt.star_chem_index,
        opt.star_chem_input_output_unit_conversion_factors, opt.star_chem_output_units,
        opt.star_chem_index_paired_calc,
        opt.star_chem_names_aperture, opt.star_chem_output_names_aperture,
        opt.star_chem_function_aperture, opt.star_chem_index_aperture,
        opt.star_chem_input_output_unit_conversion_factors_aperture, opt.star_chem_output_units_aperture
    );
    configentryname = "Star_chemistry_production_names";
    ExtraFieldCheck(configentryname,
        opt.star_chemproduction_names, opt.star_chemproduction_output_names,
        opt.star_chemproduction_function, opt.star_chemproduction_index,
        opt.star_chemproduction_input_output_unit_conversion_factors, opt.star_chemproduction_output_units,
        opt.star_chemproduction_index_paired_calc,
        opt.star_chemproduction_names_aperture, opt.star_chemproduction_output_names_aperture,
        opt.star_chemproduction_function_aperture, opt.star_chemproduction_index_aperture,
        opt.star_chemproduction_input_output_unit_conversion_factors_aperture, opt.star_chemproduction_output_units_aperture
    );
    opt.star_extraprop_aperture_calc = (opt.star_internalprop_names_aperture.size() +
        opt.star_chem_names_aperture.size()+opt.star_chemproduction_names_aperture.size()>0);

    configentryname = "BH_internal_property_names";
    ExtraFieldCheck(configentryname,
        opt.bh_internalprop_names, opt.bh_internalprop_output_names,
        opt.bh_internalprop_function, opt.bh_internalprop_index,
        opt.bh_internalprop_input_output_unit_conversion_factors, opt.bh_internalprop_output_units,
        opt.bh_internalprop_index_paired_calc,
        opt.bh_internalprop_names_aperture, opt.bh_internalprop_output_names_aperture,
        opt.bh_internalprop_function_aperture, opt.bh_internalprop_index_aperture,
        opt.bh_internalprop_input_output_unit_conversion_factors_aperture, opt.bh_internalprop_output_units_aperture
    );
    configentryname = "BH_chemistry_names";
    ExtraFieldCheck(configentryname,
        opt.bh_chem_names, opt.bh_chem_output_names,
        opt.bh_chem_function, opt.bh_chem_index,
        opt.bh_chem_input_output_unit_conversion_factors, opt.bh_chem_output_units,
        opt.bh_chem_index_paired_calc,
        opt.bh_chem_names_aperture, opt.bh_chem_output_names_aperture,
        opt.bh_chem_function_aperture, opt.bh_chem_index_aperture,
        opt.bh_chem_input_output_unit_conversion_factors_aperture, opt.bh_chem_output_units_aperture
    );
    configentryname = "BH_chemistry_production_names";
    ExtraFieldCheck(configentryname,
        opt.bh_chemproduction_names, opt.bh_chemproduction_output_names,
        opt.bh_chemproduction_function, opt.bh_chemproduction_index,
        opt.bh_chemproduction_input_output_unit_conversion_factors, opt.bh_chemproduction_output_units,
        opt.bh_chemproduction_index_paired_calc,
        opt.bh_chemproduction_names_aperture, opt.bh_chemproduction_output_names_aperture,
        opt.bh_chemproduction_function_aperture, opt.bh_chemproduction_index_aperture,
        opt.bh_chemproduction_input_output_unit_conversion_factors_aperture, opt.bh_chemproduction_output_units_aperture
    );
    opt.bh_extraprop_aperture_calc = (opt.bh_internalprop_names_aperture.size() +
        opt.bh_chem_names_aperture.size()+opt.bh_chemproduction_names_aperture.size()>0);

    configentryname = "Extra_DM_internal_property_names";
    ExtraFieldCheck(configentryname,
        opt.extra_dm_internalprop_names, opt.extra_dm_internalprop_output_names,
        opt.extra_dm_internalprop_function, opt.extra_dm_internalprop_index,
        opt.extra_dm_internalprop_input_output_unit_conversion_factors, opt.extra_dm_internalprop_output_units,
        opt.extra_dm_internalprop_index_paired_calc,
        opt.extra_dm_internalprop_names_aperture, opt.extra_dm_internalprop_output_names_aperture,
        opt.extra_dm_internalprop_function_aperture, opt.extra_dm_internalprop_index_aperture,
        opt.extra_dm_internalprop_input_output_unit_conversion_factors_aperture, opt.extra_dm_internalprop_output_units_aperture
    );
    opt.extra_dm_extraprop_aperture_calc = (opt.extra_dm_internalprop_names_aperture.size() >0);

    if (opt.gas_extraprop_aperture_calc && opt.iaperturecalc == 0){
        LOG_RANK0(error) << "Requesting extra gas properties to be calculated in apertures but apertures not set";
        LOG_RANK0(error) << "Enable aperture calculations and provide apertures or change the calculations requested for extra properties. Check config";
        ConfigExit();
    }
    if (opt.star_extraprop_aperture_calc && opt.iaperturecalc == 0){
        LOG_RANK0(error) << "Requesting extra star properties to be calculated in apertures but apertures not set";
        LOG_RANK0(error) << "Enable aperture calculations and provide apertures or change the calculations requested for extra properties. Check config";
        ConfigExit();
    }
    if (opt.bh_extraprop_aperture_calc && opt.iaperturecalc == 0){
        LOG_RANK0(error) << "Requesting extra BH properties to be calculated in apertures but apertures not set";
        LOG_RANK0(error) << "Enable aperture calculations and provide apertures or change the calculations requested for extra properties. Check config";
        ConfigExit();
    }
    if (opt.extra_dm_extraprop_aperture_calc && opt.iaperturecalc == 0){
        LOG_RANK0(error) << "Requesting extra dm properties to be calculated in apertures but apertures not set";
        LOG_RANK0(error) << "Enable aperture calculations and provide apertures or change the calculations requested for extra properties. Check config";
        ConfigExit();
    }


    //set halo 3d fof linking length if necessary
    if (opt.ellhalo3dxfac == -1) {
        opt.ellhalo3dxfac = opt.ellhalophysfac * opt.ellphys;
    }
    else {
        opt.ellhalophysfac = opt.ellhalo3dxfac / opt.ellphys;
    }

    if (ThisTask==0) {
        log_options_summary(opt);
    }

    //store the git hash
    opt.git_sha1 = velociraptor::git_sha1();
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

///define the configuration constructor so all input options are stored.
ConfigInfo::ConfigInfo(Options &opt){
    string datastring;
    //if compiler is super old and does not have at least std 11 implementation to_string does not exist
#ifndef OLDCCOMPILER
    //Add version
    AddEntry("Git_revision", opt.git_sha1);

    //general search operations
    AddEntry("Particle_search_type", opt.partsearchtype);
    AddEntry("FoF_search_type", opt.foftype);
    AddEntry("FoF_Field_search_type", opt.fofbgtype);
    AddEntry("Search_for_substructure", opt.iSubSearch);
    AddEntry("Keep_FOF", opt.iKeepFOF);
    AddEntry("Iterative_searchflag", opt.iiterflag);
    AddEntry("Baryon_searchflag", opt.iBaryonSearch);
    AddEntry("CMrefadjustsubsearch_flag", opt.icmrefadjust);
    AddEntry("Halo_core_search", opt.iHaloCoreSearch);
    AddEntry("Use_adaptive_core_search", opt.iAdaptiveCoreLinking);
    AddEntry("Use_phase_tensor_core_growth", opt.iPhaseCoreGrowth);

    //local field parameters
    AddEntry("Local_velocity_density_approximate_calculation", opt.iLocalVelDenApproxCalcFlag);
    AddEntry("Cell_fraction", opt.Ncellfac);
    AddEntry("Grid_type", opt.gridtype);
    AddEntry("Nsearch_velocity", opt.Nvel);
    AddEntry("Nsearch_physical", opt.Nsearch);

    //substructure search parameters
    AddEntry("Outlier_threshold", opt.ellthreshold);
    AddEntry("Significance_level", opt.siglevel);
    AddEntry("Velocity_ratio", opt.Vratio);
    AddEntry("Velocity_opening_angle", opt.thetaopen);
    ///\todo this configuration option will be deprecated. Replaced by Substructure_physical_linking_length
    //AddEntry("Physical_linking_length", opt.ellphys);
    AddEntry("Substructure_physical_linking_length", opt.ellphys);
    AddEntry("Velocity_linking_length", opt.ellvel);
    AddEntry("Minimum_size", opt.MinSize);

    //field object specific searches
    AddEntry("Minimum_halo_size", opt.HaloMinSize);
    ///\todo this configuration option will be deprecated. Replaced by Halo_3D_physical_linking_length
    //AddEntry("Halo_linking_length_factor", opt.ellhalophysfac);
    AddEntry("Halo_3D_linking_length", opt.ellhalo3dxfac);
    AddEntry("Halo_velocity_linking_length_factor", opt.ellhalovelfac);

    //specific to 6DFOF field search
    AddEntry("Halo_6D_linking_length_factor", opt.ellhalo6dxfac);
    AddEntry("Halo_6D_vel_linking_length_factor", opt.ellhalo6dvfac);

    //specific search for 6d fof core searches
    AddEntry("Halo_core_ellx_fac", opt.halocorexfac);
    AddEntry("Halo_core_ellv_fac", opt.halocorevfac);
    AddEntry("Halo_core_ncellfac", opt.halocorenfac);
    AddEntry("Halo_core_adaptive_sigma_fac", opt.halocoresigmafac);
    AddEntry("Halo_core_num_loops", opt.halocorenumloops);
    AddEntry("Halo_core_loop_ellx_fac", opt.halocorexfaciter);
    AddEntry("Halo_core_loop_ellv_fac", opt.halocorevfaciter);
    AddEntry("Halo_core_loop_elln_fac", opt.halocorenumfaciter);
    AddEntry("Halo_core_phase_significance", opt.halocorephasedistsig);

    //for merging structures together
    AddEntry("Structure_phase_merge_dist", opt.coresubmergemindist);
    AddEntry("Apply_phase_merge_to_host", opt.icoresubmergewithbg);

    //for changing factors used in iterative search
    AddEntry("Iterative_threshold_factor", opt.ellfac);
    AddEntry("Iterative_linking_length_factor", opt.ellxfac);
    AddEntry("Iterative_Vratio_factor", opt.vfac);
    AddEntry("Iterative_ThetaOp_factor", opt.thetafac);

    //for changing effective resolution when rescaling linking lengh
    #ifdef HIGHRES
    AddEntry("Effective_resolution", opt.Neff);
    #endif

    //for changing effective resolution when rescaling linking lengh
    AddEntry("Singlehalo_search", opt.iSingleHalo);

    //units, cosmology
    AddEntry("Length_unit", opt.lengthinputconversion);
    AddEntry("Velocity_unit", opt.velocityinputconversion);
    AddEntry("Mass_unit", opt.massinputconversion);
    AddEntry("Length_input_unit_conversion_to_output_unit", opt.lengthinputconversion);
    AddEntry("Velocity_input_unit_conversion_to_output_unit", opt.velocityinputconversion);
    AddEntry("Mass_input_unit_conversion_to_output_unit", opt.massinputconversion);
    AddEntry("Star_formation_rate_input_unit_conversion_to_output_unit", opt.SFRinputconversion);
    AddEntry("Metallicity_input_unit_conversion_to_output_unit", opt.metallicityinputconversion);
    AddEntry("Stellar_age_input_is_cosmological_scalefactor", opt.istellaragescalefactor);
    AddEntry("Hubble_unit", opt.H);
    AddEntry("Gravity", opt.G);
    AddEntry("Mass_value", opt.MassValue);
    AddEntry("Length_unit_to_kpc", opt.lengthtokpc);
    AddEntry("Velocity_to_kms", opt.velocitytokms);
    AddEntry("Mass_to_solarmass", opt.masstosolarmass);
    AddEntry("Star_formation_rate_to_solarmassperyear", opt.SFRtosolarmassperyear);
    AddEntry("Metallicity_to_solarmetallicity", opt.metallicitytosolar);
    AddEntry("Stellar_age_to_yr", opt.stellaragetoyrs);

    // simulation/cosmology info
    AddEntry("Period", opt.p);
    AddEntry("Scale_factor", opt.a);
    AddEntry("h_val", opt.h);
    AddEntry("Omega_m", opt.Omega_m);
    AddEntry("Omega_Lambda", opt.Omega_Lambda);
    AddEntry("Critical_density", opt.rhobg);
    AddEntry("Virial_density", opt.virlevel);
    AddEntry("Omega_cdm", opt.Omega_cdm);
    AddEntry("Omega_b", opt.Omega_b);
    AddEntry("Omega_r", opt.Omega_r);
    AddEntry("Omega_nu", opt.Omega_nu);
    AddEntry("Omega_k", opt.Omega_k);
    AddEntry("Omega_DE", opt.Omega_de);
    AddEntry("w_of_DE", opt.w_de);

    //unbinding
    AddEntry("Unbind_flag", opt.uinfo.unbindflag);
    AddEntry("Unbind_With_Only_Baryons_flag", opt.uinfo.iunbindwithbaryons);
    AddEntry("Unbinding_type", opt.uinfo.unbindtype);
    AddEntry("Bound_halos", opt.iBoundHalos);
    AddEntry("Allowed_kinetic_potential_ratio", opt.uinfo.Eratio);
    AddEntry("Min_bound_mass_frac", opt.uinfo.minEfrac);
    AddEntry("Keep_background_potential", opt.uinfo.bgpot);
    AddEntry("Kinetic_reference_frame_type", opt.uinfo.cmvelreftype);
    AddEntry("Min_npot_ref", opt.uinfo.Npotref);
    AddEntry("Frac_pot_ref", opt.uinfo.fracpotref);
    AddEntry("Unbinding_max_unbound_removal_fraction_per_iteration", opt.uinfo.maxunbindfrac);
    AddEntry("Unbinding_max_unbound_fraction", opt.uinfo.maxunboundfracforiterativeunbind);
    AddEntry("Unbinding_max_unbound_fraction_allowed", opt.uinfo.maxallowedunboundfrac);
    AddEntry("Softening_length", opt.uinfo.eps);
    AddEntry("Approximate_potential_calculation", opt.uinfo.iapproxpot);
    AddEntry("Approximate_potential_calculation_particle_number_fraction", opt.uinfo.approxpotnumfrac);
    AddEntry("Approximate_potential_calculation_min_particle", opt.uinfo.approxpotminnum);
    AddEntry("Approximate_potential_calculation_method", opt.uinfo.approxpotmethod);

    //property related
    AddEntry("Inclusive_halo_masses", opt.iInclusiveHalo);
    AddEntry("Extensive_halo_properties_output", opt.iextrahalooutput);
    AddEntry("Extensive_gas_properties_output", opt.iextragasoutput);
    AddEntry("Extensive_star_properties_output", opt.iextrastaroutput);
    AddEntry("Extensive_interloper_properties_output", opt.iextrainterloperoutput);
    AddEntry("Iterate_cm_flag", opt.iIterateCM);
    AddEntry("Sort_by_binding_energy", opt.iSortByBindingEnergy);
    AddEntry("Reference_frame_for_properties", opt.iPropertyReferencePosition);

    AddEntry("Calculate_aperture_quantities", opt.iaperturecalc);
    AddEntry("Number_of_apertures", opt.aperturenum);
    AddEntry("Aperture_values_in_kpc", opt.aperture_values_kpc);
    AddEntry("Number_of_projected_apertures", opt.apertureprojnum);
    AddEntry("Projected_aperture_values_in_kpc", opt.aperture_proj_values_kpc);
    AddEntry("Calculate_radial_profiles", opt.iprofilecalc);
    if(opt.iprofilecalc) {
        AddEntry("Radial_profile_min_FOF_size", opt.profileminFOFsize);
        AddEntry("Radial_profile_min_size", opt.profileminsize);
        AddEntry("Number_of_radial_profile_bin_edges", opt.profilenbins);
        AddEntry("Radial_profile_norm", opt.iprofilenorm);
        AddEntry("Radial_profile_bin_edges", opt.profile_bin_edges);
    }
    AddEntry("Number_of_overdensities", opt.SOnum);
    AddEntry("Overdensity_values_in_critical_density", opt.SOthresholds_values_crit);
    AddEntry("Spherical_overdenisty_calculation_limited_to_structure_types", (opt.SphericalOverdensitySeachMaxStructLevel-HALOSTYPE)/HALOCORESTYPE);

    //try removing index that is now stored in
    vector<string> name;
    name = StripIndexOffName(opt.gas_internalprop_names, opt.gas_internalprop_names_aperture,
        opt.gas_internalprop_index, opt.gas_internalprop_index_aperture);
    AddEntry("Gas_internal_property_names", name);
    name = StripIndexOffName(opt.gas_chem_names, opt.gas_chem_names_aperture,
        opt.gas_chem_index, opt.gas_chem_index_aperture);
    AddEntry("Gas_chemistry_names", name);
    name = StripIndexOffName(opt.gas_chemproduction_names, opt.gas_chemproduction_names_aperture,
        opt.gas_chemproduction_index, opt.gas_chemproduction_index_aperture);
    AddEntry("Gas_chemistry_production_names", name);
    name = StripIndexOffName(opt.star_internalprop_names, opt.star_internalprop_names_aperture,
        opt.star_internalprop_index, opt.star_internalprop_index_aperture);
    AddEntry("Star_internal_property_names", name);
    name = StripIndexOffName(opt.star_chem_names, opt.star_chem_names_aperture,
        opt.star_chem_index, opt.star_chem_index_aperture);
    AddEntry("Star_chemistry_names", name);
    name = StripIndexOffName(opt.star_chemproduction_names, opt.star_chemproduction_names_aperture,
        opt.star_chemproduction_index, opt.star_chemproduction_index_aperture);
    AddEntry("Star_chemistry_production_names", name);
    name = StripIndexOffName(opt.bh_internalprop_names, opt.bh_internalprop_names_aperture,
        opt.bh_internalprop_index, opt.bh_internalprop_index_aperture);
    AddEntry("BH_internal_property_names", name);
    name = StripIndexOffName(opt.bh_chem_names, opt.bh_chem_names_aperture,
        opt.bh_chem_index, opt.bh_chem_index_aperture);
    AddEntry("BH_chemistry_names", name);
    name = StripIndexOffName(opt.bh_chemproduction_names, opt.bh_chemproduction_names_aperture,
        opt.bh_chemproduction_index, opt.bh_chemproduction_index_aperture);
    AddEntry("BH_chemistry_production_names", name);
    name = StripIndexOffName(opt.extra_dm_internalprop_names, opt.extra_dm_internalprop_names_aperture,
        opt.extra_dm_internalprop_index, opt.extra_dm_internalprop_index_aperture);
    AddEntry("Extra_DM_internal_property_names", name);


    map<int, string> calcinttostring = ExtraFieldCalculationsIntFlagToString();
    vector<string> funcname;
    for (auto &x:opt.gas_internalprop_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.gas_internalprop_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("Gas_internal_property_calculation_type", funcname);
    funcname.clear();
    for (auto &x:opt.gas_chem_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.gas_chem_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("Gas_chemistry_calculation_type", funcname);
    funcname.clear();
    for (auto &x:opt.gas_chemproduction_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.gas_chemproduction_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("Gas_chemistry_production_calculation_type", funcname);
    funcname.clear();

    for (auto &x:opt.star_internalprop_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.star_internalprop_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("Star_internal_property_calculation_type", funcname);
    funcname.clear();
    for (auto &x:opt.star_chem_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.star_chem_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("Star_chemistry_calculation_type", funcname);
    funcname.clear();
    for (auto &x:opt.star_chemproduction_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.star_chemproduction_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("Star_chemistry_production_calculation_type", funcname);
    funcname.clear();

    for (auto &x:opt.bh_internalprop_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.bh_internalprop_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("BH_internal_property_calculation_type", funcname);
    funcname.clear();
    for (auto &x:opt.bh_chem_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.bh_chem_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("BH_chemistry_calculation_type", funcname);
    funcname.clear();
    for (auto &x:opt.bh_chemproduction_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.bh_chemproduction_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("BH_chemistry_production_calculation_type", funcname);
    funcname.clear();

    for (auto &x:opt.extra_dm_internalprop_function) funcname.push_back(calcinttostring[x]);
    for (auto &x:opt.extra_dm_internalprop_function_aperture) funcname.push_back(calcinttostring[x]);
    AddEntry("Extra_DM_internal_property_calculation_type", funcname);
    funcname.clear();

    AddEntry("Gas_internal_property_index", opt.gas_internalprop_index, opt.gas_internalprop_index_aperture);
    AddEntry("Gas_chemistry_index", opt.gas_chem_index, opt.gas_chem_index_aperture);
    AddEntry("Gas_chemistry_production_index", opt.gas_chemproduction_index, opt.gas_chemproduction_index_aperture);
    AddEntry("Star_internal_property_index", opt.star_internalprop_index, opt.star_internalprop_index_aperture);
    AddEntry("Star_chemistry_index", opt.star_chem_index, opt.star_chem_index_aperture);
    AddEntry("Star_chemistry_production_index", opt.star_chemproduction_index, opt.star_chemproduction_index_aperture);
    AddEntry("BH_internal_property_index", opt.bh_internalprop_index, opt.bh_internalprop_index_aperture);
    AddEntry("BH_chemistry_index", opt.bh_chem_index, opt.bh_chem_index_aperture);
    AddEntry("BH_chemistry_production_index", opt.bh_chemproduction_index, opt.bh_chemproduction_index_aperture);
    AddEntry("Extra_DM_internal_property_index", opt.extra_dm_internalprop_index, opt.extra_dm_internalprop_index_aperture);

    AddEntry("Gas_internal_property_output_units", opt.gas_internalprop_output_units, opt.gas_internalprop_output_units_aperture);
    AddEntry("Gas_chemistry_output_units", opt.gas_chem_output_units, opt.gas_chem_output_units_aperture);
    AddEntry("Gas_chemistry_production_output_units", opt.gas_chemproduction_output_units, opt.gas_chemproduction_output_units_aperture);
    AddEntry("Tcut_halogas", opt.temp_max_cut);
    AddEntry("Gas_temperature_input_unit_conversion_to_output_unit", opt.temp_input_output_unit_conversion_factor);
    AddEntry("Gas_hot_overdensity_normalisation", opt.hot_gas_overdensity_normalisation);
    AddEntry("Aperture_values_normalised_to_overdensity", opt.aperture_hotgas_normalised_to_overdensity);

    AddEntry("Star_internal_property_output_units", opt.star_internalprop_output_units, opt.star_internalprop_output_units_aperture);
    AddEntry("Star_chemistry_output_units", opt.star_chem_output_units, opt.star_chem_output_units_aperture);
    AddEntry("Star_chemistry_production_output_units", opt.star_chemproduction_output_units, opt.star_chemproduction_output_units_aperture);
    AddEntry("BH_internal_property_output_units", opt.bh_internalprop_output_units, opt.bh_internalprop_output_units_aperture);
    AddEntry("BH_chemistry_output_units", opt.bh_chem_output_units, opt.bh_chem_output_units_aperture);
    AddEntry("BH_chemistry_production_output_units", opt.bh_chemproduction_output_units, opt.bh_chemproduction_output_units_aperture);
    AddEntry("Extra_DM_internal_property_output_units", opt.extra_dm_internalprop_output_units, opt.extra_dm_internalprop_output_units_aperture);

    AddEntry("Gas_internal_property_input_output_unit_conversion_factors", opt.gas_internalprop_input_output_unit_conversion_factors, opt.gas_internalprop_input_output_unit_conversion_factors_aperture);
    AddEntry("Gas_chemistry_input_output_unit_conversion_factors", opt.gas_chem_input_output_unit_conversion_factors, opt.gas_chem_input_output_unit_conversion_factors_aperture);
    AddEntry("Gas_chemistry_production_input_output_unit_conversion_factors", opt.gas_chemproduction_input_output_unit_conversion_factors, opt.gas_chemproduction_input_output_unit_conversion_factors_aperture);
    AddEntry("Star_internal_property_input_output_unit_conversion_factors", opt.star_internalprop_input_output_unit_conversion_factors, opt.star_internalprop_input_output_unit_conversion_factors_aperture);
    AddEntry("Star_chemistry_input_output_unit_conversion_factors", opt.star_chem_input_output_unit_conversion_factors, opt.star_chem_input_output_unit_conversion_factors_aperture);
    AddEntry("Star_chemistry_production_input_output_unit_conversion_factors", opt.star_chemproduction_input_output_unit_conversion_factors, opt.star_chemproduction_input_output_unit_conversion_factors_aperture);
    AddEntry("BH_internal_property_input_output_unit_conversion_factors", opt.bh_internalprop_input_output_unit_conversion_factors, opt.bh_internalprop_input_output_unit_conversion_factors_aperture);
    AddEntry("BH_chemistry_input_output_unit_conversion_factors", opt.bh_chem_input_output_unit_conversion_factors, opt.bh_chem_input_output_unit_conversion_factors_aperture);
    AddEntry("BH_chemistry_production_input_output_unit_conversion_factors", opt.bh_chemproduction_input_output_unit_conversion_factors, opt.bh_chemproduction_input_output_unit_conversion_factors_aperture);
    AddEntry("Extra_DM_internal_property_input_output_unit_conversion_factors", opt.extra_dm_internalprop_input_output_unit_conversion_factors, opt.extra_dm_internalprop_input_output_unit_conversion_factors_aperture);

    //other options
    AddEntry("Verbose", opt.iverbose);
    AddEntry("Write_group_array_file",opt.iwritefof);
    AddEntry("Snapshot_value",opt.snapshotvalue);

    //io related
    AddEntry("Cosmological_input",opt.icosmologicalin);
    AddEntry("Input_chunk_size",opt.inputbufsize);
    AddEntry("MPI_particle_total_buf_size",opt.mpiparticletotbufsize);
    AddEntry("Separate_output_files", opt.iseparatefiles);
    AddEntry("Binary_output", opt.ibinaryout);
    AddEntry("Comoving_units", opt.icomoveunit);
    AddEntry("Extended_output", opt.iextendedoutput);

    //HDF io related info
    AddEntry("HDF_name_convention", opt.ihdfnameconvention);
    AddEntry("Input_includes_dm_particle", opt.iusedmparticles);
    AddEntry("Input_includes_gas_particle",opt.iusegasparticles);
    AddEntry("Input_includes_star_particle", opt.iusestarparticles);
    AddEntry("Input_includes_bh_particle", opt.iusesinkparticles);
    AddEntry("Input_includes_extradm_particle", opt.iuseextradarkparticles);
    AddEntry("Input_includes_wind_particle", opt.iusewindparticles);
    AddEntry("Input_includes_tracer_particle", opt.iusetracerparticles);

    //gadget io related to extra info for sph, stars, bhs,
    AddEntry("NSPH_extra_blocks", opt.gnsphblocks);
    AddEntry("NStar_extra_blocks", opt.gnstarblocks);
    AddEntry("NBH_extra_blocks", opt.gnbhblocks);

    //mpi related configuration
    AddEntry("MPI_part_allocation_fac", opt.mpipartfac);
#endif
    AddEntry("#Compilation Info");
#ifdef USEMPI
    AddEntry("#USEMPI");
#endif
#ifdef USEOPENMP
    AddEntry("#USEOPENMP");
#endif
#ifdef GASON
    AddEntry("#USEGAS");
#endif
#ifdef STARON
    AddEntry("#USESTAR");
#endif
#ifdef BHON
    AddEntry("#USEBH");
#endif
#ifdef EXTRADMON
    AddEntry("#USEEXTRADMPROPERTIES");
#endif
#ifdef HIGHRES
    AddEntry("#ZOOMSIM");
#endif
#ifdef SWIFTINTERFACE
    AddEntry("#SWIFTINTERFACE");
#endif

}
