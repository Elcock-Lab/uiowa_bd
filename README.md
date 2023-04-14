# *uiowa_bd*


## Overview
*uiowa_bd* is a parallelized program that performs Brownian dynamics (BD) simulations of macromolecules. It uses simple molecular mechanics models of the kind widely used in other simulation codes to model the internal degrees of freedom of molecules, but also has the ability to include hydrodynamic interactions (HIs) between atoms (or more generally, beads), calculated at the Rotne-Prager-Yamakawa level of theory. This makes it useful for accurately simulating the translational and rotational diffusion of macromolecules, as well as their associations, in a fundamentally implicit solvent model. 

## Example folders

A number of example directories are provided with the source code that illustrate different uses of *uiowa_bd* and that provide good starting points for understanding the formatting of the various required input files. Further examples may be added in the next few days, but for now we have the following:

1. `NATIVE_PROTEIN_DIFFUSION`

This folder contains all files necessary to run a BD-HI simulation of the small protein Chymotrypsin Inhibitor 2 (CI2) modeled at the C-alpha-only level of resolution. The energetic model used to describe the protein is essentially a Go model: equilibrium bond, angle, and dihedral angles are set to the values found in the native state structure, and favorable Lennard-Jones 12-10 potential functions are added to all pairs of beads that form a contact in the native state structure. The protein is simulated starting in its native state, and HIs are modeled at the Rotne-Prager-Yamakawa level of theory with the result that both the translational and rotational diffusion of the protein is realistically simulated. 

2. `PROTEIN_FOLDING`

This folder contains all files necessary to run a BD-HI simulation of the small protein Cold-shock Protein B (CSPB) modeled at the C-alpha-only level of resolution. Again, the energetic model used here is that of a Go model, but in this case the initial configuration of the protein is an unfolded one and the BD-HI simulation continues until the protein folds back to its native state. Using the jargon from the protein folding field, we characterize the extent of folding of each configuration of the protein using the term "Q", which represents the fraction of native contacts that are successfully formed in the given configuration. The simulation is set up to quit when the protein reaches a value of Q = 0.90.

## Referencing *uiowa_bd*
The following lists the principal papers that have marked the development of *uiowa_bd*. If you have to cite only one publication then it makes sense for this to be the most recent (Tworek & Elcock, 2023) since this paper coincides with the release of this version of the code. However, depending on what features of the code you use, you may need to cite additional publications from other people (see treecode and fixman entries below).

1. Tworek JW, Elcock AH **An orientationally averaged version of the Rotne-Prager-Yamakawa tensor provides a fast but still accurate treatment of hydrodynamic interactions in Brownian dynamics simulations of biological macromolecules.** (preprint). *bioRxiv.* 2023

2. Frembgen-Kesner T, Elcock AH (2009) **Striking effects of hydrodynamic interactions on the simulated diffusion and folding of proteins.** *J Chem Theory Comput* **5**:242-256

3. Elcock AH (2006) **Molecular simulations of cotranslational protein folding: fragment stabilities, folding cooperativity, and trapping in the ribosome.** *PLoS Comput Biol* **2**:e98

## External References and Contributions Made by Others
Most of the *uiowa_bd* code was written from the ground-up by AHE. However, certain key parts of the code were taken from other sources, and unfortunately most of those parts were incorporated so long ago that I am in some cases a bit hazy about their exact origins. I have done my best in what follows to properly attribute credit to others. Parts of the code that originate from elsewhere include:

1. code for handling bond angle and dihedral angle calculations was, if I remember correctly, adapted from code that I found online *many* years ago written by Prof. Jay Ponder (Washington University, St. Louis) and which probably formed part of his TINKER simulation package. For much more up-to-date versions of TINKER, please see Prof. Ponder's website: https://dasher.wustl.edu/tinker/

2. the treecode routine that is used for calculating long-range (Debye-Hückel) electrostatic interactions comes from the group of Prof. Robert Krasny (University of Michigan). One of the primary authors of that code – Hans Johnston – was extremely generous with his time in helping me adapt the code for use in *uiowa_bd*. If the treecode routine is used please cite:

    Li P, Johnston H, Krasny R (2009) **A Cartesian treecode for screened coulomb interactions.** *J Comput Phys* **228**:3858-3868  

3. code for writing trajectory coordinates to a compressed .xtc file came indirectly via GROMACS many years ago. A former member of my group, Dr Shun Zhu, figured out how to call that code from within *uiowa_bd*. He did this after having obtained code from two separate sources, both of whose URLs unfortunately appear now to be dead. One piece of the puzzle was the ego2xtc Fortran program which contained "writextc" and "readxtc" subroutines; that code was obtained via the following now-dead link: http://www.gromacs.org/Downloads/User_contributions/Other_software). The second piece of the puzzle was the xdrf library (written by Frans van Hoese as part of the EUROPORT project) which Shun obtained from the following now-dead link: http://hpcv100.rc.rug.nl/xdrfman.html. I cannot claim to understand how any of these routines work, but they clearly do what they are supposed to do when the final library file (`libxdrf.a`) is linked to *uiowa_bd*. Sorry, but I don't know how better to cite this part of the code. 

4. *uiowa_bd* implements two position-update algorithms. The first, invoked with the keyword `brownian` (see below) is the Brownian dynamics algorithm developed by Ermak and McCammon:

    Ermak DL, McCammon JA (1978) **Brownian dynamics with hydrodynamic interactios.** *J Chem Phys* **69**:1352.
   
The second, invoked with the keyword `langevin` is the closely-related Langevin dynamics algorithm developed by Geyer & Winter:

    Geyer T, Winter (2011) **Coarse grained simulations of a small peptide: Effects of finite damping and hydrodynamic interactions.** *J Chem Phys* **131**:104102.
    
5. the code that allows the late Prof Marshall Fixman’s Chebyshev polynomial-based method to be used to calculate correlated random displacements borrows very heavily from a corresponding C routine that was written by Tihamer Geyer when he was a faculty member at the University of Saarland and that was implemented in his BD code. If the Fixman code is used please consider citing:

    Geyer T (2011) **Many-particle Brownian and Langevin dynamics simulations with the Brownmove package.** *BMC Biophyics* **4**:7

6. while the current version of the code uses the Intel MKL routine `spotrf` to compute the Cholesky decomposition of the diffusion tensor, I want to acknowledge Dr Jonathan Hogg’s help in implementing an earlier openmp-parallelized routine for performing the same operation (HSL_MP54). It was Dr Hogg’s Cholesky decomposition code that enabled a number of our earlier studies with *uiowa_bd* to be completed.

    Hogg JD (2008) **A DAG-based parallel Cholesky Factorization for multicore systems.** Technical Report TR-RAL-2008-029

7. code for writing trajectory coordinates to movie .pdb files was mostly written by Dr Tyson Shepherd while he was rotating in my group many years ago.


## Installation and compilation
The bulk of the *uiowa_bd* source code is all contained in the single folder `SOURCE`; additional code that handles the reading and writing of .xtc trajectory files (and that I didn’t write: see above) is in the sub-folder `XTC`. A makefile is provided that “gets the job done”, but I don’t claim that this makefile is well-written: I barely understand how makefiles work, and I stopped refining the one provided as soon as it looked like it worked. All of the code is written in Fortran. I assume that the user will compile the code with Intel’s Fortran compiler (ifort) and with Intel’s Math Kernel Library (MKL) installed. The code can probably be adapted to compile with gfortran relatively easily, but care will be needed with routines that are currently handled by MKL: these include the calculation of random numbers, the spotrf routine that is used to compute the Cholesky decomposition of the diffusion tensor, and possibly some others. I am sorry to say that if you attempt to get the code working with any compiler other than ifort you will be on your own. 

**Before compiling you will need to do the following:**

1. add ifort to your `PATH` ; **for example**, with our very old installation of ifort we use the following:

`export PATH=$PATH:/opt/intel/compilers_and_libraries_2016.2.181/linux/bin/intel64`

2. add the ifort libraries to `LD_LIBRARY_PATH` ; **for example:**

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/compilers_and_libraries_2016.2.181/linux/compiler/lib/intel64`

3. add the environment variable `MKLROOT` ; **for example:**

`export MKLROOT=/opt/intel/compilers_and_libraries_2016.2.181/linux/mkl`


**Now do the actual compilation in three stages:**

1. if it doesn’t already exist, then make a library file (`libxdrf.a`) that allows *uiowa_bd* to handle .xtc files:

`cd XTC_IFORT ; make clean ; make ; cp libxdrf.a ../ ; cd ../`

2. if it doesn’t already exist, then copy the appropriate include file for MKL so that uiowa_bd can use it for random numbers:

`cp ${MKLROOT}/include/mkl_vsl.f90 .`

3. compile modules and then make uiowa_bd (may take a minute or so):

`ifort -c mkl_vsl.f90 ; ifort -c uiowa_bd.modules.f ; make uiowa_bd`

note that the first command compiles the mkl include file copied in stage 2.

**The resulting executable will be:**

`uiowa_bd.exe`

### A few comments on the development of *uiowa_bd*
The *uiowa_bd* source code has been written over a number of years and there have been many cases where features (e.g. replica exchange) have been added to the code and then deprecated before ever appearing in publication form. I have done my best to remove “dead” code but there is still likely to be some remaining. The present version represents a near-final version that, while definitely useful for production level simulations now, is unlikely to be a focus for further serious development in my group. This is for the following reasons. First, while the *uiowa_bd* code is effective for small- to medium-sized systems it is not well suited to simulating very large-scale systems that are starting to become of interest to my lab. Second, some of the decisions about code-structure that were made early in the code’s development would not be made today: the code is *far* more complicated than it needs to be in places, and this hampers efforts to build in fundamentally new approaches. Third, while the code typically achieves nice speedups using openmp, its raw (single-core) speed is not what it could be: a good chunk of this is likely attributable to the non-optimal way in which much of the data is stored in memory and accessed during calculations. Fourth, the existing code is entirely unaware of the possibility of using GPUs for compute-intensive calculations, which makes it seem something of a dinosaur now. While all of the above issues are, in principle, quite fixable, the underlying code-structure of *uiowa_bd* is sufficiently byzantine that a better approach is likely to be starting from scratch with an entirely new code.

A brief note on the language used in *uiowa_bd* is also in order: While the code makes use of a number of Fortran90 constructs (especially allocatable arrays plus the odd pointer and derived type here and there), it is written in the style of Fortran77. This is why all of the source file extensions are “.f” and not “.f90”. A professional programmer would almost certainly laugh at the way *uiowa_bd* is written. While some readers of this might be tempted to undertake reformatting and rewriting of the code to make it fully Fortran90-compliant, I don’t recommend this given that the code is unlikely to undergo any further serious development by me (see above).


## Using *uiowa_bd*: running the code
As with all simulation codes, a significant amount of time and effort will be spent setting up all of the input files necessary for performing a simulation (see below). Once that is done, the code is typically run in the following fashion:

`mkdir MOVIE ; mkdir RESTARTS`

`cp restart_file_that_you_prepared_earlier restart.file.001`

`./uiowa_bd.exe 1234 < uiowa_bd.inp > uiowa_bd.out &`

`tail -f uiowa_bd.out`


**Line 1** makes two sub-folders into which .pdb and restart.files will be written - if these folders do not exist at run-time then the code will basically hang, continually stating that it is trying to write a .pdb file.

**Line 2** copies a restart.file into the hardwired filename that *uiowa_bd* expects and that contains the initial coordinates of all atoms/beads in the system – trust me on this: this file **must** be called restart.file.001

**Line 3** runs the code: the first argument is a random seed. To obtain different replicate simulations of the same system just change this number (I typically use numbers like 1234, 2345 etc). While two otherwise identical runs that use the same random seed should, in principle, produce identical output, in practice, this is only likely to occur when running with a single openmp thread: when more than one thread is used, slight numerical differences will accumulate and cause the trajectories to begin to diverge.

**Line 4** just pipes all of the text output from *uiowa_bd* to the screen as it is written – this is useful when first running a simulation to make sure that it is working as expected, i.e. not producing insane energies or blowing up.

## Using *uiowa_bd*: summary of input

All input files and their required formats are described in detail later; probably the best way to get a "feel" for these inputs is to look in the various EXAMPLE folders that are provided. At the most basic level, however, the following files are required as inputs to *uiowa_bd*:

1. An **input file** (e.g. uiowa_bd.inp) that lists all of the details of the system to be simulated. Conceptually, this file combines features that are found in GROMACS’ grompp.mdp and topol.top files, but it uses a very different format.

2. A **parameter file** (e.g. parameter.file.01) that lists nonbonded van der Waals parameters for all atom types in the system. These essentially take the form of sigma and epsilon values.

3. A **restart file** (*always* named restart.file.001) that lists the initial coordinates of all beads in the system; this file is continually overwritten as a *uiowa_bd* simulation progresses so users must remember to keep a copy of this file prior to running a simulation.

4. For each molecule type simulated in the system, two additional files need to be provided (together, these two files are broadly equivalent to a GROMACS .itp file):

	a. a **charge parameters** file – this is a .pdb-like file that lists all atoms/beads and provides their net charges and their hydrodynamic radii – the file also contains coordinates but these are not used so can be given garbage values if desired
	
	b. an **internal parameters** file – this file identifies all bonds, bond angles, dihedral angles, and improper dihedral angles, and provides their parameters.

5. (optional) a **go parameters** file – this file lists all pairs of atoms whose contacts should be energetically rewarded during a simulation. This file is only needed when one or more molecules in the system are maintained in folded states by non-covalent Lennard-Jones (LJ) interactions that are specific to listed pairs of contacting atoms. Favorable contacts are most usually included between atom pairs that are both part of the same molecule type; this makes it possible to model the folding of molecules using a so-called Go model. Favorable contacts can also, however, be included between  atoms in different molecule types, making it possible to simulate the formation of hetero-oligomeric complexes. If none of the molecules in the system are modeled in this way – either because they are all intended to be unfolded/disordered or because they are all maintained in their folded states by covalent bonds (e.g. in an elastic network model (ENM)) – then this file can be empty.

## Using *uiowa_bd*: summary of output

The code generates four main outputs:

1. a text log file to which all messages issued by uiowa_bd are written and which writes the system energy and its components at user-specified intervals.

2. a .xtc file that stores the coordinates of all atoms, compressed, and measured in nm, written at user-specified intervals. In most applications, this will be the file that will be of most interest.

3. .pdb files that store the coordinates of all atoms, uncompressed and measured in Ångstroms, written at user-specified intervals to a sub-folder called MOVIE

4. restart files that can be used as starting points for additional uiowa_bd simulations, written at user-specified intervals to a sub-folder called RESTARTS

In typical usage, when I wish to watch a movie of a simulation (usually to make sure that it is not going crazy) I will do the following. First, I take the very first .pdb file that uiowa_bd writes to the MOVIE sub-folder, and I read it into VMD for viewing. This .pdb file has the advantage of listing all covalent bonds in the system with CONECT statements: this makes visualization of very coarse-grained systems easy as VMD’s default representaton mode (“lines”) automatically shows who is bonded to who. Second, I take the testout.001.xtc file that is generated by uiowa_bd, and I read it into VMD using “Load Data Into Molecule” tab: this adds the .xtc coordinates to the .pdb that was previously read in.

## Using *uiowa_bd*: some idiosyncracies, features, and workarounds

### How periodic boundary conditions work (or don't) in *uiowa_bd*
I am sorry to say that the implementation of periodic boundary conditions in *uiowa_bd* is complicated. First of all, let's be clear about what might happen when a simulation is run without periodic boundary conditions (i_pbc=0; see below). If molecules are unconstrained by walls, positions restraints, or capsule restraints, then there is the possibility that, given enough simulation time, they will diffuse outside of the box limits specified by xmin, xmax, ymin, ymax, zmin, zmax (see below). If that happens then the simulation will definitely crash the next time that the list of nonbonded bead pairs is updated. You *could* decide to make the box very large to try to avoid this happening, but bear in mind that the grid used to determine all nonbonded neighbors will also become quite large and zeroing it out at the beginning of each nonbonded update might then become expensive (this is a good example of non-optimal programming in *uiowa_bd*; see above).

An alternative is to use periodic boundary conditions (i_pbc=1). In understanding how these are implemented in *uiowa_bd*, it's important to note the following. **Internally**, *uiowa_bd* keeps track of the "true" coordinates of beads, i.e. the coordinates that properly describe the diffusive trajectory that beads make away from the initial coordinates passed to the program. Keeping these "true" coordinates makes it possible to easily calculate translational diffusion coefficients from trajectories later. **Externally** (by which we mean how the coordinates appear in .pdb and .xtc files), the behavior is controlled by the **wrap_molecules** flag (see below): users have the option to write out the true coordinates, or coordinates of the periodic image that lies within the central simulation box. 

With that preamble in mind, let's first deal with how periodic boundary conditions are treated in simulations that do not include HIs. For BD simulations without HI, periodic boundary conditions work fine: the "true" coordinates of beads can continue to diffuse away from the central box but all nonbonded interactions are determined using the minimum image convention so beads interact with the nearest image of all other beads in the system (so long as they are within the cutoff distance). For simulations with HIs, the situation is more complicated. The correct way to handle periodic boundary conditions in BD-HI simulations is to use an Ewald summation similar to what is used to describe electrostatic interactions. Once upon a time I programmed the Ewald summation of the RPY diffusion tensor that was derived by Beenakker in 1986. That code never ended up in a publication, even though it worked, and since I didn't have immediate use for it, I eventually stopped maintaining it. So now the treatment of HIs in periodic boundary condition simulations is, shall we say, non-optimal. For BD-HI simulations that include only a single molecule, periodic boundary conditions can be used as long as the box is large enough that we never have some atoms whose closest neighbors are in one box, and other atoms whose closet neighbors are in another box. If that happens, then the RPY diffusion tensor will for sure end up non-positive definite and it will kill the simulation. If, however, the box is large enough so that all atoms of the molecule interaction within the original version of the molecule then the RPY diffusion tensor will remain positive definite. 
 
### Walls, position restraints, and confinement potentials
As is detailed in the section describing the input options, *uiowa_bd* offers a variety of different restraint potentials that can be used to restrict the movement of beads of interest. Each of these functional forms has advantages and disadvantages. Fixed walls allowed by the code can be used to restricted selected subsets of beads to different spatial locations, and could in principle be used to determine osmotic pressures (the code did this at one time but that never reached publication). Fixed position restraints allowed by the code can cover a wide range of different scenarios: beads can be restricted to 1D, 2D or 3D positions, and can be restricted to a variety of locations relative to a capsule whose shape can be either a sphere or a bacterial-cell shape. Specifically, beads can be restricted to the surface of the capsule, the inside of a capsule or the outside of a capsule; the latter two potentials can also be combined to restrict beads to a shell of any desired thickness. Finally, *moveable* confinement potentials are also allowed by the code: while these are more functionally limited in the sense that all they do is restrict beads to remain within a sphere or a capsule, they have the advantage that they allow the radius of the sphere or capsule to vary during a simulation. This can, for example, be used to progressively shrink a large molecule so that it fits within a cell-like volume, or could be used to confine a viral genome within a sphere commensurate with its capsid.


## Overview of file formats
With only one or two exceptions (see below) all of the inputs that are provided to the code are expected to be in fixed format. This means that you should *never* mess around with the alignment of columns in files, or skip lines etc. as I cannot vouch for the behavior that will result. In what follows I will use Fortran’s description of integer and float types to specify the format – e.g. `a3` means a word (string) of 3 characters, `f15.5` means a real number that has a total of 15 characters, 5 of which come after the decimal point; `i8` means an integer that has a total of 8 characters.

Molecule-specific file formats: 1. charge.parameters file

Molecule-specific file formats: 2. internal.parameters file

### parameter file format:

This contains one line per atom type present in the system. The atom type is dictated by the atom name provided in the charge.parameters files

### uiowa_bd input file format:

Input files for a number of example situations are provided in the EXAMPLES folder. In what follows, all of the parameters listed in these input files are grouped together by the line that they appear on in the input file

---

**teprint** (ps) : time interval at which energies are written to the uiowa_bd output file

**ttprint** (ps) : time interval at which system coordinates are written to testout.001.xtc file

**tmprint** (ps) : time interval at which system coordinates are written to MOVIE/movie….pdb files

**num_lst_stp** : # of timesteps between update of nonbonded list

**num_fmd_stp** : # of timesteps between update of medium-range forces

**num_hyd_stp** : # of timesteps between update of diffusion tensor 

---

**num_threads** : # of openmp threads

**bond_dev_quit** (A) : quit if a bond deviates from its equilibrium value by more than this value

**i_do_lincs?** : if “yes” then will use lincs to constrain bonds – do not use this with HI (see below)

**i_do_YHL?** : if “yes” then use Newton’s third law to skip j:i force calculation if i:j is calculated ; this can help speed up simulations with low numbers of openmp threads 

---

**f_typs** : total number of molecule types in the system

**f_mols** : total number of molecules in the system

**i_debug?** : if “yes” then write out lots of information to the uiowa_bd output file ; don’t use this normally

**q_desired** : if simulating a folding or association event using Go potentials, then quit when this Q-value is reached – since Q must, by definition, be between 0 and 1 setting a value of 1.1 means that the simulation will never quit due to folding.

**mol_Q1** : # of the first molecule whose Go contact pairs we are monitoring

**mol_Q2** : # of the second molecule whose Go contact pairs we are monitoring

if mol_Q1 = mol_Q2 then we are monitoring the formation of intramolecular contacts, i.e. a folding event
if mol_Q1 <> mol_Q2 then we are monitoring the formation of intermolecular contacts, i.e. an association event

**go_eps_low** : minimum value of epsilon allowed for Go contact pairs that contribute to the Q-value calculation. This can be used to focus the Q-value on contacts that are assigned more favorable potential well-depths.

---

**xmin, xmax, ymin, ymax, zmin, zmax** (A) : min and max extents of the simulation box – I always set these so that they are symmetric about the origin so please do the same

**i_pbc** : =1 to use periodic boundary conditions ; =0 to not use periodic boundary conditions

**i_look_for_crashes?** : if “yes” then monitor bonds with bond_dev_quit

**i_limit_verbosity?** : if “yes” then don’t write out too much garbage to the uiowa_bd output file

**i_use_v_typ** : if “1” Lennard-Jones 12-10 potential functions are used throughout…

**steepest_descent?** : if “yes” then do a steepest-descent energy minimization instead of BD; by definition, no HIs are included and the timestep (specified below) serves as a multiplier to convert force into a displacement

---

**temperature** (K)

**ionic_strength** (mM)

**r_ion** (A) : radius of ion used in Debye-Hckel calculations ; since the treecode electrostatic routine assumes that this is zero then I would stick to using r_ion = 0.0

**dielectric** : relative dielectric constant (78.40 for water at 298K)

**viscosity** (cP)

**r_f_st**: if set to a positive number, allows the 1/r12 repulsion to be replaced by a much shallower harmonic repulsion but only when the distance between two atoms is less than the distance at which their LJ potential is most favorable ; r_f_st functions as the force constant for the harmonic repulsion. If r_f_st is set <= 0 then the usual 12-10 LJ potential is calculated.

---

**parameter_file** : name of parameter file containing Lennard-Jones parameters

**no_elec?** : if “yes” then skip electrostatic calculations entirely

**wrap_molecules** : controls behavior of molecules’ coordinates as written to .pdb and .xtc files when periodic boundary conditions are used

if “0” then coordinates are not wrapped - i.e. the "true" coordinates are written
if "1" then coordinates are wrapped on an atom-by-atom basis: for each atom, the coordinates written are those of the periodic image that lies within the box
if "2" the coordinates are wrapped on a molecule-by-molecule basis: for the first atom in each molecule, the coordinates written are again those of the periodic image that lies within the box ; for all other atoms in the molecule the coordinates written are those of the periodic image that is closest to that written for the first atom. In other words, this option endeavors to keep all molecules "whole" even as they might appear to jump from one side of the box to the other during the course of the simulation trajectory.
  
**HI_mode** : if “none” then no HI ; if “RPY” then use Rotne-Prager-Yamakawa HI ; if “OARPY” then use orientationally averaged RPY

**scale_nb** : DEPRECATED NEED TO DELETE LOTS OF HSL STUFF

**BD/LD** : if “brownian” use Ermak-McCammon algorithm ; if “langevin” use Geyer & Winter’s Langevin dynamics algorithm

---

**go_potentials_file** : name of file storing all Go contact pairs and their parameters – if you are not using Go contact pairs then just provide a junk file name here

**i_use_go_pairs?** : if “yes” then read the go_potentials_file and calculate Go pairs between appropriate atoms during the simulations; if “no” then go_potentials_file is NOT read, and all nonbonded interactions are treated as regular atoms, i.e. with LJ 12-10 potentials and electrostatic interactions

**go_primacy** : controls how Go contact potentials are combined with electrostatic interactions when atoms in a contact also bear charges ; if =1 then we only use LJ 12-10 potential function ; if =2 then we use LJ 12-10 potential function and electrostatic interaction. 

**i_skip_intermol_go?** : if “yes” then only consider Go pairs between atoms in same molecule ; all intermolecular interactions are treated as regular Lennard-Jones

**nomove_file** : name of file that lists all atoms that don’t move – they are listed one on each line – each says moltype and atomnumber

**i_dont_move_some?** : if “yes” then we will read the nomove_file to find atoms that are static. We assume that all moving atoms are listed before all static atoms.

---

for each of the “f_typs” molecule type in the system, provide:

1. name of **charge.parameters** file for this molecule type
2. name of **internal.parameters** file for this molecule type
3. number of copies of this molecule type

note that total number of all copies of all molecule types must equal **f_mols**

---

**time_step** (ps) : simulation time step

**totsimtime** (ps) : total length of the simulation

**vdw_s** (A) : short-range Lennard-Jones cutoff – interactions recalculated every step

**vdw_m** (A) : medium-range Lennard-Jones cutoff – interactions with distance > vdw_s but < vdw_m are recalculated every num_fmd_stp steps and kept constant for all intervening steps

**goe_s** (A) : short-range cutoff for Go contact pairs – interactions recalculated every step

**goe_m** (A) : medium-range cutoff for Go contact pairs – interactions with distance > goe_s but < goe_m are recalculated every num_fmd_stp steps and kept constant for all intervening steps

**ele_s** (A) : short-range electrostatic cutoff – interactions recalculated every step

**ele_m** (A) : medium-range electrostatic cutoff – interactions with distance > ele_s but < ele_m are recalculated every num_fmd_stp steps and kept constant for all intervening steps

**f_f_cell_size** (A) : size of cell used to accelerate construction of nonbonded pair list ; atom/bead pairs in all neighboring cells are examined for possible interactions – f_f_cell_size must be equal to or greater than the largest of the six cutoffs identified above *and* the box length in each of the x, y, and z dimensions should be an integer (>2) multiple of f_f_cell_size. For example, if xmin is -100A and xmax is 100A, such that the box length in x is 200 A, then possible values of f_f_cell_size would be 50A, 40A, 25A, but not 100A.

---

**position_restraint_file** : file listing any position restraints applied to atoms – if you are not using any position restraints then just provide a junk file name here

**i_do_pos_restraints?** : if “yes” then read and use the position_restraint_file

**mission_creep?** : if “yes” then the reference positions used for position_restraints are allowed to “creep” DEPRECATED REMOVE

---

**r_size** (A) : radius of encapsulating sphere centered on the origin ; if <0 then there is no encapsulating sphere

**l_size** (A) : length of encapsulating cylinder, which is assumed to be capped by hemispheres of radius r_size at either end ; if <0 then there is no encapsulating cylinder

**r_size_fac** ; factor by which to scale r_size and l_size by – only applicable if n_size > 0

**n_size** : number of steps over which to linearly scale radius back from its initial value to r_size – this works for a shrinking sphere, a shrinking capsule, and a shrinking 3D box

**f_size** (kcal/mol/A2) : force constant for the capsule-wall interaction

---

**fixman?** : if “yes” then use Fixman’s Chebyshev polynomial method instead of Cholesky decomposition to calculate correlated random displacements during BD-HI simulations

**fixman_tol** : tolerance on error estimate (Ek) used to determine convergence of the random displacements ; smaller values will require more steps to converge

**fixman_order** : # of terms to include in the Chebyshev polynomial method

**fixman_override?** : if “yes” then override the true minimum and maximum eigenvalues (which will otherwise be calculated by uiowa_bd) with the following estimates (**lmin** & **lmax**)

---

**treecode?** : if “yes” then use the Krasny group’s electrostatic treecode to recalculate all long-range electrostatic interactions at intervals of num_lst_stp steps

**theta** : multiple acceptance criterion used to decide whether to use a multipole expansion in place of the exact calculation

**order** : order of multipole expansion

**shrink** : 1

**maxatm** : max # of atoms to do something with

---

**walls?** : if “yes” then apply walls in one or more directions

**num_walls** : number of walls to apply

**wall_file** : file listing, for all atoms in the system, which of the walls they experience

(for each of the num_walls walls read:
)

## Formats of optional files:
### position_restraint_file (fixed format):
if i_do_position_restraints = “yes” then we need to provide a file that lists all of these restraints. There is one line for every restraint: if a bead is subject to no position restraints then there is no need to list it in the file; if a bead is subject to three different types of position restraint then three lines will be required to specify each of the restraints. Note also that if there are multiple copies of the same molecule type present in the system then you will need to provide separate lines for every copy (sorry about that but it is easily scripted).

**restraint type 1**: To harmonically restrain a potential to a 1D line, a 2D plane, or a 3D point list the following:

“1”,atom number, x-coord of restraint, y-coord of restraint, z-coord of restraint, force constant in x, force constant in y, force constant in z.

Format: `(a1,i8,3f20.5,3f15.5)`

To restrain an atom in 1D, e.g. along the y-axis, then set fy = 0 and fx and fz to non-zero values 

To restrain an atom in 2D, e.g. in the y-z plane, then set fy = 0, fz = 0, and fx to a non-zero value

To restrain an atom in 3D set fx, fy, fz to non-zero values.

**restraint type 2**: To harmonically restrain an atom to a specified radial distance from a point specified in 3D space, list the following:

“2”,atom number, x-coord of restraint, y-coord of restraint, z-coord of restraint, desired radial distance, radial force constant.

Format: `(a1,i8,3f20.5,2f15.5)`

**restraint type 3**: To harmonically restrain an atom to a specified radial distance around the edge of a capsule whose long-axis is aligned with the z-axis, list the following:

“3”,atom number, xyz coordinates for the center of the left-hand hemisphere, xyz coordinates for the center of the right-hand hemisphere, desired radial distance, radial force constant.

Format: `(a1,i8,6f20.5,2f15.5)`

**restraint type 4**: To restrain an atom so that it remains within a specified radial distance around the edge of a capsule whose long-axis is aligned with the z-axis, list the following:

“4”,atom number, xyz coordinates for the center of the left-hand hemisphere, xyz coordinates for the center of the right-hand hemisphere, desired radial restraint distance, radial force constant.

Format: `(a1,i8,6f20.5,2f15.5)`

**restraint type 5**: To restrain an atom so that it remains outside a specified radial distance around the edge of a capsule whose long-axis is aligned with the z-axis, list the following:

“5”,atom number, xyz coordinates for the center of the left-hand hemisphere, xyz coordinates for the center of the right-hand hemisphere, desired radial restraint distance, radial force constant.

Format: `(a1,i8,6f20.5,2f15.5)`

### wall file (free format):
If num_walls > 0 then a wall_file needs to be provided. This file MUST contain one line for every type of unique atom in the system; each line identifies which of the walls in the system are “seen” by the bead in question. Consider an example system that contains two types of molecules. The first type of molecule contains 3 atoms; the second type of molecule contains 2 atoms. Both types of molecule can be present in many copies, but we only list their atoms once. Now let’s also imagine that there are four walls specified in the system, so num_walls = 4 (see above). The wall file therefore needs to contain 5 rows and 6 columns: 5 rows because there are 5 unique types of atom in the system (3 from molecule type 1 ; 2 from molecule type 2) and 6 columns because we need to specify: the molecule type, the atom number, and a 0 or 1 flag for each of the 4 walls.

`1  1  0  0  1  0`

`1  2  0  0  0  1`

`1  3  0  0  1  1`

`2  1  1  1  0  0`

`2  2  0  0  0  0`

In this example, atom #1 of molecule type #1 “sees” wall #3 ; atom #2 of molecule type #2 “sees” wall #4 ; atom #3 of molecule type #1 “sees” walls #3 and #4 ; atom #1 of molecule type #2 “sees” walls #1 and #2 ; atom #2 of molecule type #2 “sees” no walls and so is free to diffuse freely throughout the entire system. 

### go_parameter_file (fixed format):

If i_use_go – “yes” then we will be reading the specified go_parameter_file. This is a fixed format file that specifies which atom types have favorable contacts with other atom types; atom types that form no such contacts do not need to be entered. Here is an example:
1  1  0  0  1  0 
1  2  0  0  0  1
1  3  0  0  1  1
2  1  1  1  0  0
2  2  0  0  0  0
In this example, atom #1 of molecule type #1 “sees” wall #3 ; atom #2 of molecule type #2 “sees” wall #4 ; atom #3 of molecule type #1 “sees” walls #3 and #4 ; atom #1 of molecule type #2 “sees” walls #1 and #2 ; atom #2 of molecule type #2 “sees” no walls and so is free to diffuse freely throughout the entire system. 











 
