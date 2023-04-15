                                                                 
      MODULE allocatable_arrays

C variables to remember:

C num_i_as ; ibeg_i_a ; iend_i_a : atoms that engage in nonbonded interactions

C i_res_no_nb(1:num_flx_atms)  =  0 means evaluate nonbond inters as normal
C                          = +i skip evaluation of nonbond inters for this atm
C                            'interacting' with another atm in same group
C                             note that it only makes sense for this to be done
C                             together with position restraints
C                          = -i only do evaluation of nonbond inters for this atm
C                            'interacting' with another atm in same group

C 0D arrays

      data pi/3.14159265/
      data twopi/6.28318531/

C --------------------------------------------------------
        logical      :: arbitrary_intra ! allow arb funcs on intra terms
C                       not yet implemented for high_mem...
        logical      :: periodic_bonds_okay !
        logical      :: replica_exchange
        logical      :: i_do_YHL_nonbondeds
        logical      :: i_do_rods
        logical      :: i_do_pH       
        logical      :: i_do_lincs    
        logical      :: i_do_direct_list
        logical      :: i_can_proceed
        logical      :: i_compare_go_with_others
        logical      :: i_continue_after_problem
        integer      :: sum_of_hist
        integer      :: i_append_movie  
        logical      :: i_limit_verbosity
        integer      :: i_need_to_die_now ! use with i_look_for_crashes
        logical      :: i_have_go_modes      
        logical      :: i_do_NAM ! northrup-allison-mccammon
        logical      :: i_read_bond_functions
        logical      :: i_read_nonbond_functions
        logical      :: i_use_alpha_energy    
        logical      :: i_write_user_histogram
        logical      :: i_write_bond_histogram
        logical      :: i_write_nonbond_histogram
        logical      :: i_read_ref_hist            
        logical      :: i_skip_intermol_go         
        logical      :: implicit_ions   
        logical      :: cutoff_elec     
        logical      :: ewald_elec
        logical      :: ewald_elec_grid
        logical      :: ewald_hydro
        logical      :: ewald_hydro_grid
        logical      :: treecode_elec
        logical      :: fixman     
        logical      :: cholesky   
        logical      :: fixman_override   ! overrides calc of eigenvals
C --------------------------------------------------------
        logical      :: i_do_growth ! determines whether mols grow
        logical      :: i_do_slide  ! determines whether mols slide 
        logical      :: i_do_protease  ! determines whether mols slide 
        logical      :: umbrella    ! use umbrella sampling
        logical      :: B22         ! use B22 sampling
        logical      :: B22_density_calc ! calculate a 3D dx file during B22 sampling
        logical      :: brownian
        logical      :: langevin
        logical      :: dpd
        logical      :: walls       ! if this is on then we have wall(s)
        logical      :: go_spline   ! read go from file and use spline  
        logical      :: afm         ! implement an afm tip              
        logical      :: cylinder
        logical      :: sphere      
        logical      :: shrinking_box
        logical      :: i_use_high_mem               ! determines whether we're using a more complete hydro treaent    
        logical      :: i_am_not_go                  ! determines whether we're using a more complete hydro treaent    
        logical      :: no_hydrodynamics           ! determines whether we're using a more complete hydro treaent    
        logical      :: full_hydrodynamics           ! determines whether we're using a more complete hydro treaent    
        logical      :: rpy_hydrodynamics           ! determines whether we're using a more complete hydro treaent    
        logical      :: oarpy_hydrodynamics         ! orientationally averaged HI
        logical      :: nodiv_hydrodynamics         ! same as above but without divergence term added
        logical      :: mixed_hydrodynamics           ! determines whether we're using a more complete hydro treaent    
        logical      :: tree_hydrodynamics           ! determines whether we're using a more complete hydro treaent    
        logical      :: intra_hydrodynamics           ! determines whether we're using a more complete hydro treaent    
        logical      :: multi_hydrodynamics           ! determines whether we're using a more complete hydro treaent    
        logical      :: diagonal_hydrodynamics       ! determines whether we're using a more complete hydro treaent    
        logical      :: integer_arrays               ! determines whether we're using a more complete hydro treaent    
        logical      :: hydro_cutoff                 ! determines whether we're using a more complete hydro treaent    
        logical      :: i_use_hydro                  ! determines whether we're using a more complete hydro treaent    
        logical      :: update_nonbond_list   ! forces a nonbonded list update when new as are grown in             
        logical      :: i_debug                          ! switches on various write-statements                                      
        logical      :: i_do_pos_restraint               ! denotes whether we include extra position restraints                      
        logical      :: mission_creep                    ! denotes whether we include extra position restraints                      
        logical      :: i_do_reactions                   ! denotes whether we're monitoring 'reactions' between sected as
        logical      :: i_grow_rs                     ! denotes whether we're growing rigid molecules
        logical      :: i_use_go_pairs                   ! denotes whether we're using Go potentials
        logical      :: i_use_exclusive_go ! go-pairs are exclusive for mol pairs and domain types
        logical      :: i_do_pka_calcs                   ! denotes whether we're recomputing pKas and ionization states
        logical      :: i_use_no_elec                    ! denotes whether we're omitting ectrostatic forcess
        logical      :: i_dont_move_some                 ! denotes whether we're skipping forces between certain molecules
        logical      :: i_omit_some_forces               ! denotes whether we're skipping forces between certain molecules
        logical      :: uniform_moves                    ! denotes whether we use a unifrom prob distribution for random moves
        logical      :: steepest_descent                 ! denotes whether we do a steeepest-descent-typ minimization                           
        logical      :: flex                             ! denotes whether there are flexible molecules in the system
        logical      :: rigd                             ! denotes whether there are rigid molecules in the system
        logical      :: linkers                          ! denotes whether there are linkers in the system
        integer      :: nref_oversample
        integer      :: wrap_molecules                   ! denotes whether we wrap coords in movie files:
C                       =0 no wrapping
C                       =1 atom-based wrapping
C                       =2 molecule-based wrapping

        integer      :: go_primacy ! =1 go only; =2 go + elec ; =3 all

C 2023 v1.1 - add fold_mode = 1 for folding ; =-1 for unfolding

        integer      :: fold_mode

C for mixed_hydrodynamics

        integer,allocatable :: atm_to_hyd_atm(:)! converts tot atm to hyd atm
        integer,allocatable :: atm_to_hyd_sys(:)! converts tot atm to hyd sys
        integer,allocatable :: hyd_atm_to_atm(:,:)! converts hyd atm to tot atm

        integer             :: mov_tot    ! (glob) # of moves on all threads
        integer,allocatable :: i_am_hyd(:)! (glob) =1 if atom has HI
        integer,allocatable :: mov_nod(:) ! (glob) node # (for MPI)
        integer,allocatable :: mov_thr(:) ! (glob) thread #
        integer,allocatable :: mov_typ(:) ! (glob) flex or rigd?
        integer,allocatable :: mov_beg(:) ! (glob) first atom/mol
        integer,allocatable :: mov_end(:) ! (glob) last atom/mol
        integer,allocatable :: mov_dyn(:) ! (glob) BD or LD?
        integer,allocatable :: mov_hyd(:) ! (glob) HI or noHI?
        integer,allocatable :: mov_lnc(:) ! (glob) lincs or not?
        integer,allocatable :: mov_rep(:) ! (glob) replica #
        integer,allocatable :: mythrds_rep_orig(:) ! (glob) replica #
        integer,allocatable :: mythrds_rep_curr(:) ! (glob) replica #
        integer,allocatable :: mythrds_rep_temp(:) ! (locl) replica #

        real*8              :: ene_ave        ! <E>
        real*8              :: ene_var        ! <E.E>
        real*8              :: ene_den        ! denominator
        real*8              :: qqq_ave        ! <E>
        real*8              :: qqq_var        ! <E.E>
        real*8              :: qqq_den        ! denominator
        real*8              :: rand_ave        ! <R>
        real*8              :: rand_msd        ! <R.R>
        real*8              :: rand_den        ! denominator

        integer             :: iaccept_rex     ! replica exchange
        real                :: rand_rex(1)     ! replica exchange
        real                :: boltz_rex       ! replica exchange
        real                :: deltaE_rex      ! replica exchange
        integer             :: ireptogl      ! replica exchange
        integer             :: ireptoglmax   ! replica exchange
        integer,allocatable :: irepswap(:,:)   ! replica exchange
        integer,allocatable :: irep_go(:,:)    ! replica exchange
        real,   allocatable :: rep_go_scale(:) ! replica exchange
        real,   allocatable :: ene_g_rx_old(:) ! replica exchange
        real,   allocatable :: ene_g_rx_new(:) ! replica exchange
        real,   allocatable :: f_g_rx_old(:,:,:) ! replica exchange
        real,   allocatable :: f_g_rx_new(:,:,:) ! replica exchange
        real,   allocatable :: c_a_x_rex(:,:) ! replica exchange                             
        real,   allocatable :: f_a_x_rex(:,:) ! replica exchange                             
        real                :: ene_g_rx_old_loc
        real                :: ene_g_rx_new_loc
        integer,allocatable :: ida_ff_g_rx(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_g_rx(:)             ! a pair f_f_list  array
        integer,allocatable :: kda_ff_g_rx(:)             ! m #    f_f_list  array
        integer,allocatable :: lda_ff_g_rx(:)             ! m #    f_f_list  array
        integer,allocatable :: mda_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: s11_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: s22_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: s33_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: s44_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: d11_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: d22_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: d12_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: ddd_ff_g_rx(:)             ! mode # f_f_list  array
        real,   allocatable :: e00_ff_g_rx(:)             ! mode # f_f_list  array

        integer             :: mov_num_loc    ! (priv) # of moves on this thread
        integer             :: mov_tot_loc    ! (priv) # of moving atms on this thread
        integer             :: mov_tot_loc3   ! (priv) the above * 3
        integer,allocatable :: mov_typ_loc(:) ! (priv) flex or rigd?
        integer,allocatable :: mov_beg_loc(:) ! (priv) first atom/mol
        integer,allocatable :: mov_end_loc(:) ! (priv) last atom/mol
        integer,allocatable :: mov_dyn_loc(:) ! (priv) BD or LD?
        integer,allocatable :: mov_hyd_loc(:) ! (priv) HI or noHI?
        integer,allocatable :: mov_lnc_loc(:) ! (priv) lincs or not?
        integer,allocatable :: mov_rep_loc(:) ! (priv) replica #
        integer,allocatable :: mov_beq_loc(:) ! (priv) first charge
        integer,allocatable :: mov_enq_loc(:) ! (priv) last charge

        integer      :: n_bond_spln   
        integer      :: n_angl_spln   
        integer      :: n_dihe_spln   
        integer      :: n_nbnd_spln   
        real         :: d_bond_spline
        real         :: d_angl_spline
        real         :: d_dihe_spline
        real         :: d_nbnd_spline
        real         :: h_bond_spline
        real         :: h_angl_spline
        real         :: h_dihe_spline
        real         :: h_nbnd_spline
        integer      :: num_bond_funcs   
        integer      :: num_angl_funcs   
        integer      :: num_dihe_funcs   
        integer      :: num_nbnd_funcs   
        integer      :: num_bnd_thrd     
        integer      :: num_ang_thrd     
        integer      :: num_dih_thrd     
        integer      :: num_rgd_dom_typs
        integer      :: num_ang_pairs 
        integer      :: movieframenumber 
        integer      :: NAM_runtype
        integer      :: NAM_num_q_as
        integer      :: NAM_num_runs
        integer      :: NAM_num_done
        integer      :: NAM_num_movie
        integer      :: NAM_terminate

C i_do_protease

        integer      :: num_tot_protease_actvs ! # of active sites
        integer      :: num_tot_protease_conts ! # act sites * subs
        integer      :: num_tot_protease_sites
        integer      :: num_tot_protease_bonds ! # tot 'bnds'

C B22

        integer      :: B22_bin_max ! # distance bins in the B22     
        integer      :: B22_bin     ! # distance bins in the B22     
        integer      :: B22_sample1 ! # conformers of 1st mol to try 
        integer      :: B22_sample2 ! # conformers of 2nd mol to try
C                                   !   for each of 1st mol
        integer      :: B22_sample3 ! # diagonals of 2nd mol to try
C                                   !   for each of 2nd mol
C                                   !   for each of 1st mol
        integer      :: num_B22_cnt ! counter for how often to write out the results
        integer      :: num_B22_stp ! how often to write out the results
        integer      :: B22_done1   ! # conformers done of molecule 1 
        integer      :: B22_done2   ! # conformers done of molecule 2 
        integer      :: B22_done3   ! # diagonal paths done for molecule 2 
        integer      :: B22_conf1   ! current conformer of molecule 1 
        integer      :: B22_conf2   ! current conformer of molecule 2 
        integer      :: B22_confs_1 ! # conformers of molecule 1 
        integer      :: B22_confs_2 ! # conformers of molecule 2 
        integer      :: umb_mol1    ! id of first mol in umbrella
        integer      :: umb_mol2    ! id of first mol in umbrella
        integer      :: umb_den1    ! id of first mol in umbrella
        integer      :: umb_den2    ! id of first mol in umbrella
        integer      :: num_exc_stp ! #steps b/w updating exclusive go lists
        integer      :: num_uniqs        ! #uniq types of atom
        integer      :: prinfo           ! treecode             
        integer      :: i_will_skip      ! intra_hydrodynamics
        integer      :: num_diff_rounds  ! intra_hydrodynamics
        integer      :: num_chol_rounds  ! intra_hydrodynamics
        integer      :: num_move_rounds  ! intra_hydrodynamics
        integer      :: num_move_remainder  ! intra_hydrodynamics
        integer      :: num_grow_mols    !
        integer      :: num_slide_conts  !
        integer      :: num_rdfs         ! rdfs       
        integer      :: num_rdf_bins     ! rdfs       
        integer      :: kmax             ! ewald_elec 
        integer      :: ksqmax           ! ewald_elec 
        integer      :: ksq              ! ewald_elec 
        integer      :: kx               ! ewald_elec 
        integer      :: ky               ! ewald_elec 
        integer      :: kz               ! ewald_elec 
        integer      :: totk_elec        ! ewald_elec 
        integer      :: totk_hydro       ! ewald_hydro
        integer      :: num_tot_go_pairs
        integer      :: num_cen1         ! 1st  atm for c.o.m.record  
        integer      :: num_cen2         ! last atm for c.o.m.record  
        integer      :: num_cens         ! sum of atm for c.o.m.record  
        integer      :: ivec1a           ! 1st atm of rot vector 1  
        integer      :: ivec1b           ! 1st atm of rot vector 1  
        integer      :: ivec2a           ! 1st atm of rot vector 1  
        integer      :: ivec2b           ! 1st atm of rot vector 1  
        integer      :: ivec3a           ! 1st atm of rot vector 1  
        integer      :: ivec3b           ! 1st atm of rot vector 1  
        integer      :: im_past_first    ! 1st mol used to calc Q
        integer      :: mol_Q1           ! 1st mol used to calc Q
        integer      :: mol_Q2           ! 2nd mol used to calc Q
        integer*8    :: ru               ! 2019 make *8 does this help?
        integer      :: int_scl                               
        integer      :: int_scl2                              
        integer*8    :: num_limit1                            
        integer      :: num_pos_res_f                         
        integer      :: num_pos_res_r                         
        integer      :: num_threads                    
        integer      :: num_reps  ! replica                   
        integer      :: myrep     ! replica         
        logical      :: myrep_master ! master thread for this replica         
        integer      :: num_f_f_v              
        integer      :: num_f_f_g              
        integer      :: num_uniq_m_pairs               ! number of total r as in the system                                 
        integer      :: max_r_as                    ! number of total r as in the system                                 
        integer      :: num_tot_atms                    ! total number of as in the system
        integer      :: num_mov_atms     ! 2019 total number of moving atoms in the system
        integer      :: num_flx_atms                    ! number of flexible atoms in the system                                 
        integer,allocatable      :: num_hyd_atms(:)  ! number of hydrodynamic atoms in the system                                 
        integer,allocatable      :: num_hyd_atms3(:) ! 3*number of hydrodynamic atoms in the system                                 
        integer      :: num_hyd_syss   ! number of hydrodynamic systems modeled                                  
        integer      :: max_hyd_atms   ! max number of hydrodynamic atoms in any system                                 
        integer      :: max_hyd_atms3  ! 3*max number of hydrodynamic atoms in any system                                 
        integer      :: num_flx_atms3                   ! 3*number of f as in the system                                 
        integer      :: num_q_as                    ! number of charged atoms in the system                                 
        integer      :: num_r_as                    ! number of rigid atoms in the system                                 
        integer      :: num_t_bs                    ! total number of blobs                
        integer      :: num_t_bs3                   ! 3*total number of blobs                
        integer      :: num_g_ms                    ! total number of molecules in go contacts 
        integer      :: num_t_ms                    ! total number of molecules                
        integer      :: num_t_ms3                   ! 3*total number of molecules                
        integer      :: num_f_ms                    ! number of f ms
        integer      :: num_r_ms                    ! number of r ms                   
        integer      :: num_t_atm_typs                ! number of total atom types in the system                                 
        integer      :: num_sys_bonds              ! number of total f bonds in the system                                 
        integer      :: num_sys_angls             ! number of total f angles in the system                                 
        integer      :: num_sys_dihes          ! number of total f dihedrals in the system                                 
        integer      :: num_tot_bonds              ! number of total f bonds in the system                                 
        integer      :: num_tot_angls             ! number of total f angles in the system                                 
        integer      :: num_tot_dihes          ! number of total f dihedrals in the system                                 
        integer      :: num_must_generate                ! number of r as that need regenerating at each timestep
        integer      :: num_m_pairs                    ! number of pairs of r molecules that interact
        integer      :: num_ene_scrn_cnt               ! total number of rigid molecules that grow during the simulations
        integer      :: num_tot_pka_sts                  ! total number of ionizable sites in the system
        integer      :: num_f_ts                    ! total number of fible molecule typs
        integer      :: num_t_ts                    ! total number of molecule typs
        integer      :: num_grow_r                    ! total number of rigid molecules that grow during the simulations
        integer      :: i_pbc                            ! =1 if we're using periodic boundary conditions
        integer      :: num_lst_stp                 ! number of steps required before update of nonbonded list
        integer      :: num_fmd_stp                 ! number of steps required before update of medium range forces
        integer      :: num_hyd_stp                 ! number of steps required before update of hydrodynamics
        integer      :: num_bal_stp                 ! number of steps required before update of hydrodynamics
        integer      :: num_rex_stp                 ! number of steps required before update of hydrodynamics
        integer      :: num_lst_cnt                 ! current count of steps to update of nonbonded list
        integer      :: num_fmd_cnt                 ! current count of steps to update of medium range forces
        integer      :: num_hyd_cnt                 ! current count of steps to update of hydrodynamics list
        integer      :: num_bal_cnt                 ! current count of steps to update of hydrodynamics list
        integer      :: num_rex_cnt                 ! current count of steps to update of hydrodynamics list
        integer      :: i_dum                            ! used to initialize random number generator
        integer      :: i_use_v_typ                   ! =1 if we're using 12-10 potentials,otherwise we're using 12-6 potentials
        integer      :: num_harmonics                    ! total number of harmonic interaction sites betweeen molecules 
        integer      :: num_pos_restraints               ! total number of molecules involved in restraint potentials              
        integer*8    :: nstep                            ! total number of steps taken in the simulation
        integer      :: num_x_ff             ! number of cells used for storing fible a positions in x-direction
        integer      :: num_y_ff             ! number of cells used for storing fible a positions in y-direction
        integer      :: num_z_ff             ! number of cells used for storing fible a positions in z-direction
        integer      :: num_rxn_sts                      ! total number of reactions to monitor (gs with i_do_reactions)
        integer      :: num_forces_omitted               ! total number of interactions that will be skipped in this way
        integer      :: num_clashes_allowed              ! not currently implemented
        integer      :: n_size      ! num of steps over which to shrink                        
        real         :: nbnd_func_scale ! scale all negative nbnd terms
        real         :: bal_load_fac ! exponent used to balance load in move routine for intra_hydrodynamics

        integer,allocatable :: itmp_mol_array(:) ! useful temporary array (1:num_t_ms)

        integer,allocatable :: imol_moved(:) ! does this mol move with intra HI (=1)
        integer,allocatable :: ida_rem(:) ! ids of atoms without HI  
        integer      :: num_best_match_kkk(-3:3)
        real         :: protease_cut  ! protease-substrate nbnd list cut 
        real         :: protease_cut2 
        real         :: non_protease_cut2 ! longest of the other cutoffs
        real         :: afm_tip_x_crd
        real         :: afm_tip_y_crd
        real         :: afm_tip_z_crd
        real         :: afm_tip_x_beg
        real         :: afm_tip_y_beg
        real         :: afm_tip_z_beg
        real         :: afm_tip_x_end
        real         :: afm_tip_y_end
        real         :: afm_tip_z_end
        real         :: afm_tip_x_dsp
        real         :: afm_tip_y_dsp
        real         :: afm_tip_z_dsp
        real         :: afm_tip_rad
        real         :: afm_tip_rad2
        real         :: afm_tip_frc
        real         :: NAM_radius 
        real         :: NAM_ave    
        real         :: NAM_min    
        real         :: NAM_max    
        real         :: NAM_energy
        real         :: NAM_elec1 
        real         :: NAM_elec2 
        real         :: NAM_c_d_x ! original c.o.m of non-diffusing mol
        real         :: NAM_c_d_y
        real         :: NAM_c_d_z 
        real         :: NAM_bsurf  
        real         :: NAM_qsurf  
        real         :: NAM_qsurf2 

        real         :: bond_scale
        real         :: sigma_scale
        real         :: go_scale

        integer      :: afm_tip_stp
        integer      :: afm_tip_cur
        integer      :: afm_tip_dir
C --------------------------------------------------------
        integer xd2,xd,ret,step,num_xtc
        real    prec
        character*15,allocatable :: xtc_file_out(:)     ! extra dim replica
        real,   allocatable :: x_xtc(:)
        real                :: box(9)
        integer,allocatable :: NAM_mol1(:)
        integer,allocatable :: NAM_mol2(:)
        real,   allocatable :: NAM_Qcur(:)
        real,   allocatable :: NAM_Qreq(:)
        integer,allocatable :: my_idm_MC(:) ! stores mol id of MC mover 
!                              allocated 1:num_t_ms_MC
        integer,allocatable :: my_MC_idm(:) ! stores index of MC mover for each mol
!                              allocated 1:num_t_ms
!                              =0 if mol is not a MC mover
        integer,allocatable :: my_MC_clust(:) ! stores cluster# of MC mover 
!                              allocated 1:num_t_ms_MC
        integer,allocatable :: my_MC_clust_rank(:) ! stores rank within cluster# of MC mover 
!                              allocated 1:num_t_ms_MC
        integer,allocatable :: iwant_to_add_MC(:)  ! flags to add to clstr list
        integer,allocatable :: my_MC_flg(:)  ! flags to use old or new coords
        integer,allocatable :: i_been_picked_MC(:)  ! flags whether picked
        integer,allocatable :: do_me_now_MC(:)  ! ordered list of mols to move  
        integer,allocatable :: idm_MC_pairs_tmp(:)  ! ordered list of mols to move  
        integer,allocatable :: jdm_MC_pairs_tmp(:)  ! ordered list of mols to move  
        integer,allocatable :: jdm_MC_pairs_new(:)  ! ordered list of mols to move  
        integer,allocatable :: idm_MC_pairs(:)  ! ordered list of mols to move  
        integer,allocatable :: jdm_MC_pairs(:)  ! ordered list of mols to move  
        real,   allocatable :: edm_MC_pairs_tmp(:)  ! ordered list of mols to move  
        real,   allocatable :: edm_MC_pairs(:)  ! ordered list of mols to move  
        real             :: MC1_energy  ! Etot of old position
        real             :: MC2_energy  ! Etot of new position
        real             :: MC_clstr_E  ! E crit for being in clstr
        real             :: MC_factr_E  ! E crit for being in clstr
        real             :: MC_timestep ! E crit for being in clstr
        real,allocatable :: MC_E_env(:,:) ! Eint with environment
        real,allocatable :: MC_E_mov(:,:,:,:) ! Eint with other MC movers
        real,allocatable :: dtrn_tmp(:) ! dtrans of mol type
        real,allocatable :: drot_tmp(:) ! drot   of mol type
        real,allocatable :: dtrn_MC(:)  ! dtrans of mol #
        real,allocatable :: drot_MC(:)  ! drot   of mol #
        real,allocatable :: r_d_trn_MC(:)  ! sqrt(2*dtrans*dt)
        real,allocatable :: r_d_rot_MC(:)  ! sqrt(2*drot*dt)   
        real,allocatable :: MC1_E_mol(:)  ! array of Etot of old pos 
        real,allocatable :: MC2_E_mol(:)  ! array of Etot of new pos
        real,allocatable :: MCX_E_mol(:)  ! array of Etot of accepted

C new for YHL code

        integer, allocatable :: ncc(:)    ! lincs
        real,    allocatable :: sdia(:)   ! lincs
        integer, allocatable :: icon(:,:) ! lincs
        real,    allocatable :: coef(:,:) ! lincs
        real,    allocatable :: alnc(:,:) ! lincs
        real,    allocatable :: blnc(:,:) ! lincs
        real,    allocatable :: rhsc(:,:) ! lincs
        real,    allocatable :: solc(:)   ! lincs

        integer*8            :: ftm1 ! timings of sections 2019 make *8
        integer*8            :: ftm2 ! timings of sections
        real,    allocatable :: tmtrx(:,:) ! timings of sections
        integer, allocatable :: thrd_of_atm(:,:) ! stores thread ID of each atom
C                               ! 2nd dim is replica #
        integer, allocatable :: imin_thrd_int(:,:)
        integer, allocatable :: imax_thrd_int(:,:)

C new rigid body code 

        integer, allocatable :: trd_rgd_flx_mov_tot(:)
        integer, allocatable :: typ_rgd_flx_mov_tot(:)
        integer, allocatable :: beg_rgd_flx_mov_tot(:)
        integer, allocatable :: end_rgd_flx_mov_tot(:)
        integer, allocatable :: typ_rgd_flx_mov(:)
        integer, allocatable :: beg_rgd_flx_mov(:)
        integer, allocatable :: end_rgd_flx_mov(:)
        real,    allocatable :: M11(:)
        real,    allocatable :: M12(:)
        real,    allocatable :: M13(:)
        real,    allocatable :: M21(:)
        real,    allocatable :: M22(:)
        real,    allocatable :: M23(:)
        real,    allocatable :: M31(:)
        real,    allocatable :: M32(:)
        real,    allocatable :: M33(:)
        real                 :: N11
        real                 :: N12
        real                 :: N13
        real                 :: N21
        real                 :: N22
        real                 :: N23
        real                 :: N31
        real                 :: N32
        real                 :: N33
        real,    allocatable :: fx1(:)
        real,    allocatable :: fy1(:)
        real,    allocatable :: fz1(:)
        real,    allocatable :: fx2(:)
        real,    allocatable :: fy2(:)
        real,    allocatable :: fz2(:)
        real,    allocatable :: fx3(:)
        real,    allocatable :: fy3(:)
        real,    allocatable :: fz3(:)
        real,    allocatable :: cx1(:)
        real,    allocatable :: cy1(:)
        real,    allocatable :: cz1(:)
        real,    allocatable :: cx2(:)
        real,    allocatable :: cy2(:)
        real,    allocatable :: cz2(:)
        real,    allocatable :: cx3(:)
        real,    allocatable :: cy3(:)
        real,    allocatable :: cz3(:)
        real,    allocatable :: D_11(:)
        real,    allocatable :: D_22(:)
        real,    allocatable :: D_33(:)
        real,    allocatable :: D_44(:)
        real,    allocatable :: D_55(:)
        real,    allocatable :: D_66(:)
        real,    allocatable :: E_11(:) ! sqrt(2*D_11*dt)
        real,    allocatable :: E_22(:)
        real,    allocatable :: E_33(:)
        real,    allocatable :: E_44(:)
        real,    allocatable :: E_55(:)
        real,    allocatable :: E_66(:)

C HANS
C here are the declarations relating to your treecode
C note that the treecode is invoked with the logical 'treecode_elec'  

        integer      :: treecode_order   ! order of Taylor expansion                                 
        integer      :: treecode_shrink  ! still not sure about this set to 1                        
        integer      :: treecode_maxatm  ! max # atoms in each node                                  
        real         :: treecode_theta   ! theta term for MAC in treecode                            
        integer,allocatable :: perm(:)      ! x coords of local atoms                                   
        integer,allocatable :: perm_cpy(:)      ! x coords of local atoms                                   
        integer,allocatable :: i_got_sorted(:)  ! used to random order charges                              
        integer,allocatable :: i_got_sorted2(:) ! used to random order pairs                                

C FYI here are some variables relating to a very simple attempt at a
C hydrodynamic treecode that I adapted from Josh Barnes' original code
C it's invoked with the logical 'tree_hydrodynamics'

        integer,allocatable :: idnode(:)
        integer,allocatable :: icnode(:)
        real,allocatable :: RANX(:)
        real,allocatable :: RANY(:)
        real,allocatable :: RANZ(:)
     
        real,allocatable :: c_q_x_tar(:)      ! x coords of local atoms                                   
        real,allocatable :: c_q_y_tar(:)      ! x coords of local atoms                                   
        real,allocatable :: c_q_z_tar(:)      ! x coords of local atoms                                   
        real,allocatable :: c_q_x_cpy(:)      ! x coords of local atoms                                   
        real,allocatable :: c_q_y_cpy(:)      ! x coords of local atoms                                   
        real,allocatable :: c_q_z_cpy(:)      ! x coords of local atoms                                   
        real,allocatable :: crg_of_atm_cpy(:)      ! x coords of local atoms                                   
        real,allocatable :: spengtar(:)  ! x coords of local atoms                                   
        real,allocatable :: tpengtar(:)  ! x coords of local atoms                                   
        real,allocatable :: upengtar(:)  ! x coords of local atoms                                   
        real,allocatable :: vpengtar(:)  ! x coords of local atoms                                   
        real,allocatable :: tforce(:,:)  ! x coords of local atoms                                   
        real,allocatable :: uforce(:,:)  ! x coords of local atoms                                   
        real,allocatable :: vforce(:,:)  ! x coords of local atoms                                   
C --------------------------------------------------------
        integer      :: ewald_elec_grid_data      ! total grid points in ewald elec grid                            
        integer      :: ewald_elec_grid_nx        ! grid points in x in ewald elec grid                            
        real         :: ewald_elec_grid_rx        ! grid spacing in ewald elec grid                            
        real         :: ewald_elec_grid_rx_inv    ! inv grid spacing in ewald elec grid                            
        real         :: ewald_elec_ei_with_iprime ! recip inter of i with iprime                               
        real,allocatable :: ewald_elec_grid_ff(:)    ! ewald elec grid                            
C --------------------------------------------------------
        integer      :: ewald_hydro_grid_nx        ! grid points in ewald hydro grid                            
        real         :: ewald_hydro_grid_rx        ! grid spacing in ewald hydro grid                            
        real         :: ewald_hydro_grid_rx_inv    ! inv grid spacing in ewald hydro grid                            
        real,allocatable :: ewald_hydro_grid_ff(:)    ! ewald hydro grid                            
C --------------------------------------------------------
        integer      :: ewald_elec_kmax  ! max number of k-vectors for ewald                         
        real         :: ewald_elec_self  ! self energy                               
        real         :: ewald_elec_kappa ! convergence parameter for ewald           
        real         :: ewald_elec_kappa2 ! convergence parameter for ewald           
        integer      :: ewald_hydro_kmax  ! max number of k-vectors for ewald                         
        real         :: ewald_hydro_self  ! self energy                               
        real         :: ewald_hydro_kappa ! convergence parameter for ewald           
        real         :: ewald_hydro_kappa2 ! convergence parameter for ewald           
        integer      :: fixman_order
        integer      :: fixman_steps ! counter for tree_hydro
        real*4       :: fixman_tol          ! tolerance for Fixman method
        real*4       :: fixman_y2_loc                                    
        real*4       :: fixman_y2                                        
        real*4       :: fixman_ef ! error estimate from Jendrejack
        real*4       :: fixman_l2 ! error estimate from Ando/Chous/Skolnick
        real*4       :: xme                                              
        real*4       :: yme                                              
        real*4       :: lmax                                             
        real*4       :: lmin                                             
        real*4       :: lmax_read_in                                     
        real*4       :: lmin_read_in                                     
        real*4       :: da                                               
        real*4       :: db                                             
        real*4       :: bma                                            
        real*4       :: bpa                                            
        real*4       :: dwDdw                                          
        real*4       :: dwDdw_loc                                      
        real*4       :: help                                             
        real*4       :: helpz                                            
        character*1  :: uplo                                             
C --------------------------------------------------------
        real         :: six_pi_eta  
        real         :: kT_over_8pi_eta  
        real         :: kT_over_6pi_eta  
        real         :: rootpi                                                       
        real         :: rootpi_inv                                                   
        real         :: r_size_fac  ! factor by which r_size is first scaled     
        real         :: r_size      ! radius of cell                                           
        real         :: r_size2     ! radius of cell^2                                           
        real         :: s_size      ! radius of cell (current)
        real         :: s_size2     ! radius of cell^2 (current)                                           
        real         :: l_size      ! length of cell (=-1 for cylinder)                        
        real         :: l_size2     ! length of cell^2 (=-1 for cylinder)                        
        real         :: l_sizeh     ! length of cell^0.5 (=-1 for cylinder)                        
        real         :: m_size      ! length of cell (=-1 for cylinder)                        
        real         :: m_size2     ! length of cell^2 (=-1 for cylinder)                        
        real         :: m_sizeh     ! length of cell^0.5 (=-1 for cylinder)                        
        real         :: f_size      ! force constant for cell wall         
        real         :: surface_area ! surface area of cell
        real         :: g_epsilon                                 
        real         :: beta_ii                                
        real         :: beta_ij                                
        real         :: bond_dev_quit                           
        real         :: beta_read_in                            
        real         :: Q_des    ! desired Q value for terminating                             
        real         :: go_eps_low ! only use contacts with eps>this for Qcalc                           
        real         :: inv_scl                               
        real         :: inv_scl2                              
        real         :: d_hydro                          ! rtrans term of the f a (1->num_flx_atms)
        real         :: s_hydro                          ! strans term of the f a (1->num_flx_atms)
        real         :: f_hydro                          ! strans term of the f a (1->num_flx_atms)
        real         :: df_hydro                         ! strans term of the f a (1->num_flx_atms)
        real         :: r_f_st                  
        real         :: r_f_hf                  
        real         :: time_total                 
        real         :: time_o1                    
        real         :: time_o2                    
        real         :: time_o3                    
        real         :: time_o4                    
        real         :: time_t0                    
        real         :: time_t1                    
        real         :: time_t2                    
        real         :: time_t3                    
        real         :: time_t4                    
        real         :: time_check1                
        real         :: time_check2                
        real         :: time_check3                
        real         :: time_check4                
        real         :: time_flex_list             
        real         :: time_flex_assign           
        real         :: time_get_randoms           
        real         :: time_zero_arrays           
        real         :: time_temp_2                 
        real         :: time_temp_3                 
        real         :: time_temp_4                 
        real         :: time_temp_5                 
        real         :: time_st_force
        real         :: time_md_force
        real         :: time_lg_force
        real         :: time_f_bonds
        real         :: time_f_angles 
        real         :: time_f_dihedrals
        real         :: time_f_bondeds
        real         :: time_f_f_nonbondeds
        real         :: time_assign_fs                 
        real         :: time_assign_rs                 
        real         :: time_construct_dmatrix                
        real         :: time_construct_smatrix                
        real         :: time_construct_f_f_list
        real         :: ene_tmp               
        real         :: fxtmp                 
        real         :: fytmp                 
        real         :: fztmp                 
        real         :: txtmp                 
        real         :: tytmp                 
        real         :: tztmp                 
        real         :: dist_dev_max     
        real         :: biggest_ene_jump
        real         :: tot_energy              
        real         :: ene_system_tot    
        real         :: ene_system_new    
        real         :: ene_system_old   
        real         :: rad_fit       
        real         :: rad_fit2      
        real         :: totenergy
        real         :: ene_lower_cutoff
        real         :: xmin
        real         :: ymin
        real         :: zmin
        real         :: xmax
        real         :: ymax
        real         :: zmax
        real         :: xmax_final
        real         :: ymax_final
        real         :: zmax_final
        real         :: xmin_final
        real         :: ymin_final
        real         :: zmin_final
        real         :: xlen
        real         :: ylen
        real         :: zlen
        real         :: xinv
        real         :: yinv
        real         :: zinv
        real         :: cut_h        
        real         :: cut_h2       
        real         :: cut_v_st
        real         :: cut_v_md
        real         :: cut_g_st
        real         :: cut_g_md
        real         :: cut_e_st
        real         :: cut_e_md
        real         :: cut_e_lg
        real         :: cut_v_st2
        real         :: cut_v_md2
        real         :: cut_g_st2
        real         :: cut_g_md2
        real         :: cut_e_st2
        real         :: cut_e_md2
        real         :: cut_e_lg2
        real         :: beta               
        real         :: mass       
        real         :: viscosity 
        real         :: dielectric
        real         :: dielectric_inv
        real         :: r_move_1                     
        real         :: r_move_2                     
        real         :: temperature
        real         :: r_kbt 
        real         :: ionic_strength
        real         :: kappa
        real         :: r_ion 
        real         :: r_pH
        real         :: r_probe
        real         :: one_over_one_plus_kappa_rion
        real         :: one_third
        real         :: eight_over_three
        real         :: force_limit
        real         :: force_limit2
        integer*8    :: ntimetotal ! required steps                            
        integer      :: necount 
        integer      :: neprint 
        integer      :: ntcount 
        integer      :: ntprint 
        integer      :: nmcount 
        integer      :: nmprint 
        integer      :: npcount 
        integer      :: npprint 
        real*16      :: teprint 
        real*16      :: ttprint 
        real*16      :: tmprint 
        real*16      :: tot_time
        real*16      :: tot_time_orig
        real*16      :: time_step
        real         :: x_ff
        real         :: y_ff
        real         :: z_ff
        real         :: x_ff_div2
        real         :: y_ff_div2
        real         :: z_ff_div2
        character*80 :: parameter_file
        character*80 :: goparam_file
        character*80 :: go_nonexclusive_file
        character*80 :: pos_restraint_file
        character*80 :: linker_file
        character*80 :: reaction_file
        character*80 :: noforce_file
        character*80 :: nomove_file
        character*80 :: wall_file

C 1D arrays

        integer,allocatable :: lbnd_thrd(:) 
C allocated 1:num_bnd_thrd - stores whether lincs bond done by this thread 

        integer,allocatable :: ibnd_thrd(:) 
C allocated 1:num_bnd_thrd - stores id of the bonds done by this thread 

        integer,allocatable :: iang_thrd(:) 
C allocated 1:num_ang_thrd - stores id of the angles done by this thread 

        integer,allocatable :: idih_thrd(:) 
C allocated 1:num_dih_thrd - stores id of the dihedrals done by this thread 

        real,    allocatable :: dbnd_thrd(:) 
C allocated 1:num_bnd_thrd - stores denom od the bonds done by this thread 

        real,   allocatable :: dang_thrd(:) 
C allocated 1:num_ang_thrd - stores denom of the angles done by this thread 

        real,   allocatable :: ddih_thrd(:) 
C allocated 1:num_dih_thrd - stores denom of the dihedrals done by this thread 

C i_do_protease

        integer,allocatable :: name_override(:) 
C allocated 1:num_flx_atms - flags whether name should be changed in pdb

        integer,allocatable :: my_prtase_ident(:) 
C allocated 1:num_tot_protease_conts - stores # of the active site associated 
C with this combination of protease and substrate

        integer,allocatable :: my_prtase_1_atm(:) 
C allocated 1:num_tot_protease_actvs - stores 1st atom of this active site     

        integer,allocatable :: my_prtase_inact(:) 
C allocated 1:num_tot_protease_actvs - flags whether active site is inactivated

        integer,allocatable :: my_bnd_met_crit(:) 
C allocated 1:num_tot_protease_bonds - flags whether bond met cleavage criteria

        integer,allocatable :: my_bnd_in_range(:) 
C allocated 1:num_tot_protease_bonds - flags whether bond in range to feel force

        integer,allocatable :: my_bnd_has_clvd(:) 
C allocated 1:num_tot_protease_bonds - flags whether bond has cleaved 
C = 0 if uncleaved, = 1 if cleaved, = 2 if acyl-enzyme intermediate

        real,   allocatable :: my_bnds_acy_len(:) 
C allocated 1:num_tot_protease_bonds - stores length of acyl bond 
C note that this will be undefined for most bonds - only the 1st bond
C of each cleavable site will get a value

        real,   allocatable :: my_bnds_acy_fct(:) 
C allocated 1:num_tot_protease_bonds - stores force constant of acyl bond 
C note that this will be undefined for most bonds - only the 1st bond
C of each cleavable site will get a value

        integer,allocatable :: my_bnds_acy_stp(:) 
C allocated 1:num_tot_protease_bonds - stores reqrd steps b4 acyl bond breaks
C note that this will be undefined for most bonds - only the 1st bond
C of each cleavable site will get a value

        integer,allocatable :: my_bnds_acy_cur(:) 
C allocated 1:num_tot_protease_bonds - stores curnt steps b4 acyl bond breaks
C note that this will be undefined for most bonds - only the 1st bond
C of each cleavable site will get a value

        integer,allocatable :: my_ste_met_crit(:) 
C allocated 1:num_tot_protease_sites - flags whether site met cleavage criteria

        integer,allocatable :: my_ste_has_clvd(:) 
C allocated 1:num_tot_protease_sites - flags whether site has cleaved 

        integer,allocatable :: my_stes_prtase(:) 
C allocated 1:num_tot_protease_sites - stores # of active site attached
C to this particular site

        integer,allocatable :: my_stes_1st_bnd(:) 
C allocated 1:num_tot_protease_sites - stores # of 1st bnd in each site

        integer,allocatable :: my_stes_1st_atm(:) 
C allocated 1:num_tot_protease_sites - stores # of 1st atm of site's
C                                      actual cleavable bond

        integer,allocatable :: my_stes_2nd_atm(:) 
C allocated 1:num_tot_protease_sites - stores # of 1st atm of site's
C                                      actual cleavable bond

        integer,allocatable :: my_stes_min_atm(:) 
C allocated 1:num_tot_protease_sites - stores # of min atm of site's
C                                      actual cleavable bond

        integer,allocatable :: my_stes_max_atm(:) 
C allocated 1:num_tot_protease_sites - stores # of max atm of site's
C                                      actual cleavable bond

        integer,allocatable :: my_stes_num_bnd(:) 
C allocated 1:num_tot_protease_sites - stores # of bnds in each site

        real,   allocatable :: my_stes_clv_prb(:) 
C allocated 1:num_tot_protease_sites - probability of cleaving site    

        integer,allocatable :: my_stes_prot_id(:) 
C allocated 1:num_tot_protease_sites - stores mol # of protease for this site  

        integer,allocatable :: ida_ff_p_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_p_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jdb_ff_p_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jds_ff_p_st(:)             ! a pair f_f_list  array
        real,   allocatable :: crt_ff_p_st(:)             ! a pair f_f_list  array
        real,   allocatable :: rng_ff_p_st(:)             ! a pair f_f_list  array
        real,   allocatable :: fct_ff_p_st(:)             ! a pair f_f_list  array

        real,   allocatable :: r_d_trn(:) 
        real,   allocatable :: r_d_rot(:) 
C allocated 1:num_rgd_doms - stores the random conversion factor on each domain

        real,   allocatable :: s_d_trn(:) 
        real,   allocatable :: s_d_rot(:) 
C allocated 1:num_rgd_doms - stores the diffusion conversion factor on each domain

        integer,allocatable :: my_doms_mol(:) 
C allocated 1:num_rgd_doms - stores the molecule # associated with each domain

        integer,allocatable :: my_doms_upd(:) 
C allocated 1:num_rgd_doms - stores the update frequency associated with each domain

        integer,allocatable :: my_doms_dme(:) 
C allocated 1:num_rgd_doms - =1 if doing rigid motion of this each domain on this step
C                               the dme is meant to mean 'do me'

        integer,allocatable :: my_doms_cur(:) 
C allocated 1:num_rgd_doms - stores the counter for frequency associated with each domain

        integer,allocatable :: my_doms_rgd_dom_typ(:) 
C allocated 1:num_rgd_doms - stores the domain type associated with each domain

        integer,allocatable :: my_atms_rgd_dom(:) 
C allocated 1:num_flx_atms - stores the domain number associated with each atom

        integer,allocatable :: my_atms_rgd_dom_typ(:)
C allocated 1:num_flx_atms - stores the domain type associated with each atom

        integer,allocatable :: my_atms_rgd_dom_upd(:)
C allocated 1:num_flx_atms - stores the update frequency of flex moves of this atom

        integer,allocatable :: my_atms_rgd_dom_cur(:)
C allocated 1:num_flx_atms - stores the counter for frequency of flex moves of this atom
C                        when my_atms_rgd_dom_cur=my_atms_rgd_dom_upd time for flex move

        integer,allocatable :: ths_rgd_dom_fnd_in_ths_mol(:)
C temporary array used to help assign domain numbers to each atom

        integer :: num_i_res_no_nb_grps
        integer,allocatable :: i_res_no_nb(:)
C allocated 1:num_flx_atms - =1 if atom is excused from nonbonded interactions
C (note that only interacting pairs in which *both* atoms are excused are ignored)

        integer,allocatable :: iyz(:) ! for every atom in
        integer,allocatable :: my_uniq_atm_num(:) ! for every atom in
C                 the system, this stores its # in the uniq atom list
        real,   allocatable :: x_points(:)     
        real,   allocatable :: y_points(:)     
        real,   allocatable :: e_points(:)     
        real,   allocatable :: x_spline(:)     
        real,   allocatable :: y_spline(:)     
        real,   allocatable :: u_spline(:)     
        real,   allocatable :: y2_spline(:)     
        integer,allocatable :: crd_wall(:)     
        real,   allocatable :: pos_wall(:)     
        real,   allocatable :: fct_wall(:)     

        integer,allocatable :: num_in_MC_clust(:)

        integer,allocatable :: imol_in_cont(:)
        integer,allocatable :: jmol_in_cont(:)
        integer,allocatable :: num_in_clust(:)

        integer,allocatable :: i_restraint(:)
        integer,allocatable :: ibeg_tmp_b_typ(:)
        integer,allocatable :: iend_tmp_b_typ(:)
        integer,allocatable :: itp_tot_b_typ(:)
        integer,allocatable :: ida_tot_b_typ(:)
        integer,allocatable :: jda_tot_b_typ(:)
        integer,allocatable :: ind_tot_b_typ(:)
        integer,allocatable :: jnd_tot_b_typ(:)
        integer,allocatable :: itp_tot_go_pair(:)
        integer,allocatable :: ida_tot_go_pair(:)
        integer,allocatable :: jtp_tot_go_pair(:)
        integer,allocatable :: jda_tot_go_pair(:)
        integer,allocatable :: ind_tot_go_pair(:)
        integer,allocatable :: jnd_tot_go_pair(:)
        integer,allocatable :: ibeg_m_a(:)
        integer,allocatable :: iend_m_a(:)
        integer,allocatable :: ibeg_t_a(:)
        integer,allocatable :: iend_t_a(:)
        integer,allocatable :: ibeg_bms(:)
        integer,allocatable :: ibeg_ams(:)
        integer,allocatable :: ibeg_dms(:)
        integer,allocatable :: iend_bms(:)
        integer,allocatable :: iend_ams(:)
        integer,allocatable :: iend_dms(:)
        integer,allocatable :: num_h_d_ngb(:)                                 
        integer,allocatable :: num_loc_f_as(:)                              
        integer,allocatable :: num_loc_r_as(:)                              
        integer,allocatable :: num_loc_h_as(:)                              
        integer,allocatable :: mywdth(:)                              
        integer,allocatable :: mywdth1(:)                              
        integer,allocatable :: mywdth2(:)                              
        integer,allocatable :: mywdth3(:)                              
        integer,allocatable :: mywdth4(:)                              
        integer,allocatable :: mywdth5(:)                              
        integer,allocatable :: i1_loc(:)
        integer,allocatable :: i2_loc(:)
        integer,allocatable :: i3_loc(:)
        integer,allocatable :: j1_loc(:)
        integer,allocatable :: j2_loc(:)
        integer,allocatable :: j3_loc(:)
        integer,allocatable :: k1_loc(:)
        integer,allocatable :: k2_loc(:)
        integer,allocatable :: k3_loc(:)
        integer,allocatable :: l1_loc(:)
        integer,allocatable :: l2_loc(:)
        integer,allocatable :: l3_loc(:)
        integer,allocatable :: num_m_pairs_omp(:)             ! number of pairs of r molecules that interact
        integer,allocatable :: num_a_pairs_omp(:)             ! number of pairs of r molecules that interact
        integer,allocatable :: num_qqq_pairs_omp(:)             ! number of pairs of r molecules that interact
        integer,allocatable :: num_uniq_m_pairs_omp(:)        ! number of pairs of r molecules that interact
C --------------------------------------------------------
        character*3,allocatable :: atm_tmp(:)
        character*3,allocatable :: res_tmp(:)
        character*3,allocatable :: atm_nam_typ(:)
        character*3,allocatable :: res_nam_typ(:)
        integer,    allocatable :: atm_typ(:)
        integer,    allocatable :: atm_atm_typ(:,:,:)
C --------------------------------------------------------
        character*3,allocatable :: atm_nam_atm(:)
        character*3,allocatable :: res_nam_atm(:)
        integer,    allocatable :: res_num_atm(:)
C --------------------------------------------------------
      character*16,allocatable :: name1_nbnd(:) ! resname
      character*16,allocatable :: name2_nbnd(:) ! atmname
      character*16,allocatable :: name3_nbnd(:) ! resnum
      character*16,allocatable :: name1_bond(:) ! resname
      character*16,allocatable :: name2_bond(:) ! atmname
      character*16,allocatable :: name3_bond(:) ! resnum
      character*16,allocatable :: name1_angl(:) ! resname
      character*16,allocatable :: name2_angl(:) ! atmname
      character*16,allocatable :: name3_angl(:) ! resnum
      character*16,allocatable :: name1_dihe(:) ! resname
      character*16,allocatable :: name2_dihe(:) ! atmname
      character*16,allocatable :: name3_dihe(:) ! resnum

      integer, allocatable :: hist_2d_real(:,:,:) ! hist
      real,    allocatable :: rist_2d_real(:,:,:) ! prob
      real,    allocatable :: sist_2d_real(:,:,:) ! free energy

C user_hists arrays
 
      character*54         :: new_string
      character*23         :: userfile   
      character*25         :: userfile2
      character*34         :: userfile3
      real,    allocatable :: rjac(:,:)
      character*4          :: usercrap
      integer              :: user_bin(1:6)   ! private!
      integer              :: user_bin_edge(1:6,1:2)   ! private!
      integer              :: rigid_bin(1:6)   ! private!
      integer              :: rigid_bin_edge(1:6,1:2)   ! private!
      real                 :: valu_bin(1:6)   ! private!
      real                 :: frac_bin(1:6)   ! private!
      real                 :: frac_bin_edge(1:6,1:2)   ! private!
      real                 :: tmp_v(1:6)      ! private!
      real                 :: tmp_f(1:6)      ! private! force
      real                 :: tmp_e(1:6)      ! private! energy
      real                 :: tmp_u(1:6,1:100)! private!
      integer              :: i_met_condition ! private!

      integer, allocatable :: user_energy_dim_real(:)
      integer, allocatable :: user_energy_box_real(:,:)
      integer, allocatable :: user_energy_typ_real(:,:)
      integer, allocatable :: user_energy_nbn_real(:,:)
      integer, allocatable :: user_energy_per_real(:,:)
      integer, allocatable :: user_energy_atm_real(:,:,:)
      real,    allocatable :: user_energy_min_real(:,:)
      real,    allocatable :: user_energy_hlf_real(:,:)
      real,    allocatable :: user_energy_max_real(:,:)
      real,    allocatable :: user_energy_bin_real(:,:)

      integer, allocatable :: i_get_mapped_to_grid(:)
      integer, allocatable :: rigid_energy_dim_real(:)
      integer, allocatable :: rigid_energy_box_real(:,:)
      integer, allocatable :: rigid_energy_typ_real(:,:)
      integer, allocatable :: rigid_energy_nbn_real(:,:)
      integer, allocatable :: rigid_energy_per_real(:,:)
      integer, allocatable :: rigid_energy_atm_real(:,:,:)
      real,    allocatable :: rigid_energy_min_real(:,:)
      real,    allocatable :: rigid_energy_hlf_real(:,:)
      real,    allocatable :: rigid_energy_max_real(:,:)
      real,    allocatable :: rigid_energy_bin_real(:,:)

      integer, allocatable :: rigid_energy_typ(:,:)
      real,    allocatable :: rigid_energy_cut(:,:)

      integer, allocatable :: user_hist_dim_real(:)
      integer, allocatable :: user_hist_typ_real(:,:)
      integer, allocatable :: user_hist_nbn_real(:,:)
      integer, allocatable :: user_hist_atm_real(:,:,:)
      real,    allocatable :: user_hist_min_real(:,:)
      real,    allocatable :: user_hist_max_real(:,:)
      real,    allocatable :: user_hist_bin_real(:,:)
      integer, allocatable :: user_dim_cond(:) 
      integer, allocatable :: user_typ_cond(:,:)
      integer, allocatable :: user_atm_cond(:,:,:)
      real,    allocatable :: user_hist_min_cond(:,:)
      real,    allocatable :: user_max_cond(:,:)

      real,    allocatable :: curr_angl(:) ! for 2D hist
      real,    allocatable :: curr_dihe(:) ! for 2D hist
      real,    allocatable :: sum_in_ref(:)
      real,    allocatable :: sum_in_act(:)
      real,    allocatable :: scale_fac(:)
      integer, allocatable :: hist_nbnd_ref(:)
      integer, allocatable :: hist_nbnd(:,:)
      integer, allocatable :: hist_bond(:,:)
      integer, allocatable :: hist_angl(:,:)
      integer, allocatable :: hist_dihe(:,:)

      real,    allocatable :: rist_nbnd(:,:) ! probability
      real,    allocatable :: rist_bond(:,:) ! probability
      real,    allocatable :: rist_angl(:,:)
      real,    allocatable :: rist_dihe(:,:)

      real,    allocatable :: sist_nbnd(:,:) ! free energy
      real,    allocatable :: sist_bond(:,:) ! free energy
      real,    allocatable :: sist_angl(:,:)
      real,    allocatable :: sist_dihe(:,:)
C --------------------------------------------------------
        integer,    allocatable :: loc1_nbnd(:)
        integer,    allocatable :: loc2_nbnd(:)
        integer,    allocatable :: nmax_nbnd(:)
        character*3,allocatable :: res1_nbnd(:)
        character*3,allocatable :: res2_nbnd(:)
        character*3,allocatable :: atm1_nbnd(:)
        character*3,allocatable :: atm2_nbnd(:)
        real,       allocatable :: dist_nbnd(:,:) ! read-in
        real,       allocatable :: ener_nbnd(:,:) ! read-in
        real,       allocatable :: forc_nbnd(:,:)
        real,       allocatable :: y2_nbnd(:)
        real,       allocatable :: u_nbnd(:)
        real,       allocatable :: d_nbnd(:,:)    ! splined
        real,       allocatable :: e_nbnd(:,:)    ! splined
        real,       allocatable :: f_nbnd(:,:)    ! splined
C --------------------------------------------------------
        integer,    allocatable :: loc1_bond(:)
        integer,    allocatable :: loc2_bond(:)
        integer,    allocatable :: nmax_bond(:)
        character*3,allocatable :: res1_bond(:)
        character*3,allocatable :: res2_bond(:)
        character*3,allocatable :: atm1_bond(:)
        character*3,allocatable :: atm2_bond(:)
        real,       allocatable :: dist_bond(:,:) ! read-in
        real,       allocatable :: ener_bond(:,:) ! read-in
        real,       allocatable :: forc_bond(:,:)
        real,       allocatable :: y2_bond(:)
        real,       allocatable :: u_bond(:)
        real,       allocatable :: d_bond(:,:)    ! splined
        real,       allocatable :: e_bond(:,:)    ! splined
        real,       allocatable :: f_bond(:,:)    ! splined
C --------------------------------------------------------
        integer,    allocatable :: loc1_angl(:)
        integer,    allocatable :: loc2_angl(:)
        integer,    allocatable :: loc3_angl(:)
        integer,    allocatable :: nmax_angl(:)
        character*3,allocatable :: res1_angl(:)
        character*3,allocatable :: res2_angl(:)
        character*3,allocatable :: res3_angl(:)
        character*3,allocatable :: atm1_angl(:)
        character*3,allocatable :: atm2_angl(:)
        character*3,allocatable :: atm3_angl(:)
        real,       allocatable :: dist_angl(:,:)
        real,       allocatable :: ener_angl(:,:)
        real,       allocatable :: forc_angl(:,:)
        real,       allocatable :: y2_angl(:)
        real,       allocatable :: u_angl(:)
        real,       allocatable :: d_angl(:,:)
        real,       allocatable :: e_angl(:,:)
        real,       allocatable :: f_angl(:,:)
C --------------------------------------------------------
        integer,    allocatable :: loc1_dihe(:)
        integer,    allocatable :: loc2_dihe(:)
        integer,    allocatable :: loc3_dihe(:)
        integer,    allocatable :: loc4_dihe(:)
        integer,    allocatable :: nmax_dihe(:)
        character*3,allocatable :: res1_dihe(:)
        character*3,allocatable :: res2_dihe(:)
        character*3,allocatable :: res3_dihe(:)
        character*3,allocatable :: res4_dihe(:)
        character*3,allocatable :: atm1_dihe(:)
        character*3,allocatable :: atm2_dihe(:)
        character*3,allocatable :: atm3_dihe(:)
        character*3,allocatable :: atm4_dihe(:)
        real,       allocatable :: dist_dihe(:,:)
        real,       allocatable :: ener_dihe(:,:)
        real,       allocatable :: forc_dihe(:,:)
        real,       allocatable :: y2_dihe(:)
        real,       allocatable :: u_dihe(:)
        real,       allocatable :: d_dihe(:,:)
        real,       allocatable :: e_dihe(:,:)
        real,       allocatable :: f_dihe(:,:)
C --------------------------------------------------------
        real,allocatable    :: c_o_m_x(:)     ! center of mass x
        real,allocatable    :: c_o_m_y(:)     ! x coords of center of
        real,allocatable    :: c_o_m_z(:)     ! x coords of center of
        real,   allocatable :: frecipx(:)            ! ewald_elec 
        real,   allocatable :: frecipy(:)            ! ewald_elec 
        real,   allocatable :: frecipz(:)            ! ewald_elec 
        real,   allocatable :: rkx_elec(:)           ! ewald_elec 
        real,   allocatable :: rky_elec(:)           ! ewald_elec 
        real,   allocatable :: rkz_elec(:)           ! ewald_elec 
        real,   allocatable :: akvec_elec(:)         ! ewald_elec 
        real,   allocatable :: cossum_elec(:)        ! ewald_elec
        real,   allocatable :: sinsum_elec(:)        ! ewald_elec
        real,   allocatable :: costerm_elec(:,:)     ! ewald_elec 
        real,   allocatable :: sinterm_elec(:,:)     ! ewald_elec 
        real,   allocatable :: cossum_elec_tmp(:,:)  ! ewald_elec
        real,   allocatable :: sinsum_elec_tmp(:,:)  ! ewald_elec
        real,   allocatable :: rkx_hydro(:)          ! ewald_hydro 
        real,   allocatable :: rky_hydro(:)          ! ewald_hydro 
        real,   allocatable :: rkz_hydro(:)          ! ewald_hydro 
        real,   allocatable :: ekx_hydro(:)          ! ewald_hydro 
        real,   allocatable :: eky_hydro(:)          ! ewald_hydro 
        real,   allocatable :: ekz_hydro(:)          ! ewald_hydro 
        real,   allocatable :: tm1_hydro(:)          ! ewald_hydro 
        real,   allocatable :: tm2_hydro(:)          ! ewald_hydro 
        real,   allocatable :: akvec_hydro(:)        ! ewald_hydro 
        real*4, allocatable :: ame(:,:)
        real*4, allocatable :: dme(:)
        real*4, allocatable :: eme(:)
        real*4, allocatable :: zme(:)
        real*4              :: tmp_trm(1:3)        ! tolerance for Fixman method
        real*4, allocatable :: fixman_coeff(:)     ! tolerance for Fixman method
        real*4, allocatable :: fixman_ran(:)        ! tolerance for Fixman method
        real*4, allocatable :: fixman_ran_tree(:,:)      ! tolerance for Fixman method
        real*4, allocatable :: fixman_tensor(:,:)  ! tolerance for Fixman method
        real*4, allocatable :: fixman_a_tensor(:,:,:)  ! tolerance for Fixman method
        real*4, allocatable :: fixman_m_tensor(:,:)  ! tolerance for Fixman method
        real*4, allocatable :: fixman_xx(:,:)      ! tolerance for Fixman method
        real*4, allocatable :: fixman_xx_loc(:,:)      ! tolerance for Fixman method
C following not used with multi_hydrodynamics
        real,   allocatable :: amkl(:)   ! mkl routine for eigenvalues
        real,   allocatable :: dmkl(:)   ! mkl routine for eigenvalues
        real,   allocatable :: emkl(:)   ! mkl routine for eigenvalues
        real,   allocatable :: tmkl(:)   ! mkl routine for eigenvalues
        real,   allocatable :: wmkl(:)   ! mkl routine for eigenvalues
        real,   allocatable :: work(:)   ! mkl routine for eigenvalues
        integer,allocatable :: iblock(:) ! mkl routine for eigenvalues
        integer,allocatable :: isplit(:) ! mkl routine for eigenvalues
        integer,allocatable :: iwork(:)  ! mkl routine for eigenvalues
C --------------------------------------------------------
        real,   allocatable :: cx_restraint(:)
        real,   allocatable :: cy_restraint(:)
        real,   allocatable :: cz_restraint(:)
        real,   allocatable :: dx_restraint(:)
        real,   allocatable :: dy_restraint(:)
        real,   allocatable :: dz_restraint(:)
        real,   allocatable :: x_t_restraint(:) ! trial positions
        real,   allocatable :: y_t_restraint(:)
        real,   allocatable :: z_t_restraint(:)
        real,   allocatable :: r_restraint(:)
        real,   allocatable :: fx_restraint(:)
        real,   allocatable :: fy_restraint(:)
        real,   allocatable :: fz_restraint(:)
        real,   allocatable :: fr_restraint(:)
C following all used for langevin approach
        real,   allocatable :: v_e_t(:) ! MC velocity external forces
        real,   allocatable :: v_r_t(:) ! MC velocity random   forces
        real,   allocatable :: ld_fac1(:)
        real,   allocatable :: ld_fac2(:)
        real,   allocatable :: ld_fac3(:)
        real,   allocatable :: ld_fac4(:)
C end of langevin
c       real,   allocatable :: s(:) 
        real,   allocatable :: go_t11(:)
        real,   allocatable :: go_t22(:)
        real,   allocatable :: go_t33(:)
        real,   allocatable :: go_t44(:)
        real,   allocatable :: go_dis(:)
        real,   allocatable :: go_di2(:)
        real,   allocatable :: go_eps(:)
        integer,allocatable :: go_mod(:) ! mode of binding for exclusive go stuff
        integer,allocatable :: imol_to_sort(:) ! use with exclusive go
        integer,allocatable :: jmol_to_sort(:)
        integer,allocatable :: kmod_to_sort(:)
        integer,allocatable :: lcnt_to_sort(:)
        integer,allocatable :: stil_in_list(:)
        integer             :: go_mod_max !max # of modes for any mol pair interaction
        real,   allocatable :: go_elo(:)
        real,   allocatable :: go_d12(:)
        real,   allocatable :: g_epsilon_tmp(:)
        real,   allocatable :: sm1_x(:) 
        real,   allocatable :: sm1_y(:) 
        real,   allocatable :: sm1_z(:) 
        real,   allocatable :: ran_mkl(:) 
        real,   allocatable :: xloc(:) 
        real,   allocatable :: tloc(:,:) 
        real,   allocatable :: dxp(:,:) 
        real,   allocatable :: dxm(:,:) 
        real,   allocatable :: xf(:,:) 
        real,   allocatable :: xm(:,:) 
        real,   allocatable :: dist2_moved(:)           ! max squared distance moved by a r m on each of the MPI processors 
        real,   allocatable :: rkbb_new(:)              ! 10-08-08                                             
        real,   allocatable :: r0bb_new(:)              ! equilibrium bond lengths (1-> num_tot_bonds)
        real,   allocatable :: rkba_new(:)
        real,   allocatable :: r0ba_new(:)
        real,   allocatable :: r0bd_new(:)
        real,   allocatable :: rkbb(:)                  ! bond force constants (1->num_tot_bonds)
        real,   allocatable :: r0bb(:)                  ! equilibrium bond lengths (1-> num_tot_bonds)
        real,   allocatable :: rkba(:)
        real,   allocatable :: r0ba(:)
        real,   allocatable :: r0bd(:)
        real,   allocatable :: scbeg_grow_r(:)
        real,   allocatable :: scend_grow_r(:)
        real,   allocatable :: ene_pos(:)  
        real,   allocatable :: ene_bond(:)  
        real,   allocatable :: ene_angl(:)  
        real,   allocatable :: ene_dihe(:)  
        real,   allocatable :: ene_m(:)  
        real,   allocatable :: ene_m_old(:)  
        real,   allocatable :: ene_m_new(:)  
        real,   allocatable :: rxn_dis_1(:)  
        real,   allocatable :: rxn_dis_2(:)  
        real,   allocatable :: r_t_m(:)
        real,   allocatable :: d_t_m(:)
        real,   allocatable :: s_t_m(:)
        real,   allocatable :: r_a_r(:)
        real,   allocatable :: random_trms(:) ! new random
        real,   allocatable :: r_a_r_new(:,:)
        real,   allocatable :: r_a_r_tree(:,:)
        real,   allocatable :: r_a_r_mtrx(:,:)
        real,   allocatable :: fixman_ran_mtrx(:,:) 
        real,   allocatable :: fixman_xx_mtrx(:,:,:) 
        real,   allocatable :: fixman_l2_num_mtrx(:,:) 
        real,   allocatable :: fixman_l2_den_mtrx(:,:) 
        real,   allocatable :: r_a_x(:)
        real,   allocatable :: r_a_y(:)
        real,   allocatable :: r_a_z(:)
        real,   allocatable :: div_x(:) !divergence contribution to move
        real,   allocatable :: div_y(:)
        real,   allocatable :: div_z(:)
        real,   allocatable :: pka0(:)    
        real,   allocatable :: actualq(:)                              
        real,   allocatable :: time_h_gen(:)  ! time to make hyd mat
        real,   allocatable :: time_h_mov(:)  ! time to make hyd mov 
        real,   allocatable :: time_t_acc(:,:)  ! time accumulators    
        real,   allocatable :: tmp_local(:)                              
        real,   allocatable :: c_n_x(:) ! box-normalized for ewald_elec                             
        real,   allocatable :: c_n_y(:)                              
        real,   allocatable :: c_n_z(:)                              

C these have had extra dim added for num_reps

        real,   allocatable :: c_a_x(:,:,:)                             
        real,   allocatable :: c_a_x_fail(:,:,:)
        real,   allocatable :: c_a_x_old(:,:,:) 
        real,   allocatable :: t_a_x(:,:,:)                             
        real,   allocatable :: c_c_x(:,:,:,:) ! previous frames for crashes                              
        real,   allocatable :: f_a_x(:,:,:)                             
        real,   allocatable :: f_u_x(:,:,:)                              
        real,   allocatable :: f_a_r(:,:) 
        real,   allocatable :: r_c_r(:,:,:) ! new random
        real,   allocatable :: r_u_r(:,:,:) ! new random
        real,   allocatable :: h_a_r(:,:,:) 
        real,   allocatable :: h_e_r(:,:,:) ! f-tilda  external forces
        real,   allocatable :: h_c_r(:,:,:,:) ! new random
        real,   allocatable :: h_u_r(:,:,:,:) ! new random
        real,   allocatable :: d_a_r(:,:,:,:)  
        real,   allocatable :: s_a_r(:,:,:,:)
        real,   allocatable :: v_e_r(:,:) ! velocity external forces
        real,   allocatable :: v_r_r(:,:) ! velocity random   forces
        real,   allocatable :: f_e_r(:,:) ! f-tilda  external forces
        real,   allocatable :: f_r_r(:,:) ! f-tilda  random   forces
        real,   allocatable :: e_a_x_st(:,:,:)                              
        real,   allocatable :: e_a_x_md(:,:,:)                              
        real,   allocatable :: e_a_x_lg(:,:,:)                              
        real,   allocatable :: e_e_x_st(:,:,:) ! ewald-elec and treecode    
        real,   allocatable :: e_e_x_md(:,:,:) ! treecode                   
        real,   allocatable :: g_a_x(:,:,:,:) ! new for YHL code           

        real,   allocatable :: NAM_c_a_x(:,:)                              
        real,   allocatable :: c_t_x(:,:) ! trial locations of atoms for MC
        real,   allocatable :: c_t_y(:,:)                              
        real,   allocatable :: c_t_z(:,:)                              
        real,   allocatable :: cx_tmp(:) ! used to get cluster center for MC
        real,   allocatable :: cy_tmp(:) ! used to get cluster center for MC
        real,   allocatable :: cz_tmp(:) ! used to get cluster center for MC
        real,   allocatable :: dx_tmp(:) ! used to get cluster center for MC
        real,   allocatable :: dy_tmp(:) ! used to get cluster center for MC
        real,   allocatable :: dz_tmp(:) ! used to get cluster center for MC

C stuff for blob models ---------------------------------------------
        real,   allocatable :: c_b_x(:) ! coords for blobs                           
        real,   allocatable :: c_b_y(:)                              
        real,   allocatable :: c_b_z(:)                              
        real,   allocatable :: r_t_b(:)                              
        real,   allocatable :: d_t_b(:)                              
        real,   allocatable :: s_t_b(:)                              
        real,   allocatable :: r_b_r(:)                              
        real,   allocatable :: s_b_r(:)                              
        real,   allocatable :: f_b_r(:)                              
        real,   allocatable :: d_b_m(:,:)                              
        real,   allocatable :: s_b_m(:,:)                              
        integer,allocatable :: ic1map(:)                              
        integer,allocatable :: ic2map(:)                              
        integer,allocatable :: ic3map(:)                              
        integer,allocatable :: nmems(:)                              
        integer,allocatable :: id_blob(:)                            
        integer,allocatable :: id_beg_blobs(:)
        integer,allocatable :: id_end_blobs(:)
C -------------------------------------------------------------------  
        real,   allocatable :: c_r_x(:)                              
        real,   allocatable :: c_r_y(:)                              
        real,   allocatable :: c_r_z(:)                              
        real,   allocatable :: c_m_x(:)                 ! dim 1->num_t_ms
        real,   allocatable :: c_m_y(:)                 ! dim 1->num_t_ms
        real,   allocatable :: c_m_z(:)                 ! dim 1->num_t_ms
        real,   allocatable :: d_m_m(:,:) ! multi-scale HI mol tensor
        real,   allocatable :: s_m_m(:,:) ! multi-scale HI mol tensor                             
        real,   allocatable :: d_a_g(:,:)                              
        real,   allocatable :: d_a_g_x(:)                              
        real,   allocatable :: d_a_g_y(:)                              
        real,   allocatable :: d_a_g_z(:)                              
        real,   allocatable :: d_a_x(:,:)                              
        real,   allocatable :: d_a_y(:,:)                              
        real,   allocatable :: d_a_z(:,:)                              
        real,   allocatable :: s_a_x(:,:)                              
        real,   allocatable :: s_a_y(:,:)                              
        real,   allocatable :: s_a_z(:,:)                              
        real,   allocatable :: f_m_r(:)                              
        real,   allocatable :: r_m_r(:)                              
        real,   allocatable :: s_m_r(:)                              
        real,   allocatable :: f_e_x_st(:,:) ! ewald-elec and treecode    
        real,   allocatable :: f_e_x_md(:,:) ! treecode                   
        real,   allocatable :: f_a_x_pt(:,:)   ! force on protease     
        real,   allocatable :: f_a_x_st(:,:)                              
        real,   allocatable :: f_a_x_md(:,:)                              
        real,   allocatable :: f_a_x_lg(:,:)                              
        real,   allocatable :: f_a_x_bd(:,:)                              
        real,   allocatable :: f_a_x_ang(:,:)                              
        real,   allocatable :: f_a_x_dih(:,:)                              

        real,   allocatable :: pki_of_atm(:) ! intrinsic DG due to pKa 44-07-12
        real,   allocatable :: pka_of_atm(:) ! pKa of each atom - 44-07-12
        real,   allocatable :: qnf_a(:) ! current charge of each atom - 44-07-12
        real,   allocatable :: q_f_old_acc(:) ! old charge of accepted MC site
        real,   allocatable :: q_f_new_acc(:) ! old charge of accepted MC site
        integer             :: num_p_as ! stores # of protonatable sites
        integer,allocatable :: i_p_a(:) ! stores the protonatable site # of this atom   

C note that both the follwing had extra dimension 1:num_reps

        integer,allocatable :: n_f_a(:,:) ! =1 if protonated, =0 if not - 44-07-12
        integer*8,allocatable :: npf_a(:,:) ! =#steps this group protonated - 44-07-12

        integer,allocatable :: iip_acc(:) ! stores the site # where proton state changed
        integer,allocatable :: iii_acc(:) ! stores the corresponding atom #

        integer             :: num_osm
        real, allocatable   :: surface(:) ! SA for pressure calc
        real, allocatable   :: osm_cur(:)
        real, allocatable   :: osm_cur_loc(:)
        real, allocatable   :: osm_ave(:)
        real, allocatable   :: osm_var(:)

        real,   allocatable :: crg_of_atm(:) ! stores the partial charge of each flexible atom
        integer,allocatable :: stt_of_atm(:) ! stores the dynamic state of each flexible atom
        real,   allocatable :: r_t_a(:) ! stores the hydro radius of each hydrodynamic atom
        real,   allocatable :: d_t_a(:) ! stores the diff coeff of each hydrodynamic atom
        real,   allocatable :: d_t_a_3frm(:) ! same but done as a 1:3*atoms copy
        real,   allocatable :: d_t_a_root(:) ! stores the diff coeff of each hydrodynamic atom
        real,   allocatable :: d_t_a_root_inv(:) ! stores the diff coeff of each hydrodynamic atom
        integer,allocatable :: i_f_a(:) ! stores whether the atom is charged
        integer,allocatable :: i_q_a(:) ! (treecode_elec) stores the charge # of this atom   
        real,   allocatable :: q_r_a(:) ! stores the partial charge of each rigid atom
        real,   allocatable :: can1_x_rdf(:)      ! stores x of c.o.g. of each mol          
        real,   allocatable :: can1_y_rdf(:)      ! stores x of c.o.g. of each mol          
        real,   allocatable :: can1_z_rdf(:)      ! stores x of c.o.g. of each mol          
        integer(kind=1),allocatable :: i_a_g(:,:) ! integer version of d_a_g
        integer(kind=1),allocatable :: i_a_g_x(:) ! integer version of d_a_g
        integer(kind=1),allocatable :: i_a_g_y(:) ! integer version of d_a_g
        integer(kind=1),allocatable :: i_a_g_z(:) ! integer version of d_a_g
        integer,allocatable :: i_r_a(:) ! stores whether the atom is charged
        integer,allocatable :: h_r_a(:) ! stores the id of the nrst hydrodynamic atom
        logical,allocatable :: imv_unif_harm(:)                       ! logical - true if this particular harm restraint moves in uniform increments                
        integer,allocatable :: b_f_a(:)     ! stores the blob number of each fible a                                            
        integer,allocatable :: idm_molmol(:) ! thread-local id 
        integer,allocatable :: jdm_molmol(:) ! thread-local id 
        integer,allocatable :: idb_blbblb(:) ! thread-local id 
        integer,allocatable :: jdb_blbblb(:) ! thread-local id 
        integer,allocatable :: ida_atmatm(:) ! thread-local id 
        integer,allocatable :: jda_atmatm(:) ! thread-local id 
        integer,allocatable :: i_molmol(:) ! thread-local flag
        integer,allocatable :: i_blbblb(:,:) ! thread-local flag
        integer,allocatable :: i_blb(:,:)   ! stores the ids of atoms assigned to each blob                                     
        integer,allocatable :: ifree(:)     ! stores whether atom is free (default=1)                                           
        integer,allocatable :: jfree(:)     ! stores last free atom on same mol as current atom                                   
        integer,allocatable :: c_f_a(:)     ! when i_do_MC_moves stores the cluster number of each fible a                                            
        integer,allocatable :: mol_of_atm(:)                          ! stores the molecule number of each fible a                                            
        integer,allocatable :: m_r_a(:)              ! stores the molecule number of each rigd atom                                            
        integer,allocatable :: n_r_a(:)              ! stores the atom number of each rigd atom (e.g. 3 for 3rd atom of mol 3)     
        integer,allocatable :: my_num_t_typ(:)                    ! stores i_f_tp of this t typ if f,stores i_r_tp otherwise                    
        integer,allocatable :: my_typ_t_typ(:)                   ! =1 if this t typ of m is f,=0 if it's rigid                      
        integer,allocatable :: idm_uniq_m_pair(:)                   ! id of 1st m in each unique pair of mecular interactions                
        integer,allocatable :: jdm_uniq_m_pair(:)                   ! id of 2nd m in each unique pair of mecular interactions                
        integer,allocatable :: idm_f_tmp(:)                        ! used to temporarily assign to each f a its m num                    
        integer,allocatable :: ida_f_tmp(:)                        ! used to temporarily assign to each f a its a num                    
        integer,allocatable :: idm_temp_generate(:)                   ! m num of r a that needs to be generated each timestep               
        integer,allocatable :: ida_temp_generate(:)                   ! a num of r a that needs to be generated each timestep               
        integer,allocatable :: idm_must_generate(:)                   ! m num of r a that needs to be generated each timestep               
        integer,allocatable :: ida_must_generate(:)                   ! a num of r a that needs to be generated each timestep               
        integer,allocatable :: my_cell_1_ref_x(:)                     ! temp holds the x-cell ref of each a or m - speeds up cell assignment
        integer,allocatable :: my_cell_1_ref_y(:)
        integer,allocatable :: my_cell_1_ref_z(:)
        integer,allocatable :: my_cell_2_ref_x(:)
        integer,allocatable :: my_cell_2_ref_y(:)
        integer,allocatable :: my_cell_2_ref_z(:)
        integer,allocatable :: idm_grow_r(:)
        integer,allocatable :: nstep_grow_r(:)
        integer,allocatable :: irdf(:)          ! stores id of 1st mol type in rdf pair
        integer,allocatable :: jrdf(:)          ! stores id of 2nd mol type in rdf pair
        integer,allocatable :: id_beg(:)  ! number of 1st a minus one in this molecule in 1D list of all as in system (excludes dyn as)
        integer,allocatable :: id_end(:)
        integer,allocatable :: id_beg_bonds(:)
        integer,allocatable :: id_beg_angls(:)
        integer,allocatable :: id_beg_dihes(:)
        integer,allocatable :: num_bonds(:)
        integer,allocatable :: num_angls(:)
        integer,allocatable :: num_dihes(:)
        integer,allocatable :: nbp1(:)! for i_write_nonbonded_histogram
        integer,allocatable :: nbp2(:)! for i_write_nonbonded_histogram
        integer,allocatable :: iiib_new(:)    ! new for 10-08-08
        integer,allocatable :: nbb0_new(:)                           
        integer,allocatable :: nbb1_new(:)    ! new for 10-08-08
        integer,allocatable :: nbb2_new(:)
        integer,allocatable :: iiia_new(:)    ! new for 10-08-08
        integer,allocatable :: nba0_new(:)
        integer,allocatable :: nba1_new(:)
        integer,allocatable :: nba2_new(:)
        integer,allocatable :: nba3_new(:)
        integer,allocatable :: iiid_new(:)    ! new for 10-08-08
        integer,allocatable :: nbd0_new(:)
        integer,allocatable :: nbd1_new(:)
        integer,allocatable :: nbd2_new(:)
        integer,allocatable :: nbd3_new(:)
        integer,allocatable :: nbd4_new(:)
        integer,allocatable :: bnd_dnm(:)
        integer,allocatable :: nbb1(:)
        integer,allocatable :: nbb2(:)
        integer,allocatable :: nba1(:)
        integer,allocatable :: nba2(:)
        integer,allocatable :: nba3(:)
        integer,allocatable :: nbd0(:)
        integer,allocatable :: nbd1(:)
        integer,allocatable :: nbd2(:)
        integer,allocatable :: nbd3(:)
        integer,allocatable :: nbd4(:)
        integer,allocatable :: nxxmx(:)
        integer,allocatable :: nyxmx(:)
        integer,allocatable :: nzxmx(:)
        integer,allocatable :: nxxmn(:)
        integer,allocatable :: nyxmn(:)
        integer,allocatable :: nzxmn(:)
        integer,allocatable :: ityp(:)                     ! =1 if f,=0 if r,        dim: 1->num_t_ms
        integer,allocatable :: i_r_tp(:)                ! stores r typ of current m,dim: 1->num_t_ms
        integer,allocatable :: i_f_tp(:)                ! stores f typ of current m,dim: 1->num_t_ms
        integer,allocatable :: i_t_tp(:)                ! stores t typ of current m,dim: 1->num_t_ms
        integer,allocatable :: num_ths_f_typ(:)      ! stores number of this f typ of molecule dim: 1->num_t_ts
        integer,allocatable :: num_ths_r_typ(:)      ! stores number of this r typ of molecule dim: 1->num_t_ts
        integer,allocatable :: num_ths_t_typ(:)      ! stores number of this t typ of molecule dim: 1->num_t_ts
        integer,allocatable :: my_num_ths_t_typ(:)   ! stores the rank-number of this t typ of molecule dim: 1->num_t_ms
        integer,allocatable :: id1tmp(:)
        integer,allocatable :: id2tmp(:)
        integer,allocatable :: num_r_hs_pt(:) ! #hyd atms in rigd mol   
        integer,allocatable :: num_r_qs_pt(:) ! #crg atms in rigd mol
        integer,allocatable :: num_r_vs_pt(:) ! #vdw atms in rigd mol
        integer,allocatable :: num_r_ds_pt(:) ! #dyn atms in rigd mol
        integer,allocatable :: num_r_ss_pt(:) ! #all atms in rigd mol: sum
        integer,allocatable :: num_f_bs_pt(:) ! #blb atms in flex mol: all
        integer,allocatable :: num_atms_pt(:) ! #atms in mol type
        integer,allocatable :: num_dyns_pt(:) ! #atms in mol type (wrt)
        integer,allocatable :: np1xA(:)
        integer,allocatable :: np1yA(:)
        integer,allocatable :: np1zA(:)
        integer,allocatable :: np1xB(:)
        integer,allocatable :: np1yB(:)
        integer,allocatable :: np1zB(:)
        integer,allocatable :: inophiA(:)
        integer,allocatable :: inophiB(:)
        integer,allocatable :: i_done_bonds(:)
        integer,allocatable :: i_done_angls(:)
        integer,allocatable :: i_done_dihes(:)
        integer,allocatable :: link_f_m(:)
        integer,allocatable :: link_f_a(:)
        integer,allocatable :: link_r_m(:)
        integer,allocatable :: link_r_a(:)
        integer,allocatable :: idm_f_harm(:)              
        integer,allocatable :: idm_r_harm(:)              
        integer,allocatable :: fconst_harmonic(:)   
        integer,allocatable :: nstep_harm_req(:)   
        integer,allocatable :: istep_complete(:)  
        integer,allocatable :: istep_harmonic(:) 
        integer,allocatable :: ioffset_harmonic(:)                   
        integer,allocatable :: num_f_harm(:)                     
        integer,allocatable :: num_r_harm(:)                    
        integer,allocatable :: n_rxn_m_1(:)                     
        integer,allocatable :: n_rxn_m_2(:)                    
        integer,allocatable :: n_rxn_a_1(:)                   
        integer,allocatable :: n_rxn_a_2(:)                  
        integer,allocatable :: i_rxn_pbc(:)                   
        integer,allocatable :: i_omit_force_m_1(:)    
        integer,allocatable :: i_omit_force_m_2(:)   
        integer,allocatable :: i_omit_move_mtyp(:)    
        integer,allocatable :: i_omit_move_a(:)   
        integer,allocatable :: i_clash(:)                  
        integer,allocatable :: i_need_recompute(:)                  
        integer,allocatable :: z(:)                        
        integer,allocatable :: idq1(:)                    
        integer,allocatable :: idq2(:)                   

C stuff for f-f nonbonded interactions

        real,   allocatable :: s11_ff_rods(:)             ! a pair f_f_list  array
        real,   allocatable :: s22_ff_rods(:)             ! a pair f_f_list  array
        real,   allocatable :: s33_ff_rods(:)             ! a pair f_f_list  array
        real,   allocatable :: s44_ff_rods(:)             ! a pair f_f_list  array
        real,   allocatable :: d11_ff_rods(:)             ! a pair f_f_list  array
        real,   allocatable :: d22_ff_rods(:)             ! a pair f_f_list  array
        real,   allocatable :: e00_ff_rods(:)             ! a pair f_f_list  array
        real,   allocatable :: s11_ff_v_st(:)             ! a pair f_f_list  array
        real,   allocatable :: s22_ff_v_st(:)             ! a pair f_f_list  array
        real,   allocatable :: s33_ff_v_st(:)             ! a pair f_f_list  array
        real,   allocatable :: s44_ff_v_st(:)             ! a pair f_f_list  array
        real,   allocatable :: d11_ff_v_st(:)             ! a pair f_f_list  array
        real,   allocatable :: d22_ff_v_st(:)             ! a pair f_f_list  array
        real,   allocatable :: e00_ff_v_st(:)             ! a pair f_f_list  array
        real,   allocatable :: s11_ff_v_md(:)             ! a pair f_f_list  array
        real,   allocatable :: s22_ff_v_md(:)             ! a pair f_f_list  array
        real,   allocatable :: s33_ff_v_md(:)             ! a pair f_f_list  array
        real,   allocatable :: s44_ff_v_md(:)             ! a pair f_f_list  array
        real,   allocatable :: d11_ff_v_md(:)             ! a pair f_f_list  array
        real,   allocatable :: d22_ff_v_md(:)             ! a pair f_f_list  array
        real,   allocatable :: e00_ff_v_md(:)             ! a pair f_f_list  array
        real,   allocatable :: s11_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: s22_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: s33_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: s44_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: d11_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: d22_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: d12_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: ddd_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: e00_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: rep_ff_g_st(:)             ! a pair f_f_list  array
        real,   allocatable :: s11_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: s22_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: s33_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: s44_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: d11_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: d22_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: d12_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: ddd_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: e00_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: rep_ff_g_md(:)             ! a pair f_f_list  array
        real,   allocatable :: qqq_ff_e_st(:)             ! a pair f_f_list  array
        real,   allocatable :: qqq_ff_e_md(:)             ! a pair f_f_list  array
        real,   allocatable :: qqq_ff_e_lg(:)             ! a pair f_f_list  array
        real,   allocatable :: qqq_ff_e_st_t(:)             ! a pair f_f_list  array
        real,   allocatable :: qqq_ff_e_md_t(:)             ! a pair f_f_list  array
        real,   allocatable :: qqq_ff_e_lg_t(:)             ! a pair f_f_list  array
        real,   allocatable :: ene_pos_loc(:)  
        real,   allocatable :: ene_bond_loc(:)  
        real,   allocatable :: ene_angl_loc(:)  
        real,   allocatable :: ene_dihe_loc(:)  
        integer,allocatable :: ida_ff_a_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_a_st(:)             ! a pair f_f_list  array
        integer,allocatable :: kda_ff_a_st(:)             ! a pair f_f_list  array
        integer,allocatable :: lda_ff_a_st(:)             ! a pair f_f_list  array
        integer,allocatable :: ida_ff_a_md(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_a_md(:)             ! a pair f_f_list  array
        integer,allocatable :: kda_ff_a_md(:)             ! a pair f_f_list  array
        integer,allocatable :: lda_ff_a_md(:)             ! a pair f_f_list  array
        integer,allocatable :: ida_ff_b_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_b_st(:)             ! a pair f_f_list  array
        integer,allocatable :: kda_ff_b_st(:)             ! a pair f_f_list  array
        integer,allocatable :: lda_ff_b_st(:)             ! a pair f_f_list  array
        integer,allocatable :: ida_ff_b_md(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_b_md(:)             ! a pair f_f_list  array
        integer,allocatable :: kda_ff_b_md(:)             ! a pair f_f_list  array
        integer,allocatable :: lda_ff_b_md(:)             ! a pair f_f_list  array
        integer,allocatable :: ida_ff_u_st(:)             ! 42-07-12 f_list  array
        integer,allocatable :: jda_ff_u_st(:)             ! 42-07-12 f_list  array
        integer,allocatable :: kda_ff_u_st(:)             ! 42-07-12 f_list  array
        integer,allocatable :: ida_ff_i_st(:)             ! 44-07-12 f_list  array
        integer,allocatable :: jda_ff_i_st(:)             ! 44-07-12_f_list  array

        integer,allocatable :: ida_ff_r_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_r_st(:)             ! a pair f_f_list  array
        integer,allocatable :: kda_ff_r_st(:)             ! a pair f_f_list  array

        integer,allocatable :: ida_ff_rods(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_rods(:)             ! a pair f_f_list  array

        integer,allocatable :: ida_ff_v_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_v_st(:)             ! a pair f_f_list  array

        integer,allocatable :: ida_ff_v_md(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_v_md(:)             ! a pair f_f_list  array

        integer,allocatable :: ida_ff_g_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_g_st(:)             ! a pair f_f_list  array
        integer,allocatable :: kda_ff_g_st(:)             ! m #    f_f_list  array
        integer,allocatable :: lda_ff_g_st(:)             ! m #    f_f_list  array
        integer,allocatable :: mda_ff_g_st(:)             ! mode # f_f_list  array
        integer,allocatable :: nda_ff_g_st(:)             ! mode # f_f_list  array
        integer,allocatable :: ida_ff_g_md(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_g_md(:)             ! a pair f_f_list  array
        integer,allocatable :: kda_ff_g_md(:)             ! m #    f_f_list  array
        integer,allocatable :: lda_ff_g_md(:)             ! m #    f_f_list  array
        integer,allocatable :: mda_ff_g_md(:)             ! mode # f_f_list  array
        integer,allocatable :: nda_ff_g_md(:)             ! mode # f_f_list  array
        integer,allocatable :: ida_ff_e_st(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_e_st(:)             ! a pair f_f_list  array
        integer,allocatable :: ida_ff_e_md(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_e_md(:)             ! a pair f_f_list  array
        integer,allocatable :: ida_ff_e_lg(:)             ! a pair f_f_list  array
        integer,allocatable :: jda_ff_e_lg(:)             ! a pair f_f_list  array

        integer,allocatable :: idm_tmp_nbnd(:)                  
        character*80,allocatable :: qef_file(:)
        character*80,allocatable :: dyn_file(:)
        character*80,allocatable :: bnd_file(:)
        character*80,allocatable :: phi_file_A(:)
        character*80,allocatable :: phi_file_B(:)

C 2D arrays

        integer,allocatable :: idm_in_MC_clust(:,:)  ! for cluster
        integer,allocatable :: imol_in_clust(:,:)  ! for cluster
        integer,allocatable :: go_mod_arr(:,:,:)   ! for i_use_exclusive_go
        integer,allocatable :: go_mod_exc(:,:,:,:) ! for i_use_exclusive_go
        integer(kind=4),   allocatable :: d_status(:,:)! for i_write_nonbonded_histogram
        integer(kind=1),   allocatable :: i_status(:,:)                                 
        integer(kind=4),   allocatable :: g_status(:,:)                                 
        integer,allocatable :: first_go_entry(:,:)
        integer,allocatable :: last_go_entry(:,:)
        integer,allocatable :: hist_of_Q(:,:,:,:)
        integer,allocatable :: num_cur_Q(:,:,:) ! 3rd dim is replica #
        integer,allocatable :: num_req_Q(:,:)                                 
        integer,allocatable :: ida_h_d_ngb(:,:)                                 
        integer,allocatable :: ida_h_d_ngb_tot(:)                                 
        integer,allocatable :: jda_h_d_ngb_tot(:)                                 
        integer,allocatable :: idm_uniq_m_pair_omp(:,:)                   ! id of 1st m in each unique pair of mecular interactions                
        integer,allocatable :: jdm_uniq_m_pair_omp(:,:)                   ! id of 2nd m in each unique pair of mecular interactions                
        integer(kind=2),   allocatable :: id_a_r(:,:)                              
        integer(kind=2),   allocatable :: id_a_x(:,:)                              
        integer(kind=2),   allocatable :: id_a_y(:,:)                              
        integer(kind=2),   allocatable :: id_a_z(:,:)                              
        integer(kind=2),   allocatable :: is_a_r(:,:)                              
        integer(kind=2),   allocatable :: is_a_x(:,:)                              
        integer(kind=2),   allocatable :: is_a_y(:,:)                              
        integer(kind=2),   allocatable :: is_a_z(:,:)                              
c       real,   allocatable :: cut_g_st2(:,:)    
c       real,   allocatable :: cut_g_md2(:,:)    
        real,   allocatable :: rij_h_d_ngb_tot(:)                                 
        real,   allocatable :: B22_numer(:)                                        
        real,   allocatable :: B22_denom(:)                                        
        real,   allocatable :: B22_elec1(:)                                        
        real,   allocatable :: B22_elec2(:)                                        
        real,   allocatable :: B22_x1(:,:)                                        
        real,   allocatable :: B22_y1(:,:)                                        
        real,   allocatable :: B22_z1(:,:)                                        
        real,   allocatable :: B22_x2(:,:)                                        
        real,   allocatable :: B22_y2(:,:)                                        
        real,   allocatable :: B22_z2(:,:)                                        
        real,   allocatable :: tors1_new(:,:)
        real,   allocatable :: tors3_new(:,:)
        real,   allocatable :: tors1(:,:)
        real,   allocatable :: tors3(:,:)
        real,   allocatable :: rcutoff(:,:)
        real,   allocatable :: rcutoff2(:,:)
c       real,   allocatable :: d_a_prf(:,:)                !                                                         
        real,   allocatable :: dtrans(:,:)
        real,   allocatable :: drot(:,:)
        real,   allocatable :: d11(:,:) ! sigma for the following interactions
        real,   allocatable :: d22(:,:) ! squared sigma for the following interactions
        real,   allocatable :: e00(:,:) ! lowest epsilon for the following interactions
        real,   allocatable :: s11(:,:)
        real,   allocatable :: s22(:,:)
        real,   allocatable :: s33(:,:)
        real,   allocatable :: s44(:,:)
        real,   allocatable :: e11(:,:)
        real,   allocatable :: e22(:,:)
        real,   allocatable :: g(:,:)    
        real,   allocatable :: rhistdist(:,:)     ! distance histogram for each mol type pair
        integer,allocatable :: this_pair_can_go(:,:)        ! =1 for all m typ pairs that have go pairs
        integer,allocatable :: num_cur_go_pairs(:,:,:)  ! for all m typ pairs stores the # of go pairs
C                            ! note that 3rd dimension is the replica #
        integer,allocatable :: num_des_go_pairs(:,:)            ! for all m typ pairs stores the # of go pairs
        integer,allocatable :: jjj_corresponds_to_jdm(:,:)  ! for each m,records its # in list of iii's interactions
        integer,allocatable :: iflag_interaction(:,:)       ! notes whether two ms interact
        integer,allocatable :: maxhist(:,:)
        integer,allocatable :: nbnd_f_f(:,:)          ! =1 if this pair of f-f molecules are interacting
        integer,allocatable :: ihistdist(:,:)     ! distance histogram for each mol type pair
        integer,allocatable :: mrdf(:,:)          ! stores rdf-index of each mol type pair

C 3D arrays

        real,   allocatable :: d_a_m(:,:,:)   ! multi-scale HI atm tensor
        real,   allocatable :: s_a_m(:,:,:)   ! multi-scale HI atm tensor  
        real,   allocatable :: ene_v_tbl(:,:,:)
        real,   allocatable :: frc_v_tbl(:,:,:)
        integer,allocatable :: num_cur_Q_mod(:,:,:)                               
        integer,allocatable :: num_req_Q_mod(:,:,:)                               
        integer,allocatable :: num_des_go_mod_pairs(:,:,:)            ! for all m typ pairs stores the # of go pairs
        integer,allocatable :: num_cur_Q_loc(:,:,:)                               
        integer,allocatable :: num_cur_Q_mod_loc(:,:,:,:)                               

C user_hists will store 6D histograms of user_defined terms
C it will be allocated 1:num_user_hists

        type user_hist_arrays        
          integer,allocatable         :: hist(:,:,:,:,:,:)
        end type user_hist_arrays        

        type (user_hist_arrays),allocatable :: user_hists(:)        

C rigid_energy will store up to 6D user-defined energy functions 
C it will be allocated 1:num_rigid_energy

        logical                  :: i_rigid_energy   
        integer                  :: num_rigid_energy
        real,allocatable         :: rigid_energy_scal(:)
        real,allocatable         :: rigid_energy_ceil(:)
        character*80             :: rigid_energy_matchup_file
        character*80,allocatable :: rigid_energy_file(:)
      
        type rigid_energy_arrays        
          real,   allocatable         :: ene(:,:,:,:,:,:)
          real,   allocatable         :: frc(:,:,:,:,:,:,:)
        end type rigid_energy_arrays        

        type (rigid_energy_arrays),allocatable :: rigid_energy(:)        

C user_energy will store up to 6D user-defined energy functions 
C it will be allocated 1:num_user_funcs

        logical                  :: i_user_energy
        integer                  :: num_user_energy
        integer                  :: ltmp(6)
        real,allocatable         :: user_energy_scal(:)
        real,allocatable         :: user_energy_ceil(:)
        character*48             :: junkstring
        character*80             :: user_energy_matchup_file
        character*80,allocatable :: user_energy_file(:)
      
        type user_energy_arrays        
          real,   allocatable         :: ene(:,:,:,:,:,:)
          real,   allocatable         :: frc(:,:,:,:,:,:,:)
        end type user_energy_arrays        

        type (user_energy_arrays),allocatable :: user_energy(:)        

C user_match(1:num_uniqs) will store all possible user_energy
C functions that would be used for this unique atom

C e.g. for atom 11, involved in an alpha-helical interaction with both
C atoms 7 and 15, and a possible anti-parallel beta-sheet with atom 51,
C we might have the following setup:

C   user_match(11)%num=3
C   user_match(11)%typ(1)=1      ! for alpha-helical backwards
C   user_match(11)%typ(2)=1      ! for alpha-helical forwards
C   user_match(11)%typ(3)=2      ! for anti-parallel beta-sheet
C   user_match(11)%jda(1)=7      ! for alpha-helical backwards
C   user_match(11)%jda(2)=11     ! for alpha-helical backwards
C   user_match(11)%jda(3)=51     ! for anti-parallel beta-sheet
C   user_match(11)%scl(1)=1.0    ! extra scaling for propensity 
C   user_match(11)%scl(2)=1.0    ! extra scaling for propensity 
C   user_match(11)%scl(3)=1.0    ! extra scaling for propensity 
C   user_match(11)%npr(1)=5      ! 5 atoms define this user_energy
C   user_match(11)%npr(2)=5      ! 5 atoms define this user_energy
C   user_match(11)%npr(3)=6      ! 6 atoms define this user_energy
C   user_match(11)%kda(1,1) = 7
C   user_match(11)%kda(1,2) = 8
C   user_match(11)%kda(1,3) = 9
C   user_match(11)%kda(1,4) = 10
C   user_match(11)%kda(1,5) = 11
C   user_match(11)%kda(2,1) = 11
C   user_match(11)%kda(2,2) = 12
C   user_match(11)%kda(2,3) = 13
C   user_match(11)%kda(2,4) = 14
C   user_match(11)%kda(2,5) = 15
C   user_match(11)%kda(3,1) = 10
C   user_match(11)%kda(3,2) = 11
C   user_match(11)%kda(3,3) = 12
C   user_match(11)%kda(3,4) = 52
C   user_match(11)%kda(3,5) = 51
C   user_match(11)%kda(3,6) = 50

C all of these things would be defined when reading
C user_energy_function_matchup.txt
C then the idea would be every time we see a short-range interaction
C between 11 and any other atom, we loop from i=1,user_match(11)%num
C and see if user_match(11)%jda(i)= the other atom...

        type user_func_matchup
          integer                     :: num
          integer,allocatable         :: typ(:)
          integer,allocatable         :: jda(:)
          integer,allocatable         :: npr(:)
          real,   allocatable         :: scl(:)
          integer,allocatable         :: kda(:,:)
        end type user_func_matchup

        type (user_func_matchup),allocatable :: user_match(:)

C this derived type will be allocated(1:num_flx_atms)
C wall array - stores, for each unique atom type, the # of walls that
C it sees and their identity - can be up to 100...
C note that elsewhere we have crd_wall(:),pos_wall(:),fct_wall(:)

        type wall_arrays        
          integer                                          :: num    
c         integer,dimension(:),pointer                     :: iwl    
          integer,allocatable                              :: iwl(:)
        end type wall_arrays        

        type (wall_arrays),allocatable :: wall_ptr(:)        

C this derived type will be allocated(1:num_flx_atms)
C %nbterm  : # bodies interacting with this atom 
C %ncterm  : # cells  interacting with this atom 
C %ndterm  : # cells  interacting with this atom 
C %jda(:)  : body/cell id of the interacting body
C %d1(:)   : diff tensor term of the interacting cell/body
C %d9(:)   : diff tensor term of the interacting cell/body

        type hydro_tree_arrays    
          integer                                          :: nbterm 
          integer                                          :: ncterm 
          integer                                          :: ndterm 
          integer,allocatable                              :: jda(:)
          real,   allocatable                              :: d1(:)
          real,   allocatable                              :: d2(:)
          real,   allocatable                              :: d3(:)
          real,   allocatable                              :: d4(:)
          real,   allocatable                              :: d5(:)
          real,   allocatable                              :: d6(:)
          real,   allocatable                              :: d7(:)
          real,   allocatable                              :: d8(:)
          real,   allocatable                              :: d9(:)
        end type hydro_tree_arrays  

        type (hydro_tree_arrays),allocatable :: hi_tree(:)     

C i_do_protease :

C note that this derived type will be allocated (1:num_flx_atms)
C this derived type will be allocated(1:num_q_as)
C %ida     : atom # corresponding to this charge #
C %num_bnd : # charges bonded to this charge
C %jda(:)  : charge id of the bonded charge

        type charge_arrays    
          integer                                          :: ida    
          integer                                          :: num_bnd
          integer,allocatable                              :: jda(:)
        end type charge_arrays  

        type (charge_arrays),allocatable :: q_info(:)     

C note that this derived type will be allocated (1:num_p_as)
C %ida     : atom # corresponding to this site   (1<%jda<=num_flx_atms)
C %idq     : charge # corresponding to this site (1<%jda<=num_q_as)
C           (remember that there will probably be more charges than protonatable sites)

        type protonatable_arrays    
          integer                                          :: ida    
          integer                                          :: idq
        end type protonatable_arrays  

        type (protonatable_arrays),allocatable :: p_info(:) ! protonatable sites

C note that this derived type will be allocated (1:num_p_as2)
C this one is for the pairs of protonatable sites that are close by

C %ida     : atom # corresponding to this site   (1<%jda<=num_flx_atms)
C %idq     : charge # corresponding to this site (1<%jda<=num_q_as)
C           (remember that there will probably be more charges than protonatable sites)
C %jda     : atom # corresponding to paired site   (1<%jda<=num_flx_atms)
C %jdq     : charge # corresponding to paired site (1<%jda<=num_q_as)

        type protonatable_pair_arrays    
          integer                                          :: ida    
          integer                                          :: idq
          integer                                          :: jda
          integer                                          :: jdq
        end type protonatable_pair_arrays  

        type (protonatable_pair_arrays),allocatable :: p_info2(:) ! protonatable sites

C i_do_protease :

C note that this derived type will be allocated (1:num_flx_atms)

C %num     : # of cleavage bonds associated with this atom
C %jda(:)  : id of the substrate atom associated with this bond
C %jdb(:)  : # of the bond in the full 1D list  (poss. value from 1 to num_tot_protease_bonds)
C %jds(:)  : # of the site in the full 1D list  (poss. value from 1 to num_tot_protease_sites)
C %jdx(:)  : # of the active site in the full 1D list  (poss. value from 1 to num_tot_protease_actvs)
C %crt(:)  : dist criterion of the bond in the full 1D list 
C %rng(:)  : dist at which force kicks in between protease & substrate
C %fct(:)  : force constant that acts between atoms of protease & substrate

        type protease_arrays
          integer                                          :: num    
          integer,allocatable                              :: jda(:)
          integer,allocatable                              :: jdb(:)
          integer,allocatable                              :: jds(:)
          integer,allocatable                              :: jdx(:)
          real,   allocatable                              :: crt(:)
          real,   allocatable                              :: rng(:)
          real,   allocatable                              :: fct(:)
        end type protease_arrays

        type (protease_arrays),allocatable :: protease_info(:)     


C coord arrays

C real_coord arrays

        type real_coord_arrays
          integer                                         :: num    
c         integer,dimension(:),pointer                     :: id1    
          integer,allocatable                              :: id1(:)
        end type real_coord_arrays  

        type (real_coord_arrays),allocatable :: pkref(:)

C f coord arrays
C first set depend only on # of molecules

        type c_f_m_arrays
          integer                                  :: num    
c         real,   dimension(:),pointer             :: x
c         real,   dimension(:),pointer             :: y
c         real,   dimension(:),pointer             :: z
c         integer,dimension(:),pointer             :: free
c         integer,dimension(:),pointer             :: harm      ! denotes whether involved in a harmonic bond
c         integer,dimension(:),pointer             :: harm_gone ! denotes whether just released from a harmonic bond
c         integer,dimension(:),pointer             :: just_harm ! denotes whether just placed in a harmonic bond
c         integer,dimension(:),pointer             :: just_free ! used to protect newly-free as from force-limit
c         integer,dimension(:),pointer             :: near_free ! used to denote next free as
          real,   allocatable                   :: x(:)
          real,   allocatable                   :: y(:)
          real,   allocatable                   :: z(:)
          integer,allocatable                   :: free(:)
          integer,allocatable                   :: harm(:)      ! denotes whether involved in a harmonic bond
          integer,allocatable                   :: harm_gone(:) ! denotes whether just released from a harmonic bond
          integer,allocatable                   :: just_harm(:) ! denotes whether just placed in a harmonic bond
          integer,allocatable                   :: just_free(:) ! used to protect newly-free as from force-limit
          integer,allocatable                   :: near_free(:) ! used to denote next free as
        end type c_f_m_arrays                                   
        type (c_f_m_arrays),allocatable :: c_f_m(:)


C second set depends on # of typs of molecules - note that this applies
C to both flexible and rigid molecules

C qnoproton is the charge on the group when unprotonated
C pKa is the intrinsic pKa of the group

        type c_typ_arrays
          integer                             :: num    
          real,   allocatable                 :: x(:)
          real,   allocatable                 :: y(:)
          real,   allocatable                 :: z(:)
          real,   allocatable                 :: q(:)  ! charge
          real,   allocatable                 :: m(:)  ! mass (for dpd)
          real,   allocatable                 :: r(:)  ! hyd rad
          real,   allocatable                 :: qnoproton(:)  ! charge
          real,   allocatable                 :: pKa(:)  ! hyd rad
          integer,allocatable                 :: a(:) ! atom typ
          integer,allocatable                 :: s(:) ! =1 if no move
                                                      ! 2019
        end type c_typ_arrays                                   

        type (c_typ_arrays),allocatable :: c_typ(:)

C character arrays

        type character_arrays
          integer                             :: num    
          character*6,allocatable             :: ch1(:)    
          character*4,allocatable             :: ch2(:)    
          character*3,allocatable             :: ch3(:)    
          character*4,allocatable             :: ch4(:)    
          character*2,allocatable             :: ch5(:)    
        end type character_arrays                                            

        type (character_arrays),allocatable :: n_typ(:)             

C energy arrays

        type energy_arrays
          integer                                    :: num    
c         integer,dimension(:),pointer             :: jdm  ! id of interacting partner    
c         real,dimension(:),pointer                :: e    
c         real,dimension(:),pointer                :: v    
c         real,dimension(:),pointer                :: g    
c         real,dimension(:),pointer                :: dis    
c         real,dimension(:),pointer                :: tot    
          integer,allocatable                 :: jdm(:)  ! id of interacting partner    
          real,allocatable                    :: e(:)    
          real,allocatable                    :: v(:)    
          real,allocatable                    :: g(:)    
          real,allocatable                    :: dis(:)    
          real,allocatable                    :: tot(:)    
        end type energy_arrays                                            

        type (energy_arrays),allocatable :: ene(:)                

C integer_coord arrays

        type integer_coord_arrays
          integer                                    :: num    
c         integer,dimension(:),pointer             :: id1    
          integer,allocatable                 :: id1(:)
        end type integer_coord_arrays               

        type (integer_coord_arrays),allocatable :: nptyp(:)  

C rigid domain type arrays

        type rigid_domain_arrays 
          real                                :: dtrn    ! translational diff coeff
          real                                :: drot    ! rotational diff coeff
          integer                             :: num     ! num of atoms in this domain type
          integer                             :: mol_typ ! moltype associated with this domain
          integer                             :: req_stp ! how often to flexibly move
          integer                             :: cur_stp ! current counter for flex move 
          integer,allocatable                 :: ida(:)  ! atom # (within mol typ)
        end type rigid_domain_arrays 
                  
        type (rigid_domain_arrays),allocatable :: rgd_dom_typ(:)

C position restraint arrays

        type pos_restraint_arrays
          real                                       :: fconst
          integer                                    :: idm       ! id of m that's restrained 
          integer                                    :: num       ! num of as that are restrained
c         integer,dimension(:),pointer             :: ida      
c         real,dimension(:),pointer                :: x               
c         real,dimension(:),pointer                :: y               
c         real,dimension(:),pointer                :: z               
          integer,allocatable                 :: ida(:)   
          real,allocatable                    :: x(:)            
          real,allocatable                    :: y(:)            
          real,allocatable                    :: z(:)            
        end type pos_restraint_arrays
                  
        type (pos_restraint_arrays),allocatable :: rest_info_f(:)
        type (pos_restraint_arrays),allocatable :: rest_info_r(:)

C intra_hydrodynamics arrays - each allocated as 1->num_mols

        type intra_hydro_int_arrays
          integer                                       :: num    !
          integer,allocatable                           :: iii(:) !
        end type intra_hydro_int_arrays                                   

        type (intra_hydro_int_arrays),allocatable :: imp_int(:)
        type (intra_hydro_int_arrays),allocatable :: jmp_int(:)

        type intra_hydro_real_arrays
          integer                                       :: num    !
          real,allocatable                              :: iii(:) !
        end type intra_hydro_real_arrays                                   

        type (intra_hydro_real_arrays),allocatable :: amp_int(:)

C this derived type holds the schedule for doing intra_hydrodynamics
C note that not all parts are used for diff_scdl, chol_scdl, move_scdl
C note also that we assume that in move_scdl each thread acts on no more
C than *one* molecule per round...
C
C diff_scdl is difftens schedule: gets allocated 1->num_diff_rounds
C chol_scdl is cholesky schedule: gets allocated 1->num_chol_rounds
C move_scdl is move     schedule: gets allocated 1->num_move_rounds
C
C %num : number of molecules done per round
C %upd : how often to update the molecules treated in this round
C %cnt : keeps count of the updates for molecules treated in this round
C %mst : =1 if this thread is a master thread this round
C %idm : id of molecules to be done this round
C %jdm : id of molecules attached to this thread this round
C %ida : id of first atom in this mol to be done this round
C %jda : id of last atom in this mol to be done this round
C %nda : (for domain stuff) total # of atoms in this domain
C %trd : #threads to be used for each molecule this round
C %scl : scale_nb value for each molecule this round

        type intra_hydro_schedule_arrays
          integer                                       :: num  !
          integer                                       :: upd  !
          integer                                       :: cnt  !
c         integer,dimension(:),pointer                  :: mst  !
c         integer,dimension(:),pointer                  :: idm  !
c         integer,dimension(:),pointer                  :: jdm  !
c         integer,dimension(:),pointer                  :: ida  !
c         integer,dimension(:),pointer                  :: jda  !
c         integer,dimension(:),pointer                  :: trd  !
c         real,   dimension(:),pointer                  :: scl  !
          integer,allocatable                      :: mst(:)  !
          integer,allocatable                      :: idd(:)  ! domain
          integer,allocatable                      :: jdd(:)  ! domain
          integer,allocatable                      :: idm(:)  ! mol
          integer,allocatable                      :: jdm(:)  ! mol
          integer,allocatable                      :: ida(:)  ! atm
          integer,allocatable                      :: jda(:)  ! atm
          integer,allocatable                      :: nda(:)  !ats in do
          integer,allocatable                      :: trd(:)  !
          real,   allocatable                      :: scl(:)  !
        end type intra_hydro_schedule_arrays                                   

        type (intra_hydro_schedule_arrays),allocatable :: diff_scdl(:)
        type (intra_hydro_schedule_arrays),allocatable :: chol_scdl(:)
        type (intra_hydro_schedule_arrays),allocatable :: move_scdl(:)

C this derived type lets us store intramolecular diffusion tensors nicely

        type intra_hydro_diffusion_tensor                                                               
          integer                                       :: idm                                       
          integer                                       :: jdm                                       
          real,allocatable                              :: val(:,:)                                               
        end type intra_hydro_diffusion_tensor
                                                           
        type (intra_hydro_diffusion_tensor),allocatable   ::  d_i_m(:)                                   
        type (intra_hydro_diffusion_tensor),allocatable   ::  s_i_m(:)                                   

C this derived type controls the schedule for growing atoms
 
C %req_stp : number of steps of growth that this molecule undergoes
C           (i.e. number of residues to be grown/disappeared)
C %cur_stp : which growth step this particular molecule is currently on
C %cur_mini_stp : which mini-step of growth this molecule is on
C %idm : id of the molecule 
C %ida : id of first atom in this mol to be free this step  
C %jda : id of last  atom in this mol to be free this step  
C %req_mini_stp : number of mini-steps before this molecule grows again

        type growth_schedule_arrays
          integer                                 :: req_stp  !
          integer                                 :: cur_stp  ! 
          integer                                 :: cur_mini_stp  !
          integer                                 :: idm  !
c         integer,dimension(:),pointer            :: ida  !
c         integer,dimension(:),pointer            :: jda  !
c         integer,dimension(:),pointer            :: req_mini_stp  !
          integer,allocatable                :: ida(:)  !
          integer,allocatable                :: jda(:)  !
          integer,allocatable                :: req_mini_stp(:)  !
        end type growth_schedule_arrays                                   

        type (growth_schedule_arrays),allocatable :: grow_scdl(:)

C this derived type controls the schedule for sliding atoms
  
C %num_bnd : number of bonds that make up this particular contact      
C %req_stp : number of steps of sliding that this molecule undergoes
C           (i.e. number of residues to be grown/disappeared)
C %cur_stp : which growth step this particular molecule is currently on
C %cur_mini_stp : which mini-step of growth this molecule is on
C %idm : id of the first (non-sliding) molecule 
C %jdm : id of the second (sliding) molecule 
C %ida : ids of the atoms on the first molecule in this contact
C %jda : ids of the atoms on the second molecule in this contact
C        note that there are req_stp copies of these
C %bnd : bond lengths of the bonds in this contact
C %fct : force constants of the bonds in this contact
C %req_mini_stp : number of mini-steps before this molecule grows again
C %req_trns_stp : number of mini-steps over which the movement occurs 

        type slide_schedule_arrays
          integer                                 :: num_bnd  !
          integer                                 :: req_stp  !
          integer                                 :: cur_stp  ! 
          integer                                 :: cur_mini_stp  !
          integer                                 :: idm  !
          integer                                 :: jdm  !
c         integer,dimension(:),pointer            :: ida  !
c         integer,dimension(:),pointer            :: req_mini_stp  !
c         integer,dimension(:),pointer            :: req_trns_stp  !
c         real,   dimension(:),pointer            :: bnd  !
c         real,   dimension(:),pointer            :: fct  !
c         integer,dimension(:,:),pointer          :: jda  !
          integer,allocatable                :: ida(:)  !
          integer,allocatable                :: req_mini_stp(:)  !
          integer,allocatable                :: req_trns_stp(:)  !
          real,   allocatable                :: bnd(:)  !
          real,   allocatable                :: fct(:)  !
          integer,allocatable                :: jda(:,:)  !
        end type slide_schedule_arrays                                   

        type (slide_schedule_arrays),allocatable :: slide_scdl(:)

C this derived type controls protease activity
  
C %num_bnd : number of bonds that make up this particular contact      
C %num_sts : number of sites where cutting can occur for this combo   
C %num_lim : number of steps before limited proteolysis is tried      
C %num_cur : counter for if we're doing limited proteolysis this step (unused)
C %acy_len : length of acyl bond b/w atoms of the 1st bond 
C %acy_fct : force constant of acyl bond b/w atoms of the 1st bond 
C %acy_stp : # of steps before the acyl bond decays 
C %idm     : id of the protease molecule 
C %jdm     : id of the substrate molecule 
C %clv     : has the bond been cleaved? =0 if no
C %prb     : probability that cut will happen if all bonds satisfied -
C            note that we allow different cleavage probs for diff bonds
C %ida     : ids of the atoms on the protease 
C %jda     : ids of the atoms on the substrate
C            note that there are num_sts copies of these
C %crt     : required distances for the bonds in this contact (distance criteria)
C %rng     : range criterion - force acts within this distance
C %fct     : force constants of the bonds in this contact 

        type protease_schedule_arrays
          integer                                 :: num_bnd  !
          integer                                 :: num_sts  !
          integer                                 :: num_lim  !
          integer                                 :: num_cur  !
          real                                    :: acy_len  !
          real                                    :: acy_fct  !
          real                                    :: acy_stp  !
          integer                                 :: idm  !
          integer                                 :: jdm  !
c         integer,dimension(:),pointer            :: clv  ! (1:num_sts)
c         real,   dimension(:),pointer            :: prb  ! (1:num_sts)
c         integer,dimension(:),pointer            :: ida  ! (1:num_bnd)
c         real,   dimension(:),pointer            :: crt  ! .
c         real,   dimension(:),pointer            :: rng  ! .
c         real,   dimension(:),pointer            :: fct  ! .
c         integer,dimension(:,:),pointer          :: jda  !
          integer,allocatable                :: clv(:)  ! (1:num_sts)
          real,   allocatable                :: prb(:)  ! (1:num_sts)
          integer,allocatable                :: ida(:)  ! (1:num_bnd)
          real,   allocatable                :: crt(:)  ! .
          real,   allocatable                :: rng(:)  ! .
          real,   allocatable                :: fct(:)  ! .
          integer,allocatable                :: jda(:,:)  !
        end type protease_schedule_arrays                                   

        type (protease_schedule_arrays),allocatable :: protease_scdl(:)

C force_exclusion arrays

        type force_exclusion_arrays
          integer                                         :: num    ! 1 -> num_f_ts
c         integer,dimension(:),pointer                  :: id1    ! 1 -> num_flx_atms
c         integer,dimension(:,:),pointer                :: id2    ! 1 -> number of neighbor bonded as
          integer,allocatable                           :: id1(:)    ! 1 -> num_flx_atms
          integer,allocatable                           :: id2(:,:)    ! 1 -> number of neighbor bonded as
        end type force_exclusion_arrays                                   

        type (force_exclusion_arrays),allocatable :: b_typ(:)

C new_nonbonded arrays

        type f_f_nonbonded_arrays
          real                              :: qqq           ! product of charges                                   
          real                              :: s11           ! v term                                                
          real                              :: s22           ! v term                                                
          real                              :: s33           ! v term                                                
          real                              :: s44           ! v term                                                
          real                              :: d11           ! v term                                                
          real                              :: d22           ! v term                                                
          real                              :: e00           ! v term                                                
          integer                           :: idm           ! id of iii m in the list
          integer                           :: ida           ! id of iii a in the list of f as in that m
          integer                           :: iat           ! id of iii a in the list of all f as
          integer                           :: jdm           ! id of jjj m in the list
          integer                           :: jda           ! id of jjj a in the list of f as in that m
          integer                           :: jat           ! id of iii a in the list of all f as
          integer                           :: igo           ! flag = 1 if go pair = 0 otherwise
          integer                           :: noq           ! flag = 1 if no ectrostatic component necessary
        end type f_f_nonbonded_arrays                                 
        type (f_f_nonbonded_arrays),
     &        allocatable::f_f_list(:) ! dimension here is num_a_pairs    


C cell arrays used for storing location of all as of fible molecules or molecule-centers of rigid molecules  

        type cll_f_array                                                                    
          integer                                  :: num     ! total number of flexible atoms contained in each cell 
c         integer,dimension(:),pointer             :: ida     ! a number of the i'th member of the cell
c         integer,dimension(:),pointer             :: idm     ! m number of the i'th member of the cell
          integer,allocatable                      :: ida(:)  ! a number of the i'th member of the cell
          integer,allocatable                      :: idm(:)  ! m number (not yet used) 
          integer,allocatable                      :: flg(:)  ! use in MC: =1 for old, =2 for new
        end type cll_f_array                                                        

        type (cll_f_array),allocatable  ::   cll_f_1(:,:,:,:) ! first 3D array containing cell members
        type (cll_f_array),allocatable  ::   cll_f_2(:,:,:,:) ! as above but for static atoms only
        type (cll_f_array),allocatable  ::   cll_f_r(:,:,:,:) ! same as above but for rods
                                        ! extra dimension is num_reps


      END MODULE allocatable_arrays


