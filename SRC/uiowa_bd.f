
      include 'mkl_vsl.fi'
      use allocatable_arrays    
      USE MKL_VSL_TYPE
      USE MKL_VSL
           
      implicit real(a-h,o-z) 



      interface
      subroutine finalize_long_range_forces(
     &             miv_num_loc,miv_beq_loc,miv_enq_loc,
     &             num_q_zs,num_tot_atmz,
     &             fff_ll_e_st,fff_ll_e_md,fff_ll_e_lg,
     &             zpengtar,zforce,perm_cpz,
     &             eee_e_st_loc,eee_e_md_loc,eee_e_lg_loc)

      integer, intent(in)     :: miv_num_loc
      integer, intent(in)     :: miv_beq_loc(miv_num_loc)
      integer, intent(in)     :: miv_enq_loc(miv_num_loc)
      integer, intent(in)     :: num_q_zs                        ! num_q_as
      integer, intent(in)     :: perm_cpz(1:num_q_zs)
      real,    intent(inout)  :: fff_ll_e_st(1:3,1:num_tot_atmz) ! fff_ll_e_st
      real,    intent(inout)  :: fff_ll_e_md(1:3,1:num_tot_atmz) ! fff_ll_e_md
      real,    intent(inout)  :: fff_ll_e_lg(1:3,1:num_tot_atmz) ! fff_ll_e_lg
      real,    intent(inout)  :: zpengtar(1:num_tot_atmz)        ! zpengtar
      real,    intent(inout)  :: zforce(1:3,1:num_tot_atmz)      ! zforce
      real,    intent(in)     :: eee_e_st_loc                    !
      real,    intent(in)     :: eee_e_md_loc                    !
      real,    intent(inout)  :: eee_e_lg_loc                    !

      end subroutine
      end interface


      interface
      subroutine assign_atoms_to_grid(ndum1,mithrd,mirep,
     &           i1_lic,i2_lic,i3_lic)

      integer ndum1 ! how many atoms
      integer mirep ! myrep
      integer i1_lic(ndum1),i2_lic(ndum1),i3_lic(ndum1)

      end subroutine
      end interface


C 2019 add a second grid for storing static atoms

      interface
      subroutine assign_static_atoms_to_grid(ndum1,mithrd,mirep,
     &           i1_lic,i2_lic,i3_lic)

      integer ndum1 ! how many atoms
      integer mirep ! myrep
      integer i1_lic(ndum1),i2_lic(ndum1),i3_lic(ndum1)

      end subroutine
      end interface


      interface
      subroutine calculate_treecode_elec(mythrd_loc,num_q_zs,
     &             miv_num_loc,miv_beq_loc,miv_enq_loc,
     &             crg_of_atm_cpz,perm_cpz,
     &             zpengtar,zforce,
     &             c_q_x_cpz,c_q_y_cpz,c_q_z_cpz,
     &             c_q_x_taz,c_q_y_taz,c_q_z_taz)

      integer, intent(in)    :: mythrd_loc
      integer, intent(in)    :: num_q_zs
      integer, intent(in)    :: miv_num_loc
      integer, intent(in)    :: miv_beq_loc(miv_num_loc)
      integer, intent(in)    :: miv_enq_loc(miv_num_loc)
      integer, intent(inout) :: perm_cpz(1:num_q_zs)
      real,    intent(inout) :: crg_of_atm_cpz(1:num_q_zs)
      real,    intent(inout) :: zpengtar(1:num_q_zs)
      real,    intent(inout) :: zforce(1:3,1:num_q_zs)
      real,    intent(inout) :: c_q_x_cpz(1:num_q_zs)
      real,    intent(inout) :: c_q_y_cpz(1:num_q_zs)
      real,    intent(inout) :: c_q_z_cpz(1:num_q_zs)
      real,    intent(inout) :: c_q_x_taz(1:num_q_zs)
      real,    intent(inout) :: c_q_y_taz(1:num_q_zs)
      real,    intent(inout) :: c_q_z_taz(1:num_q_zs)

      end subroutine
      end interface


      interface
      subroutine move_atoms(num_cnt,mirep,
     &                      miv_num_loc,miv_typ_loc,miv_beg_loc,
     &                      miv_end_loc,miv_dyn_loc,miv_hyd_loc)


      integer,              intent(in) :: num_cnt
      integer,              intent(in) :: mirep
      integer,              intent(in) :: miv_num_loc
      integer,              intent(in) :: miv_typ_loc(miv_num_loc)
      integer,              intent(in) :: miv_beg_loc(miv_num_loc)
      integer,              intent(in) :: miv_end_loc(miv_num_loc)
      integer,              intent(in) :: miv_dyn_loc(miv_num_loc)
      integer,              intent(in) :: miv_hyd_loc(miv_num_loc)

      end subroutine
      end interface


      interface 
      subroutine make_nonbond_list(mythrd,num_threads_loc,mirep,
     &           miv_num_loc,miv_beg_loc,miv_end_loc,
     &           n_ll_r_st,ida_ll_r_st,jda_ll_r_st,kda_ll_r_st,
     &           n_ll_u_st,ida_ll_u_st,jda_ll_u_st,kda_ll_u_st,
     &           n_ll_a_st,ida_ll_a_st,jda_ll_a_st,kda_ll_a_st,
     &           n_ll_a_md,ida_ll_a_md,jda_ll_a_md,kda_ll_a_md,
     &           n_ll_b_st,ida_ll_b_st,jda_ll_b_st,kda_ll_b_st,
     &                     lda_ll_b_st,
     &           n_ll_b_md,ida_ll_b_md,jda_ll_b_md,kda_ll_b_md,
     &                     lda_ll_b_md,
     &           n_ll_v_st,ida_ll_v_st,jda_ll_v_st,s11_ll_v_st,
     &                     s22_ll_v_st,s33_ll_v_st,s44_ll_v_st,
     &                     d11_ll_v_st,d22_ll_v_st,e00_ll_v_st,
     &           n_ll_v_md,ida_ll_v_md,jda_ll_v_md,s11_ll_v_md,
     &                     s22_ll_v_md,s33_ll_v_md,s44_ll_v_md,
     &                     d11_ll_v_md,d22_ll_v_md,e00_ll_v_md,
     &           n_ll_g_st,ida_ll_g_st,jda_ll_g_st,kda_ll_g_st,
     &         lda_ll_g_st,mda_ll_g_st,nda_ll_g_st,s11_ll_g_st,
     &                     s22_ll_g_st,s33_ll_g_st,s44_ll_g_st,
     &                     d11_ll_g_st,d22_ll_g_st,d12_ll_g_st,
     &                     ddd_ll_g_st,e00_ll_g_st,rep_ll_g_st,
     &           n_ll_g_md,ida_ll_g_md,jda_ll_g_md,kda_ll_g_md,
     &         lda_ll_g_md,mda_ll_g_md,nda_ll_g_md,s11_ll_g_md,
     &                     s22_ll_g_md,s33_ll_g_md,s44_ll_g_md,
     &                     d11_ll_g_md,d22_ll_g_md,d12_ll_g_md,
     &                     ddd_ll_g_md,e00_ll_g_md,rep_ll_g_md,
     &           n_ll_g_rx,ida_ll_g_rx,jda_ll_g_rx,kda_ll_g_rx,
     &                     lda_ll_g_rx,mda_ll_g_rx,s11_ll_g_rx,
     &                     s22_ll_g_rx,s33_ll_g_rx,s44_ll_g_rx,
     &                     d11_ll_g_rx,d22_ll_g_rx,d12_ll_g_rx,
     &                     ddd_ll_g_rx,e00_ll_g_rx,
     &           n_ll_e_st,ida_ll_e_st,jda_ll_e_st,qqq_ll_e_st,
     &           n_ll_e_md,ida_ll_e_md,jda_ll_e_md,qqq_ll_e_md,
     &     n_ll_r_st_old,n_ll_u_st_old,n_ll_a_st_old,n_ll_a_md_old,
     &     n_ll_b_st_old,n_ll_b_md_old,n_ll_v_st_old,n_ll_v_md_old,
     &     n_ll_g_st_old,n_ll_g_md_old,n_ll_e_st_old,n_ll_e_md_old,
     &     n_ll_g_rx_old,goscale)

      integer,              intent(in)    :: mirep
      integer,              intent(in)    :: miv_num_loc
      integer,              intent(in)    :: miv_beg_loc(miv_num_loc)
      integer,              intent(in)    :: miv_end_loc(miv_num_loc)
      integer, allocatable, intent(inout) :: ida_ll_r_st(:)
      integer, allocatable, intent(inout) :: jda_ll_r_st(:)
      integer, allocatable, intent(inout) :: kda_ll_r_st(:)
      integer, allocatable, intent(inout) :: ida_ll_u_st(:)
      integer, allocatable, intent(inout) :: jda_ll_u_st(:)
      integer, allocatable, intent(inout) :: kda_ll_u_st(:)
      integer, allocatable, intent(inout) :: ida_ll_a_st(:)
      integer, allocatable, intent(inout) :: jda_ll_a_st(:)
      integer, allocatable, intent(inout) :: kda_ll_a_st(:)
      integer, allocatable, intent(inout) :: ida_ll_a_md(:)
      integer, allocatable, intent(inout) :: jda_ll_a_md(:)
      integer, allocatable, intent(inout) :: kda_ll_a_md(:)
      integer, allocatable, intent(inout) :: ida_ll_b_st(:)
      integer, allocatable, intent(inout) :: jda_ll_b_st(:)
      integer, allocatable, intent(inout) :: kda_ll_b_st(:)
      integer, allocatable, intent(inout) :: lda_ll_b_st(:)
      integer, allocatable, intent(inout) :: ida_ll_b_md(:)
      integer, allocatable, intent(inout) :: jda_ll_b_md(:)
      integer, allocatable, intent(inout) :: kda_ll_b_md(:)
      integer, allocatable, intent(inout) :: lda_ll_b_md(:)
      integer, allocatable, intent(inout) :: ida_ll_v_st(:)
      integer, allocatable, intent(inout) :: jda_ll_v_st(:)
      real,    allocatable, intent(inout) :: s11_ll_v_st(:)
      real,    allocatable, intent(inout) :: s22_ll_v_st(:)
      real,    allocatable, intent(inout) :: s33_ll_v_st(:)
      real,    allocatable, intent(inout) :: s44_ll_v_st(:)
      real,    allocatable, intent(inout) :: d11_ll_v_st(:)
      real,    allocatable, intent(inout) :: d22_ll_v_st(:)
      real,    allocatable, intent(inout) :: e00_ll_v_st(:)
      integer, allocatable, intent(inout) :: ida_ll_v_md(:)
      integer, allocatable, intent(inout) :: jda_ll_v_md(:)
      real,    allocatable, intent(inout) :: s11_ll_v_md(:)
      real,    allocatable, intent(inout) :: s22_ll_v_md(:)
      real,    allocatable, intent(inout) :: s33_ll_v_md(:)
      real,    allocatable, intent(inout) :: s44_ll_v_md(:)
      real,    allocatable, intent(inout) :: d11_ll_v_md(:)
      real,    allocatable, intent(inout) :: d22_ll_v_md(:)
      real,    allocatable, intent(inout) :: e00_ll_v_md(:)
      integer, allocatable, intent(inout) :: ida_ll_g_st(:)
      integer, allocatable, intent(inout) :: jda_ll_g_st(:)
      integer, allocatable, intent(inout) :: kda_ll_g_st(:)
      integer, allocatable, intent(inout) :: lda_ll_g_st(:)
      integer, allocatable, intent(inout) :: mda_ll_g_st(:)
      integer, allocatable, intent(inout) :: nda_ll_g_st(:)
      real,    allocatable, intent(inout) :: s11_ll_g_st(:)
      real,    allocatable, intent(inout) :: s22_ll_g_st(:)
      real,    allocatable, intent(inout) :: s33_ll_g_st(:)
      real,    allocatable, intent(inout) :: s44_ll_g_st(:)
      real,    allocatable, intent(inout) :: d11_ll_g_st(:)
      real,    allocatable, intent(inout) :: d22_ll_g_st(:)
      real,    allocatable, intent(inout) :: d12_ll_g_st(:)
      real,    allocatable, intent(inout) :: ddd_ll_g_st(:)
      real,    allocatable, intent(inout) :: e00_ll_g_st(:)
      real,    allocatable, intent(inout) :: rep_ll_g_st(:)
      integer, allocatable, intent(inout) :: ida_ll_g_md(:)
      integer, allocatable, intent(inout) :: jda_ll_g_md(:)
      integer, allocatable, intent(inout) :: kda_ll_g_md(:)
      integer, allocatable, intent(inout) :: lda_ll_g_md(:)
      integer, allocatable, intent(inout) :: mda_ll_g_md(:)
      integer, allocatable, intent(inout) :: nda_ll_g_md(:)
      real,    allocatable, intent(inout) :: s11_ll_g_md(:)
      real,    allocatable, intent(inout) :: s22_ll_g_md(:)
      real,    allocatable, intent(inout) :: s33_ll_g_md(:)
      real,    allocatable, intent(inout) :: s44_ll_g_md(:)
      real,    allocatable, intent(inout) :: d11_ll_g_md(:)
      real,    allocatable, intent(inout) :: d22_ll_g_md(:)
      real,    allocatable, intent(inout) :: d12_ll_g_md(:)
      real,    allocatable, intent(inout) :: ddd_ll_g_md(:)
      real,    allocatable, intent(inout) :: e00_ll_g_md(:)
      real,    allocatable, intent(inout) :: rep_ll_g_md(:)
      integer, allocatable, intent(inout) :: ida_ll_g_rx(:)
      integer, allocatable, intent(inout) :: jda_ll_g_rx(:)
      integer, allocatable, intent(inout) :: kda_ll_g_rx(:)
      integer, allocatable, intent(inout) :: lda_ll_g_rx(:)
      integer, allocatable, intent(inout) :: mda_ll_g_rx(:)
      real,    allocatable, intent(inout) :: s11_ll_g_rx(:)
      real,    allocatable, intent(inout) :: s22_ll_g_rx(:)
      real,    allocatable, intent(inout) :: s33_ll_g_rx(:)
      real,    allocatable, intent(inout) :: s44_ll_g_rx(:)
      real,    allocatable, intent(inout) :: d11_ll_g_rx(:)
      real,    allocatable, intent(inout) :: d22_ll_g_rx(:)
      real,    allocatable, intent(inout) :: d12_ll_g_rx(:)
      real,    allocatable, intent(inout) :: ddd_ll_g_rx(:)
      real,    allocatable, intent(inout) :: e00_ll_g_rx(:)
      integer, allocatable, intent(inout) :: ida_ll_e_st(:)
      integer, allocatable, intent(inout) :: jda_ll_e_st(:)
      real,    allocatable, intent(inout) :: qqq_ll_e_st(:)
      integer, allocatable, intent(inout) :: ida_ll_e_md(:)
      integer, allocatable, intent(inout) :: jda_ll_e_md(:)
      real,    allocatable, intent(inout) :: qqq_ll_e_md(:)
      integer,              intent(inout) :: n_ll_r_st_old
      integer,              intent(inout) :: n_ll_u_st_old
      integer,              intent(inout) :: n_ll_a_st_old
      integer,              intent(inout) :: n_ll_a_md_old
      integer,              intent(inout) :: n_ll_b_st_old
      integer,              intent(inout) :: n_ll_b_md_old
      integer,              intent(inout) :: n_ll_v_st_old
      integer,              intent(inout) :: n_ll_v_md_old
      integer,              intent(inout) :: n_ll_g_st_old
      integer,              intent(inout) :: n_ll_g_md_old
      integer,              intent(inout) :: n_ll_e_st_old
      integer,              intent(inout) :: n_ll_e_md_old
      integer,              intent(inout) :: n_ll_g_rx_old
        end subroutine
      end interface

C                                                                                                     
C welcome to uiowa_bd                                                                                                     
C                                                                                                     
      real*16                  :: tecount       
      real*16                  :: tmcount      
      real*16                  :: ttcount     
      integer*8                :: ktime1 ! used to calc Time_total
      integer*8                :: ktime2 ! used to calc Time_total

C                                                                                                     
C parameter-related stuff                                                                             
C                                                                                                     
      character*80             :: string                      
      character*80             :: crapster                    
      character*80             :: fname1                     
      character*80             :: fname2                    
      character*80             :: fname3                   
      character*80             :: fname4                  
      character*80             :: fname2a                
      character*80             :: fname2b               
      character*80             :: growr_file   
      character*80             :: ewald_elec_grid_file   
      character*80             :: ewald_hydro_grid_file   
      character*3              :: atemp               
      character*3              :: btemp              
      character*3              :: ctemp             
      character*3              :: dtemp             
      character*3              :: xtemp             
      character*3              :: ytemp             
      character*5              :: etemp             
      character*7              :: ftemp             
      character*8              :: htemp             
      character*55             :: krop              
      character*1              :: krup              
      integer                  :: ierror(200)
      integer      :: omp_get_thread_num
      integer      :: omp_get_num_threads
      integer      :: omp_get_max_threads
      integer      :: omp_get_nested
      character*30             :: junk
      character*19             :: umbrella_output
      character*6              :: crap6            
      character*7              :: defunct          
      character*80             :: B22_restart_1
      character*80             :: B22_restart_2
      character*80             :: B22_density_pdb
      character*80             :: bond_func_file
      character*80             :: angl_func_file
      character*80             :: dihe_func_file
      character*80             :: nbnd_func_file
      character*80             :: rigid_domain_file
      character*80             :: monte_carlo_file
      character*80             :: NAM_mol_file
      character*80             :: mont_schedule_file
      character*9              :: temp_string   
      character*1610           :: long_string   
      character*180            :: long_string2   
      character*23             :: fname1h
      character*23             :: fname2h
      character*3              :: plop
      character*16 restart_file_name
      integer,allocatable      :: ifound(:,:) 

C variables need for HSL_MP54 cholesky routine

      integer,parameter :: wp = kind(0e0)
      integer,parameter :: long = selected_int_kind(18)
c     type(mp54_keep),allocatable :: keep(:)
c     type(mp54_control),allocatable :: ctrl(:)
      integer,allocatable :: keep(:)
      integer,allocatable :: ctrl(:)
      integer :: info,jid
      integer,allocatable :: nh3(:)
      integer,allocatable :: nb(:)
      integer,allocatable :: nrblk(:)
      integer(long) :: lrhs,ma
      integer(long),allocatable :: aoft(:)
      integer(long),allocatable :: la(:)
      integer(long) :: num_tmp_long
      integer,allocatable :: nthrd(:)
      real,allocatable    :: scale_nb(:)
      real(wp),dimension(:,:,:),allocatable :: amp
      integer,dimension(:,:),allocatable  :: imp
      integer,dimension(:,:),allocatable  :: jmp

C variables needed for MKL random numbers

      integer(kind=4) errcode,errcode_rex
      real(kind=4) a,sigma
      integer brng,method,seed,seed_rex
      TYPE (VSL_STREAM_STATE) :: stream
      TYPE (VSL_STREAM_STATE) :: stream_rex
      parameter (zero=0.0) ! eigenvalue stuff

C variables needed for reseting rotation matrices

      dimension dm(4),vm(4,4),cmov(3),cfix(3) 
      dimension tjac(3,3),qjac(4,4),qcac(4,4) 
      dimension sumf(3),summ(3),stmp(3)      
      dimension xunit(3,3)                  
      data xunit/1.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,1.000/

      real, allocatable :: amkl_new(:)

      call getarg(1,junk)
      read(junk,*)iseed_new   ! random seed

      i_do_rods=.false.

C set some constants/variables                                                                        
 
      r_kbt=1.987e-3                                                    
      rootpi=1.77245385
      rootpi_inv=1.0/rootpi

C set some MKL related stuff

c     brng=VSL_BRNG_MT19937
      brng=VSL_BRNG_MT2203  ! this recently mentioned

      method=0
      a_for_mkl=0.0
      sigma=1.0

      do iii=1,200      
        ierror(iii)=0     
      enddo    

C read the input file                                                                             

      read(5,*)                                                        
      read(5,*)                                                        
      read(5,*)teprint,ttprint,tmprint,                         
     &         num_lst_stp,num_fmd_stp,num_hyd_stp
      num_bal_stp=-1 ! set for compatibility with old code

      read(5,*)
      read(5,*)
      read(5,*)num_threads,bond_dev_quit,atemp,btemp

      write(*,*)
      write(*,*)'NOTE NOTE NOTE NOTE NOTE NOTE NOTE'
      write(*,*)'this version of the code will update the box'
      write(*,*)'dimensions at every step - this has been done'
      write(*,*)'as a simple fix to allow us to skip nonbonded'
      write(*,*)'interactions and updates entirely...'
      write(*,*)

      if(mod(num_hyd_stp,num_lst_stp).ne.0) then
        write(*,*)
        write(*,*)'sorry - the number of steps between updates'
        write(*,*)'of the diffusion tensor (num_hyd_stp) should be'
        write(*,*)'a multiple of the number of steps between'
        write(*,*)'updates of the nonbonded list (num_lst_stp)'
        write(*,*)
        write(*,*)'NOTE NOTE NOTE NOTE we used to quit here :('
        write(*,*)
c       stop
      endif
     
      if(atemp.eq.'yes') then
        i_do_lincs=.true.
        write(*,*)
        write(*,*)'we will be using lincs to constrain bonds'
        write(*,*)'this can be dangerous when using HI'
        write(*,*)
      else
        i_do_lincs=.false.
        write(*,*)
        write(*,*)'we will not be using lincs to constrain bonds'
        write(*,*)'note that it is possible we will change our minds'
        write(*,*)'on this if we are using mixed_hydrodynamics and if '
        write(*,*)'we find noHI atoms + lincs in mixed_schedule_file'
        write(*,*)
      endif

      if(btemp.eq.'yes') then
        i_do_YHL_nonbondeds=.true.
        write(*,*)
        write(*,*)"we will use Newton's second law for nonbondeds"
        write(*,*)"this means that we exploit the fact that Fij=-Fji"
        write(*,*)"this generally works best for low thread counts"
        write(*,*)"- should definitely be the choice for one thread"
        write(*,*)
      else
        i_do_YHL_nonbondeds=.false.
        write(*,*)
        write(*,*)"we will not use Newton's second law for nonbondeds"
        write(*,*)"this is generally the best for high thread counts"
        write(*,*)
      endif

C initialize openmp 

      call omp_set_num_threads(num_threads)
      write(*,*)'setting number of threads to ',num_threads,
     &          omp_get_max_threads()

C MUST set mkl_set_dynamic to zero if we want MKL routines to
C automatically use multiple threads when called by a single thread in
C openmp parallel region...

      call mkl_set_dynamic(0)
c     call mkl_set_num_threads(num_threads)
      call omp_set_nested(2)
      write(*,*)'check openmp nested? ',omp_get_nested()

C for now,read in number of total molecule typs - avoid later                                       
C note that we read in num_t_ms here,but this gets recomputed below                                  

      write(*,*)
      write(*,*)'about to read #moltypes #mols and other stuff'
      write(*,*)'note that if you crash here you may need to add'
      write(*,*)'the term go_eps_low which is the epsilon value'
      write(*,*)'used to determine which go contacts to track for'
      write(*,*)'calculating Q: if eps>go_eps_low then we use it '
      write(*,*)'in the calculation of Q, otherwise we use it only'
      write(*,*)'for calculating the forces on go contacts'
      write(*,*)
      write(*,*)'to use *all* go contacts in the calculation of Q'
      write(*,*)'simply set go_eps_low to a negative number, e.g. -1.0'
      write(*,*)
      read(5,*) 
      read(5,*)
      read(5,*)num_f_ts,num_f_ms,atemp,
     &         Q_des,mol_Q1,mol_Q2,go_eps_low

C 01-02-08 - include a debug logical

      i_debug=.false.
      if(atemp.eq.'yes') i_debug=.true.

      num_t_ts=num_f_ts          ! totl # mol types                             
      num_t_ms=num_f_ms          ! totl # molecules                            
      num_t_ms3=num_t_ms*3
      allocate(itmp_mol_array(1:num_t_ms))
      allocate(c_o_m_x(1:num_t_ms))
      allocate(c_o_m_y(1:num_t_ms))
      allocate(c_o_m_z(1:num_t_ms))
      allocate(i_f_tp(1:num_t_ms))                  
      allocate(i_r_tp(1:num_t_ms))                   
      allocate(i_t_tp(1:num_t_ms))                  
      allocate(my_num_ths_t_typ(1:num_t_ms))
      allocate(ityp(1:num_t_ms))                      
      allocate(dtrans(1:3,1:num_t_ms))                 
      allocate(drot(1:3,1:num_t_ms))                   
      allocate(ene(1:num_t_ms))                        
      allocate(num_ths_f_typ(1:num_f_ts))   
      allocate(num_ths_t_typ(1:num_t_ts))   
      allocate(qef_file(1:num_t_ts))                    
      allocate(dyn_file(1:num_t_ts))                    
      allocate(bnd_file(1:num_t_ts))                    
      allocate(phi_file_A(1:num_t_ts))                  
      allocate(phi_file_B(1:num_t_ts))                  
      allocate(my_num_t_typ(1:num_t_ts))       
      allocate(my_typ_t_typ(1:num_t_ts))
      allocate(ene_m(1:num_t_ms))
      allocate(ene_m_old(1:num_t_ms))
      allocate(ene_m_new(1:num_t_ms))
      allocate(g_epsilon_tmp(1:num_threads*16))
      allocate(num_m_pairs_omp(1:num_threads))
      allocate(num_a_pairs_omp(1:num_threads))
      allocate(num_qqq_pairs_omp(1:num_threads))
      allocate(num_uniq_m_pairs_omp(1:num_threads))
      do iii=1,num_t_ms
        allocate(ene(iii)%jdm(0))
        allocate(ene(iii)%tot(0))
        allocate(ene(iii)%e(0))
        allocate(ene(iii)%v(0))
        allocate(ene(iii)%g(0))
        allocate(ene(iii)%dis(0))
        ene(iii)%num=0
      enddo

C teprint,ttprint,tmprint are for output of energies,c_m_xs(+rdfs)                                 
C and movies respectively                                                                             

      read(5,*)
      read(5,*)
      read(5,*)xmin,xmax,ymin,ymax,zmin,zmax,i_pbc
      read(5,*)
      read(5,*)
      read(5,*)atemp,btemp,ctemp
      i_do_pka_calcs=.false.

      i_look_for_crashes=.false.
      if(atemp.eq.'yes') i_look_for_crashes=.true.
      periodic_bonds_okay=.false.
      if(btemp.eq.'yes') periodic_bonds_okay=.true.
      i_use_high_mem=.false.
      if(ctemp.eq.'yes') i_use_high_mem=.true.

      xlen=xmax-xmin                                                    
      ylen=ymax-ymin                                                    
      zlen=zmax-zmin                                                    
      xinv=1.000/xlen                                                   
      yinv=1.000/ylen                                                   
      zinv=1.000/zlen                                                  
      write(*,*)'dimensions of simulation system are ',xlen,ylen,zlen 

C 31-07-12 additions
C i_append_movie is used to determine whether we want to add to movies
C rather than start again at zero - it looks for the movie with the
C higher number and either starts there or adds one depending on the
C value of i_append_movie
C i_append_movie = 1 --> overwrite the last movie
C                = 2 --> add to the list
C note that the starting movieframenumber can also be read from the 
C restart.file if provided there

C i_limit_verbosity limits how much is written out when we use the
C arbitrary functions of our force field if 'yes' then don't write out
C much

      read(5,*)
      read(5,*)
      read(5,*)i_append_movie,ctemp,
     &         i_use_v_typ,btemp             

      replica_exchange=.false.
      i_limit_verbosity=.false.
      if(ctemp.eq.'yes') i_limit_verbosity=.true.
      movieframenumber=0                                         
      if(i_append_movie.gt.0) then
        call system('ls -lsr MOVIE/mo*pdb* | head -1 > temp_movie')
        open(unit=18,file='temp_movie',status='unknown')
        read(18,11711,end=11713)long_string
11711   format(a120)
        do i=1,120
          if(long_string(i:i+4).eq.'movie') then
            temp_string=long_string(i+6:i+14)
            read(temp_string,11712) movieframenumber
11712       format(i9)
            goto 11713
          endif
        enddo
11713   continue
        close(18)
        write(*,*)'the most recent movie is # ',movieframenumber
        if(i_append_movie.eq.2) movieframenumber=movieframenumber+1
        write(*,*)'i will start writing at movie # ',movieframenumber
      endif

      uniform_moves=.false.
      steepest_descent=.false.
      if(btemp.eq.'yes') steepest_descent=.true.

      read(5,*)
      read(5,*)
      read(5,*)temperature,ionic_strength,r_ion,dielectric,
     &         viscosity,     ! 051-05-10 viscosity :)
     &         r_f_st         ! 22-03-08 this is the force constant for short-range v

      r_f_hf=r_f_st*0.5
      beta=1.000/(r_kbt*temperature)    
      dielectric_inv=1.0/dielectric
      implicit_ions=.false.

      i_do_pH=.false.
      write(*,*)'protonation states will not change'

C 64-11-11 change here - previously we only defined kappa etc when
C ionic_strength.gt.0 - i don't see any reason *not* to define them even
C when ionic_strength=0 - it means the B22/NAM calculations go easier

C 32-07-12 - since we're using the yukawa debye-huckel code now, which
C assumes r_ion = 0, we should probably enforce that here...

      if(r_ion.gt.0.0) then
        write(*,*)'sorry but this is a bit shaky to use',
     &            ', better off setting r_ion=0.0'
        write(*,*)'(the issue is the treecode elec routine)'
        write(*,*)
        write(*,*)'quitting :('
        write(*,*)
        stop
      endif
      if(ionic_strength.gt.0.0) then
        implicit_ions=.true.
      endif
      kappa=0.010375*sqrt(ionic_strength)                     
      kappa=kappa*sqrt(298.150)/sqrt(temperature)      
      kappa=kappa*sqrt(77.96521)/sqrt(dielectric)      
      one_over_one_plus_kappa_rion=1.0/(1.0+kappa*r_ion)
      one_third=1.0/3.0
      eight_over_three=8.0/3.0

C 2023 v1.1 read fold_mode as integer (used to be scale_nb)

      read(5,*)
      read(5,*)
      read(5,*)parameter_file,atemp,wrap_molecules,
     &         etemp,fold_mode,htemp

      i_use_no_elec=.false.
      if(atemp.eq.'yes') i_use_no_elec=.true.

      write(*,*)
      if(wrap_molecules.eq.0) then
        write(*,*)'will not wrap coordinates in movie files'
      elseif(wrap_molecules.eq.1) then
        write(*,*)'will wrap coordinates in movie files by atom'
      elseif(wrap_molecules.eq.2) then
        write(*,*)'will wrap coordinates in movie files by molecule'
      endif
      write(*,*)

C 2023 v1.1 write whether we are in folding or unfolding mode

      if(fold_mode.eq.1) then
        write(*,*)
        write(*,*)'this is nominally a folding simulation'
        write(*,*)
      elseif(fold_mode.eq.-1) then
        write(*,*)
        write(*,*)'this is nominally an unfolding simulation'
        write(*,*)
      else
        write(*,*)
        write(*,*)'fold_mode is set to ',fold_mode
        write(*,*)'this is a fatal error: must be +1 or -1'
        write(*,*)'quitting :('
        write(*,*)
        stop
      endif

C we default to saying i_use_hydro=.true.

      brownian=.false.
      langevin=.false.
      i_use_hydro=.true.
      write(*,*)
      write(*,*)'we default to setting i_use_hydro = .true.'
      write(*,*)

C set num_reps to 1 as default

      num_reps=1

      no_hydrodynamics=.false.
      full_hydrodynamics=.false.
      rpy_hydrodynamics=.false.
      oarpy_hydrodynamics=.false.
      nodiv_hydrodynamics=.false.
      mixed_hydrodynamics=.false.

      if(etemp(1:4).eq.'none') then
        no_hydrodynamics=.true.
      elseif(etemp(1:3).eq.'RPY') then
        rpy_hydrodynamics=.true.
      elseif(etemp(1:5).eq.'OARPY') then
        oarpy_hydrodynamics=.true.
      elseif(etemp(1:5).eq.'NODIV') then
        nodiv_hydrodynamics=.true.
      else
        write(*,*)'sorry, only options are: none/RPY/OARPY/CUT'
        write(*,*)'all other options now defunct'
        write(*,*)
        write(*,*)'quitting :('
        write(*,*)
      endif

C now set full_hydrodynamics to true if *either* RPY or OARPY...

      if(rpy_hydrodynamics)   full_hydrodynamics=.true.
      if(oarpy_hydrodynamics) full_hydrodynamics=.true.
      if(nodiv_hydrodynamics) full_hydrodynamics=.true.

C note that the following were formerly in the mixed_hydrodynamics
C section but we need them to be done for all kinds of sims

C allocate memory for storing - for each thread - which replica it works
C on - note that no thread can work on more than one replica...

C we also have mythrds_rep_curr - this starts as the same thing as
C mythrds_rep_orig but is used to keep track of what the current replica's
C Hamiltonian really is (i.e. what its current value of rep_go_scale is)
C - it's not used for calculating any forces, moves etc - it's just used
C for accounting work on the Q values

      allocate(mythrds_rep_orig(1:num_threads))
      allocate(mythrds_rep_curr(1:num_threads))
      mythrds_rep_orig=1 ! defaults make it work w/o replica_exchange
      mythrds_rep_curr=1

      if(i_do_lincs) then
        if(full_hydrodynamics) then
          write(*,*)
          write(*,*)'you cannot use lincs with full_hydrodynamics'
          write(*,*)
          write(*,*)'quitting :('
          write(*,*)
          stop
        endif
      endif

C note that multi-scale hydrodynamics is currently only implemented for
C cutoff and Ewald electrostatics - not for Ewald grid electrostatics

      integer_arrays=.false.
      hydro_cutoff=.false.
      dpd=.false.
      i_grow_rs=.false.

C keep scaling factor for 'nb' entry when using Cholesky
C 2023 v1.1 just set to 1.0 here - we will use gtemp for fold_mode

      scale_nb_temp=1.0

C read langevin from last column (htemp)

      if(htemp(1:8).eq.'brownian') brownian=.true.
      if(brownian) write(*,*)'using brownian dynamics'
      if(htemp(1:8).eq.'langevin') langevin=.true.
      if(langevin) write(*,*)'using langevin dynamics'
      if(.not.brownian.and..not.langevin) then
        write(*,*)
        write(*,*)'need to specify a valid dynamics algorithm  - your'
        write(*,*)'choices are brownian or langevin'
        write(*,*)
        write(*,*)'quitting :('
        write(*,*)
        stop
      endif

      read(5,*)
      read(5,*)
      write(*,*)
      write(*,*)'about to read go-related stuff - if we crash here you '
      write(*,*)'may need to add an integer for go_primacy as the '
      write(*,*)'3rd argument on the line of the input file: '
      write(*,*)
      write(*,*)'if = 1 then go contacts use only go potential'
      write(*,*)'if = 2 then go contacts use go potential + elec'
      write(*,*)'if = 3 then go contacts use go potential + elec + vdw'
      write(*,*)
      write(*,*)'note that 1 & 2 will use standard 10-12 LJ potential'
      write(*,*)'3 will use a special capped harmonic well term'
      write(*,*)
      write(*,*)'if in doubt, set it to 2...'
      write(*,*)
      write(*,*)'you probably also need to set i_skip_intermol_go'
      write(*,*)'as the 4th argument on the line - if set to yes'
      write(*,*)'this means that all intermoleculare go contacts'
      write(*,*)'are ignored: useful if you want them modeled using'
      write(*,*)'arbitrary potentials etc'
      write(*,*)

      read(5,*)goparam_file,atemp,go_primacy,dtemp

      if(atemp.eq.'yes') then
        i_use_go_pairs=.true.
      else
        i_use_go_pairs=.false.
      endif

      if(dtemp.eq.'yes') then
        i_skip_intermol_go=.true.
      else
        i_skip_intermol_go=.false.
      endif

      i_use_exclusive_go=.false.
      i_compare_go_with_others=.false.
      i_have_go_modes=.false.
      i_do_reactions=.false.

      read(5,*)
      read(5,*)
      read(5,*)nomove_file,atemp                                        
      i_dont_move_some=.false.
      if(atemp.eq.'yes') i_dont_move_some=.true.
      if(i_dont_move_some) then
        open(unit=96,file=nomove_file,status='unknown')
        iii=0
635     format(1x)
636     read(96,635,end=637)
        iii=iii+1
        goto 636
637     continue
        num_moves_omitted=iii
        allocate(i_omit_move_mtyp(1:num_moves_omitted))                 
        allocate(i_omit_move_a(1:num_moves_omitted))                    
        rewind(96)
        iii=0
606     read(96,*,end=677)i1,i2
        iii=iii+1 
        i_omit_move_mtyp(iii)=i1
        i_omit_move_a(iii)=i2
        goto 606
677     continue
        close(96)                                                       
        write(*,*)
        write(*,*)'# of static atoms = ',num_moves_omitted
        write(*,*)
      endif                                                             
      read(5,*)                                                        
      read(5,*)                                                        
C                                                                                                     
C set all molecule numbers and typs to zero and accumulate                                                
C                                                                                                     
      num_t_ms=0                                                      
      num_f_ms=0                                                      
      num_r_ms=0                                                      
      num_t_tmp=0                                                      
      num_f_tmp=0                                                      
      num_r_tmp=0                                                      
      icrashed=0
7788  continue             

      read(5,*)fname1
      read(5,*)fname2A
      read(5,*)nofthistyp
      read(5,*)

C read diffusion coefficients                                                                         

      num_t_tmp=num_t_tmp+1                                       
      qef_file(num_t_tmp)=fname1                                      
      num_f_tmp=num_f_tmp+1                                      
      my_typ_t_typ(num_t_tmp)=1
      my_num_t_typ(num_t_tmp)=num_f_tmp
      bnd_file(num_t_tmp)=fname2A                                   
      num_ths_t_typ(num_t_tmp)=nofthistyp
      num_ths_f_typ(num_f_tmp)=nofthistyp
      
C note that each of the following variables runs over num_t_ms...                                   
                                                                        
      do ijk=1,nofthistyp                                             
        num_t_ms=num_t_ms+1                                  
        i_t_tp(num_t_ms)=num_t_tmp                       
        ityp(num_t_ms)=1                                          
        num_f_ms=num_f_ms+1                                   
        i_f_tp(num_f_ms)=num_f_tmp                        
        flex=.true.                                                 

        dtrans(1,num_t_ms)   =ftmp1                              
        dtrans(2,num_t_ms)   =ftmp1                               
        dtrans(3,num_t_ms)   =ftmp1                              
        drot(1,num_t_ms)     =ftmp2                              
        drot(2,num_t_ms)     =ftmp2                              
        drot(3,num_t_ms)     =ftmp2                              

      enddo 
      if(num_t_tmp.eq.num_t_ts) goto 7799                  
      goto 7788                                                         
7799  do i=1,num_t_ts
        write(*,912)i,num_ths_t_typ(i)
912     format('for total typ ',i8,' there are ',i8,' copies')
      enddo
      read(5,*)
      read(5,*) time_step,totsimtime,
     &          cut_v_st,cut_v_md,
     &          cut_g_st,cut_g_md,
     &          cut_e_st,cut_e_md,f_f_cell_size

      neprint=teprint/time_step
      ntprint=ttprint/time_step
      nmprint=tmprint/time_step

      write(*,*)'NOTE - we will use time_step = ',time_step

C 32-07-12 - if we set tmprint to a negative number we never write movie

      if(tmprint.lt.0.0) nmprint=-1
      ntimetotal=totsimtime/time_step
      write(*,*)'will write out energy every ',neprint,' steps'
      write(*,*)'will write out xtc    every ',ntprint,' steps'
      write(*,*)'will write out movie  every ',nmprint,' steps'
      write(*,*)'will carry out a total of   ',ntimetotal,' steps'
      cut_v_st2=cut_v_st**2
      cut_v_md2=cut_v_md**2
      cut_g_st2=cut_g_st**2
      cut_g_md2=cut_g_md**2
      cut_e_st2=cut_e_st**2
      cut_e_md2=cut_e_md**2
 501  format(/,a)

C read position restraint information

      write(*,*)'about to read lines about position restraints'
      read(5,*)                                                        
      read(5,*)                                                        
      read(5,*) pos_restraint_file,atemp,btemp                    
      i_do_pos_restraint=.false.
      if(atemp.eq.'yes') i_do_pos_restraint=.true.
      mission_creep=.false.
      if(btemp.eq.'yes') mission_creep=.true.
      if (mission_creep) then
        write(*,*)'note that mission_creep is on!'
        write(*,*)'position restraints will drift every num_hyd_cnt ',
     &            'steps'
      endif

C read sphere/cylinder confinement information

      write(*,*)
      write(*,*)'about to read lines about sphere/cylinder containers'
      write(*,*)'r_size is radius of sphere/spherical caps'
      write(*,*)'l_size is the length of the cylinder'
      write(*,*)
      read(5,*)                                                        
      read(5,*)                                                        
      read(5,*) r_size,l_size,r_size_fac,n_size,f_size
      sphere=.false.
      cylinder=.false.
c     if(r_size.eq.0.0.and.l_size.eq.0.0) shrinking_box=.true.
      if(n_size.gt.0) shrinking_box=.true.
      if(r_size.gt.0.0.and.l_size.lt.0.0) sphere=.true.
      if(r_size.gt.0.0.and.l_size.gt.0.0) cylinder=.true.
      if(sphere) write(*,*)'will impose sphere conditions'
      if(cylinder) write(*,*)'will impose cylinder conditions'
      if(shrinking_box) then
        write(*,*)'will impose shrinking box conditions'
        xmax_final=xmax
        ymax_final=ymax
        zmax_final=zmax
        xmin_final=xmin
        ymin_final=ymin
        zmin_final=zmin
      endif
      r_size2=r_size**2
      l_size2=l_size**2
      l_sizeh=l_size*0.5
      s_size=r_size
      m_size=l_size

      cutoff_elec=.true.

      ewald_hydro=.false.
      ewald_hydro_grid=.false.

C --------------------------------------------------------
C read whether we're doing using Fixman method

      write(*,*)'about to read lines about fixman random terms' 
      read(5,*)                                                        
      read(5,*)                                                        
      read(5,*) atemp,fixman_tol,fixman_order,btemp,               
     &          lmin_read_in,lmax_read_in
      cholesky=.true. ! default to assuming this is on
      fixman=.false.
      if(atemp.eq.'yes') then
        fixman=.true.
        cholesky=.false.
      endif
      fixman_override=.false.
      if(btemp.eq.'yes') fixman_override=.true.
      if(fixman_override.and..not.fixman) then
        write(*,*)'must use fixman with fixman_override'
        write(*,*)'quitting :('
        stop
      endif
      if(no_hydrodynamics) cholesky=.false.
      if(cholesky) write(*,*)'using cholesky method'
      if(fixman)   write(*,*)'using fixman method'  
C --------------------------------------------------------
C read whether we're using a Treecode method   

C here's where your treecode-related stuff is read from file

      write(*,*)'about to read lines about treecode electrostatics'
      read(5,*)                                                        
      read(5,*)                                                        
      read(5,*) atemp,treecode_theta,treecode_order,
     &                treecode_shrink,treecode_maxatm
      treecode_elec=.false.

      if(atemp.eq.'yes') treecode_elec=.true.

      if(treecode_elec.and.i_do_YHL_nonbondeds.and.
     &   num_threads.gt.1) then
        write(*,*)
        write(*,*)'if using more than one thread: '
        write(*,*)'there are unresolved race conditions with using '
        write(*,*)"treecode and newton's law for nonbondeds"
        write(*,*)'set i_do_YHL_nonbondeds to "no" instead'
        write(*,*)'or do not use treecode...'
        write(*,*)
        write(*,*)'quitting :('
        write(*,*)
        stop
      endif

      prinfo=0
      if(treecode_elec.and.i_debug) prinfo=1

      i_do_growth=.false.
      i_do_slide =.false.

C --------------------------------------------------------
C look for walls

      write(*,*)'about to read lines about planar walls'     
      read(5,*)                                                        
      read(5,*)                                                        
      read(5,*)atemp,num_walls,wall_file
      walls=.false.       
      if(atemp.eq.'yes') then
        walls=.true.
        write(*,*)
        write(*,*)'will write out pressure on each wall '
        write(*,*)'to osmotic_pressure.txt'
        write(*,*)
        allocate(crd_wall(1:num_walls)) ! +/-x, y or z if +/-1,2 or 3
        allocate(pos_wall(1:num_walls)) ! position of wall
        allocate(fct_wall(1:num_walls)) ! force constant
        do i=1,num_walls
          read(5,*)crd_wall(i),pos_wall(i),fct_wall(i)
        enddo
      endif                              
 
C allocate space for osmotic pressure - note that we allocate for 7
C dimensions +x,+y,+z,-x,-y,-z,radial (for sphere and cylinder)

      if(walls.or.sphere.or.cylinder) then
        open(unit=25,file='osmotic_pressure.txt',status='unknown')
        allocate(osm_cur(1:7)) ! current osmotic pressure
        allocate(osm_ave(1:7)) ! cumulative osmotic pressure
        allocate(osm_var(1:7)) ! variance osmotic pressure
        allocate(surface(1:7)) ! surface aread of wall    
        num_osm=7              ! 7 poss dimensions for osm pressure
        osm_ave=0.0
        osm_var=0.0
        surface=0.0
        if(sphere) then
          surface(7)=4.0*pi*s_size*s_size
          write(*,*)
          write(*,*)'surface area of sphere = ',surface
          write(*,*)
        elseif(cylinder) then
          surface(7)=2.0*pi*s_size*m_size + 4.0*pi*s_size2
          write(*,*)
          write(*,*)'surface area of cylinder = ',surface
          write(*,*)
        else
          surface(1)=ylen*zlen
          surface(2)=xlen*zlen
          surface(3)=xlen*ylen
          surface(4)=ylen*zlen
          surface(5)=xlen*zlen
          surface(6)=xlen*ylen
        endif
      else
        num_osm=0
      endif                              

      go_spline=.false.       
      afm=.false.
      umbrella=.false.
      B22=.false.
      i_do_protease=.false.
      i_read_bond_functions=.false.
      i_read_nonbond_functions=.false.
      i_write_bond_histogram=.false.
      i_write_nonbond_histogram=.false.
      i_write_user_histogram=.false.
      i_read_ref_hist=.false.
      arbitrary_intra=.false.
      i_rigid_energy=.false.
      i_do_NAM=.false.
      i_user_energy=.false.

C redetermine the types of the molecules - whether flexible or rigid
C note that if we set a molecule to have a non-zero flex typ (i_f_tp)
C then we also set to have a zero rigid typ (i_r_tp) and vice versa

      num_f_ts=0
      num_r_ts=0
      do iii=1,num_t_ms
        itmp=i_t_tp(iii)
        do jjj=iii-1,1,-1 ! check backwards - it's quicker
          if(itmp.eq.i_t_tp(jjj)) then
            if(ityp(iii).eq.1) then
              i_f_tp(iii)=i_f_tp(jjj)
              i_r_tp(iii)=i_r_tp(jjj)
            else
              i_r_tp(iii)=i_r_tp(jjj)
              i_f_tp(iii)=i_f_tp(jjj)
            endif
            goto 82128
          endif
        enddo
        if(ityp(iii).eq.1) then
          num_f_ts=num_f_ts+1
          i_f_tp(iii)=num_f_ts
          i_r_tp(iii)=0
        else
         num_r_ts=num_r_ts+1
          i_r_tp(iii)=num_r_ts
          i_f_tp(iii)=0
        endif
82128   continue
        if(i_debug) then
          write(*,8890)iii,i_t_tp(iii),i_f_tp(iii),i_r_tp(iii)
8890      format(' molecule ',i10,' is of t typ ',i5,
     &           ' of f typ ',i5,' of r typ ',i5)
        endif
      enddo
 
      write(*,*)
      write(*,*)'after double-checking, there are ',
     &          num_t_ts, ' total types '
      write(*,*)'after double-checking, there are ',
     &          num_f_ts, ' flexible types '
      write(*,*)'after double-checking, there are ',
     &          num_r_ts, ' rigid types '
      write(*,*)

C allocate arrays for generalized diffusion tensors

      allocate(D_11(1:num_r_ts))
      allocate(D_22(1:num_r_ts))
      allocate(D_33(1:num_r_ts))
      allocate(D_44(1:num_r_ts))
      allocate(D_55(1:num_r_ts))
      allocate(D_66(1:num_r_ts))
      allocate(E_11(1:num_r_ts))
      allocate(E_22(1:num_r_ts))
      allocate(E_33(1:num_r_ts))
      allocate(E_44(1:num_r_ts))
      allocate(E_55(1:num_r_ts))
      allocate(E_66(1:num_r_ts))

C allocate arrays for generating r coords on the fly

      allocate(f_f_list(1:1),stat=ierror(1))                         
      num_grow_r=0                                                  
                                                                                                     
C allocate arrays for growing/scaling of coordinates                                                         
                                                                                                     
      allocate(idm_grow_r(1:num_r_ms),stat=ierror(1))      
      allocate(scbeg_grow_r(1:num_r_ms),stat=ierror(2))      
      allocate(scend_grow_r(1:num_r_ms),stat=ierror(3))      
      allocate(nstep_grow_r(1:num_r_ms),stat=ierror(4))      
                                                                                                     
C 01-04-13 allocate some arrays at this point                                                           

      allocate(id_beg(1:num_t_ms),stat=ierror(1))                   
      allocate(id_end(1:num_t_ms),stat=ierror(2))                  
      allocate(id_beg_bonds(1:num_f_ts),stat=ierror(3))      
      allocate(id_beg_angls(1:num_f_ts),stat=ierror(4))    
      allocate(id_beg_dihes(1:num_f_ts),stat=ierror(5))

C allocate atom info for flexible and rigid molecules
C note that num_atms_pt will contain all atoms for flexible mols
C and 3 axis atoms for rigid molecules - importantly, we will
C all dynamic atoms will be held in num_dyns_pt

      allocate(num_atms_pt(1:num_t_ts),stat=ierror(36))     
      allocate(num_dyns_pt(1:num_t_ts),stat=ierror(36))     

      allocate(num_r_qs_pt(1:num_r_ts),stat=ierror(35))
      allocate(num_r_ds_pt(1:num_r_ts),stat=ierror(36))

C AHE TEMP TEMP TEMP we need a proper zero'ing of all these arrays

      do itp=1,num_t_ts
        num_atms_pt(itp)=0
        num_dyns_pt(itp)=0
      enddo

      allocate(num_bonds(1:num_t_ts),stat=ierror(6))                 
      allocate(num_angls(1:num_t_ts),stat=ierror(7))                
      allocate(num_dihes(1:num_t_ts),stat=ierror(8))             
      allocate(i_done_bonds(1:num_t_ts),stat=ierror(114))   
      allocate(i_done_angls(1:num_t_ts),stat=ierror(115))     
      allocate(i_done_dihes(1:num_t_ts),stat=ierror(116))  

C n_typ used to be name_f
C c_typ used to be c_f_tp
C b_typ used to be bnd_exc

      allocate(n_typ(1:num_t_ts),stat=ierror(81))              
      allocate(c_typ(1:num_t_ts),stat=ierror(132))      
      allocate(b_typ(1:num_t_ts),stat=ierror(133))                 

      allocate(c_m_x(1:num_t_ms),stat=ierror(83))       
      allocate(c_m_y(1:num_t_ms),stat=ierror(83))       
      allocate(c_m_z(1:num_t_ms),stat=ierror(83))       
      allocate(ene_pos(1:num_f_ms),stat=ierror(124))               
      allocate(ene_bond(1:num_f_ms),stat=ierror(124))               
      allocate(ene_angl(1:num_f_ms),stat=ierror(125))              
      allocate(ene_dihe(1:num_f_ms),stat=ierror(126))       
      allocate(i_clash(1:num_t_ms),stat=ierror(128))          ! should eventually be dimd to 1->num_r_ms
      allocate(i_need_recompute(1:num_t_ms),stat=ierror(129)) ! should eventually be dimd to 1->num_r_ms
      allocate(c_f_m(1:num_f_ms),stat=ierror(134))      
      allocate(idm_tmp_nbnd(1:num_t_ms),stat=ierror(137))      

C now read the charge & dynamics files to get the number of as                                

      do 666 iii=1,num_t_ms                                     

C first ask if we've seen this typ of molecule before                                                

        do mmm=iii-1,1,-1                                             
          if(i_t_tp(iii).eq.i_t_tp(mmm)) goto 666             
        enddo                                                         

C if not go ahead and get its atom numbers etc

        ntmp=0                                                          
600     format(1x)                          
                            
        open(32,file=qef_file(i_t_tp(iii)),form='formatted',
     &       status='old')                                               
601     read(32,700,end=602)krop,krup                                         
700     format(a55,a)

C make arrays storing name/type information etc - c_typ and n_typ
C note that for flexible mols we look for 'Q' atoms; for rigid mols we
C look for both 'Q' and 'D' atoms - it is intended that there will only
C be 3 of the former however...

        if(krop(1:4).eq.'ATOM'.and.krup.eq.'Q') ntmp=ntmp+1

        if(ityp(iii).eq.0) then
          if(krop(1:4).eq.'ATOM'.and.krup.eq.'D') ntmp=ntmp+1
        endif

        goto 601                                                         
602     close(32)                                                        

        write(*,*)
        write(*,*)'allocating memory for mol type ',i_t_tp(iii),ntmp
        write(*,*)
      
        if(ityp(iii).eq.1) then
          n_typ(i_t_tp(iii))%num=ntmp                          
          c_typ(i_t_tp(iii))%num=ntmp
          b_typ(i_t_tp(iii))%num=ntmp 
        else
          n_typ(i_t_tp(iii))%num=ntmp  ! poss many atoms to display                      
          c_typ(i_t_tp(iii))%num=ntmp  ! (but only 3 active 'atoms')
          b_typ(i_t_tp(iii))%num=0     ! zero bonds
        endif                                                              

        allocate(n_typ(i_t_tp(iii))%ch1(n_typ(i_t_tp(iii))%num))
        allocate(n_typ(i_t_tp(iii))%ch2(n_typ(i_t_tp(iii))%num))
        allocate(n_typ(i_t_tp(iii))%ch3(n_typ(i_t_tp(iii))%num))
        allocate(n_typ(i_t_tp(iii))%ch4(n_typ(i_t_tp(iii))%num))

        allocate(c_typ(i_t_tp(iii))%a(c_typ(i_t_tp(iii))%num))
        allocate(c_typ(i_t_tp(iii))%x(c_typ(i_t_tp(iii))%num))
        allocate(c_typ(i_t_tp(iii))%y(c_typ(i_t_tp(iii))%num))
        allocate(c_typ(i_t_tp(iii))%z(c_typ(i_t_tp(iii))%num))
        allocate(c_typ(i_t_tp(iii))%q(c_typ(i_t_tp(iii))%num))
        allocate(c_typ(i_t_tp(iii))%r(c_typ(i_t_tp(iii))%num))
        allocate(c_typ(i_t_tp(iii))%s(c_typ(i_t_tp(iii))%num)) ! 2019

        allocate(b_typ(i_t_tp(iii))%id1(b_typ(i_t_tp(iii))%num))

666   enddo                                                               
C                                                                                                     
C make molecule-dependent coords part of allocatable arrays...                                        

      do 676 iii=1,num_t_ms                                         
        if(ityp(iii).eq.1) then                                         
          c_f_m(iii)%num=c_typ(i_t_tp(iii))%num   
          allocate(c_f_m(iii)%x(c_f_m(iii)%num))
          allocate(c_f_m(iii)%y(c_f_m(iii)%num))
          allocate(c_f_m(iii)%z(c_f_m(iii)%num))
          allocate(c_f_m(iii)%free(c_f_m(iii)%num))
          allocate(c_f_m(iii)%harm(c_f_m(iii)%num)) 
          allocate(c_f_m(iii)%harm_gone(c_f_m(iii)%num))
          allocate(c_f_m(iii)%just_harm(c_f_m(iii)%num))
          allocate(c_f_m(iii)%just_free(c_f_m(iii)%num)) 
          allocate(c_f_m(iii)%near_free(c_f_m(iii)%num)) 
        endif                                                 
676   enddo                                                  
                                                                                                      
C now read the charges for each molecule                                                    

      write(*,*)'about to make call to make_charge_parameters' 
      call make_charge_parameters                              
      write(*,*)'returned from call to make_charge_parameters' 

C figure out cell sizes for f-f interactions                                     

      num_x_ff=int((xmax-xmin)/f_f_cell_size)                
      num_y_ff=int((ymax-ymin)/f_f_cell_size)               
      num_z_ff=int((zmax-zmin)/f_f_cell_size)              

      x_ff=f_f_cell_size                                    
      y_ff=f_f_cell_size                                   
      z_ff=f_f_cell_size                                  
      x_ff_div2=x_ff*0.500                          
      y_ff_div2=y_ff*0.500                         
      z_ff_div2=z_ff*0.500                        

C RELEASE put in a catch to make sure that cell sizes are big enough...

      cut_max=max(cut_v_st,cut_v_md)
      cut_max=max(cut_max,cut_g_st)
      cut_max=max(cut_max,cut_g_md)
      cut_max=max(cut_max,cut_e_st)
      cut_max=max(cut_max,cut_e_md)
      if(f_f_cell_size.lt.cut_max) then
        write(*,*)
        write(*,*)'f_f_cell_size = ',f_f_cell_size
        write(*,*)'but biggest cutoff = ',cut_max
        write(*,*)'you need to increase f_f_cell_size'
        write(*,*)'quitting :('
        write(*,*)
        stop
      endif

C RELEASE also put in a catch if system dimensions don't cleanly
C divide up using the input f_f_cell_size

      if(abs(num_x_ff*f_f_cell_size-(xmax-xmin)).gt.0.1.or.
     &   abs(num_y_ff*f_f_cell_size-(ymax-ymin)).gt.0.1.or.
     &   abs(num_z_ff*f_f_cell_size-(zmax-zmin)).gt.0.1) then
        write(*,*)
        write(*,*)'f_f_cell_size = ',f_f_cell_size
        write(*,*)'this does not fit the system dimensions'
        write(*,*)'(xmax-xmin / ymax-ymin / zmax-zmin'
        write(*,*)'make sure all dims are divisible by f_f_cell_size'
        write(*,*)'quitting :('
        write(*,*)
        stop
      endif

      write(*,*)'flexible atom cell dims ',num_x_ff,num_y_ff,num_z_ff
                                                                                                      
C now read the go contacts for each molecule (if necessary)                                           
C first assume that all pairs of molecule typs don't have Go interactions

      if(i_use_go_pairs) then

        allocate(num_des_go_pairs(1:num_t_ts,1:num_t_ts))
        allocate(num_cur_go_pairs(1:num_t_ts,1:num_t_ts,1:num_reps))
        allocate(first_go_entry(1:num_t_ts,1:num_t_ts))
        allocate(last_go_entry(1:num_t_ts,1:num_t_ts))
        do itp=1,num_t_ts
          do jtp=1,num_t_ts
            num_des_go_pairs(jtp,itp)=0
            first_go_entry(jtp,itp)=0
            last_go_entry(jtp,itp)=0
          enddo
        enddo
        write(*,*)'about to make call to make_go_parameters' 
        call make_go_parameters                         
        write(*,*)'returned from call to make_go_parameters' 

        write(*,*)'max # of binding ',
     &            'modes found for any pair = ',go_mod_max

C make a list of how many contacts are needed for each mode of each type
C of go molecule interaction - note that this is similar to
C num_des_go_pairs but with an additional dimension

C note also that it's an obvious area where we may have memory issues

        allocate(num_des_go_mod_pairs
     &          (1:go_mod_max,1:num_t_ts,1:num_t_ts))

        do itp=1,num_t_ts    
          do jtp=1,num_t_ts
            do ktp=1,go_mod_max
              num_des_go_mod_pairs(ktp,jtp,itp)=0
            enddo
          enddo
        enddo

        open(32,file=goparam_file,form='formatted',status='unknown')
        read(32,789)     
789     format(1x)
802     read(32,803,end=299)itp,iat,jtp,jat,t0,e1,id
803     format(4i10,2f15.5,i10)                                                           

C need the following line if there are no mode identifiers
C note that in such a case we shouldn't use i_use_exclusive_go

        if(id.eq.0) id=1

C note that we use go_eps_low here so that only those go contacts with
C an epsilon above a threshold value get included in the calculation of
C Q values 

        if(e1.ge.go_eps_low) 
     &    num_des_go_mod_pairs(id,itp,jtp)=
     &    num_des_go_mod_pairs(id,itp,jtp)+1

        goto 802                                                                            
299     close(32)      

C crappy way to limit the size of the go arrays by finding the last of
C the proteins involved in go contacts... poor but will work

        num_g_ms=0 ! last mol# involved in go contacts
        do iii=1,num_t_ms
          do jjj=1,num_t_ms
            itp=i_t_tp(iii)
            jtp=i_t_tp(jjj)
            if(num_des_go_pairs(itp,jtp).gt.0) then
              if(iii.gt.num_g_ms) num_g_ms=iii
              if(jjj.gt.num_g_ms) num_g_ms=jjj
            endif
          enddo
        enddo

C simple fix to make sure no seg faults due to num_cur_Q

        if(num_g_ms.eq.0) num_g_ms=1

        write(*,*)'# of molecules with Go contacts = ',num_g_ms

        allocate(hist_of_Q(1:num_g_ms,1:num_g_ms,1:40,1:num_reps))
        allocate(num_cur_Q(1:num_g_ms,1:num_g_ms,1:num_reps))                   
        allocate(num_req_Q(1:num_g_ms,1:num_g_ms),stat=ierror(1))                   

        hist_of_Q=0

        do itp=1,num_t_ts
          do jtp=1,num_t_ts
            if(num_des_go_pairs(itp,jtp).ge.1)
     &      write(*,*)'for moltypes ',itp,' & ',jtp,' #go-contacts = ',
     &        num_des_go_pairs(itp,jtp)
          enddo
        enddo

        do iii=1,num_g_ms
          do jjj=1,num_g_ms
            itp=i_t_tp(iii)
            jtp=i_t_tp(jjj)
            num_req_Q(iii,jjj)=num_des_go_pairs(itp,jtp)

            if(num_req_Q(iii,jjj).ge.1.and.
     &         iii.eq.mol_Q1.and.jjj.eq.mol_Q2) 
     &      write(*,*)'for mols ',iii,' & ',jjj,' #go-contacts = ',
     &        num_req_Q(iii,jjj)

          enddo
        enddo

        allocate(num_req_Q_mod
     &          (1:go_mod_max,1:num_g_ms,1:num_g_ms))
        do mmm=1,go_mod_max
          do iii=1,num_g_ms
            do jjj=1,num_g_ms
              itp=i_t_tp(iii)
              jtp=i_t_tp(jjj)
              num_req_Q_mod(mmm,iii,jjj)=
     &        num_des_go_mod_pairs(mmm,itp,jtp)
            enddo
          enddo
        enddo

      endif

C note that id_beg etc now refer to the number of first something MINUS 1    

C num_tot_atms is the total atoms in system (3 for every rigid)
C num_flx_atms is the total flexible atoms in system 
C num_rgd_atms is the total rigid    atoms in system
C num_hyd_atms is the hydrodynamic atoms in each system (1 atm for every rigid)
C 2019 num_mov_atoms is the total # of non-static atoms in the system

C it's questionable whether we'll use the lists other than num_tot_atms

      num_tot_atms=0                                                   
      num_flx_atms=0                                             
      num_rgd_atms=0                                             
      num_mov_atms=0                                             

      num_t_bs=0                                                   
      max_f_as=0                                              
      max_f_bs=0 ! max # of blobs in a flexible molecule
      max_r_as=0                                             
      max_h_as=0 ! max # of hydrodynamic atoms in a rigid molecule

      do iii=1,num_t_ms

        id_beg(iii)=num_tot_atms                            
        id_end(iii)=num_tot_atms

        num_tot_atms=num_tot_atms+num_atms_pt(i_t_tp(iii))

        if(ityp(iii).eq.1) then
          num_flx_atms=num_flx_atms+num_atms_pt(i_t_tp(iii))
        else
          num_rgd_atms=num_rgd_atms+num_atms_pt(i_t_tp(iii))
        endif

      enddo                     

      write(*,*)'num_tot_atms = ',num_tot_atms
      write(*,*)'num_flx_atms = ',num_flx_atms
      write(*,*)'num_rgd_atms = ',num_rgd_atms

C now we want to flag atoms of rigid molecules that will not be used for
C identifying nonbonded interactions - we'll only use the 3rd atom of
C each rigid molecule for this (it's at the origin)

      allocate(i_get_mapped_to_grid(1:num_tot_atms))

      i_get_mapped_to_grid=1

      do iii=1,num_t_ms
        if(ityp(iii).eq.0) then
          i_get_mapped_to_grid(id_beg(iii)+1)=0
          i_get_mapped_to_grid(id_beg(iii)+2)=0
        endif
      enddo

C now count up # of uniq atom and non-uniq atom types 
C (the latter should equal # atoms in system)
       
      allocate(my_uniq_atm_num(1:num_tot_atms))

      num_nonuniqs=0
      num_uniqs=0

      do iii=1,num_t_ts
        do jjj=1,num_ths_t_typ(iii)
          do kkk=1,num_atms_pt(iii)
            if(jjj.eq.1) num_uniqs=num_uniqs+1
            num_nonuniqs=num_nonuniqs+1
            nnn=0        
            do lll=1,iii-1
              nnn=nnn+num_atms_pt(lll)
            enddo
            my_uniq_atm_num(num_nonuniqs)=nnn+kkk   
          enddo
        enddo
      enddo

      if(num_nonuniqs.ne.num_flx_atms) then
        write(*,*)'problem adding up #non-unique atoms ',num_nonuniqs,
     &            ' should be ',num_flx_atms
        write(*,*)'quitting :('
        stop
      endif        

      write(*,*)'num unique atoms    = ',num_uniqs
      write(*,*)'num nonunique atoms = ',num_nonuniqs

      write(*,*)
      write(*,*)'allocating i_am_hyd with #atoms = ',num_tot_atms
      write(*,*)

      allocate(i_am_hyd(1:num_tot_atms))  ! whether HI or not
      i_am_hyd=0                          ! default

      if(full_hydrodynamics) then

        i_am_hyd=1
        num_hyd_syss=1
        allocate(num_hyd_atms(1:num_hyd_syss))
        allocate(num_hyd_atms3(1:num_hyd_syss))
        num_hyd_atms(1)=num_tot_atms
        max_hyd_atms=num_tot_atms
        allocate(atm_to_hyd_atm(1:num_tot_atms))
        allocate(atm_to_hyd_sys(1:num_tot_atms))
        allocate(hyd_atm_to_atm(1:max_hyd_atms,num_hyd_syss))
        do jjj=1,num_tot_atms
          atm_to_hyd_atm(jjj)=jjj
          atm_to_hyd_sys(jjj)=1
          hyd_atm_to_atm(jjj,1)=jjj
        enddo

C allocate memory for a bunch of cholesky-related terms:

        allocate(keep(1:num_hyd_syss))
        allocate(ctrl(1:num_hyd_syss))
        allocate(aoft(1:num_hyd_syss))
        allocate(nh3(1:num_hyd_syss))
        allocate(nb(1:num_hyd_syss))
        allocate(nrblk(1:num_hyd_syss))
        allocate(la(1:num_hyd_syss))
        allocate(scale_nb(1:num_hyd_syss))
        allocate(nthrd(1:num_hyd_syss))

C now add entries so that we can use the same cholesky code for full HI
C and for the mixed schedule HI...

        num_chol_rounds=1
        allocate(chol_scdl(1:num_chol_rounds)) 
        chol_scdl(1)%num=1
        chol_scdl(1)%cnt=1
        allocate(chol_scdl(1)%idm(1:chol_scdl(1)%num))
        allocate(chol_scdl(1)%trd(1:chol_scdl(1)%num))
        allocate(chol_scdl(1)%scl(1:chol_scdl(1)%num)) 
        allocate(chol_scdl(1)%mst(1:num_threads))
        allocate(chol_scdl(1)%jdm(1:num_threads))

C the following would normally be read from chol_schedule_file:

        chol_scdl(1)%idm(1:chol_scdl(1)%num)=1
        chol_scdl(1)%trd(1:chol_scdl(1)%num)=num_threads
        chol_scdl(1)%scl(1:chol_scdl(1)%num)=scale_nb_temp

C at this point the code is identical to that used above when we have
C multiple different HI systems specified by mixed_schedule_file

        nnn_thrd=0
        do jjj=1,chol_scdl(1)%num
          do kkk=1,chol_scdl(1)%trd(jjj)
            nnn_thrd=nnn_thrd+1
            if(kkk.eq.1) then
              chol_scdl(1)%mst(nnn_thrd)=1
            else
              chol_scdl(1)%mst(nnn_thrd)=0
            endif
            chol_scdl(1)%jdm(nnn_thrd)=chol_scdl(1)%idm(jjj)
            write(*,745)1,nnn_thrd,
     &      chol_scdl(1)%mst(nnn_thrd),
     &      chol_scdl(1)%jdm(nnn_thrd)
          enddo
          if(jjj.eq.chol_scdl(1)%num.and.
     &       nnn_thrd.ne.num_threads) then
            write(*,*)'error: nnn_thrd = ',nnn_thrd,
     &                ' num_threads = ',num_threads
            stop
          endif
        enddo
        do nnn=1,num_hyd_syss
          nthrd(nnn)=0
          do iii=1,num_chol_rounds
            do jjj=1,chol_scdl(iii)%num
              if(chol_scdl(iii)%idm(jjj).eq.nnn) then
                nthrd(nnn)   =chol_scdl(iii)%trd(jjj)
                scale_nb(nnn)=chol_scdl(iii)%scl(jjj)
                write(*,744)nnn,iii,
     &                      chol_scdl(iii)%trd(jjj),
     &                      chol_scdl(iii)%scl(jjj)
              endif
            enddo
          enddo
        enddo

745     format('In round # ',i5,' thread # ',i6,
     &         ' is a master? ',i3,' and works on HI system # ',i6)
744     format('cholesky decomposition for HI system ',i6,
     &         ' will be done in round ',i6,
     &         ' using ',i6,' threads and scale factor ',f10.3)

      elseif(no_hydrodynamics) then

        i_am_hyd=0
        num_hyd_syss=0
        allocate(atm_to_hyd_atm(1:num_tot_atms))
        allocate(atm_to_hyd_sys(1:num_tot_atms))
        do jjj=1,num_tot_atms
          atm_to_hyd_atm(jjj)=0
          atm_to_hyd_sys(jjj)=0
        enddo

      endif

      write(*,*)
      write(*,*)'after considering mixed/full/noHI, we have ',
     &          'max_hyd_atms = ',max_hyd_atms
      write(*,*)

      num_tot_atms3=num_tot_atms*3
      num_flx_atms3=num_flx_atms*3
      max_hyd_atms3=max_hyd_atms*3
      do iii=1,num_hyd_syss
        num_hyd_atms3(iii)=num_hyd_atms(iii)*3
      enddo

C allocate hydrodynamic entries - start with a bunch that need to be
C allocated for all atoms regardless...

      if(i_use_hydro) then

        allocate(r_t_a(1:num_tot_atms))  
        allocate(d_t_a(1:num_tot_atms))  

        allocate(d_t_a_3frm(1:num_tot_atms3)) ! same as d_t_a
        allocate(d_t_a_root(1:num_tot_atms3))  
        allocate(d_t_a_root_inv(1:num_tot_atms3))

        if(langevin) then
          allocate(f_e_r(1:num_tot_atms3,1:num_reps))
          allocate(f_r_r(1:num_tot_atms3,1:num_reps))
          allocate(v_e_r(1:num_tot_atms3,1:num_reps))
          allocate(v_r_r(1:num_tot_atms3,1:num_reps))
          v_e_r=0.0
          v_r_r=0.0
          allocate(ld_fac1(1:num_tot_atms3))
          allocate(ld_fac2(1:num_tot_atms3))
          allocate(ld_fac3(1:num_tot_atms3))
          allocate(ld_fac4(1:num_tot_atms3))
        endif

        if(full_hydrodynamics .and.cholesky) then

          num_tmp_long=max_hyd_atms3 ! num_tmp_long is a long-integer...
          ma=num_tmp_long*(num_tmp_long+1)/2

C note that we added a dimension for num_reps to amp since this stores
C the actual diffusion tensor; no need for it on imp and jmp

          allocate(amp(1:ma,1:num_hyd_syss,1:num_reps))
          allocate(imp(1:ma,1:num_hyd_syss))
          allocate(jmp(1:ma,1:num_hyd_syss))

          do n=1,num_hyd_syss
            nh3(n)=num_hyd_atms3(n) ! used for cleaning up calls to HSL
            write(*,*)'allocating hydrodynamic arrays for ',nh3(n)
            write(*,*)'note we include only 3rd atom of each rigid mol'
c           nb(n)=mp54_get_nb(num_hyd_atms3(n),num_hyd_atms3(n),
c    &                        nthrd(n),ctrl(n))
            nb(n)=1
            nb(n)=nb(n)*scale_nb(n)
            aoft(n)=1      
            nrblk(n)=(num_hyd_atms3(n)-1)/nb(n)+1
            ma=0
            do i=1,num_hyd_atms3(n)
              do j=i,num_hyd_atms3(n)
                ma=ma+1
                jmp(ma,n)=j
                imp(ma,n)=i
              enddo
            enddo
            la(n)=ma
          enddo ! end loop over HI systems

        endif

      endif
      write(*,*)'finished allocating arrays for hydrodynamics'

C more allocations 
C note that we will almost always be using num_tot_atms

      write(*,*)
      write(*,*)'allocating ifree with #atoms = ',num_tot_atms
      write(*,*)

      allocate(ifree(1:num_tot_atms))  ! whether free
      allocate(jfree(1:num_tot_atms))  ! whether free
      ifree=1 ! default

C new for YHL code - allocate memory for thread-thread interactions

      allocate(thrd_of_atm(1:num_tot_atms,1:num_reps))  ! thread# of each atom
      thrd_of_atm=0 ! default

      allocate(mol_of_atm(1:num_tot_atms))  ! molecule number of atom
      allocate(crg_of_atm(1:num_tot_atms))  ! charge of atom
      allocate(stt_of_atm(1:num_tot_atms))  ! 2019 =1 if static,
                                            ! =0 if dynamic

C probably easier to have these allocated as 1:num_tot_atms as they'll
C be used a lot during the dynamics

      allocate(i_f_a(1:num_tot_atms))  ! flag whether atom is charged
      i_f_a=0                          ! default to zero
      allocate(i_q_a(1:num_tot_atms))  ! charge number of this atom

C 2019 - if we have some atoms that don't move then we need to put
C entries into c_typ(:)%s(:) - do this now...

      if(i_dont_move_some) then                                           
        do iii=1,num_t_ts
          c_typ(iii)%s=0 ! assume all are dynamic for now
        enddo
        do iii=1,num_moves_omitted
          jjj=i_omit_move_mtyp(iii)
          kkk=i_omit_move_a(iii)
          c_typ(jjj)%s(kkk)=1
        enddo
      endif

C assign other parameters

      ntmp=0
      do iii=1,num_t_ms                                     
        do jjj=1,num_atms_pt(i_t_tp(iii))   
          ntmp=ntmp+1

C we can set each of these for flexible and rigid mols

          mol_of_atm(ntmp)=iii
          crg_of_atm(ntmp)=c_typ(i_t_tp(iii))%q(jjj)
          stt_of_atm(ntmp)=c_typ(i_t_tp(iii))%s(jjj) ! 2019

C 2019 increment num_mov_atms also here

          if(stt_of_atm(ntmp).ne.1) num_mov_atms=num_mov_atms+1

          r_t_a(ntmp)=c_typ(i_t_tp(iii))%r(jjj)

C need to be careful when converting hydrodynamic radius to diffusion
C coefficient (rigid mols have some zero radii); also want to correct 
C the dtrans values for non-standard viscosity

C d_t_a_root_inv etc are dimensioned 1:3N for convenience when
C computing random displacements later...

          if(r_t_a(ntmp).ne.0.0) then
            d_t_a(ntmp)=0.246897022/r_t_a(ntmp)
            d_t_a(ntmp)=d_t_a(ntmp)*0.89/viscosity
            kt=(ntmp-1)*3
            do j=1,3
              d_t_a_3frm(kt+j)=d_t_a(ntmp)
              d_t_a_root(kt+j)=sqrt(d_t_a(ntmp))
              d_t_a_root_inv(kt+j)=1.0/d_t_a_root(kt+j)
            enddo
          else
            d_t_a(ntmp)=0.0
            kt=(ntmp-1)*3
            do j=1,3
              d_t_a_3frm(kt+j)=d_t_a(ntmp)
              d_t_a_root(kt+j)=0.0
              d_t_a_root_inv(kt+j)=0.0
            enddo
          endif

          if(langevin) then
 
C first get gama (friction coefficient) in kcal.ps/mol/A^2
 
            gama=(r_kbt*temperature)/d_t_a(ntmp)

C figure out the mass - this is a mess, just use default
C RELEASE we set mass=0.001 - this used to be read-in from the input
C file but it's something we'll never mess with so hardwire it here

C           mass=rtemp ! this was read-in earlier
            mass=0.001
 
C convert mass from Daltons to kcal.ps^2/mol/A^2 such that
C mass/gama/time_step is unitless

            mass=mass*2390.057362

C _fac1 is the first factor  of eq.13 in winter & geyer
C _fac2 is the second factor of eq.13
C _fac3 is the first factor  of eq.14 in winter & geyer
C _fac4 is the second factor of eq.14

            kt=(ntmp-1)*3
            do j=1,3
              ld_fac1(kt+j)=1.0-exp(-gama*time_step/mass)
              ld_fac1(kt+j)=1.0-(mass*ld_fac1(kt+j)/(gama*time_step))
              ld_fac2(kt+j)=1.0-exp(-gama*time_step/mass)
              ld_fac2(kt+j)=ld_fac2(kt+j)* mass/time_step
              ld_fac3(kt+j)=exp(-gama*time_step/mass)
              ld_fac4(kt+j)=1.0-exp(-gama*time_step/mass)
              ld_fac4(kt+j)=ld_fac4(kt+j)/gama
            enddo
          endif

C now do some pH/charge stuff that is only for flexible mols

          if(ityp(iii).eq.1) then                                

            if(abs(crg_of_atm(ntmp)).gt.0.001) then
              i_f_a(ntmp)=1
            else
              i_f_a(ntmp)=0
            endif

          endif
        enddo
      enddo

C 2019 write out total # of moving atoms

      write(*,*)'num_mov_atms = ',num_mov_atms

C compute 6pi.eta, kt/6pi.eta etc here once...
C probably should check these at some point - this is pretty crap

      six_pi_eta=r_kbt*temperature/d_t_a(1)/r_t_a(1)
      kT_over_6pi_eta=d_t_a(1)*r_t_a(1)
      kT_over_8pi_eta=kT_over_6pi_eta*6.0/8.0  

C allocate memory for coords of all atoms - note padding by 16 entries
C note that t_a_x is just a temporary array that we use so that we can
C do a NUMA-like parallel initialization of c_a_x in openmp

      write(*,*)'allocating t_a_x,c_a_x,f_a_x to ',num_tot_atms

      allocate(t_a_x(1:3,1:num_tot_atms,1:num_reps)) 
      allocate(c_a_x(1:3,1:num_tot_atms,1:num_reps)) 
      allocate(f_a_x(1:3,1:num_tot_atms,1:num_reps))  

C 2022 GENERAL - add div_x,div_y,div_z for OA

      allocate(div_x(1:num_tot_atms))
      allocate(div_y(1:num_tot_atms))
      allocate(div_z(1:num_tot_atms))

C allocate memory for current positions if using lincs

      if(i_do_lincs) then
        allocate(c_a_x_old(1:3,1:num_tot_atms,1:num_reps)) 
      endif

C new for YHL code - need a counter for all forces

      allocate(g_a_x(1:3,1:num_tot_atms,1:num_threads,1:num_reps))                
                                                                                                      
C now read the bond parameters for each molecule                                                      
 
      write(*,*)
      write(*,*)'making call to internal parameters routine'
      call make_internal_parameters                                   
      write(*,*)'returned from internal parameters routine'
      write(*,*)

C assign some more hydrodynamic parameters

C note that our 'default' value is time_step - we only use time_step
C for terms that explicitly require the shorter time_step
C however, since all of that code is poorly supported we will require at
C this point that the two timesteps be identical...

      if(i_use_hydro) then
  
        s_hydro=time_step/(r_kbt*temperature) 
        d_hydro=sqrt(2.000*time_step)             
        f_hydro=sqrt(2.000)*r_kbt*temperature/sqrt(time_step)
        df_hydro=d_hydro/f_hydro

        write(*,*)'s_hydro  = ',s_hydro
        write(*,*)'d_hydro  = ',d_hydro
        write(*,*)'f_hydro  = ',f_hydro
        write(*,*)'df_hydro = ',df_hydro

      endif
                                                                                                      
C now add up all the bonds, angls, dihes in the system - this includes 
C all copies of each bond/angl/dihe - in contrast to num_tot_bonds  
C which is the number of total unique types of bond/angl/dihe
       
      write(*,*)'adding up all bonds, angles, dihedrals in system'

      num_sys_bonds=0                                          
      num_sys_angls=0                                        
      num_sys_dihes=0                                    
      do iii=1,num_t_ms                                      
        if(ityp(iii).eq.1) then                                 
          num_sys_bonds=num_sys_bonds+num_bonds(i_t_tp(iii))
          num_sys_angls=num_sys_angls+num_angls(i_t_tp(iii))
          num_sys_dihes=num_sys_dihes+num_dihes(i_t_tp(iii))
        endif                                                           
      enddo    

      write(*,*)'Total System Bonds,Angles & Dihedrals Before Cull ',
     &           num_sys_bonds,num_sys_angls,num_sys_dihes

      allocate(iiib_new(1:num_sys_bonds))
      allocate(nbb0_new(1:num_sys_bonds))
      allocate(nbb1_new(1:num_sys_bonds))
      allocate(nbb2_new(1:num_sys_bonds))
      allocate(iiia_new(1:num_sys_angls))
      allocate(nba0_new(1:num_sys_angls))
      allocate(nba1_new(1:num_sys_angls))
      allocate(nba2_new(1:num_sys_angls))
      allocate(nba3_new(1:num_sys_angls))
      allocate(iiid_new(1:num_sys_dihes))
      allocate(nbd0_new(1:num_sys_dihes))

      nbb0_new=0
      nba0_new=0
      nbd0_new=0

      allocate(nbd1_new(1:num_sys_dihes))
      allocate(nbd2_new(1:num_sys_dihes))
      allocate(nbd3_new(1:num_sys_dihes))
      allocate(nbd4_new(1:num_sys_dihes))
      allocate(r0bd_new(1:num_sys_dihes))
      allocate(rkbb_new(1:num_sys_bonds))
      allocate(r0bb_new(1:num_sys_bonds))
      allocate(rkba_new(1:num_sys_angls))
      allocate(r0ba_new(1:num_sys_angls))
      allocate(tors1_new(1:4,1:num_sys_dihes))
      allocate(tors3_new(1:4,1:num_sys_dihes))

      numtmp1=0
      numtmp2=0
      numtmp3=0
      do iii=1,num_t_ms                                      
        if(ityp(iii).eq.1) then                                 
          do jjj=1,num_bonds(i_t_tp(iii))     
            numtmp1=numtmp1+1
            iiib_new(numtmp1)=iii                                               
            nbb1_new(numtmp1)=nbb1(jjj+id_beg_bonds(i_t_tp(iii)))+
     &                       id_beg(iii)
            nbb2_new(numtmp1)=nbb2(jjj+id_beg_bonds(i_t_tp(iii)))+
     &                       id_beg(iii)
            r0bb_new(numtmp1)=r0bb(jjj+id_beg_bonds(i_t_tp(iii)))
            rkbb_new(numtmp1)=rkbb(jjj+id_beg_bonds(i_t_tp(iii)))
          enddo
          do jjj=1,num_angls(i_t_tp(iii))     
            numtmp2=numtmp2+1
            iiia_new(numtmp2)=iii                                               
            nba1_new(numtmp2)=nba1(jjj+id_beg_angls(i_t_tp(iii)))+
     &                       id_beg(iii)
            nba2_new(numtmp2)=nba2(jjj+id_beg_angls(i_t_tp(iii)))+
     &                       id_beg(iii)
            nba3_new(numtmp2)=nba3(jjj+id_beg_angls(i_t_tp(iii)))+
     &                       id_beg(iii)
            r0ba_new(numtmp2)=r0ba(jjj+id_beg_angls(i_t_tp(iii)))
            rkba_new(numtmp2)=rkba(jjj+id_beg_angls(i_t_tp(iii)))
          enddo
          do jjj=1,num_dihes(i_t_tp(iii))     
            numtmp3=numtmp3+1
            iiid_new(numtmp3)=iii                                               
            nbd1_new(numtmp3)=
     &        nbd1(jjj+id_beg_dihes(i_t_tp(iii)))+id_beg(iii)
            nbd2_new(numtmp3)=
     &        nbd2(jjj+id_beg_dihes(i_t_tp(iii)))+id_beg(iii)
            nbd3_new(numtmp3)=
     &        nbd3(jjj+id_beg_dihes(i_t_tp(iii)))+id_beg(iii)
            nbd4_new(numtmp3)=
     &        nbd4(jjj+id_beg_dihes(i_t_tp(iii)))+id_beg(iii)
            r0bd_new(numtmp3)=
     &        r0bd(jjj+id_beg_dihes(i_t_tp(iii)))+id_beg(iii)
            tors1_new(1,numtmp3)=
     &        tors1(1,jjj+id_beg_dihes(i_t_tp(iii)))
            tors1_new(2,numtmp3)=
     &        tors1(2,jjj+id_beg_dihes(i_t_tp(iii)))
            tors1_new(3,numtmp3)=
     &        tors1(3,jjj+id_beg_dihes(i_t_tp(iii)))
            tors1_new(4,numtmp3)=
     &        tors1(4,jjj+id_beg_dihes(i_t_tp(iii)))
            tors3_new(1,numtmp3)=
     &        tors3(1,jjj+id_beg_dihes(i_t_tp(iii)))
            tors3_new(2,numtmp3)=
     &        tors3(2,jjj+id_beg_dihes(i_t_tp(iii)))
            tors3_new(3,numtmp3)=
     &        tors3(3,jjj+id_beg_dihes(i_t_tp(iii)))
            tors3_new(4,numtmp3)=
     &        tors3(4,jjj+id_beg_dihes(i_t_tp(iii)))
          enddo
        endif
      enddo

      num_sys_bonds=numtmp1
      num_sys_angls=numtmp2
      num_sys_dihes=numtmp3

      write(*,*)'Total System Bonds,Angles & Dihedrals After  Cull ',
     &           num_sys_bonds,num_sys_angls,num_sys_dihes


C lincs additions here

      if(i_do_lincs) then

        write(*,*)
        write(*,*)'allocating memory for lincs related stuff'
        write(*,*)

        allocate(ncc(1:num_sys_bonds))
        allocate(sdia(1:num_sys_bonds))
        allocate(icon(1:10,1:num_sys_bonds))
        allocate(coef(1:10,1:num_sys_bonds))

C note that till very recently the following 4 arrays were all private
C declarations - I don't know why - I was hoping that making them global
C would help lincs work better for multithreaded cases but it doesn't
C seem to help at all - there must be something else wrong somewhere...

        allocate(alnc(1:10,num_sys_bonds))
        allocate(blnc(1:3,1:num_sys_bonds))
        allocate(rhsc(1:num_sys_bonds,1:2))
        allocate(solc(1:num_sys_bonds))

        do jjj=1,num_sys_bonds
          sdia(jjj)=1.0/(sqrt(2.0/1.0)) ! denom is generic mass of 1.000
        enddo
        do jjj=1,num_sys_bonds
          ncc(jjj)=0
          do kkk=jjj+1,num_sys_bonds
            if(nbb1_new(jjj).eq.nbb1_new(kkk).or.
     &         nbb1_new(jjj).eq.nbb2_new(kkk).or.
     &         nbb2_new(jjj).eq.nbb1_new(kkk).or.
     &         nbb2_new(jjj).eq.nbb2_new(kkk)) then
              ncc(jjj)=ncc(jjj)+1
              icon(ncc(jjj),jjj)=kkk
              rsign=1.0
              if(nbb1_new(jjj).eq.nbb1_new(kkk).or.
     &           nbb2_new(jjj).eq.nbb2_new(kkk)) rsign=-1.0
              coef(ncc(jjj),jjj)=rsign*(1.0/1.000)* ! mass of 1.000
     &           sdia(jjj)*sdia(icon(ncc(jjj),jjj))
            endif
          enddo
        enddo
      endif

C 30-11-11 

C here we can look to see if we're using special terms for 
C bonds, angles and dihedrals...
C in order to do this we need to have, for every atom in the system,
C its atom name, its residue name and its residue number...

      allocate(atm_nam_atm(1:num_tot_atms))                                                           
      allocate(res_nam_atm(1:num_tot_atms))                                                           
      allocate(res_num_atm(1:num_tot_atms))                                                           

      ntmp=0                                                                                    
      do iii=1,num_t_ms                                                                             
        do jjj=1,num_atms_pt(i_t_tp(iii))                                                          
          ntmp=ntmp+1                                                                   
          atm_nam_atm(ntmp)=n_typ(i_t_tp(iii))%ch2(jjj)(2:4)
          res_nam_atm(ntmp)=n_typ(i_t_tp(iii))%ch3(jjj)
          read(n_typ(i_t_tp(iii))%ch4(jjj),'(i4)') res_num_atm(ntmp)
        enddo
      enddo

C prepare for position restraints - allocate flag here

      allocate(i_restraint(1:num_tot_atms)) 
      i_restraint=0

C note that the next allocation is relevant to position restraints but
C must be done regardless of whether i_do_pos_restraint=.true.
C i_res_no_nb = 0 means evaluate nonbond inters as normal
C             = + skip evaluation of nonbond inters for this atm
C             = - do evaluation of nonbond inters for this atm

      allocate(i_res_no_nb(1:num_tot_atms)) ! flag =1 to skip nb inters
      i_res_no_nb=0
      num_i_res_no_nb_grps=0

      if(i_do_pos_restraint) then

C new for 32-11-11 - allow option of skipping nonbondeds for restraineds

        allocate(fx_restraint(1:num_tot_atms)) ! force constant
        allocate(fy_restraint(1:num_tot_atms)) ! force constant
        allocate(fz_restraint(1:num_tot_atms)) ! force constant
        allocate(fr_restraint(1:num_tot_atms)) ! force constant
        allocate(cx_restraint(1:num_tot_atms)) ! x ref position
        allocate(cy_restraint(1:num_tot_atms)) ! y ref position
        allocate(cz_restraint(1:num_tot_atms)) ! z ref position
        allocate(dx_restraint(1:num_tot_atms)) ! x ref position
        allocate(dy_restraint(1:num_tot_atms)) ! y ref position
        allocate(dz_restraint(1:num_tot_atms)) ! z ref position
        allocate(r_restraint(1:num_tot_atms))  ! ref radius     

C now actually read the restraints

        call read_position_restraints 

C before going on we will determine num_i_as: the total number of atoms
C that see nonbonded interactions

C upate - the following stuff is either unnecessary or is currently not
C implemented - it's the beginnings of an attempt to load-balance the
C non-interacting atoms. may not be necessary at all.

C second update - don't just delete num_i_as: it's used to provide better
C estimates of the memory requirements of the code
 
        num_i_as=0
        do i=1,num_tot_atms
          if(i_res_no_nb(i).eq.0) num_i_as=num_i_as+1
        enddo
 
        write(*,*)'out of ',num_tot_atms,' atoms in the system, only ',
     &            num_i_as,' of them engage in nonbonded  interactions'

        write(*,*)'the number of non-bonded skipping groups is ',
     &            num_i_res_no_nb_grps 

        if(num_i_as.eq.0) num_i_as=1 ! kluge to ensure no seg fault
 
C close the if statment over i_do_pos_restraint

      endif
                                                                                                      
      allocate(my_cell_1_ref_x(1:num_tot_atms),stat=ierror(27))
      allocate(my_cell_1_ref_y(1:num_tot_atms),stat=ierror(28)) 
      allocate(my_cell_1_ref_z(1:num_tot_atms),stat=ierror(29)) 
      allocate(my_cell_2_ref_x(1:num_tot_atms),stat=ierror(30)) 
      allocate(my_cell_2_ref_y(1:num_tot_atms),stat=ierror(31)) 
      allocate(my_cell_2_ref_z(1:num_tot_atms),stat=ierror(32)) 

      biggest_ene_jump=0.0
                        
C reset the number of steps,time etc.                                                                
                       
      nstep=0                                                     
      tot_time=0.0000                                            
      necount=neprint                                           
      ntcount=ntprint                                          
      if(nmprint.ne.-1) nmcount=nmprint
      num_umb_cnt=num_umb_stp                   
      num_exc_cnt=num_exc_stp                   
      num_lst_cnt=num_lst_stp                   
      num_fmd_cnt=num_fmd_stp                   
      num_hyd_cnt=num_hyd_stp                    
      num_bal_cnt=0                 
      num_ene_scrn_cnt=1000                              
                      
C now read in all coordinates from a restart.file                                                     
C note that this does not require that all molecules be written out
C in order - don't know if this flexibility is a good idea but...                        
              
      tot_time_orig=0.000         

      write(*,*)
      write(*,*)'starting reading the restart.file'
      call read_restart_file
      write(*,*)'finished reading the restart.file'
      write(*,*)

C allocate memory for running last ten frames to recreate crashes

      if(i_look_for_crashes) then
        allocate(c_c_x(1:3,1:num_tot_atms,1:10,1:num_reps))
      endif

      ene_ave=0.0d0
      ene_var=0.0d0
      ene_den=0.0d0
      qqq_ave=0.0d0
      qqq_var=0.0d0
      qqq_den=0.0d0
      rand_ave=0.0d0
      rand_msd=0.0d0
      rand_den=0.0d0

C open up the wall file and read, for every unique type of atom, which
C walls it sees - #columns must be the same as #walls..
C need to read them and then assign on the basis of #totatoms

      if(walls) call read_wall_functions(num_walls)

C 09-08-07 
C include cotranslational option here by initially setting all
C flexible atoms to be free

      do iii=1,num_t_ms                                
        do jjj=1,num_atms_pt(i_t_tp(iii))       
          c_f_m(iii)%free(jjj)=1                   
          c_f_m(iii)%harm(jjj)=0                   
          c_f_m(iii)%harm_gone(jjj)=0                   
          c_f_m(iii)%just_harm(jjj)=0                   
          c_f_m(iii)%just_free(jjj)=0                   
          c_f_m(iii)%near_free(jjj)=0                   
        enddo
      enddo

C Having assigned coordinates to every molecule,it's time to start                                   
C all of the following could go to a subroutine...

C note that we make the coordinates,total forces,dmatrix and smatrix
C global,shared arrays                                                

      allocate(time_h_gen(1:num_threads))
      allocate(time_h_mov(1:num_threads))
      allocate(tmp_local(1:num_threads))
      allocate(mywdth(1:num_threads))
      allocate(mywdth1(1:num_threads))
      allocate(mywdth2(1:num_threads))
      allocate(mywdth3(1:num_threads))
      allocate(mywdth4(1:num_threads))
      allocate(mywdth5(1:num_threads))
      allocate(num_loc_f_as(1:num_threads))
      allocate(num_loc_r_as(1:num_threads))
      allocate(num_loc_h_as(1:num_threads))

C allocate so as to minimize false sharing even at boundaries
C I read somewhere that EM64T cpus have a 64-byte cache line
C hence the padding by 16 4-byte reals
C                 
      allocate(c_n_x(1:num_tot_atms))  ! ewald elec terms                       
      allocate(c_n_y(1:num_tot_atms))
      allocate(c_n_z(1:num_tot_atms))
      allocate(c_r_x(1:num_r_as))
      allocate(c_r_y(1:num_r_as))
      allocate(c_r_z(1:num_r_as))

      allocate(ibeg_m_a(1:num_threads)) ! SR terms temporary
      allocate(iend_m_a(1:num_threads)) ! SR terms temporary
      allocate(ibeg_t_a(1:num_threads)) ! SR terms
      allocate(iend_t_a(1:num_threads)) ! SR terms
      allocate(ibeg_bms(1:num_threads)) ! bonds on master
      allocate(ibeg_ams(1:num_threads)) ! angles
      allocate(ibeg_dms(1:num_threads)) ! dihedrals
      allocate(iend_bms(1:num_threads)) ! bonds on master
      allocate(iend_ams(1:num_threads))
      allocate(iend_dms(1:num_threads))

C these allocations used to be private

      allocate(cll_f_1(1:num_x_ff,1:num_y_ff,1:num_z_ff,1:num_reps))
      cll_f_1(num_x_ff,num_y_ff,num_z_ff,num_reps)%num=1 

C 2019 make a second grid for static atoms

      allocate(cll_f_2(1:num_x_ff,1:num_y_ff,1:num_z_ff,1:num_reps))
      cll_f_2(num_x_ff,num_y_ff,num_z_ff,num_reps)%num=1 

      do i4=1,num_reps 
        do i3=1,num_z_ff
          do i2=1,num_y_ff
            do i1=1,num_x_ff
              allocate(cll_f_1(i1,i2,i3,i4)%ida(1))
              allocate(cll_f_2(i1,i2,i3,i4)%ida(1)) ! 2019
            enddo
          enddo 
        enddo  
      enddo  

      if(no_hydrodynamics) then
  
        allocate(f_a_r(1:num_tot_atms3,1:num_reps))                 
        allocate(r_u_r(1:num_tot_atms3,1:num_hyd_stp,1:num_reps))

      elseif(full_hydrodynamics) then

        if(cholesky) 
     &  allocate(s_a_r(1:max_hyd_atms3,1:max_hyd_atms3,1:num_hyd_syss,
     &                 1:num_reps))
C       s_a_r=0.0
        allocate(d_a_r(1:max_hyd_atms3,1:max_hyd_atms3,1:num_hyd_syss,
     &                 1:num_reps))

        allocate(f_a_r(1:num_tot_atms3,1:num_reps))                 
        allocate(r_u_r(1:num_tot_atms3,1:num_hyd_stp,1:num_reps))               
        allocate(r_c_r(1:num_tot_atms3,1:num_hyd_stp,1:num_reps))               

C make force arrays that hold just the hydrodynamic forces
C these are the exact analogs of f_a_r, r_u_r, r_c_r

        allocate(h_a_r(1:max_hyd_atms3,1:num_hyd_syss,1:num_reps))
        allocate(h_u_r(1:max_hyd_atms3,1:num_hyd_stp,1:num_hyd_syss,
     &                 1:num_reps))
        allocate(h_c_r(1:max_hyd_atms3,1:num_hyd_stp,1:num_hyd_syss,
     &                 1:num_reps))

        if(langevin) then
          allocate(h_e_r(1:max_hyd_atms3,1:num_hyd_syss,1:num_reps))
        endif

      endif

C need to assign first_ & last_go_entry for each atom in the system
C note that for very large systems it might be worth reading this
C from a file instead of recomputing it - it's very inefficient code

      if(i_use_go_pairs) then

        write(*,*)'indexing 1D Go pairs array '
        allocate(ind_tot_go_pair(1:num_flx_atms))           
        allocate(jnd_tot_go_pair(1:num_flx_atms))           

        do iii=1,num_flx_atms   
          if(i_debug.and.iii.eq.int(iii/1000)*iii) write(*,*)'AA ',iii
          kkk=iii-id_beg(mol_of_atm(iii))
          itp=i_f_tp(mol_of_atm(iii))
          ind_tot_go_pair(iii)=-1000 
          jnd_tot_go_pair(iii)=-1000 
          do jjj=1,num_tot_go_pairs          
            if(itp_tot_go_pair(jjj).lt.itp) cycle
            if(itp_tot_go_pair(jjj).gt.itp) goto 939
            if(ida_tot_go_pair(jjj).eq.kkk.and.
     &         ind_tot_go_pair(iii).eq.-1000)
     &         ind_tot_go_pair(iii)=jjj
            if(ida_tot_go_pair(jjj).eq.kkk)      
     &         jnd_tot_go_pair(iii)=jjj
          enddo
939       continue
          if(ind_tot_go_pair(iii).eq.-1000) ind_tot_go_pair(iii)=0
          if(jnd_tot_go_pair(iii).eq.-1000) jnd_tot_go_pair(iii)=-1
        enddo 

        if(.not.i_limit_verbosity) then
        open(unit=18,file='go_pair_check.txt',status='unknown')
        do iii=1,num_flx_atms   
          write(18,81)iii,ind_tot_go_pair(iii),jnd_tot_go_pair(iii)
81        format('go pair atom   indices ',3i10)
        enddo
221     format(4i10)                                      
        close(18)                      
        endif                    

      endif
     
C call to initialize the treecode here - note that we not do the
C assignment of i_q_a in the main code as we seem to use this even if we
C don't subsequently use the treecode...

C note that this is done before we enter the huge 'openmp parallel' region
C that contains the actual simulation routines

      num_q_as=0
      num_p_as=0 ! number of protonatable sites

      write(*,*)'computing i_q_a with num_flx_atms = ',num_flx_atms

C 2022 note that use of temporary counter ntmp here - for some reason
C the code can miscount if instead we use num_q_as within the loop -
C e.g. for a system of 5000 charges it comes back with 625 - weird...

      ntmp=0
      do iii=1,num_flx_atms
        if(i_f_a(iii).eq.1.and.ifree(iii).eq.1) then
          ntmp=ntmp+1
          i_q_a(iii)=ntmp
        else
          i_q_a(iii)=0
        endif
      enddo
      num_q_as=ntmp

      if(treecode_elec) call initialize_treecode_elec

C finally - at least for all flexible atoms - if we're using
C the high_mem option then we need to go through all atom pairs
C and produce a huge 2D array marking the status of each pair

        if(i_use_high_mem) then

          allocate(i_status(1:num_flx_atms,1:num_flx_atms))
          allocate(g_status(1:num_flx_atms,1:num_flx_atms))
          write(*,*)'assigning all atom pairs a status'

          if(i_write_nonbond_histogram) then
            allocate(d_status(1:num_flx_atms,1:num_flx_atms))
            d_status=0
          endif

          write(*,*)'i_status for each atom pair takes on the ',
     &              'following values: '
          write(*,*)'(arb. func. means use our IB-derived nonbonded ',
     &              'potential functions)'
          write(*,*)
          write(*,*)'    =0 : atoms are bonded so no nonbonded term'
          write(*,*)'    =1 : atoms interact via vdw + elec '
          write(*,*)'    =2 : atoms interact via vdw only '
          write(*,*)'    =3 : atoms interact via go  + elec '
          write(*,*)'    =4 : atoms interact via go only '
          write(*,*)'    =5 : atoms interact via vdw + elec; arb. func.'
          write(*,*)'    =6 : atoms interact via vdw only;   arb. func.'
          write(*,*)'    =7 : atoms interact via go  + elec; arb. func.'
          write(*,*)'    =8 : atoms interact via go only;    arb. func.'
          write(*,*)

          ncxy=0
          nbxy=0

C default to assuming each pair of atoms is involved in both vdw 
C interaction + electrostatic interaction - note that use of F90 here

C if i_status(iii,jjj)=0 then they're in a bonded interaction...
C                     =1 vdw + elec
C                     =2 vdw only
C                     =3 go + elec
C                     =4 go only
C                     =5 vdw + elec - arbitrary potential
C                     =6 vdw only   - arbitrary potential
C                     =7 go + elec  - arbitrary potential
C                     =8 go only    - arbitrary potential

C n_ff_v_st : all those with i_status = 1/2
C n_ff_g_st : all those with i_status = 3/4
C n_ff_a_st : all those with i_status = 5/6
C n_ff_b_st : all those with i_status = 7/8

          i_status=1
          g_status=0

          do iii=1,num_flx_atms           

            i_status(iii,iii)=0
            if(i_read_nonbond_functions) kkk=atm_typ(iii)

C look at the bonding terms 

            do m=ind_tot_b_typ(iii),
     &           jnd_tot_b_typ(iii)
              jjj=jda_tot_b_typ(m)+id_beg(mol_of_atm(iii))
              i_status(jjj,iii)=0                    
              write(*,*)'setting i_status to zero ',iii,jjj
            enddo

C look at the go terms if necessary

            if(i_use_go_pairs) then
              do m=ind_tot_go_pair(iii),
     &             jnd_tot_go_pair(iii)

C get the type of molecule that this atom forms Go contact with

                jtp=jtp_tot_go_pair(m)

C find all molecules of this type
C note that we only check the kkk,lll,0 dimension of atm_atm_typ because
C it is immaterial whether it is a regular nonbonded interaction or a
C more local nonbonded interaction - the former will always be present
C 
                do n=1,num_t_ms
                  if(i_f_tp(n).eq.jtp) then
                    jjj=id_beg(n)+jda_tot_go_pair(m)
                    if(i_read_nonbond_functions) lll=atm_typ(jjj)
                    if(i_f_a(iii).eq.1. and.i_f_a(jjj).eq.1) then
                      i_status(jjj,iii)=3
                      if(i_read_nonbond_functions) then
                        if(atm_atm_typ(kkk,lll,0).ne.0)
     &                    i_status(jjj,iii)=7
                      endif
                    else
                      i_status(jjj,iii)=4
                      if(i_read_nonbond_functions) then
                        if(atm_atm_typ(kkk,lll,0).ne.0)
     &                    i_status(jjj,iii)=8
                      endif
                    endif
                    g_status(jjj,iii)=m
                    ncxy=ncxy+1
                  endif
                enddo
              enddo
            endif

C now go through and change necessary atoms to i_status=2,5 or 6
C note that the form is unusual here - first, decide whether to change
C to i_status=5; then decide whether to from 1->2 or 5->6 depending on
C whether either of the two atoms is uncharged

            do jjj=1,num_flx_atms
              if(i_f_a(iii).ne.1.or.i_f_a(jjj).ne.1) then
                if(i_status(jjj,iii).ne.0) 
     &             i_status(jjj,iii)=i_status(jjj,iii)+1
              endif
            enddo
          enddo
          write(*,*)'finished determining atom-pair status',ncxy,nbxy
        endif

C regardless of whether i_use_high_mem is true, deallocate b_typ

        write(*,*)'done indexing, deallocating b_typ array '
        do itp=1,num_f_ts
          deallocate(b_typ(itp)%id1)
          deallocate(b_typ(itp)%id2)
        enddo
        deallocate(b_typ)

        if(i_use_go_pairs) then
          allocate(num_cur_Q_loc(1:num_g_ms,1:num_g_ms,1:num_threads))
          allocate(num_cur_Q_mod(1:go_mod_max,
     &             1:num_g_ms,1:num_g_ms))
          allocate(num_cur_Q_mod_loc(1:go_mod_max,
     &             1:num_g_ms,1:num_g_ms,1:num_threads))
        endif

C allocate arrays for Fixman method if necessary
C note that we have to make num_hyd_atms3 be (:) so here we've just
C opted for assuming there's only one hydrodynamic system...

        if(fixman) then

          allocate(fixman_coeff(0:fixman_order))
          allocate(zme(0:fixman_order))
          allocate(fixman_ran(1:num_hyd_atms3(1)))
          allocate(fixman_xx(1:num_hyd_atms3(1),0:fixman_order))
          if(full_hydrodynamics) then
          allocate(r_a_r_mtrx(1:num_hyd_atms3(1),1:num_hyd_stp))               
          allocate(fixman_tensor(1:num_hyd_atms3(1),1:num_hyd_atms3(1)))
          allocate(fixman_ran_mtrx(1:num_hyd_atms3(1),1:num_hyd_stp))
          allocate(fixman_xx_mtrx(1:num_hyd_atms3(1),1:num_hyd_stp,
     &                            0:fixman_order))
          allocate(fixman_l2_num_mtrx(1:num_hyd_stp,1:num_threads))
          allocate(fixman_l2_den_mtrx(1:num_hyd_stp,1:num_threads))
          num_tmp_long=num_hyd_atms3(1) ! num_tmp_long is a long-integer...
          la(1)=num_tmp_long*(num_tmp_long+1)/2
          allocate(amkl(1:la(1)))
          allocate(dmkl(1:num_hyd_atms3(1)))
          allocate(emkl(1:num_hyd_atms3(1)))
          allocate(tmkl(1:num_hyd_atms3(1)))
          allocate(wmkl(1:num_hyd_atms3(1)))
          allocate(work(1:5*num_hyd_atms3(1)))
          endif

        else

C to use spotrf we need to make sure to allocate amkl for cholesky

          if(full_hydrodynamics) then
          num_tmp_long=num_hyd_atms3(1) ! num_tmp_long is a long-integer...
          la(1)=num_tmp_long*(num_tmp_long+1)/2
          allocate(amkl_new(1:la(1)))
          endif

        endif

C new for YHL code - also need to allocate space for imax_thrd_int and
C imin_thrd_int for each thread - these store the min and max atoms on
C the other threads that this thread needs to add forces to

        allocate(imin_thrd_int(1:num_threads,1:num_threads))
        allocate(imax_thrd_int(1:num_threads,1:num_threads))

C allocate some stuff that used to be needed for replica_exchange
C the following lines were from "if not replica_exchange" so keep

        num_rex_cnt=-999
        allocate(irep_go(1:num_f_ts,1:num_f_ts))
        allocate(rep_go_scale(1:num_reps))
        irep_go=0
        rep_go_scale=1.0

        write(*,*)
        write(*,*)'Hold on to your hats, here we go!'
        write(*,*)

        write(*,*)'opening up the .xtc file'
        write(*,*)
        allocate(x_xtc(1:num_flx_atms*3))
        allocate(xtc_file_out(1:num_reps))

        xtc_file_out='testout.XXX.xtc'
        do i=1,num_reps
          write(atemp,'(i3)')i
          do j=1,3
            if(atemp(j:j).eq.' ') atemp(j:j)='0'
          enddo
          xtc_file_out(i)(9:11)=atemp
        enddo

        num_xtc=0

C allocate memory for timings

        allocate(tmtrx(1:50,1:num_threads))
        tmtrx=0.0

        call system_clock(count=ktime1,count_rate=ru)

C if not using newton's law we can just use global force arrays
C (if we are using newton's law then we use private force arrays)

        if(.not.i_do_YHL_nonbondeds) then
          allocate(e_a_x_st(1:3,1:num_tot_atms,1:num_reps))
          allocate(e_a_x_md(1:3,1:num_tot_atms,1:num_reps))
          allocate(e_a_x_lg(1:3,1:num_tot_atms,1:num_reps))
          allocate(e_e_x_st(1:3,1:num_tot_atms,1:num_reps))
          allocate(e_e_x_md(1:3,1:num_tot_atms,1:num_reps))
        endif

C this is the point at which we start the actual simulation - we open up
C an 'openmp parallel' section and stay within it for the entire duration
C of the simulation so that there is no need to pay the overhead for
C continually starting up and closing down threads. this makes it
C somewhat difficult for me to call other people's routines that are themselves
C threaded - e.g. I have still not really figured out how to get the
C master thread to call a threaded MKL routine and have all the other
C openmp threads stand by so that the master+MKL can make use of their
C resources. but that's beside the point...

C now start!

!$omp   parallel default (shared) 
!$omp$  private(ibeg_bnd,ibeg_ang,ibeg_dih,ibeg_f_a,mythrd,
!$omp$          iend_bnd,iend_ang,iend_dih,iend_f_a,
!$omp&          ddd_ff_g_st,ddd_ff_g_md,
!
! mixed_hydrodynamics stuff
!
!$omp&          mov_num_loc,mov_typ_loc,mov_beg_loc,
!$omp&          mov_end_loc,mov_dyn_loc,mov_hyd_loc,
!$omp&          mov_beq_loc,mov_enq_loc,mov_lnc_loc,
!$omp&          mov_tot_loc,mov_tot_loc3,
!$omp&          m3,ncol,trm,jid,info,i_will_skip,
!
! replica_exchange stuff
!
!$omp&          myrep,myrep_master,myswap,
!$omp&          stream_rex,errcode_rex,seed_rex,
!$omp&          n_ff_g_rx,n_ff_g_rx_old,
!$omp&          ida_ff_g_rx,jda_ff_g_rx,kda_ff_g_rx,lda_ff_g_rx,
!$omp&          mda_ff_g_rx,s11_ff_g_rx,s22_ff_g_rx,s33_ff_g_rx,
!$omp&          s44_ff_g_rx,d11_ff_g_rx,d22_ff_g_rx,d12_ff_g_rx,
!$omp&          ddd_ff_g_rx,e00_ff_g_rx,
!$omp&          ene_g_rx_old_loc,ene_g_rx_new_loc,
!$omp&          rand_rex,boltz_rex,deltaE_rex,iaccept_rex,
!$omp&          c_a_x_rex,f_a_x_rex,mythrds_rep_temp,
!
! osmotic pressure stuff
!
!$omp&          osm_cur_loc,
!
! new rigid molecule stuff
!
!$omp&          num_rgd_flx_mov,typ_rgd_flx_mov,n_ff_r_st,
!$omp&          beg_rgd_flx_mov,end_rgd_flx_mov,
!$omp&          ida_ff_r_st,jda_ff_r_st,kda_ff_r_st,
!$omp&          cx_lab,cy_lab,cz_lab,
!$omp&          fx_tmp,fy_tmp,fz_tmp,
!$omp&          fx_mol,fy_mol,fz_mol,
!$omp&          tx_mol,ty_mol,tz_mol,
!$omp&          q_B1,q_B2,q_B3,q_B4,q_B5,q_B6,
!$omp&          q_D1,q_D2,q_D3,q_D4,q_D5,q_D6,
!$omp&          q_T1,q_T2,q_T3,q_T4,q_T5,q_T6,
!$omp&          omg2,omg2_inv,omg,omg_inv,cos_omg,sin_omg,
!$omp&          V11,V12,V13,V21,V22,V23,V31,V32,V33,
!$omp&          N11,N12,N13,N21,N22,N23,N31,N32,N33,
!$omp&          rigid_bin,rigid_bin_edge,
!$omp&          t11,t12,t21,t22,t31,t32,t41,t42,t51,t52,t61,t62,
!$omp&          u11,u12,u21,u22,u31,u32,u41,u42,u51,u52,u61,u62,
!$omp&          mymol1,mymol2,
!
! i_do_lincs stuff
!
!$omp&          n,na1,na2,plnc,rlen,nrec,lbnd_thrd,
!
! new timing stuff
!
!$omp&          ftm1,ftm2,
!
! new random stuff
!
!$omp&          num_randoms,random_trms,
!$omp&          seed,ran_mkl,errcode,stream,
!
! new bond, angle, dihedral thread handlers
!
!$omp$          num_bnd_thrd,ibnd_thrd,dbnd_thrd,
!$omp$          num_ang_thrd,iang_thrd,dang_thrd,
!$omp$          num_dih_thrd,idih_thrd,ddih_thrd,
!$omp&          iloc,jloc,kloc,lloc,
!
! debug bonds
!
!$omp$          rij_old,
!$omp&          ccx,ccy,ccz,ddx,ddy,ddz,eex,eey,eez,ffx,ffy,ffz,
!
! terms related to atoms avoiding nonbonded inters due to pos restraints
! note that rrr1 and rrr3 should have been made private years ago!
!
!$omp&          ibeg_i_a,iend_i_a,ii1,rrr1,rrr2,rrr3,
!
! arbitrary potentials
!
!$omp&          angle_deg,
!$omp&          n_ff_a_st,n_ff_a_md,                       
!$omp&          n_ff_b_st,n_ff_b_md,                       
!$omp&          ida_ff_a_st,ida_ff_a_md,ida_ff_b_st,ida_ff_b_md,
!$omp&          jda_ff_a_st,jda_ff_a_md,jda_ff_b_st,jda_ff_b_md,
!$omp&          kda_ff_a_st,kda_ff_a_md,kda_ff_b_st,kda_ff_b_md,
!$omp&          lda_ff_b_st,lda_ff_b_md,
!
! histogram terms
!
!$omp&          ihist,jhist,ang1,ang2,num_cur_pairs,
!
! protease terms    
! 
!$omp&          jj1,jj1_match,jj2_match,i_cleave_now,unif,
!$omp&          jj1_loc,jj2_loc,factor,fraction_base,charge,
!$omp&          n_ff_p_st,ida_ff_p_st,jda_ff_p_st,
!$omp&          jdb_ff_p_st,jds_ff_p_st,crt_ff_p_st,
!$omp&          rng_ff_p_st,fct_ff_p_st,
!$omp&          jdb,jds,crt,rng,fct,
!$omp&          f_a_x_pt,
!
! i_compare_go_with_others terms
! 
!$omp&          y2,y3,xdist2,ikk,
!
! afm and confinement terms
!
! exclusive go terms
! 
!$omp&          mno,
!
! afm and confinement terms
!
!$omp&          rad,rad2,
!
! go_spline terms
!
!$omp&          klo,khi,abc,bcd,cde,
!
! brute-force N^2 elec terms
!
!$omp&          u1,u2,u3,v2,             
!
! fixman terms
!
!$omp&          help,dwDdw_loc,tmp_trm,
!
! ewald_hydro terms
!
!$omp&          aaa,aaa3,sm1,sm2,bracket,ekx,eky,ekz,
!$omp&          tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10,
!$omp&          tm11,tm12,tm13,tm14,tm15,
!$omp&          crap1m,crap2m,crap1p,crap2p,
!$omp&          dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,dx9,
!$omp&          ex1,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ex9,
!
! ewald_elec terms
!
!$omp&          rkx,rky,rkz,akvec,dotik,e_e,
!$omp&          ene_ewld_real_loc,ene_ewld_rcip_loc,
!$omp&          ibeg_kvec,iend_kvec, ! no longer used
!
! ewald_elec_grid & ewald_hydro_grid terms
!
!$omp&          ijk,jkl,klm,lmn,lmn16a,lmn16b,lmn24a,lmn24b,g1,g2,g3,
!$omp&          dx,dy,dz,dxdy,dxdz,dydz,dxdydz,
!$omp&          dy_dydz,dz_dydz,dx_dxdy,dz_dxdz,
!$omp&          u000,u010,u001,u011,u100,u110,u101,u111,
!$omp&          o1,o2,o3,e1,e2,e3,
!$omp&          u_g,u_p,f_g1,f_g2,f_g3,
!
! multi-scale HI terms 
!
!$omp&          sum_ri,sum_ri_rj,ftmp,
!$omp&          ibeg_f_m,iend_f_m,i4,i5,i6,  ! molecules treated by each thread
!$omp&          ib1,ib2,ib3,ib4,lcur,dmin,nnd,
!$omp&          ii4,jj2,jj4,
!$omp&          i_molmol,i_blbblb,
!$omp&          num_molmol,idm_molmol,jdm_molmol,
!$omp&          num_blbblb,idb_blbblb,jdb_blbblb,
!$omp&          num_atmatm,ida_atmatm,jda_atmatm,
!
! 44-07-12 i_do_pH stuff here 
!
!$omp&          n_ff_i_st,ida_ff_i_st,jda_ff_i_st,i_got_sorted2,
!$omp&          i_got_sorted,iip_acc,iii_acc,num_acc,iaccept,
!$omp&          q_f_old_acc,q_f_new_acc,qjf_old,qjf_new,
!$omp&          iiq,jjq,jk2,jk3,ene_old,ene_new,
!$omp&          nfa_iii_old,nfa_iii_new,q_f_iii_old,q_f_iii_new,
!$omp&          nfa_jjj_old,nfa_jjj_new,q_f_jjj_old,q_f_jjj_new,
!
! 41-07-12 user_hists stuff here 
!
!$omp&          i_met_condition,user_bin,user_bin_edge,frac_bin_edge,
!
! 42-07-12 user_energy terms
!
!$omp&          nii,njj,mytmp1,mytmp2,mytmp3,mytmp4,rrr,
!$omp&          n_ff_u_st,mii,mjj,mkk,tmp_f,tmp_e,tmp_u,tmp_v,
!$omp&          user_lo,user_hi,tmp_ene,tmp_lo,tmp_hi,
!$omp&          valu_bin,frac_bin,
!$omp&          ida_ff_u_st,jda_ff_u_st,kda_ff_u_st,
!$omp&          k1,k2,k3,k4,k5,k6,k7,
!
! deallocate code
!
!$omp&          n_ff_p_st_old,n_ff_u_st_old,n_ff_a_st_old,
!$omp&          n_ff_a_md_old,n_ff_b_st_old,n_ff_b_md_old,
!$omp&          n_ff_v_st_old,n_ff_v_md_old,n_ff_g_st_old,
!$omp&          n_ff_g_md_old,n_ff_e_st_old,n_ff_e_md_old,
!$omp&          n_ff_i_st_old,n_ff_r_st_old,n_ff_e_lg_old,
!
! rods stuff
!
!$omp&          n_ff_rods,n_ff_rods_old,
!$omp&          ida_ff_rods,jda_ff_rods,
!$omp&          s11_ff_rods,s22_ff_rods,s33_ff_rods,s44_ff_rods,
!$omp&          d11_ff_rods,d22_ff_rods,e00_ff_rods,
!
! treecode_elec terms
! HANS
! here are openmp 'private' terms relating to your treecode
!
!$omp&          ibeg_q_a,iend_q_a,
!$omp&          c_q_x_cpy,c_q_y_cpy,c_q_z_cpy,
!$omp&          crg_of_atm_cpy,perm_cpy,
!
! 44-07-12 removing private declaration of these force arrays
! YHL code - reimposing private conditions 
!
!$omp&          f_e_x_st,f_e_x_md,
!
!$omp&          ibeg_f_a3,iend_f_a3,ii3,jj3,kk3,ik3,
!$omp&          ibeg_r_a,iend_r_a,ibeg_r_m,iend_r_m,mtp,ntp,idis,
!$omp&          ibeg_h_a,iend_h_a,ibeg_h_a3,iend_h_a3,
!$omp&          nloc_h_a,nloc_h_a3,
!$omp&          ibeg_l_a,iend_l_a, ! used for cholesky decomposition
!$omp&          dxp,dxm,xf,xm,dm,vm,cmov,cfix,tjac,qjac,sumf,summ,
!$omp&          qcac,xloc,tloc,imov,nmrot,idiedinjacobi,stmp,
!$omp&          iii,ii2,jjj,xij1,xij2,xij3,rij,fbond,
!$omp&          mynum1,mynum2,mynum3,mynum4,i_status_tmp,
!$omp&          myid,myjd,
!$omp&          i1_loc,i2_loc,i3_loc,j1_loc,j2_loc,j3_loc,
!$omp&          i1,i2,i3,j1,j2,j3,q1,q2,l1,l2,
!$omp&          yij1,yij2,yij3,
!$omp&          it1,it2,it3,ic1,ic2,ic3,kkk,lll,mmm,nnn,
!$omp&          i_am_not_go,itp,jtp,ktp,
!$omp&          s1a,s2a,s3a,s4a,d1a,d2a,d3a,e0a,f_v,u_v,
!$omp&          dist,dist_inv,dist2,dist2_inv,qqq,temp,
!$omp&          n_ff_v_st,n_ff_v_md,                       
!$omp&          n_ff_g_st,n_ff_g_md,                       
!$omp&          n_ff_e_st,n_ff_e_md,n_ff_e_lg,
!$omp&          ida_ff_v_st,ida_ff_v_md,ida_ff_g_st,ida_ff_g_md,
!$omp&          ida_ff_e_st,ida_ff_e_md,ida_ff_e_lg,
!$omp&          jda_ff_v_st,jda_ff_v_md,jda_ff_g_st,jda_ff_g_md,
!$omp&          jda_ff_e_st,jda_ff_e_md,jda_ff_e_lg,
!$omp&          mda_ff_v_st,mda_ff_v_md,                         
!$omp&          mda_ff_e_st,mda_ff_e_md,mda_ff_e_lg,
!$omp&          nda_ff_v_st,nda_ff_v_md,                         
!$omp&          nda_ff_e_st,nda_ff_e_md,nda_ff_e_lg,
!$omp&          mda_ff_g_st,mda_ff_g_md,
!$omp&          nda_ff_g_st,nda_ff_g_md,
!$omp&          kda_ff_g_st,lda_ff_g_st,kda_ff_g_md,lda_ff_g_md,
!$omp&          s11_ff_v_st,s11_ff_v_md,s11_ff_g_st,s11_ff_g_md,
!$omp&          s22_ff_v_st,s22_ff_v_md,s22_ff_g_st,s22_ff_g_md,
!$omp&          s33_ff_v_st,s33_ff_v_md,s33_ff_g_st,s33_ff_g_md,
!$omp&          s44_ff_v_st,s44_ff_v_md,s44_ff_g_st,s44_ff_g_md,
!$omp&          d11_ff_v_st,d11_ff_v_md,d11_ff_g_st,d11_ff_g_md,
!$omp&          e00_ff_v_st,e00_ff_v_md,e00_ff_g_st,e00_ff_g_md,
!$omp&          d22_ff_v_st,d22_ff_v_md,rep_ff_g_st,rep_ff_g_md,
!$omp&          d22_ff_g_st,d22_ff_g_md,d12_ff_g_st,d12_ff_g_md,
!$omp&          qqq_ff_e_st,qqq_ff_e_md,qqq_ff_e_lg,
!$omp&          time_f_move_loc,
!$omp&          x1,y1,z1,
!$omp&          f_a_x_st,f_a_x_md,f_a_x_lg,f_a_x_bd,
!$omp&          xia,yia,zia,xib,yib,zib,xic,yic,zic,xid,yid,zid,
!$omp&          xab,yab,zab,xcb,ycb,zcb,rab2,rcb2,xp1,yp1,zp1,rp,
!$omp&          dot,cosine,angle,dt,dt2,deddt,terma,termc,                   
!$omp&          dedxia,dedyia,dedzia,
!$omp&          dedxib,dedyib,dedzib,
!$omp&          dedxic,dedyic,dedzic,
!$omp&          dedxid,dedyid,dedzid,
!$omp&          xba,yba,zba,xdc,ydc,zdc,force,
!$omp&          xt,yt,zt,xu,yu,zu,xtu,ytu,ztu,rt2,ru2,rtru,rcb,
!$omp&          cosine1,sine1,v1,c1,s1,v3,c3,s3,
!$omp&          cosine2,sine2,cosine3,sine3,phi1,phi3,dphi1,dphi3,
!$omp&          dedphi,xxca,yyca,zzca,xdb,ydb,zdb,
!$omp&          dedxt,dedyt,dedzt,dedxu,dedyu,dedzu,
!$omp&          ene_bond_loc,ene_angl_loc,ene_dihe_loc,
!$omp&          ene_pos_loc,
!$omp&          ene_nbnd_loc,ene_gofav_loc,ene_gounfav_loc,
!
! 32-07-12 - local versions of energies for replica-exchange sims
!
!$omp&          ene_v_st_loc,ene_g_st_loc,ene_e_st_loc,
!$omp&          ene_v_md_loc,ene_g_md_loc,ene_e_md_loc,
!$omp&          ene_e_lg_loc,
!$omp&          i,j,k,l,m,w1,w2,w3,d1,d2,d3,rrijsq,
!$omp&          crap0,crap1,crap2,crap3,
!$omp&          gaus1,gaus2,gaus3,sumtx,sumty,sumtz,
!$omp&          sn1_x,sn1_y,sn1_z,
!$omp&          ida,jda,tmp1,tmp2,tmp3,
!$omp&          ida_h_d_ngb_tot,jda_h_d_ngb_tot,
!$omp&          rij_h_d_ngb_tot,beta_tp,
!$omp&          d_a_g_x,d_a_g_y,d_a_g_z,
!$omp&          i_a_g_x,i_a_g_y,i_a_g_z,
!$omp&          so1_x,so1_y,so1_z,si1_x,si1_y,si1_z)

C xxxxx Section
C xxxxx Evaluate all thread-specific counters

        mythrd=omp_get_thread_num()+1
        m=mythrd

C identify which replica this thread works on - note that we only need
C to do this if we're using mixed_hydrodynamics - we also, at the same
C time set myrep_master according to whether this is the first thread to
C deal with a particular replica

        myrep_master=.false.
        myrep=1
        if(mythrd.eq.1) myrep_master=.true.

        if(walls.or.sphere.or.cylinder) allocate(osm_cur_loc(1:7))

C for backwards compatibility get ibeg_f_a and iend_f_a
C this will probably get deleted eventually but DON'T DELETE IT YET!

C 2019 - since ibeg_f_a and iend_f_a appear to be almost completely
C deprecated (they still get used for umbrella calculations but that
C seems like old code), we can probably reuse them here to help divide
C up the moving atoms amongst the threads properly - all we have to do
C is change num_tot_atoms to num_mov_atoms below...

C ******************************************************************
C note that for this all to work properly we probably should require
C that all moving atoms appear before all static atoms...
C ******************************************************************

c       rrr1=real(num_tot_atms)/real(num_threads)
        rrr1=real(num_mov_atms)/real(num_threads)
        ibeg_f_a=1+int((mythrd-1)*rrr1)
        if(mythrd.eq.num_threads) then
c         iend_f_a=num_tot_atms
          iend_f_a=num_mov_atms
        else
          iend_f_a=int(mythrd)*rrr1
        endif
        ibeg_f_a3=ibeg_f_a*3-2                       
        iend_f_a3=iend_f_a*3                       

C now set up other thread local arrays for mixed_hydrodynamics
C note that we only accept a move if it's the correct node, thread and
C replica number (mov_nod, mov_thr, mov_rep)

C   num_rgd_flx_mov_tot is the total number of entries in the file
C   listing all various moves schedule for each thread

C   mov_typ_loc is the type of move (flex or rigid etc)
C   mov_beg_loc is the first atom or molecule of the move
C   mov_end_loc is the last atom or molecule of the move
C   mov_dyn_loc is the type of dynamics (BD or LD)
C   mov_hyd_loc is the type of hydrodynamics (full or none)
C   mov_lnc_loc is whether lincs is being used
C   mov_beq_loc is the first charge of the move
C   mov_enq_loc is the last charge of the move

        if(full_hydrodynamics) then

          mov_num_loc=1

          allocate(mov_typ_loc(1:mov_num_loc))
          allocate(mov_beg_loc(1:mov_num_loc))
          allocate(mov_end_loc(1:mov_num_loc))
          allocate(mov_dyn_loc(1:mov_num_loc))
          allocate(mov_hyd_loc(1:mov_num_loc))
          allocate(mov_lnc_loc(1:mov_num_loc))
          allocate(mov_beq_loc(1:mov_num_loc))
          allocate(mov_enq_loc(1:mov_num_loc))

          mov_typ_loc(1)=1
          mov_beg_loc(1)=ibeg_f_a
          mov_end_loc(1)=iend_f_a
          mov_hyd_loc(1)=1
          mov_lnc_loc(1)=0 ! no lincs if full HI
          if(brownian) mov_dyn_loc(1)=1
          if(langevin) mov_dyn_loc(1)=2

        elseif(no_hydrodynamics) then

          mov_num_loc=1

          allocate(mov_typ_loc(1:mov_num_loc))
          allocate(mov_beg_loc(1:mov_num_loc))
          allocate(mov_end_loc(1:mov_num_loc))
          allocate(mov_dyn_loc(1:mov_num_loc))
          allocate(mov_hyd_loc(1:mov_num_loc))
          allocate(mov_lnc_loc(1:mov_num_loc))
          allocate(mov_beq_loc(1:mov_num_loc))
          allocate(mov_enq_loc(1:mov_num_loc))

          mov_typ_loc(1)=1
          mov_beg_loc(1)=ibeg_f_a
          mov_end_loc(1)=iend_f_a
          mov_hyd_loc(1)=0
          if(i_do_lincs) then
            mov_lnc_loc(1)=1
          else
            mov_lnc_loc(1)=0
          endif
          if(brownian) mov_dyn_loc(1)=1
          if(langevin) mov_dyn_loc(1)=2

        endif

        write(*,*)'for thread # ',mythrd,' mov_num_loc = ',mov_num_loc
!$omp barrier ! purely cosmetic

!$omp single ! purely cosmetic
        write(*,580)
580     format(36x,'    node    thrd   round  movtyp  movbeg  movend ',
     &             ' movdyn  movhyd ')
!$omp end single ! purely cosmetic

        do nnn=1,mov_num_loc
          write(*,581)mynode,mythrd,nnn,
     &          mov_typ_loc(nnn),mov_beg_loc(nnn),mov_end_loc(nnn),
     &          mov_dyn_loc(nnn),mov_hyd_loc(nnn),mov_lnc_loc(nnn)
581       format('check types of moves on each thread ',9i8)
        enddo
!$omp barrier ! purely cosmetic

C now figure out how many random terms we need - since this is quite
C complicated we will, for now, just get random terms for *every* atom
C that is held by this thread - mov_tot_loc replaces nloc_h_a

        mov_tot_loc=0
        do nnn=1,mov_num_loc
          do iii=mov_beg_loc(nnn),mov_end_loc(nnn)
            mov_tot_loc=mov_tot_loc+1
          enddo
        enddo
        mov_tot_loc3=mov_tot_loc*3
        num_randoms=mov_tot_loc3*num_hyd_stp

        write(*,*)
        write(*,*)'thread # ',mythrd,' mov_tot_loc = ',mov_tot_loc
        write(*,*)'thread # ',mythrd,' generates ',num_randoms,' rnds'
        write(*,*)

        allocate(random_trms(1:num_randoms))

C allocate some private array stuff for hydrodynamic treecode
C idnode stores (temporarily) the
C identities of the bodies and cells that interact with the current atom
C - they are immediately transferred into the derived type hi_tree

C now calculate all the nonbonded terms again
C 
C note that we could save memory by making thread-local arrays smaller

        allocate(xloc(1:3))
        allocate(tloc(1:3,1:3))
        allocate(dxp(1:max_h_as,1:3))
        allocate(dxm(1:max_h_as,1:3))
        allocate(xf(1:max_h_as,1:3))
        allocate(xm(1:max_h_as,1:3))
        allocate(ran_mkl(1:3*num_tot_atms))
        allocate(i1_loc(1:num_tot_atms))
        allocate(i2_loc(1:num_tot_atms))
        allocate(i3_loc(1:num_tot_atms))
        allocate(j1_loc(1:num_tot_atms))
        allocate(j2_loc(1:num_tot_atms))
        allocate(j3_loc(1:num_tot_atms))

        write(*,*)'Allocating local nonbonded list arrays ',num_ff_lim1

C deallocate code

        n_ff_i_st_old=0
        n_ff_p_st_old=0
        n_ff_u_st_old=0
        n_ff_r_st_old=0
        n_ff_a_st_old=0
        n_ff_a_md_old=0
        n_ff_b_st_old=0
        n_ff_b_md_old=0
        n_ff_v_st_old=0
        n_ff_v_md_old=0
        n_ff_g_st_old=0
        n_ff_g_md_old=0
        n_ff_e_st_old=0
        n_ff_e_md_old=0
        n_ff_e_lg_old=0
        n_ff_g_rx_old=0
        n_ff_rods_old=0

        allocate(ida_ff_r_st(1:1))
        allocate(jda_ff_r_st(1:1))
        allocate(kda_ff_r_st(1:1))

C 42-07-12

        allocate(ida_ff_u_st(1:1))
        allocate(jda_ff_u_st(1:1))
        allocate(kda_ff_u_st(1:1))

        allocate(ida_ff_a_st(1:1)) ! atm i
        allocate(jda_ff_a_st(1:1)) ! atm j
        allocate(kda_ff_a_st(1:1)) ! pair type
        allocate(ida_ff_a_md(1:1)) ! atm i
        allocate(jda_ff_a_md(1:1)) ! atm j
        allocate(kda_ff_a_md(1:1)) ! pair type
        allocate(ida_ff_b_st(1:1)) ! atm i
        allocate(jda_ff_b_st(1:1)) ! atm j
        allocate(kda_ff_b_st(1:1)) ! pair type
        allocate(lda_ff_b_st(1:1)) ! go dist
        allocate(ida_ff_b_md(1:1)) ! atm i
        allocate(jda_ff_b_md(1:1)) ! atm j
        allocate(kda_ff_b_md(1:1)) ! pair type
        allocate(lda_ff_b_md(1:1)) ! go dist  

        allocate(ida_ff_i_st(1:1))
        allocate(jda_ff_i_st(1:1))
        allocate(i_got_sorted2(1:1))

        allocate(ida_ff_v_st(1:1))
        allocate(jda_ff_v_st(1:1))
        allocate(s11_ff_v_st(1:1))
        allocate(s22_ff_v_st(1:1))
        allocate(s33_ff_v_st(1:1))
        allocate(s44_ff_v_st(1:1))
        allocate(d11_ff_v_st(1:1))
        allocate(d22_ff_v_st(1:1))
        allocate(e00_ff_v_st(1:1))

        allocate(ida_ff_g_rx(1:1))
        allocate(jda_ff_g_rx(1:1))
        allocate(kda_ff_g_rx(1:1))
        allocate(lda_ff_g_rx(1:1))
        allocate(mda_ff_g_rx(1:1))
        allocate(s11_ff_g_rx(1:1))
        allocate(s22_ff_g_rx(1:1))
        allocate(s33_ff_g_rx(1:1))
        allocate(s44_ff_g_rx(1:1))
        allocate(d11_ff_g_rx(1:1))
        allocate(d22_ff_g_rx(1:1))
        allocate(d12_ff_g_rx(1:1))
        allocate(ddd_ff_g_rx(1:1))
        allocate(e00_ff_g_rx(1:1))

        allocate(ida_ff_v_md(1:1))
        allocate(ida_ff_g_st(1:1))
        allocate(ida_ff_g_md(1:1))
        allocate(ida_ff_e_st(1:1))
        allocate(ida_ff_e_md(1:1))
        allocate(jda_ff_v_md(1:1))
        allocate(jda_ff_g_st(1:1))
        allocate(jda_ff_g_md(1:1))
        allocate(jda_ff_e_st(1:1))
        allocate(jda_ff_e_md(1:1))
        allocate(kda_ff_g_st(1:1))
        allocate(kda_ff_g_md(1:1))
        allocate(lda_ff_g_st(1:1))
        allocate(lda_ff_g_md(1:1))
        allocate(mda_ff_g_st(1:1))
        allocate(mda_ff_g_md(1:1))
        allocate(nda_ff_g_st(1:1))
        allocate(nda_ff_g_md(1:1))
        allocate(s11_ff_v_md(1:1))
        allocate(s22_ff_v_md(1:1))
        allocate(s33_ff_v_md(1:1))
        allocate(s44_ff_v_md(1:1))
        allocate(d11_ff_v_md(1:1))
        allocate(d22_ff_v_md(1:1))
        allocate(e00_ff_v_md(1:1))
        allocate(s11_ff_g_st(1:1))
        allocate(s22_ff_g_st(1:1))
        allocate(s33_ff_g_st(1:1))
        allocate(s44_ff_g_st(1:1))
        allocate(d11_ff_g_st(1:1))
        allocate(d22_ff_g_st(1:1))
        allocate(d12_ff_g_st(1:1))
        allocate(ddd_ff_g_st(1:1))
        allocate(e00_ff_g_st(1:1))
        allocate(rep_ff_g_st(1:1))
        allocate(s11_ff_g_md(1:1))
        allocate(s22_ff_g_md(1:1))
        allocate(s33_ff_g_md(1:1))
        allocate(s44_ff_g_md(1:1))
        allocate(d11_ff_g_md(1:1))
        allocate(d22_ff_g_md(1:1))
        allocate(d12_ff_g_md(1:1))
        allocate(ddd_ff_g_md(1:1))
        allocate(e00_ff_g_md(1:1))
        allocate(rep_ff_g_md(1:1))
        allocate(qqq_ff_e_st(1:1))
        allocate(qqq_ff_e_md(1:1))

C we use the treecode_elec to calculate all electrostatic interactions
C at intervals of 'num_lst_stp' simulation steps. we use a separate loop to
C calculate short-range electrostatic interactions *every* step, and a
C separate loop to calculate medium-range electrostatic interactions
C every 'num_fmd_stp' simulation steps - the short-range and medium-range
C forces are computed directly (i.e. N^2 dependence) and substracted
C from the treecode forces to give long-range forces that are updated
C only intermittently (i.e. every 'num_lst_stp' steps). In other words,
C the treecode_elec is used to provide only the long-range
C electrostatics...

C replica exchange comment - note that since all of the following are
C private arrays I don't think I need to add the num_reps dimension:
C this would only be needed if each thread worked on >1 replica...

        if(i_do_YHL_nonbondeds) then
          allocate(f_a_x_st(1:3,1:num_tot_atms))
          allocate(f_a_x_md(1:3,1:num_tot_atms))
          allocate(f_a_x_lg(1:3,1:num_tot_atms))
          if(treecode_elec) then
            allocate(f_e_x_st(1:3,1:num_tot_atms))
            allocate(f_e_x_md(1:3,1:num_tot_atms))
          endif
        endif

C note that f_a_x_bd is a private array allocated for all atoms but we
C only ever use those entries that correspond to thread-local atoms

        allocate(f_a_x_bd(1:3,1:num_tot_atms))

        allocate(ene_pos_loc(1:num_f_ms))
        allocate(ene_bond_loc(1:num_f_ms))
        allocate(ene_angl_loc(1:num_f_ms))
        allocate(ene_dihe_loc(1:num_f_ms))

        time_f_move_loc=0.0
 
C initialize the random number generator separately for each thread
C note that we currently try to get the same numbers on each thread
C - a better way to do this is make everything NOT private - i.e.
C allocate ran_mkl before the big openmp loop
 
c       seed=iseed_new
        seed=iseed_new+mythrd
        write(*,*)
        write(*,*)'CHECK THREAD + SEED ',mythrd,seed
        write(*,*)
        errcode=vslnewstream(stream,brng,seed)

C if treecode_elec, find the first and last charged atoms that are being
C treated by this particular thread: these are 'ibeg_q_a' and 'iend_q_a'
C 'i_f_a(iii)' is a flag for each atom 'iii' that = 1 if the atom has a
C charge, and = 0 if it is uncharged (we omit the latter from the tree)

C note that this must use all atoms - even ones that are not feeling
C nonbonded interactions, otherwise the ones that are feeling nonbonded
C interactions will not see the others' charges
C note however that we do *not* include non-free atoms here

C AHE is not currently sure why we have commented out the if statement
C for treecode_elec - we'll continue with it as is for now but it would
C be nice to get to the bottom of why this is used for all cases

c       if(treecode_elec) then

C used to find 'ibeg_q_a' and 'iend_q_a' for this thread but now we're
C replacing them with mov_beq_loc and mov_enq_loc...

C first set initial values for them so that we don't
C get crazy loops from 0 to 0 later in the program... also note that we
C loop over all mov_num_loc rounds done on this thread

          do nnn=1,mov_num_loc
            mov_beq_loc(nnn)=1
            do iii=mov_beg_loc(nnn),mov_end_loc(nnn)
              if(i_f_a(iii).eq.1.and.ifree(iii).eq.1) then
                mov_beq_loc(nnn)=i_q_a(iii)
                exit
              endif
            enddo
          enddo
          do nnn=1,mov_num_loc
            mov_enq_loc(nnn)=0
            do iii=mov_end_loc(nnn),mov_beg_loc(nnn),-1
              if(i_f_a(iii).eq.1.and.ifree(iii).eq.1) then
                mov_enq_loc(nnn)=i_q_a(iii)
                exit
              endif
            enddo
          enddo

c       endif
  
!$omp barrier

        write(*,*)
        write(*,*)'atoms on threads:'
        write(*,*)

        write(*,176)
176     format('    node    thrd    mov     ')
        do nnn=1,mov_num_loc
          write(*,177)mynode,mythrd,nnn,
     &                mov_beg_loc(nnn),mov_end_loc(nnn),
     &                mov_beq_loc(nnn),mov_enq_loc(nnn)
        enddo
177     format(8i8)
        write(*,*)

C new for YHL code - need to record which thread each atom belongs to:
C also need to quit if an atom remains unassigned - this is a global

        do nnn=1,mov_num_loc
          do iii=mov_beg_loc(nnn),mov_end_loc(nnn)
            thrd_of_atm(iii,myrep)=mythrd
          enddo
        enddo
!$omp barrier

C 2019 skip this check as it'll definitely fail for static atoms -
C these, by definition, will not have been assigned to a thread above -
C this means (at least for now) that we can't use i_do_YHL since this
C seems to rely on having thrd_of_atm assigned in make_nonbonded_list.f

      if(.not.i_dont_move_some) then 
!$omp single 
        do jjj=1,num_reps
          do iii=1,num_tot_atms
            if(thrd_of_atm(iii,jjj).eq.0) then
              write(*,*)'atom ',iii,' not assigned to thread'
              write(*,*)' quitting :('
              stop
            endif
          enddo
        enddo
!$omp end single
      endif

C treecode_elec private allocations here

        if(treecode_elec) then
          allocate(c_q_x_cpy(1:num_q_as))  
          allocate(c_q_y_cpy(1:num_q_as))  
          allocate(c_q_z_cpy(1:num_q_as))  
          allocate(crg_of_atm_cpy(1:num_q_as))  
          allocate(perm_cpy(1:num_q_as))  
        endif

C check all bonds looking for ones that are relevant to the
C thread-local list of atoms - keep note of the first & last bond that
C is relevant to the thread-local list

        num_bnd_thrd=0
        do nnn=1,mov_num_loc
          do iii=1,num_sys_bonds
            iloc=0
            jloc=0
            if(nbb1_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbb1_new(iii).le.mov_end_loc(nnn)) iloc=1
            if(nbb2_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbb2_new(iii).le.mov_end_loc(nnn)) jloc=1
            if(iloc.eq.1.or.jloc.eq.1) then
              num_bnd_thrd=num_bnd_thrd+1
            endif
          enddo
        enddo

C allocate memory for the bonds on this thread

C ibnd_thrd : gives the id of this bond from 1:num_sys_bonds
C lbnd_thrd : denotes whether this is a lincs bond - this will be used
C             to skip non-lincs bonds - eventually we'll replace this
C dbnd_thrd : is the factor used to scale energy by so that we don't add
C             up the energy wrongly when doing same bond on multiple threads

        allocate(ibnd_thrd(1:num_bnd_thrd))
        allocate(lbnd_thrd(1:num_bnd_thrd))
        allocate(dbnd_thrd(1:num_bnd_thrd))
        write(*,*)'# of bonds on this local thread = ',
     &             mythrd,num_bnd_thrd

        lbnd_thrd=0 ! default to assuming not a lincs bond

C now put in all the information for this thread

        num_bnd_thrd=0
        do nnn=1,mov_num_loc
          do iii=1,num_sys_bonds
            iloc=0
            jloc=0
            if(nbb1_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbb1_new(iii).le.mov_end_loc(nnn)) iloc=1
            if(nbb2_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbb2_new(iii).le.mov_end_loc(nnn)) jloc=1
            if(iloc.eq.1.and.jloc.eq.1) then
              num_bnd_thrd=num_bnd_thrd+1
              ibnd_thrd(num_bnd_thrd)=iii
              dbnd_thrd(num_bnd_thrd)=1.0
              if(mov_hyd_loc(nnn).eq.0.and.mov_lnc_loc(nnn).eq.1)
     &          lbnd_thrd(num_bnd_thrd)=1
            elseif(iloc.eq.1.and.jloc.eq.0.or.
     &             iloc.eq.0.and.jloc.eq.1) then
              num_bnd_thrd=num_bnd_thrd+1
              ibnd_thrd(num_bnd_thrd)=iii
              dbnd_thrd(num_bnd_thrd)=1.0/2.0

C if it's a bond that's spread over 2 threads then we just arbitrarily
C say that it's done by the thread for which iloc=1 and jloc=0

              if(mov_hyd_loc(nnn).eq.0.and.mov_lnc_loc(nnn).eq.1
     &          .and.iloc.eq.1.and.jloc.eq.0)
     &          lbnd_thrd(num_bnd_thrd)=2
            endif
          enddo
        enddo

C do the same thing for angles 

        num_ang_thrd=0
        do nnn=1,mov_num_loc
          do iii=1,num_sys_angls
            iloc=0
            jloc=0
            kloc=0
            if(nba1_new(iii).ge.mov_beg_loc(nnn).and.
     &         nba1_new(iii).le.mov_end_loc(nnn)) iloc=1
            if(nba2_new(iii).ge.mov_beg_loc(nnn).and.
     &         nba2_new(iii).le.mov_end_loc(nnn)) jloc=1
            if(nba3_new(iii).ge.mov_beg_loc(nnn).and.
     &         nba3_new(iii).le.mov_end_loc(nnn)) kloc=1
            if(iloc.eq.1.or.jloc.eq.1.or.kloc.eq.1) then
              num_ang_thrd=num_ang_thrd+1
            endif
          enddo
        enddo

C allocate memory for the bonds on this thread
C iang_thrd : gives the id of this angle from 1:num_sys_angls
C dang_thrd : is the factor used to scale energy by so that we don't add
C up the energy incorrectly when doing same angle on multiple threads

        allocate(iang_thrd(1:num_ang_thrd))
        allocate(dang_thrd(1:num_ang_thrd))
        write(*,*)'# of angls on this local thread = ',
     &             mythrd,num_ang_thrd

C now put in all the information for this thread

        num_ang_thrd=0
        do nnn=1,mov_num_loc
          do iii=1,num_sys_angls
            iloc=0
            jloc=0
            kloc=0
            if(nba1_new(iii).ge.mov_beg_loc(nnn).and.
     &         nba1_new(iii).le.mov_end_loc(nnn)) iloc=1
            if(nba2_new(iii).ge.mov_beg_loc(nnn).and.
     &         nba2_new(iii).le.mov_end_loc(nnn)) jloc=1
            if(nba3_new(iii).ge.mov_beg_loc(nnn).and.
     &         nba3_new(iii).le.mov_end_loc(nnn)) kloc=1
            if(iloc.eq.1.and.jloc.eq.1.and.kloc.eq.1) then
              num_ang_thrd=num_ang_thrd+1
              iang_thrd(num_ang_thrd)=iii
              dang_thrd(num_ang_thrd)=1.0
            elseif((iloc.eq.1.and.jloc.eq.0.and.kloc.eq.0).or.
     &             (iloc.eq.0.and.jloc.eq.1.and.kloc.eq.0).or.
     &             (iloc.eq.0.and.jloc.eq.0.and.kloc.eq.1)) then
              num_ang_thrd=num_ang_thrd+1
              iang_thrd(num_ang_thrd)=iii
              dang_thrd(num_ang_thrd)=1.0/3.0
            elseif((iloc.eq.1.and.jloc.eq.1.and.kloc.eq.0).or.
     &             (iloc.eq.1.and.jloc.eq.0.and.kloc.eq.1).or.
     &             (iloc.eq.0.and.jloc.eq.1.and.kloc.eq.1)) then
              num_ang_thrd=num_ang_thrd+1
              iang_thrd(num_ang_thrd)=iii
              dang_thrd(num_ang_thrd)=2.0/3.0
            endif
          enddo
        enddo

C do the same thing for dihedrals

        num_dih_thrd=0
        do nnn=1,mov_num_loc
          do iii=1,num_sys_dihes
            iloc=0
            jloc=0
            kloc=0
            lloc=0
            if(nbd1_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbd1_new(iii).le.mov_end_loc(nnn)) iloc=1
            if(nbd2_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbd2_new(iii).le.mov_end_loc(nnn)) jloc=1
            if(nbd3_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbd3_new(iii).le.mov_end_loc(nnn)) kloc=1
            if(nbd4_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbd4_new(iii).le.mov_end_loc(nnn)) lloc=1
            if(iloc.eq.1.or.jloc.eq.1.or.kloc.eq.1.or.lloc.eq.1) then
              num_dih_thrd=num_dih_thrd+1
            endif
          enddo
        enddo

C allocate memory for the dihedrals on this thread
C idih_thrd : gives the id of this dihedral from 1:num_sys_dihes
C ddih_thrd : is the factor used to scale energy by so that we don't add
C up the energy incorrectly when doing same dihedral on multiple threads

        allocate(idih_thrd(1:num_dih_thrd))
        allocate(ddih_thrd(1:num_dih_thrd))
        write(*,*)'# of dihes on this local thread = ',
     &             mythrd,num_dih_thrd

C now put in all the information for this thread

        num_dih_thrd=0
        do nnn=1,mov_num_loc
          do iii=1,num_sys_dihes
            iloc=0
            jloc=0
            kloc=0
            lloc=0
            if(nbd1_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbd1_new(iii).le.mov_end_loc(nnn)) iloc=1
            if(nbd2_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbd2_new(iii).le.mov_end_loc(nnn)) jloc=1
            if(nbd3_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbd3_new(iii).le.mov_end_loc(nnn)) kloc=1
            if(nbd4_new(iii).ge.mov_beg_loc(nnn).and.
     &         nbd4_new(iii).le.mov_end_loc(nnn)) lloc=1
            if(iloc.eq.1.and.jloc.eq.1.and.
     &         kloc.eq.1.and.lloc.eq.1) then
              num_dih_thrd=num_dih_thrd+1
              idih_thrd(num_dih_thrd)=iii
              ddih_thrd(num_dih_thrd)=1.0
            elseif((iloc.eq.1.and.jloc.eq.0.and.
     &              kloc.eq.0.and.lloc.eq.0).or.
     &             (iloc.eq.0.and.jloc.eq.1.and.
     &              kloc.eq.0.and.lloc.eq.0).or.
     &             (iloc.eq.0.and.jloc.eq.0.and.
     &              kloc.eq.1.and.lloc.eq.0).or.
     &             (iloc.eq.0.and.jloc.eq.0.and.
     &              kloc.eq.0.and.lloc.eq.1)) then
              num_dih_thrd=num_dih_thrd+1
              idih_thrd(num_dih_thrd)=iii
              ddih_thrd(num_dih_thrd)=1.0/4.0
            elseif((iloc.eq.1.and.jloc.eq.1.and.
     &              kloc.eq.0.and.lloc.eq.0).or.
     &             (iloc.eq.1.and.jloc.eq.0.and.
     &              kloc.eq.1.and.lloc.eq.0).or.
     &             (iloc.eq.1.and.jloc.eq.0.and.
     &              kloc.eq.0.and.lloc.eq.1).or.
     &             (iloc.eq.0.and.jloc.eq.1.and.
     &              kloc.eq.1.and.lloc.eq.0).or.
     &             (iloc.eq.0.and.jloc.eq.1.and.
     &              kloc.eq.0.and.lloc.eq.1).or.
     &             (iloc.eq.0.and.jloc.eq.0.and.
     &              kloc.eq.1.and.lloc.eq.1)) then
              num_dih_thrd=num_dih_thrd+1
              idih_thrd(num_dih_thrd)=iii
              ddih_thrd(num_dih_thrd)=2.0/4.0
            elseif((iloc.eq.1.and.jloc.eq.1.and.
     &              kloc.eq.1.and.lloc.eq.0).or.
     &             (iloc.eq.1.and.jloc.eq.1.and.
     &              kloc.eq.0.and.lloc.eq.1).or.
     &             (iloc.eq.1.and.jloc.eq.0.and.
     &              kloc.eq.1.and.lloc.eq.1).or.
     &             (iloc.eq.0.and.jloc.eq.1.and.
     &              kloc.eq.1.and.lloc.eq.1)) then
              num_dih_thrd=num_dih_thrd+1
              idih_thrd(num_dih_thrd)=iii
              ddih_thrd(num_dih_thrd)=3.0/4.0
            endif
          enddo
        enddo

C try a parallel initialization of everything that might be an
C issue in the move loop - this is probably overkill, but I can't think
C of a good reason not to do it right now

!$omp barrier ! temporary
        do nnn=1,mov_num_loc
          do iii=mov_beg_loc(nnn),mov_end_loc(nnn)
            f_a_x(1:3,iii,myrep)=0.0
            c_a_x(1:3,iii,myrep)=t_a_x(1:3,iii,myrep)
            if(i_do_lincs) c_a_x_old(1:3,iii,myrep)=t_a_x(1:3,iii,myrep)
          enddo
        enddo

C 2019 - the above will not give us coordinates for all the static atoms
C in the system so let's do that now - note that the "move_atoms"
C routine will not mess with these once they've been assigned...

!$omp single
        do iii=num_mov_atms+1,num_tot_atms
          c_a_x(1:3,iii,myrep)=t_a_x(1:3,iii,myrep)
        enddo
!$omp end single
        
C set long-range elec energy to zero so that if we *don't* use
C treecode_elec there is still a reasonable value posted

        ene_e_lg_loc=0.0

!$omp barrier ! temporary
        do nnn=1,mov_num_loc
          if(mov_hyd_loc(nnn).ge.1) then
            kkk=mov_hyd_loc(nnn)
            do ijk=mov_beg_loc(nnn),mov_end_loc(nnn)
              iii=atm_to_hyd_atm(ijk)
              ii3=iii*3-2
              d_a_r(1:num_hyd_atms3(kkk),ii3  ,kkk,myrep)=0.0
              d_a_r(1:num_hyd_atms3(kkk),ii3+1,kkk,myrep)=0.0
              d_a_r(1:num_hyd_atms3(kkk),ii3+2,kkk,myrep)=0.0
              s_a_r(1:num_hyd_atms3(kkk),ii3  ,kkk,myrep)=0.0
              s_a_r(1:num_hyd_atms3(kkk),ii3+1,kkk,myrep)=0.0
              s_a_r(1:num_hyd_atms3(kkk),ii3+2,kkk,myrep)=0.0
            enddo
          endif
        enddo 

C put in entries for self-self interactions (j=i,i) once

!$omp barrier ! necessary & not worth worrying about

C if using fixman and override need to calculate Chebyshev coefficients 

       if(fixman.and.fixman_override) then

!$omp master

         lmax=lmax_read_in
         lmin=lmin_read_in
         da=2.0/(lmax-lmin)
         db=(lmax+lmin)/(lmin-lmax)
         write(*,881)lmin,lmax,da,db
         do i1=1,fixman_order
           fixman_coeff(i1-1)=0.0
         enddo
         bma=0.5*(lmax-lmin)
         bpa=0.5*(lmax+lmin)
         do i1=1,fixman_order
           yme=cos(3.14159265d0*(i1-0.5d0)/fixman_order)
           zme(i1)=sqrt(yme*bma+bpa)
         enddo
         do i1=1,fixman_order
           xme=0.000
           do i2=1,fixman_order
             xme=xme+zme(i2)*
     &            cos((3.14159265d0*(i1-1))*((i2-0.5d0)/
     &            (fixman_order)))
            enddo
            fixman_coeff(i1-1)=xme*2.000/(fixman_order)
            write(*,*)'fixman coeffs ',i1-1,fixman_coeff(i1-1)
          enddo

!$omp end master
!$omp barrier

       endif

!$omp barrier ! 99% necessary
            
C xxxxx Section
C xxxxx Main loop of the program starts here

        write(*,223)num_tot_atms,num_flx_atms,num_hyd_atms
223     format('FINAL SANITY CHECK ',3i10)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCC main loop begins here CCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

222     continue                                       

        i_need_to_die_now=0

!$omp barrier ! TEMPORARY FOR TIMINGS

        call system_clock(count=ftm1,count_rate=ru)

C zero energies only when we're going to record energy

        if(i_debug) write(*,*)'new step  ',mythrd
        if (necount.eq.neprint) then                     

          ene_pos_loc=0.0
          ene_bond_loc=0.0
          ene_angl_loc=0.0
          ene_dihe_loc=0.0
          ene_nbnd_loc=0.0
          ene_gofav_loc=0.0
          ene_gounfav_loc=0.0

!$omp single
          ene_m=0.0
          ene_pos=0.0
          ene_bond=0.0
          ene_angl=0.0
          ene_dihe=0.0
          etvst=0.0
          etgst=0.0
          etest=0.0
          etvmd=0.0
          etgmd=0.0
          etemd=0.0
          etelg=0.0
          etewldreal=0.0
          etewldrcip=0.0
          ene_nbnd_tot=0.0
          ene_gofav_tot=0.0
          ene_gounfav_tot=0.0
!$omp end single 

        endif
        if(i_debug) write(*,*)'new2 step  ',mythrd

C now zero forces on each thread at each timestep

        if(walls.or.sphere.or.cylinder) osm_cur_loc=0.0

        if(walls.or.sphere.or.cylinder) then
          if(mythrd.eq.1) osm_cur=0.0
        endif

C for bonded forces - which are private - we only need to zero those
C entries that correspond to thread-local atoms

        do nnn=1,mov_num_loc
          f_a_x_bd(1:3,mov_beg_loc(nnn):mov_end_loc(nnn))=0.0   
        enddo
        if(i_debug) write(*,*)'new3 step  ',mythrd

C for nonbonded forces we can go one of two ways
C if using newton's law then we zero our entire private force arrays
C if not using newton's law then we zero only those entries of the
C global force array that correspond to thread-local atoms - note that
C since e_a_x_st,e_e_x_st are global they need myrep (replica) dimension

        if(i_do_YHL_nonbondeds) then
          f_a_x_st=0.0 
          if(treecode_elec) f_e_x_st=0.0
        else
          do nnn=1,mov_num_loc
            e_a_x_st(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)=0.0
          enddo
          if(treecode_elec) then
            do nnn=1,mov_num_loc
              e_e_x_st(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)=0.0
            enddo
          endif
        endif
         if(i_debug) write(*,*)'new4 step  ',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(1,mythrd)=tmtrx(1,mythrd)+real(ftm2-ftm1)/real(ru)
 
C xxxxx Section
C xxxxx Update the various nonbonded lists

C allocate atoms to grid cells here

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

!$omp single
        if(shrinking_box) call shrink_box
!$omp end single

        if(i_debug) write(*,*)'new5 step  ',mythrd
        if(num_lst_cnt.eq.num_lst_stp.or.update_nonbond_list) then

!$omp barrier    ! DEFO NEEDED - this serves as the barrier before assigning
                 ! atoms to the grid or for diffusion tensor calculations
                 ! and ensures that all forces are zero'd correctly

C note that only the master thread of each replica assigns atoms to grid

          if(i_debug) write(*,*)'assign     ',mythrd

          if(myrep_master)     
     &      call assign_atoms_to_grid(num_tot_atms,mythrd,myrep,
     &                                i1_loc,i2_loc,i3_loc)

C 2019 also put static atoms on a grid - note that in the call above and
C the call below we pass "num_tot_atms" as we temporarily use arrays
C called i1_lic,i2_lic,i3_lic to store the grid locations of the atoms -
C I don't know why we pass these back as i1_loc etc as they're never
C used again - the code below will probably mess up the passed-back
C values but since they're not used again I don't think we need worry

          if(i_debug) write(*,*)'assign0    ',mythrd
          if(myrep_master.and.nstep.eq.0)     
     &      call assign_static_atoms_to_grid(num_tot_atms,mythrd,myrep,
     &                                       i1_loc,i2_loc,i3_loc)

          if(i_debug) write(*,*)'done assign     ',mythrd

        endif

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(2,mythrd)=tmtrx(2,mythrd)+real(ftm2-ftm1)/real(ru)
 
C now if we're using exclusive go-pairs then we need to ask if it's time 
C to update the list of assigned exclusive interactions
C next do a zip through the nonbonded pair counting - should basically be a copy
C of what we do immediately afterwards
C note that this section likely needs a lot of work especially now that
C we have replica added as an option - the c_a_x entries should all work
C but the go_mod_arr stuff may/may not need an extra replica dimension

!$omp barrier ! TEMPORARY FOR TIMINGS
 
C now, having checked the exclusive go stuff we can return to putting together
C a nonbonded list (we jumped out of the 'num_lst_cnt' if loop previously)
C note that we now use the mixed_schedule mov_num_loc entries instead of
C ibeg_f_a,iend_f_a in order to build the lists - allows for better
C load-balancing with systems like chromosome surrounded by rigid bodies

        if(num_lst_cnt.eq.num_lst_stp.or.update_nonbond_list) then

!$omp barrier ! TEMPORARY FOR TIMINGS
          call system_clock(count=ftm1,count_rate=ru)

          imin_thrd_int=num_tot_atms+1 ! new for YHL code
          imax_thrd_int=0              ! new for YHL code
          if(i_debug) write(*,*)'start nonbond   ',mythrd

          call make_nonbond_list(mythrd,num_threads,myrep,
     &           mov_num_loc,mov_beg_loc,mov_end_loc,
     &           n_ff_r_st,ida_ff_r_st,jda_ff_r_st,kda_ff_r_st,
     &           n_ff_u_st,ida_ff_u_st,jda_ff_u_st,kda_ff_u_st,
     &           n_ff_a_st,ida_ff_a_st,jda_ff_a_st,kda_ff_a_st,
     &           n_ff_a_md,ida_ff_a_md,jda_ff_a_md,kda_ff_a_md,
     &           n_ff_b_st,ida_ff_b_st,jda_ff_b_st,kda_ff_b_st,
     &                     lda_ff_b_st,
     &           n_ff_b_md,ida_ff_b_md,jda_ff_b_md,kda_ff_b_md,
     &                     lda_ff_b_md,
     &           n_ff_v_st,ida_ff_v_st,jda_ff_v_st,s11_ff_v_st,
     &                     s22_ff_v_st,s33_ff_v_st,s44_ff_v_st,
     &                     d11_ff_v_st,d22_ff_v_st,e00_ff_v_st,
     &           n_ff_v_md,ida_ff_v_md,jda_ff_v_md,s11_ff_v_md,
     &                     s22_ff_v_md,s33_ff_v_md,s44_ff_v_md,
     &                     d11_ff_v_md,d22_ff_v_md,e00_ff_v_md,
     &           n_ff_g_st,ida_ff_g_st,jda_ff_g_st,kda_ff_g_st,
     &         lda_ff_g_st,mda_ff_g_st,nda_ff_g_st,s11_ff_g_st,
     &                     s22_ff_g_st,s33_ff_g_st,s44_ff_g_st,
     &                     d11_ff_g_st,d22_ff_g_st,d12_ff_g_st,
     &                     ddd_ff_g_st,e00_ff_g_st,rep_ff_g_st,
     &           n_ff_g_md,ida_ff_g_md,jda_ff_g_md,kda_ff_g_md,
     &         lda_ff_g_md,mda_ff_g_md,nda_ff_g_md,s11_ff_g_md,
     &                     s22_ff_g_md,s33_ff_g_md,s44_ff_g_md,
     &                     d11_ff_g_md,d22_ff_g_md,d12_ff_g_md,
     &                     ddd_ff_g_md,e00_ff_g_md,rep_ff_g_md,
     &           n_ff_g_rx,ida_ff_g_rx,jda_ff_g_rx,kda_ff_g_rx,
     &                     lda_ff_g_rx,mda_ff_g_rx,s11_ff_g_rx,
     &                     s22_ff_g_rx,s33_ff_g_rx,s44_ff_g_rx,
     &                     d11_ff_g_rx,d22_ff_g_rx,d12_ff_g_rx,
     &                     ddd_ff_g_rx,e00_ff_g_rx,
     &           n_ff_e_st,ida_ff_e_st,jda_ff_e_st,qqq_ff_e_st,
     &           n_ff_e_md,ida_ff_e_md,jda_ff_e_md,qqq_ff_e_md,
     &     n_ff_r_st_old,n_ff_u_st_old,n_ff_a_st_old,n_ff_a_md_old,
     &     n_ff_b_st_old,n_ff_b_md_old,n_ff_v_st_old,n_ff_v_md_old,
     &     n_ff_g_st_old,n_ff_g_md_old,n_ff_e_st_old,n_ff_e_md_old,
     &     n_ff_g_rx_old,rep_go_scale(mythrds_rep_curr(mythrd)))
          n_ff_rods=0
          if(i_debug) write(*,*)'done  nonbond   ',mythrd
        endif                                           

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(4,mythrd)=tmtrx(4,mythrd)+real(ftm2-ftm1)/real(ru)
 
C zero out all the go modes and contact numbers - note that we only
C count up go contacts when we update the nonbonded list and we only add
C them all up when we want to write out the energy - this was horribly
C slow code previously because of the critical statements but these seem
C completely unnecessary given that the arrays are dimensioned by thread

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

        if (necount.eq.neprint) then                     
          if(i_use_go_pairs) then

          if(i_debug) write(*,*)'start ene zero  ',mythrd

            num_cur_Q_loc(:,:,mythrd)=0       ! legal?
            num_cur_Q_mod_loc(:,:,:,mythrd)=0 ! legal?

C note that here we use nda_ff_g_st and nda_ff_g_md as a flag so that we
C only include this go contact if its epsilon is about the thresthold
C go_eps_low read in at run time

            if(i_do_YHL_nonbondeds) then

              do mmm=1,n_ff_g_st ! add up st contacts
                if(ddd_ff_g_st(mmm).le.d12_ff_g_st(mmm).and.
     &             nda_ff_g_st(mmm).eq.1) then
                  kkk=kda_ff_g_st(mmm) ! go mol 1
                  lll=lda_ff_g_st(mmm) ! go mol 2
                  mno=mda_ff_g_st(mmm) ! go mode
                  num_cur_Q_loc(kkk,lll,mythrd)=
     &            num_cur_Q_loc(kkk,lll,mythrd)+1
                  num_cur_Q_mod_loc(mno,kkk,lll,mythrd)=
     &            num_cur_Q_mod_loc(mno,kkk,lll,mythrd)+1
                  num_cur_Q_loc(lll,kkk,mythrd)=
     &            num_cur_Q_loc(lll,kkk,mythrd)+1
                  num_cur_Q_mod_loc(mno,lll,kkk,mythrd)=
     &            num_cur_Q_mod_loc(mno,lll,kkk,mythrd)+1
                endif
              enddo

              do mmm=1,n_ff_g_md ! add up md contacts
                if(ddd_ff_g_md(mmm).le.d12_ff_g_md(mmm).and.
     &             nda_ff_g_md(mmm).eq.1) then
                  kkk=kda_ff_g_md(mmm) ! go mol 1
                  lll=lda_ff_g_md(mmm) ! go mol 2
                  mno=mda_ff_g_md(mmm) ! go mode
                  num_cur_Q_loc(kkk,lll,mythrd)=
     &            num_cur_Q_loc(kkk,lll,mythrd)+1
                  num_cur_Q_mod_loc(mno,kkk,lll,mythrd)=
     &            num_cur_Q_mod_loc(mno,kkk,lll,mythrd)+1
                  num_cur_Q_loc(lll,kkk,mythrd)=
     &            num_cur_Q_loc(lll,kkk,mythrd)+1
                  num_cur_Q_mod_loc(mno,lll,kkk,mythrd)=
     &            num_cur_Q_mod_loc(mno,lll,kkk,mythrd)+1
                endif
              enddo

            else ! do this if .not.i_do_YHL_nonbondeds

              do mmm=1,n_ff_g_st ! add up st contacts
                if(ddd_ff_g_st(mmm).le.d12_ff_g_st(mmm).and.
     &             nda_ff_g_st(mmm).eq.1) then
                  kkk=kda_ff_g_st(mmm) ! go mol 1
                  lll=lda_ff_g_st(mmm) ! go mol 2
                  mno=mda_ff_g_st(mmm) ! go mode
                  num_cur_Q_loc(kkk,lll,mythrd)=
     &            num_cur_Q_loc(kkk,lll,mythrd)+1
                  num_cur_Q_mod_loc(mno,kkk,lll,mythrd)=
     &            num_cur_Q_mod_loc(mno,kkk,lll,mythrd)+1
                endif
              enddo

              do mmm=1,n_ff_g_md ! add up md contacts
                if(ddd_ff_g_md(mmm).le.d12_ff_g_md(mmm).and.
     &             nda_ff_g_md(mmm).eq.1) then
                  kkk=kda_ff_g_md(mmm) ! go mol 1
                  lll=lda_ff_g_md(mmm) ! go mol 2
                  mno=mda_ff_g_md(mmm) ! go mode
                  num_cur_Q_loc(kkk,lll,mythrd)=
     &            num_cur_Q_loc(kkk,lll,mythrd)+1
                  num_cur_Q_mod_loc(mno,kkk,lll,mythrd)=
     &            num_cur_Q_mod_loc(mno,kkk,lll,mythrd)+1
                endif
              enddo
              if(i_debug) write(*,*)'don e ene zero  ',mythrd

            endif
          endif

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(5,mythrd)=tmtrx(5,mythrd)+real(ftm2-ftm1)/real(ru)
 
        endif                                           

C xxxxx Section
C xxxxx Hydrodynamics-related stuff here

        if(num_hyd_cnt.eq.num_hyd_stp.or.update_nonbond_list) then

!$omp barrier ! TEMPORARY FOR TIMINGS
          call system_clock(count=ftm1,count_rate=ru)

          if(i_debug) write(*,*)'start diff calc ',mythrd
          call calculate_diffusion_tensor(mythrd,myrep,
     &           mov_num_loc,mov_typ_loc,mov_beg_loc,
     &           mov_end_loc,mov_dyn_loc,mov_hyd_loc)

          if(i_debug) write(*,*)'done  diff calc ',mythrd
          call system_clock(count=ftm2,count_rate=ru)
          tmtrx(6,mythrd)=tmtrx(6,mythrd)+real(ftm2-ftm1)/real(ru)
 
C obtain a bunch of random vectors - easier to do this in main program
C and we do it here because it avoids us explicitly adding a barrier
C between this routine and the next, which makes the correlated disps

!$omp barrier ! TEMPORARY FOR TIMINGS

          call system_clock(count=ftm1,count_rate=ru)

          if(i_debug) write(*,*)'call randoms lc ',mythrd
          errcode=vsrnggaussian(method,stream,
     &             num_randoms,random_trms,a_for_mkl,sigma)
          if(i_debug) write(*,*)'done randoms lc ',mythrd

          call system_clock(count=ftm2,count_rate=ru)
          tmtrx(7,mythrd)=tmtrx(7,mythrd)+real(ftm2-ftm1)/real(ru)
 
C convert these random vectors into uncorrelated random displacements

!$omp barrier ! TEMPORARY FOR TIMINGS
          call system_clock(count=ftm1,count_rate=ru)

          if(i_debug) write(*,*)'call Uncorr   c ',mythrd
          call calculate_uncorrelated_randoms(mythrd,myrep,
     &           mov_num_loc,mov_typ_loc,mov_beg_loc,
     &           mov_end_loc,mov_dyn_loc,mov_hyd_loc,
     &           num_randoms,random_trms)
          if(i_debug) write(*,*)'done Uncorr   c ',mythrd

          call system_clock(count=ftm2,count_rate=ru)
          tmtrx(8,mythrd)=tmtrx(8,mythrd)+real(ftm2-ftm1)/real(ru)
 
C the HSL code is so complex it's probably not worth trying to get it to
C work in a subroutine, so keep the code here in the main routine 
C note that use of col_from_hybrid required a minor alteration to
C the HSL_MP54_CODE.f90 to make this routine 'public'

C some ideas to test:
C (1) loop loading amp could be written vectorizable NO I THINK NOT
C (2) loop loading s_a_r might be reorderable using 'N' in sgemm DONE!
C (3) barriers might be removable if parts dispersed through code
C (4) 'master' might be replaceable by 'single' SOUNDS DANGEROUS
C (5) col_from_hybrid_block might be threadable PROBABLY NOT

!$omp barrier ! TEMPORARY FOR TIMINGS
          call system_clock(count=ftm1,count_rate=ru)

          if(cholesky) then

!$omp single 

C first, copy the diffusion tensor (d_a_r) into amkl

            if(i_debug) write(*,*)'copying d_a_r to amkl'
            s_a_r=d_a_r

C now get the cholesky decomposition

            call omp_set_num_threads(num_threads)
CCC         call omp_set_num_threads(12)

            if(i_debug) write(*,*)'calling spotrf'
            call spotrf('L',num_hyd_atms3(1),s_a_r,
     &                      num_hyd_atms3(1),info)
            if(i_debug) write(*,*)'back from spotrf'

C now zero out the upper triangle in s_a_r

            do i1=1,num_hyd_atms3(1)
              do j1=i1+1,num_hyd_atms3(1)
                s_a_r(i1,j1,1,1)=0.0
              enddo
            enddo
            if(i_debug) write(*,*)'done copying amkl to s_a_r'

C quit if the cholesky decomposition failed...

            if(info.ne.0) then
              write(*,*)
              write(*,*)'failure in spotrf ',info
              write(*,*)'quitting :('
              write(*,*)
              stop
            endif

C let's see if this helps...

C           call omp_set_num_threads(1)

!$omp end single

C let's see if this helps...

            call omp_set_num_threads(1)

C fixman and not fixman_override here...

          elseif(full_hydrodynamics.and.fixman.and.
     &           .not.fixman_override) then

!$omp single 

C first, copy the diffusion tensor ('d_a_r') into UPLO form required
C again, note here that we are only implementing for 1 hydrodynamic
C system and for 1 replica currently (hence ...1,1) below)

            k1=0
            do i1=1,num_hyd_atms3(1)
              do j1=1,i1
                k1=k1+1
                amkl(k1)=d_a_r(j1,i1,1,1)
              enddo
            enddo

C now convert matrix into tridiagonal form

            call ssptrd('U',num_hyd_atms3(1),amkl,dmkl,emkl,tmkl,info)

C now get the eigenvalues of this matrix - dmkl stores eigenvalues

            call ssteqr('N',num_hyd_atms3(1),dmkl,emkl,zmkl,
     &                  num_hyd_atms3(1),work,info)

C now find maxima and minima of eigenvalues

            lmin=1000000000.0
            lmax=-1000000000.0
            do i=1,num_hyd_atms3(1) ! 1 hydro sys
              if(dmkl(i).gt.lmax) lmax=dmkl(i)
              if(dmkl(i).lt.lmin) lmin=dmkl(i)
            enddo
            da=2.0/(lmax-lmin)
            db=(lmax+lmin)/(lmin-lmax)
            write(*,881)lmin,lmax,da,db
881         format('lmin, lmax, da, db = ',4f12.8)

            stop

C calculate Chebyshev coefficients - there are fixman_order coeffs but
C they're stored in entries (0:fixman_order-1) of the array fixman_coeff

            do i1=1,fixman_order
              fixman_coeff(i1-1)=0.0
            enddo
            bma=0.5*(lmax-lmin)
            bpa=0.5*(lmax+lmin)
            do i1=1,fixman_order
              yme=cos(3.14159265d0*(i1-0.5d0)/fixman_order)
              zme(i1)=sqrt(yme*bma+bpa)
            enddo
            do i1=1,fixman_order
              xme=0.000
              do i2=1,fixman_order
                xme=xme+zme(i2)*
     &              cos((3.14159265d0*(i1-1))*((i2-0.5d0)/
     &              (fixman_order)))
              enddo
              fixman_coeff(i1-1)=xme*2.000/(fixman_order)
            enddo
!$omp end single
!$omp barrier

          endif

          call system_clock(count=ftm2,count_rate=ru)
          tmtrx(9,mythrd)=tmtrx(9,mythrd)+real(ftm2-ftm1)/real(ru)
 
C convert the uncorrelated randoms into correlated randoms
C note that a barrier is needed prior to this...

!$omp barrier ! TEMPORARY FOR TIMINGS
          call system_clock(count=ftm1,count_rate=ru)

!$omp barrier ! DEFO NEEDED - we need randoms/D/S fully defined
              ! INEXPENSIVE because only performed infrequently

C note that this is only currently implemented for flex molecules
C also note that it's much cleaner to call this subroutine multiple
C times than it is to have the nnn=1,mov_num_loc within the routine

          if(cholesky) then

            if(i_debug) write(*,*)'calling cholesky',mythrd
            do nnn=1,mov_num_loc
              if(mov_typ_loc(nnn).eq.1) then
                if(mov_hyd_loc(nnn).ge.1) then
                  jjj=atm_to_hyd_atm(mov_beg_loc(nnn))
                  jj3=jjj*3-2
                  kk3=(mov_end_loc(nnn)-mov_beg_loc(nnn)+1)*3
                  kkk=mov_hyd_loc(nnn) ! pass the identity of HI system 
                  call calculate_cholesky_randoms(jj3,kk3,kkk,
     &                                            mythrd,myrep)
                endif
              elseif(mov_typ_loc(nnn).eq.0) then
                write(*,*)'calculate_cholesky_randoms not yet'
                write(*,*)'implemented for rigd molecules'
                write(*,*)'quitting :('
                stop
              endif
            enddo

          elseif(fixman) then

C copy and scale the diffusion tensor

            do i1=mov_beg_loc(1)*3-2,mov_end_loc(1)*3
              do i2=1,num_hyd_atms3(1)
                fixman_tensor(i2,i1)=d_a_r(i2,i1,1,1)*da
              enddo
            enddo
            do i1=mov_beg_loc(1)*3-2,mov_end_loc(1)*3
              fixman_tensor(i1,i1)=fixman_tensor(i1,i1)+db
            enddo

C store the zero'th fixman tensor term - again note that the first '1'
C is the hydrodynamic system and the second '1' is the replica number

            do i1=1,num_hyd_stp
              do i2=mov_beg_loc(1)*3-2,mov_end_loc(1)*3
                fixman_xx_mtrx(i2,i1,0)=h_u_r(i2,i1,1,1)
              enddo
            enddo

!$omp barrier ! needed?

            call sgemm('T', ! needed?
     &                 'N',
     &                 mov_tot_loc3,  ! assumes 1 system etc...
     &                 num_hyd_stp,
     &                 num_hyd_atms3(1), ! ditto...
     &                 1.0,
     &                 fixman_tensor(1,mov_beg_loc(1)*3-2),
     &                 num_hyd_atms3(1),
     &                 h_u_r(1,1,1,1),
     &                 num_hyd_atms3(1),
     &                 0.0,
     &                 fixman_xx_mtrx(mov_beg_loc(1)*3-2,1,1),
     &                 num_hyd_atms3(1))

!$omp barrier ! needed?

C calculate the Chebyshev total for the zero'th and first terms
C these are the approximation to the correlated random displacements
C that we have after adding up only 2 terms in the expansion
C (typically we'll need 10-30 terms depending on the system and on the
C error tolerance that we put up with)

            do i1=1,num_hyd_stp
              do i2=mov_beg_loc(1)*3-2,mov_end_loc(1)*3
                fixman_ran_mtrx(i2,i1)=
     &      0.5*fixman_coeff(0)*fixman_xx_mtrx(i2,i1,0)+
     &          fixman_coeff(1)*fixman_xx_mtrx(i2,i1,1)
              enddo
            enddo
 
C continue adding on terms till we're converged or we
C exceed the maximum number of iteratons (fixman_order)

C note that we first calculate the fixman_tensor term and then
C multiply it by its Chebyshev coefficient & check the error too...

            do i1=2,fixman_order-1

              call sgemm('T',
     &                   'N',
     &                   mov_tot_loc3,  ! assumes 1 system etc...
     &                   num_hyd_stp,
     &                   num_hyd_atms3(1), ! ditto...
     &                   2.0, ! note alpha=2 here
     &                   fixman_tensor(1,mov_beg_loc(1)*3-2),
     &                   num_hyd_atms3(1),
     &                   fixman_xx_mtrx(1,1,i1-1),
     &                   num_hyd_atms3(1),
     &                   0.0,
     &                   fixman_xx_mtrx(mov_beg_loc(1)*3-2,1,i1),
     &                   num_hyd_atms3(1))

!$omp barrier ! needed?

C having added on a new term we can update the recursive terms
C in the Chebyshev expansion and start to compute the error

              do i2=1,num_hyd_stp
                fixman_l2_num_mtrx(i2,mythrd)=0.0
                fixman_l2_den_mtrx(i2,mythrd)=0.0
                do i3=mov_beg_loc(1)*3-2,mov_end_loc(1)*3
                  fixman_xx_mtrx(i3,i2,i1)=
     &            fixman_xx_mtrx(i3,i2,i1)-
     &            fixman_xx_mtrx(i3,i2,i1-2)

C calculate components for the l2-norm - need the numerator and
C denominator updating simultaneously

                  fixman_l2_den_mtrx(i2,mythrd)=
     &            fixman_l2_den_mtrx(i2,mythrd)+
     &            fixman_ran_mtrx(i3,i2)**2
                  fixman_l2_num_mtrx(i2,mythrd)=
     &            fixman_l2_num_mtrx(i2,mythrd)+
     &           (fixman_coeff(i1)*
     &            fixman_xx_mtrx(i3,i2,i1))**2

C here's where we actually update our approximation to the correlated
C random displacements...

                  fixman_ran_mtrx(i3,i2)=
     &            fixman_ran_mtrx(i3,i2)+
     &            fixman_coeff(i1)*
     &            fixman_xx_mtrx(i3,i2,i1)
                enddo
              enddo

!$omp barrier ! needed?

C now calculate the l2norm of the difference between i'th iteration and
C the i-1'th iteration of the correlated random displacements and
C normalize by that of the i-1'th iteration
C note that this involves adding up components obtained from each thread
C and computing the error for all 'num_hyd_stp' sets of 3N correlated
C displacements. only if *all* sets are below the specified error do we
C jump out of the Chebyshev expansion

!$omp single
              fixman_l2_mtrx=-999.9
              do i2=1,num_hyd_stp
                fixman_l2_num=0.0
                fixman_l2_den=0.0
                do i3=1,num_threads
                  fixman_l2_num=fixman_l2_num+
     &                          fixman_l2_num_mtrx(i2,i3)
                  fixman_l2_den=fixman_l2_den+
     &                          fixman_l2_den_mtrx(i2,i3)
                enddo
                temp_l2=sqrt(fixman_l2_num/fixman_l2_den)
                if(temp_l2.gt.fixman_l2_mtrx) then
                  fixman_l2_mtrx=temp_l2
                endif
              enddo
!$omp end single

C if the worst error that we found is below the required tolerance then
C we are done and we can skip ahead

              if(fixman_l2_mtrx.lt.fixman_tol) goto 7887

            enddo

C come here if the method converged

7887        continue

C finish by copying the random terms into the h_c_r array
C again, this currently assumes only one hydrodynamic system

            do kkk=1,num_hyd_stp
              do jjj=mov_beg_loc(1)*3-2,mov_end_loc(1)*3
                h_c_r(jjj,kkk,1,1)=fixman_ran_mtrx(jjj,kkk)
              enddo
            enddo

!$omp barrier ! defo needed if we thread the h_c_r update above

          endif

          if(i_debug) write(*,*)'done calling cholesky',mythrd

          call system_clock(count=ftm2,count_rate=ru)
          tmtrx(10,mythrd)=tmtrx(10,mythrd)+real(ftm2-ftm1)/real(ru)

        endif

C xxxxx Section
C xxxxx Write out movie and restart files

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

        if(ntcount.eq.ntprint) then
          if(mythrd.eq.1) call write_xtc(mythrd,myrep)
        endif ! in 01-04-13 this was prior to barrier...

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(11,mythrd)=tmtrx(11,mythrd)+real(ftm2-ftm1)/real(ru)

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

        if(nmcount.eq.nmprint) then                        
          if(mythrd.eq.1) then
            write(*,544)myrep,tot_time
544         format('replica # ',i5,' asks write restart.file @ ',
     &             'time = ',f15.5)   
          endif
          if(myrep_master) then
            write(*,*)'replica # ',myrep,' writes movie.pdb'
            call write_movie_file(movieframenumber,icrashed,myrep)       
            write(*,*)'replica # ',myrep,' writes restart.file at ',
     &                'time = ',tot_time
            call write_restart(movieframenumber,myrep)
            write(*,*)'replica # ',myrep,' done with writing'  
          endif
        endif
                                                                                                      
!$omp barrier ! need this before resetting nmcount

        if(nmcount.eq.nmprint.and.mythrd.eq.1) then
          movieframenumber=movieframenumber+1                      
          nmcount=0
        endif

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(12,mythrd)=tmtrx(12,mythrd)+real(ftm2-ftm1)/real(ru)
 
C now call all of the force routines                                                  
                           
C!$omp barrier ! NOT NEEDED COS WE USE A SINGLE TO ZERO ENERGY/FORCE

C xxxxx Section
C xxxxx Calculate bond forces 

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

         if(i_debug) write(*,*)'calling bonded forces',mythrd
        call calculate_bonded_forces
     &        (num_tot_atms,num_f_ms,ibeg_f_a,iend_f_a,
     &         mov_num_loc,mov_beg_loc,mov_end_loc,
     &         num_bnd_thrd,ibnd_thrd,dbnd_thrd,ene_bond_loc,
     &         num_ang_thrd,iang_thrd,dang_thrd,ene_angl_loc,
     &         num_dih_thrd,idih_thrd,ddih_thrd,ene_dihe_loc,
     &         ene_pos_loc,ene_wll_loc,ene_umb_loc,
     &         i_need_to_die_now,i_am_master,
     &         n_ff_rods,ida_ff_rods,jda_ff_rods,
     &         s11_ff_rods,s22_ff_rods,s33_ff_rods,s44_ff_rods,
     &         d11_ff_rods,d22_ff_rods,e00_ff_rods,
     &         umb_fx1,umb_fy1,umb_fz1,
     &         umb_fx2,umb_fy2,umb_fz2,
     &         f_a_x_bd,num_osm,osm_cur_loc,myrep,mythrd)
         if(i_debug) write(*,*)'done calling bonded forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(13,mythrd)=tmtrx(13,mythrd)+real(ftm2-ftm1)/real(ru)
 
C!$omp barrier ! NOT NEEDED SURELY - SO LONG AS NONBOND LIST LOCAL

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(14,mythrd)=tmtrx(14,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

        ene_v_st_loc=0.0
        ene_g_st_loc=0.0
        ene_e_st_loc=0.0

C note that since we update a thread-dependent part of e_a_x_st 
C in the no_newton routines we *may* need to have an explicit 
C interface for each routine - not yet sure but seems likely

         if(i_debug) write(*,*)'calling vdw forces',mythrd
        if(i_do_YHL_nonbondeds) then
          call calculate_vdw_forces_newton(n_ff_v_st,
     &           ida_ff_v_st,jda_ff_v_st,s11_ff_v_st,
     &           s22_ff_v_st,s33_ff_v_st,s44_ff_v_st,
     &           d11_ff_v_st,d22_ff_v_st,e00_ff_v_st,
     &           ene_v_st_loc,f_a_x_st)
        else
          call calculate_vdw_forces_no_newton(n_ff_v_st,
     &           ida_ff_v_st,jda_ff_v_st,s11_ff_v_st,
     &           s22_ff_v_st,s33_ff_v_st,s44_ff_v_st,
     &           d11_ff_v_st,d22_ff_v_st,e00_ff_v_st,
     &           ene_v_st_loc,e_a_x_st(:,:,myrep),myrep)
        endif
         if(i_debug) write(*,*)'done calling vdw forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(15,mythrd)=tmtrx(15,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

         if(i_debug) write(*,*)'calling go forces',mythrd
        if(i_do_YHL_nonbondeds) then
          call calculate_go_forces_newton(n_ff_g_st,
     &           ida_ff_g_st,jda_ff_g_st,kda_ff_g_st,
     &           lda_ff_g_st,mda_ff_g_st,s11_ff_g_st,
     &           s22_ff_g_st,s33_ff_g_st,s44_ff_g_st,
     &           d11_ff_g_st,d22_ff_g_st,d12_ff_g_st,
     &           e00_ff_g_st,rep_ff_g_st,ene_g_st_loc,
     &           f_a_x_st,myrep)
        else
          call calculate_go_forces_no_newton(n_ff_g_st,
     &           ida_ff_g_st,jda_ff_g_st,kda_ff_g_st,
     &           lda_ff_g_st,mda_ff_g_st,s11_ff_g_st,
     &           s22_ff_g_st,s33_ff_g_st,s44_ff_g_st,
     &           d11_ff_g_st,d22_ff_g_st,d12_ff_g_st,
     &           e00_ff_g_st,rep_ff_g_st,ene_g_st_loc,
     &           e_a_x_st(:,:,myrep),myrep)
        endif
         if(i_debug) write(*,*)'done calling go forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(17,mythrd)=tmtrx(17,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

         if(i_debug) write(*,*)'calling elec forces',mythrd
        if(i_do_YHL_nonbondeds) then
          call calculate_elec_forces_st_newton(n_ff_e_st,
     &           ida_ff_e_st,jda_ff_e_st,qqq_ff_e_st,
     &           ene_e_st_loc,f_a_x_st,f_e_x_st,myrep)
        else
          call calculate_elec_forces_st_no_newton(n_ff_e_st,
     &           ida_ff_e_st,jda_ff_e_st,qqq_ff_e_st,
     &           ene_e_st_loc,e_a_x_st(:,:,myrep),e_e_x_st(:,:,myrep),
     &           myrep)
        endif
         if(i_debug) write(*,*)'done calling elec forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(18,mythrd)=tmtrx(18,mythrd)+real(ftm2-ftm1)/real(ru)

C xxxxx Section
C xxxxx Calculate medium-range forces (at every num_fmd_stp)

        if(num_fmd_cnt.eq.num_fmd_stp) then

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

C first zero out the force and energy terms

         if(i_debug) write(*,*)'zeroing elec md forces',mythrd
          if(i_do_YHL_nonbondeds) then
            f_a_x_md=0.0 
            if(treecode_elec) f_e_x_md=0.0
          else
            do nnn=1,mov_num_loc
              e_a_x_md(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)=0.0
            enddo
            if(treecode_elec) then
              do nnn=1,mov_num_loc
              e_e_x_md(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)=0.0
              enddo
            endif
          endif

          ene_v_md_loc=0.0
          ene_g_md_loc=0.0
          ene_e_md_loc=0.0

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(19,mythrd)=tmtrx(19,mythrd)+real(ftm2-ftm1)/real(ru)

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

         if(i_debug) write(*,*)'calling vdw md forces',mythrd
          if(i_do_YHL_nonbondeds) then
            call calculate_vdw_forces_newton(n_ff_v_md,
     &           ida_ff_v_md,jda_ff_v_md,s11_ff_v_md,
     &           s22_ff_v_md,s33_ff_v_md,s44_ff_v_md,
     &           d11_ff_v_md,d22_ff_v_md,e00_ff_v_md,
     &           ene_v_md_loc,f_a_x_md)
          else
            call calculate_vdw_forces_no_newton(n_ff_v_md,
     &           ida_ff_v_md,jda_ff_v_md,s11_ff_v_md,
     &           s22_ff_v_md,s33_ff_v_md,s44_ff_v_md,
     &           d11_ff_v_md,d22_ff_v_md,e00_ff_v_md,
     &           ene_v_md_loc,e_a_x_md(:,:,myrep),myrep)
          endif
         if(i_debug) write(*,*)'done calling vdw md forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(20,mythrd)=tmtrx(20,mythrd)+real(ftm2-ftm1)/real(ru)

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

         if(i_debug) write(*,*)'calling go md forces',mythrd
          if(i_do_YHL_nonbondeds) then
            call calculate_go_forces_newton(n_ff_g_md,
     &           ida_ff_g_md,jda_ff_g_md,kda_ff_g_md,
     &           lda_ff_g_md,mda_ff_g_md,s11_ff_g_md,
     &           s22_ff_g_md,s33_ff_g_md,s44_ff_g_md,
     &           d11_ff_g_md,d22_ff_g_md,d12_ff_g_md,
     &           e00_ff_g_md,rep_ff_g_md,ene_g_md_loc,
     &           f_a_x_md,myrep)
          else
            call calculate_go_forces_no_newton(n_ff_g_md,
     &           ida_ff_g_md,jda_ff_g_md,kda_ff_g_md,
     &           lda_ff_g_md,mda_ff_g_md,s11_ff_g_md,
     &           s22_ff_g_md,s33_ff_g_md,s44_ff_g_md,
     &           d11_ff_g_md,d22_ff_g_md,d12_ff_g_md,
     &           e00_ff_g_md,rep_ff_g_md,ene_g_md_loc,
     &           e_a_x_md(:,:,myrep),myrep)
          endif
         if(i_debug) write(*,*)'done calling go md forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(22,mythrd)=tmtrx(22,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

         if(i_debug) write(*,*)'calling elec md forces',mythrd
          if(i_do_YHL_nonbondeds) then
            call calculate_elec_forces_md_newton(n_ff_e_md,n_ff_e_st,
     &           ida_ff_e_md,jda_ff_e_md,qqq_ff_e_md,
     &           ida_ff_e_st,jda_ff_e_st,qqq_ff_e_st,
     &           ene_e_md_loc,f_a_x_md,f_e_x_md,myrep)
          else
            call calculate_elec_forces_md_no_newton(n_ff_e_md,n_ff_e_st,
     &           ida_ff_e_md,jda_ff_e_md,qqq_ff_e_md,
     &           ida_ff_e_st,jda_ff_e_st,qqq_ff_e_st,
     &           ene_e_md_loc,e_a_x_md(:,:,myrep),e_e_x_md(:,:,myrep),
     &           myrep)
          endif
         if(i_debug) write(*,*)'done calling elec md forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(23,mythrd)=tmtrx(23,mythrd)+real(ftm2-ftm1)/real(ru)
 
        endif

C xxxxx Section
C xxxxx Calculate long-range forces (at every num_hyd_stp)
C note that is new - we used to do it every num_lst_stp...
C update - we've changed it back to num_lst_stp...

        if((num_lst_cnt.eq.num_lst_stp.or.update_nonbond_list).and.
     &      treecode_elec) then

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

          if(i_do_YHL_nonbondeds) then
            f_a_x_lg=0.0 
          else
            do nnn=1,mov_num_loc
              e_a_x_lg(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)=0.0
            enddo
          endif

          ene_e_lg_loc=0.0

         if(i_debug) write(*,*)'calling treecode forces',mythrd

C you need to get this code from a previous directory
C ctt
          call calculate_treecode_elec(mythrd,num_q_as,
     &           mov_num_loc,mov_beq_loc,mov_enq_loc,
     &           crg_of_atm_cpy,perm_cpy,
     &           tpengtar,tforce,
     &           c_q_x_cpy,c_q_y_cpy,c_q_z_cpy,
     &           c_q_x_tar,c_q_y_tar,c_q_z_tar)
         if(i_debug) write(*,*)'done calling treecode forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(24,mythrd)=tmtrx(24,mythrd)+real(ftm2-ftm1)/real(ru)
 
C TEMP MAY NEED TO REEXAMINE THIS CODE NOW FOR STRIPED NONBONDED LIST
C NOTE THAT THIS ENTIRE ROUTINE IS CURRENTLY COMMENTED OUT DUE TO US
C CHANGING THE FORCE ARRAYS TO ADD THE REPLICA DIMENSIONS

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

C note that for this particular routine we don't need separate copies
C for newton and no-newton - call same routine but with different args
C ctt
          if(i_do_YHL_nonbondeds) then
            call finalize_long_range_forces(
     &             mov_num_loc,mov_beq_loc,mov_enq_loc,
     &             num_q_as,num_tot_atms,
     &             f_e_x_st,f_e_x_md,f_a_x_lg,
     &             tpengtar,tforce,perm_cpy,
     &             ene_e_st_loc,ene_e_md_loc,ene_e_lg_loc)
          else
            call finalize_long_range_forces(
     &             mov_num_loc,mov_beq_loc,mov_enq_loc,
     &             num_q_as,num_tot_atms,
     &             e_e_x_st,e_e_x_md,e_a_x_lg,
     &             tpengtar,tforce,perm_cpy,
     &             ene_e_st_loc,ene_e_md_loc,ene_e_lg_loc)
          endif

c!$omp barrier ! f_e_x_st
c        if(mythrd.eq.1) then
c          do i=1,num_tot_atms
c            write(89,8998)i,f_e_x_st(1,i)
c8998        format(i10,f15.5)
c          enddo
c          stop
c        endif

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(27,mythrd)=tmtrx(27,mythrd)+real(ftm2-ftm1)/real(ru)
 
C!$omp barrier ! NOT NEEDED SINCE WE HAVE ONE LATER 

        endif 

C here we have to switch off the logical forcing us to update
C nonbonded lists and treecode every time there's a cleavage event

        if(update_nonbond_list) update_nonbond_list=.false.

C xxxxx Section
C xxxxx Add up all local forces                                       

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

c         write(*,*)'adding up forces',mythrd

C again, for this particular routine we don't need separate copies
C for newton and no-newton - call same routine but with different args

         if(i_debug) write(*,*)'calling add up forces',mythrd
        if(i_do_YHL_nonbondeds) then
          call add_up_forces(
     &         mov_num_loc,mov_beg_loc,mov_end_loc,
     &         f_a_x_bd,f_a_x_st,f_a_x_md,f_a_x_lg,f_u_x,
     &         f_a_x(:,:,myrep),myrep)
        else
          call add_up_forces(
     &         mov_num_loc,mov_beg_loc,mov_end_loc,
     &         f_a_x_bd,e_a_x_st(:,:,myrep),e_a_x_md(:,:,myrep),
     &         e_a_x_lg,f_u_x,
     &         f_a_x(:,:,myrep),myrep)
        endif
         if(i_debug) write(*,*)'done add up forces',mythrd

C add up osmotic pressure contributions at every step
C this isn't elegant code but it should get the job done

        if(walls.or.sphere.or.cylinder) then
!$omp critical
          do i=1,7
            osm_cur(i)=osm_cur(i)+osm_cur_loc(i)
          enddo
!$omp end critical
        endif 
         if(i_debug) write(*,*)'done add up osmotic ',mythrd

c         write(*,*)'done adding up forces',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(28,mythrd)=tmtrx(28,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

C only if using newton's law - at this point each thread has all forces
C calculated on this thread for its thread-local atoms. now, 
C each thread needs to copy its nonbonded forces into a 3D global
C force array so that other threads can access their entries; then,
C after a barrier, we collect the entries from the other threads

c         write(*,*)'adding up gax ',mythrd
         if(i_debug) write(*,*)'add up gax ',mythrd
        if(i_do_YHL_nonbondeds) then
          g_a_x(1:3,1:num_tot_atms,mythrd,myrep)=
     &      f_a_x_st(1:3,1:num_tot_atms)+
     &      f_a_x_md(1:3,1:num_tot_atms)+
     &      f_a_x_lg(1:3,1:num_tot_atms)
        endif
         if(i_debug) write(*,*)'done add up gax ',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(29,mythrd)=tmtrx(29,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp barrier ! DEFO NEEDED 

C calculate the cumulative osmotic pressure at every step
 
         if(i_debug) write(*,*)'add up osm_cur ',mythrd
          if(walls.or.sphere.or.cylinder) then
            if(mythrd.eq.1) then
              do i=1,7
                osm_cur(i)=osm_cur(i)/surface(i) ! convert F to P
                osm_cur(i)=osm_cur(i)*69478.58 ! convert to bar
                ee=osm_cur(i)-osm_ave(i)
                osm_ave(i)=osm_ave(i)+ee/real(nstep+1)
                osm_var(i)=osm_var(i)+ee*(osm_cur(i)-osm_ave(i))
c               if(i.eq.7) write(*,*)'check osmotic pres = ',osm_cur(i),
c    &                     osm_ave(i),osm_var(i),ee,real(nstep+1)
              enddo
            endif
          endif
         if(i_debug) write(*,*)'done add up osm_cur ',mythrd

C note that we don't need the following if we're not using newton's law

C not sure if the following statements are legal...
C note also that a barrier after this point is only needed if we're
C using hydrodynamics as other threads need to know all forces...

C note that the following could probably be speeded up if we limited the
C loops to only those thread-local atoms that have had forces calculated
C for them on other nodes - I have imin_thrd_int and imax_thrd_int
C calculated for that purpose but they don't work well when using
C mixed_hydrodynamics (we need versions for each of the mov_num_locs)

C note that the expense of this step does not diminish with added
C threads - it's effectively serial in cost and this is the main reason
C why the newton approach doesn't scale as well as no-newton - if we
C could figure out an easy way to only add on those terms that are
C non-zero it might be possible to speed this approach up - also note
C however that we'll still have the cost of zeroing the arrays at the
C beginning of the timestep (tmtrx(1,...): this is also serial in form

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

c       write(*,*)'adding up fax ',mythrd
         if(i_debug) write(*,*)'add up fax ',mythrd
        if(i_do_YHL_nonbondeds) then
          do nnn=1,mov_num_loc
            do jjj=1,mythrd-1
              f_a_x(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)=
     &        f_a_x(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)+
     &        g_a_x(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),jjj,myrep)
            enddo
            do jjj=mythrd+1,num_threads
              f_a_x(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)=
     &        f_a_x(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),myrep)+
     &        g_a_x(1:3,mov_beg_loc(nnn):mov_end_loc(nnn),jjj,myrep)
            enddo
          enddo
        endif
         if(i_debug) write(*,*)'done add up fax ',mythrd

!$omp barrier ! TEST DOES THIS FIX PROBLEMS - NO IT DOESN'T

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(30,mythrd)=tmtrx(30,mythrd)+real(ftm2-ftm1)/real(ru)
 
C transfer from f_a_x(:,:) to f_a_r(:) - may be better ways to do this

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

        if(i_debug) write(*,*)'add up far ',mythrd
        do nnn=1,mov_num_loc
          do iii=mov_beg_loc(nnn),mov_end_loc(nnn)
            ii3=iii*3-2
            f_a_r(ii3  ,myrep)=f_a_x(1,iii,myrep)
            f_a_r(ii3+1,myrep)=f_a_x(2,iii,myrep)
            f_a_r(ii3+2,myrep)=f_a_x(3,iii,myrep)
          enddo
        enddo

        if(i_debug) write(*,*)'add up fer ',mythrd
        do nnn=1,mov_num_loc
          if(mov_dyn_loc(nnn).eq.2) then ! if langevin
            do iii=mov_beg_loc(nnn),mov_end_loc(nnn)
              ii3=iii*3-2
              f_e_r(ii3  ,myrep)=ld_fac1(ii3  )*f_a_r(ii3  ,myrep)+
     &                           ld_fac2(ii3  )*v_e_r(ii3  ,myrep)
              f_e_r(ii3+1,myrep)=ld_fac1(ii3+1)*f_a_r(ii3+1,myrep)+
     &                           ld_fac2(ii3+1)*v_e_r(ii3+1,myrep)
              f_e_r(ii3+2,myrep)=ld_fac1(ii3+2)*f_a_r(ii3+2,myrep)+
     &                           ld_fac2(ii3+2)*v_e_r(ii3+2,myrep)
              v_e_r(ii3  ,myrep)=ld_fac3(ii3  )*v_e_r(ii3  ,myrep)+
     &                           ld_fac4(ii3  )*f_a_r(ii3  ,myrep)
              v_e_r(ii3+1,myrep)=ld_fac3(ii3+1)*v_e_r(ii3+1,myrep)+
     &                           ld_fac4(ii3+1)*f_a_r(ii3+1,myrep)
              v_e_r(ii3+2,myrep)=ld_fac3(ii3+2)*v_e_r(ii3+2,myrep)+
     &                           ld_fac4(ii3+2)*f_a_r(ii3+2,myrep)
            enddo
          endif
        enddo

C we now copy these forces into special arrays for hydrodynamic forces

        if(i_debug) write(*,*)'add up har ',mythrd
        do nnn=1,mov_num_loc
          if(mov_hyd_loc(nnn).eq.0) cycle
          if(mov_dyn_loc(nnn).eq.1) then
            ikk=atm_to_hyd_atm(mov_beg_loc(nnn))
            ik3=ikk*3
            ik3=ik3-mov_beg_loc(nnn)*3
            kkk=mov_hyd_loc(nnn)
            do ii3=mov_beg_loc(nnn)*3-2,mov_end_loc(nnn)*3
              h_a_r(ii3+ik3,kkk,myrep)=f_a_r(ii3,myrep)
            enddo
          elseif(mov_dyn_loc(nnn).eq.2) then
            ikk=atm_to_hyd_atm(mov_beg_loc(nnn))
            ik3=ikk*3
            ik3=ik3-mov_beg_loc(nnn)*3
            kkk=mov_hyd_loc(nnn)
            do ii3=mov_beg_loc(nnn)*3-2,mov_end_loc(nnn)*3
              h_a_r(ii3+ik3,kkk,myrep)=f_a_r(ii3,myrep)
              h_e_r(ii3+ik3,kkk,myrep)=f_e_r(ii3,myrep)
            enddo
          endif
        enddo
        if(i_debug) write(*,*)'done with har ',mythrd

!$omp barrier ! DEFO NEEDED with full_hydrodynamics

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(31,mythrd)=tmtrx(31,mythrd)+real(ftm2-ftm1)/real(ru)
 
C if we had a problem with a bond earlier we need to die now

        if(i_need_to_die_now.gt.0) then

!$omp single
          call clean_up_crash(myrep)
!$omp end single

C if i_continue_after_problem then we'll reinstate old coordinates and
C velocities, reset various counters and hope to continue :)

          if(.not.i_continue_after_problem) then
            write(*,*)'quitting :('
            stop
          endif
        endif

C we calculate go_pairs here cos we need to pass a BARRIER after the
C earlier operations performed this round

        if(i_debug) write(*,*)'calc go pairs ',mythrd
        if (necount.eq.neprint) then                     
          if(i_use_go_pairs) then
!$omp barrier ! we just added this to see if stops us getting Q=0
!$omp single
            call calculate_go_pairs
            if(rat.ge.Q_des) then
              write(*,87)mol_Q1,mol_Q2,rat,tot_time
87            format('for mols ',i5,' & ',i5,' Q = ',f10.5,
     &               ' at time = ',f20.5)
              call write_movie_file(movieframenumber,icrashed)
              write(*,*)'quitting but in the right way  :)'
              stop
            endif
!$omp end single
          endif
        endif
        if(i_debug) write(*,*)'done calc go pairs ',mythrd

C if appropriate, add up all of the energies - yet again there is room
C for improvement here - there's criticals and barriers here...

        if(necount.eq.neprint) then                     

          if(i_debug) write(*,*)'strt calc ener ',mythrd
!$omp barrier ! TEMPORARY FOR TIMINGS
          call system_clock(count=ftm1,count_rate=ru)

          if(i_debug) write(*,*)'cont calc ener ',mythrd
!$omp critical   
          do iii=1,num_f_ms
            ene_m(iii)=ene_m(iii)+
     &                 ene_pos_loc(iii)+
     &                 ene_bond_loc(iii)+
     &                 ene_angl_loc(iii)+
     &                 ene_dihe_loc(iii)
            ene_pos(iii)= ene_pos(iii)+ ene_pos_loc(iii)
            ene_bond(iii)=ene_bond(iii)+ene_bond_loc(iii)
            ene_angl(iii)=ene_angl(iii)+ene_angl_loc(iii)
            ene_dihe(iii)=ene_dihe(iii)+ene_dihe_loc(iii)
          enddo
          ene_nbnd_tot=ene_nbnd_tot+ene_nbnd_loc
          ene_gofav_tot=ene_gofav_tot+ene_gofav_loc
          ene_gounfav_tot=ene_gounfav_tot+ene_gounfav_loc
          etvst=etvst+ene_v_st_loc
          etgst=etgst+ene_g_st_loc
          etest=etest+ene_e_st_loc
          etvmd=etvmd+ene_v_md_loc
          etgmd=etgmd+ene_g_md_loc
          etemd=etemd+ene_e_md_loc
          etelg=etelg+ene_e_lg_loc
!$omp end critical

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(32,mythrd)=tmtrx(32,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp barrier ! NEEDED - INEXPENSIVE


!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)
         if(i_debug) write(*,*)'cont2 calc ener ',mythrd

!$omp master

          ep=ene_pos(1)
          eb=ene_bond(1)
          ea=ene_angl(1)
          ed=ene_dihe(1)
          do iii=2,num_f_ms
            ep=ep+ene_pos(iii)
            eb=eb+ene_bond(iii)
            ea=ea+ene_angl(iii)
            ed=ed+ene_dihe(iii)
          enddo
          et=eb+ea+ed+etvst+etgst+etest+
     &       etvmd+etgmd+etemd+etelg+ep+eu

C write out the cumulative osmotic pressure

          if(walls.or.sphere.or.cylinder) then
            write(25,399)tot_time,(osm_ave(i),sqrt(osm_var(i))/
     &                   real(nstep+1),i=1,7)
399         format('time = ',f15.5,' osmotic press. & stdev = ',14f15.8)
          endif

          if(shrinking_box) write(*,42)tot_time,m_size,s_size
42        format('box at time ',f15.3,' length = ',f10.3,
     &           ' radius = ',f10.3)

          write(*,39)
39        format('                         time    total_ene',
     &           '    bond_ene    angl_ene    dihe_ene   ',
     &           'vdw_short    go_short  elec_short  vdw_medium',
     &           '   go_medium elec_medium   elec_long    elec_tot  ',
     &           'restraints sphere_size')

          write(*,41)tot_time,et,eb,ea,ed,
     &      etvst,etgst,etest,etvmd,etgmd,etemd,etelg,
     &      etest+etemd+etelg,ep,s_size

41        format('time & E  ',f20.3,14f12.3)

C calculate average and variance of the energy

c         ene_den=ene_den+1.0d0
c         ee=et-ene_ave
c         ene_ave=ene_ave+ee/real(ene_den)
c         ene_var=ene_var+ee*(et-ene_ave)

c         write(*,'(a17,6f15.5)')'rand/E ave/msd = ',rand_ave/rand_den,
c    &                                               rand_msd/rand_den,
c    &                                               ene_ave,
c    &                                          sqrt(ene_var/ene_den),
c    &                                               qqq_ave,
c    &                                          sqrt(qqq_var/qqq_den)

C write out the Go contact numbers

          necount=0

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(33,mythrd)=tmtrx(33,mythrd)+real(ftm2-ftm1)/real(ru)
 

!$omp end master
!$omp barrier ! NEEDED - INEXPENSIVE

        endif ! for writing out energy
         if(i_debug) write(*,*)'done calc ener ',mythrd

C at this point we have the option to quit - note that we have
C up-to-date coordinates here *and* for what it's worth we also have the
C forces and energies appropriate to these coordinates

        if(nstep.ge.ntimetotal) then
          if(mythrd.eq.1) then
            write(*,*)'quitting at time = ',tot_time
          endif
          goto 555
        endif

C if we got here then we're continuing with the sim, so make the move


!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

        if(i_debug) write(*,*)'strt move atom ',mythrd
        call move_atoms(num_hyd_cnt,myrep,
     &                  mov_num_loc,mov_typ_loc,mov_beg_loc,
     &                  mov_end_loc,mov_dyn_loc,mov_hyd_loc)

        if(i_debug) write(*,*)'done move atom ',mythrd
        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(34,mythrd)=tmtrx(34,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp barrier ! NEEDED BUT MAYBE WE CAN PUT LATER? - WE MUST HAVE ALL
              ! THINGS THAT DEPEND ON THESE COUNTERS FINISHED...

C do lincs stuff here ! we assume that we're post-barrier here

        call system_clock(count=ftm1,count_rate=ru)
        if(i_do_lincs) then

C note that lbnd_thrd = 0 indicates that one of the atoms in this bond
C is held by another thread so the bond is done by that thread instead
C lbnd_thrd=2 indicates that this is a bond that is done by this thread
C but which also appears on another thread
C lbnd_thrd=1 is a bond that is not shared by another thread

          if(i_debug) write(*,*)'starting lincs 1st stage'
          do ijk=1,num_bnd_thrd
            iii=ibnd_thrd(ijk)
            blnc(1,iii)=c_a_x_old(1,nbb1_new(iii),myrep)-
     &                  c_a_x_old(1,nbb2_new(iii),myrep)
            blnc(2,iii)=c_a_x_old(2,nbb1_new(iii),myrep)-
     &                  c_a_x_old(2,nbb2_new(iii),myrep)
            blnc(3,iii)=c_a_x_old(3,nbb1_new(iii),myrep)-
     &                  c_a_x_old(3,nbb2_new(iii),myrep)
            rlen=sqrt(blnc(1,iii)**2+blnc(2,iii)**2+blnc(3,iii)**2)
            blnc(1,iii)=blnc(1,iii)/rlen
            blnc(2,iii)=blnc(2,iii)/rlen
            blnc(3,iii)=blnc(3,iii)/rlen
          enddo

!$omp barrier

          if(i_debug) write(*,*)'starting lincs 2nd stage'
          do ijk=1,num_bnd_thrd
            if(lbnd_thrd(ijk).eq.0) cycle
            iii=ibnd_thrd(ijk)
            do nnn=1,ncc(iii)
              kkk=icon(nnn,iii)
              alnc(nnn,iii)=coef(nnn,iii)*
     &         (blnc(1,iii)*blnc(1,kkk)+
     &          blnc(2,iii)*blnc(2,kkk)+ 
     &          blnc(3,iii)*blnc(3,kkk))
            enddo
            na1=nbb1_new(iii)
            na2=nbb2_new(iii)
            rhsc(iii,1)=sdia(iii)*
     &                 (blnc(1,iii)*
     &                 (c_a_x(1,na1,myrep)-c_a_x(1,na2,myrep))+
     &                  blnc(2,iii)*
     &                 (c_a_x(2,na1,myrep)-c_a_x(2,na2,myrep))+
     &                  blnc(3,iii)*
     &                 (c_a_x(3,na1,myrep)-c_a_x(3,na2,myrep))-
     &                  r0bb_new(iii))
            solc(iii)=rhsc(iii,1)
          enddo

!omp barrier ! probably needed

          if(i_debug) write(*,*)'starting lincs 3rd stage'
          jjj=2
          do nrec=1,8 ! 4 seems good enough
            do ijk=1,num_bnd_thrd
              if(lbnd_thrd(ijk).eq.0) cycle
              iii=ibnd_thrd(ijk)
              rhsc(iii,jjj)=0.0
              do nnn=1,ncc(iii)
                rhsc(iii,jjj)=rhsc(iii,jjj)+
     &                        alnc(nnn,iii)*rhsc(icon(nnn,iii),3-jjj)
              enddo
              solc(iii)=solc(iii)+rhsc(iii,jjj)
            enddo
            jjj=3-jjj
          enddo

!omp barrier ! added to see if it helps

          do ijk=1,num_bnd_thrd
            if(lbnd_thrd(ijk).eq.0) cycle
            iii=ibnd_thrd(ijk)
            na1=nbb1_new(iii)        
            na2=nbb2_new(iii)        
            c_a_x(1,na1,myrep)=c_a_x(1,na1,myrep)-
     &                         blnc(1,iii)*sdia(iii)*solc(iii)
            c_a_x(2,na1,myrep)=c_a_x(2,na1,myrep)-
     &                         blnc(2,iii)*sdia(iii)*solc(iii)
            c_a_x(3,na1,myrep)=c_a_x(3,na1,myrep)-
     &                         blnc(3,iii)*sdia(iii)*solc(iii)
            c_a_x(1,na2,myrep)=c_a_x(1,na2,myrep)+
     &                         blnc(1,iii)*sdia(iii)*solc(iii)
            c_a_x(2,na2,myrep)=c_a_x(2,na2,myrep)+
     &                         blnc(2,iii)*sdia(iii)*solc(iii)
            c_a_x(3,na2,myrep)=c_a_x(3,na2,myrep)+
     &                         blnc(3,iii)*sdia(iii)*solc(iii)
          enddo

!omp barrier ! may not be needed - looks entirely thread-local

          if(i_debug) write(*,*)'starting lincs 4th stage'
          do ijk=1,num_bnd_thrd
            if(lbnd_thrd(ijk).eq.0) cycle
            iii=ibnd_thrd(ijk)
            na1=nbb1_new(iii)
            na2=nbb2_new(iii)
            plnc=sqrt(2.0*r0bb_new(iii)**2-
     &          (c_a_x(1,na1,myrep)-c_a_x(1,na2,myrep))**2-
     &          (c_a_x(2,na1,myrep)-c_a_x(2,na2,myrep))**2-
     &          (c_a_x(3,na1,myrep)-c_a_x(3,na2,myrep))**2)
            rhsc(iii,1)=sdia(iii)*(r0bb_new(iii)-plnc)
            solc(iii)=rhsc(iii,1)
          enddo

!omp barrier ! probably needed - uses rhs for connected bonds

          jjj=2
          do nrec=1,8 ! 4 seems good enough
            do ijk=1,num_bnd_thrd
              if(lbnd_thrd(ijk).eq.0) cycle
              iii=ibnd_thrd(ijk)
              rhsc(iii,jjj)=0.0
              do nnn=1,ncc(iii)
                rhsc(iii,jjj)=rhsc(iii,jjj)+
     &                        alnc(nnn,iii)*rhsc(icon(nnn,iii),3-jjj)
              enddo
              solc(iii)=solc(iii)+rhsc(iii,jjj)
            enddo
            jjj=3-jjj
          enddo
          
!omp barrier ! added to see if it helps

          if(i_debug) write(*,*)'starting lincs 5th stage'
          do ijk=1,num_bnd_thrd
            if(lbnd_thrd(ijk).eq.0) cycle
            iii=ibnd_thrd(ijk)
            na1=nbb1_new(iii)
            na2=nbb2_new(iii)
            c_a_x(1,na1,myrep)=c_a_x(1,na1,myrep)-
     &                         blnc(1,iii)*sdia(iii)*solc(iii)
            c_a_x(2,na1,myrep)=c_a_x(2,na1,myrep)-
     &                         blnc(2,iii)*sdia(iii)*solc(iii)
            c_a_x(3,na1,myrep)=c_a_x(3,na1,myrep)-
     &                         blnc(3,iii)*sdia(iii)*solc(iii)
            c_a_x(1,na2,myrep)=c_a_x(1,na2,myrep)+
     &                         blnc(1,iii)*sdia(iii)*solc(iii)
            c_a_x(2,na2,myrep)=c_a_x(2,na2,myrep)+
     &                         blnc(2,iii)*sdia(iii)*solc(iii)
            c_a_x(3,na2,myrep)=c_a_x(3,na2,myrep)+
     &                         blnc(3,iii)*sdia(iii)*solc(iii)
          enddo
          if(i_debug) write(*,*)'finishng lincs'

!$omp barrier ! PROBABLY NEED ONE AFTER LINCS IS DONE?

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(39,mythrd)=tmtrx(39,mythrd)+real(ftm2-ftm1)/real(ru)
 
!$omp single ! TEMPORARY for checking
           tot_dev=0.0
           do iii=1,num_sys_bonds
             na1=nbb1_new(iii)
             na2=nbb2_new(iii)
             dist=(c_a_x(1,na1,myrep)-c_a_x(1,na2,myrep))**2+
     &            (c_a_x(2,na1,myrep)-c_a_x(2,na2,myrep))**2+
     &            (c_a_x(3,na1,myrep)-c_a_x(3,na2,myrep))**2
             dist=sqrt(dist)
             dist=abs(dist-r0bb_new(iii))
             tot_dev=tot_dev+dist
           enddo
           ave_dev=tot_dev/real(num_sys_bonds)
          write(*,*)'deviation = ',ave_dev,' on ',num_sys_bonds,' bonds'
!$omp end single

        endif

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

        if(i_debug) write(*,*)'strt update  ',mythrd
!$omp master 
        call update_counters
!$omp end master
        if(i_debug) write(*,*)'done update  ',mythrd

        call system_clock(count=ftm2,count_rate=ru)
        tmtrx(35,mythrd)=tmtrx(35,mythrd)+real(ftm2-ftm1)/real(ru)
 
C!$omp barrier ! APPARENTLY NEEDED FOR NEXT STEP BUT WHY?
               ! WILL COMMENT IT OUT FOR NOW

!$omp barrier ! TEMPORARY FOR TIMINGS
        call system_clock(count=ftm1,count_rate=ru)

C go back and do another step of the simulation

        if(i_debug) write(*,*)'strt new stp ',mythrd
        goto 222

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCC main loop ends here CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

555     continue


!$omp barrier
!$omp master

        time_max=0.0
        do i=1,num_threads
          time_cur=0.0
          do j=1,40
            time_cur=time_cur+tmtrx(j,i)
          enddo
          if(time_cur.gt.time_max) time_max=time_cur
        enddo

        open(unit=28,file='timings_percents.txt',status='unknown')
        open(unit=29,file='timings_absolute.txt',status='unknown')

        do j=1,40 
          write(28,546)j,(100*tmtrx(j,i)/time_max,i=1,num_threads)
546       format(i5,100f6.1)
          write(29,547)j,(tmtrx(j,i),i=1,num_threads)
547       format(i5,100f10.2)
        enddo

        close(28)
        close(29)

!$omp end master
!$omp barrier    

!$omp master

        write(*,*)'now writing restart.files using thread # ',mythrd
        write(*,*)'tot_time = ',tot_time

        do iii=1,num_reps

          restart_file_name='restart.file.XXX'
          write(atemp,'(i3)')iii
          do j=1,3
            if(atemp(j:j).eq.' ') atemp(j:j)='0'
          enddo
          restart_file_name(14:16)=atemp

          open(unit=21,file=restart_file_name,status='unknown')
          write(21,10)tot_time                             
10        format('time = ',f20.5,' ps')          
          do kkk=1,num_t_ms                               
            write(21,20)kkk,0.0,0.0,0.0
            do lll=1,num_atms_pt(i_f_tp(kkk))                    
              write(21,731)kkk,lll,
     &        c_a_x(1,id_beg(kkk)+lll,iii),
     &        c_a_x(2,id_beg(kkk)+lll,iii),
     &        c_a_x(3,id_beg(kkk)+lll,iii)
            enddo                                    
          enddo                                     
          close(21)                                
731       format(2i8,3f20.5)                      
734       format(2i8,3f20.5,i8)
20        format(i8,'       0',3f20.5)          

        enddo

      call xdrfclose(39,ret)

!$omp end master
!$omp barrier   
       
!$omp barrier ! conditional on finishing
      write(*,*)'thread # ',mythrd,' done'

C                                                                                                     
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                
C XX                                                 XX                                                
C XX here's where the main loop over time_steps ends XX                                                
C XX                                                 XX                                                
C XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                                
C                                                                                                     
!     ***** Deinitialize *****
!     errcode=vsldeletestream( stream )
!     call CheckVslError(errcode)
!$omp end parallel 
        write(*,*)'thread # ',mythrd,' completely done'

        call system_clock(count=ktime2,count_rate=ru)
        time_total=real(ktime2-ktime1)/real(ru)

        write(*,*)                                                  
        write(*,*)'Time_total                 = ',Time_total
        write(*,*)          

C                                                                                                     
C exit and end                                                                                        
C               
      stop                                                                                           
      end                                                                 
                                                                         
                                                                                                      
                                                      
      SUBROUTINE WrtRstrtBk(movieframenumber,mirep)         
C ----------------------------------------------------------------------
C Tyson 011807 -                                                             
C  Subroutine that reads in the current restart file and copies the     
C  data into a backup restart file in the RESTARTS directory. The file  
C  corresponds to the movie file number in the MOVIE directory.              
C ----------------------------------------------------------------------
      IMPLICIT NONE                               
      integer :: mirep,j   ! (added by AHE)                          
      character*3 :: atemp ! (added by AHE)
      CHARACTER(LEN=26) :: restart_file_name ! (added by AHE)               
      INTEGER :: movieframenumber,ios=0,i                                 
      CHARACTER(LEN=26) :: tmpFileName                                      
      CHARACTER(LEN=35) :: filename                                         
      CHARACTER(LEN=80) :: tmpStr                                           
      CHARACTER(LEN=9) :: tmpStrTwo,FileNmStr                              

      restart_file_name="restart.file.XXX"
      tmpFileName = "RESTARTS/restart.file.XXX."                            
      write(atemp,'(i3)')mirep
      do j=1,3
        if(atemp(j:j).eq.' ') atemp(j:j)='0'
      enddo
      restart_file_name(14:16)=atemp
      tmpFileName(23:25)=atemp

      tmpStrTwo = FileNmStr(movieframenumber)                               
      filename=tmpFileName//tmpStrTwo   
                                    
301   continue
      OPEN(UNIT=201+mirep,FILE=restart_file_name,
     &     STATUS='unknown',err=301)            
302   continue
      OPEN(UNIT=301+mirep,FILE=filename,STATUS='unknown',err=302)                  
      DO                                                                    
        READ(201+mirep,'(80a)',IOSTAT=ios) tmpStr                                  
        IF (ios.gt.0) THEN                                                  
          WRITE(*,*) "Error ",ios," reading restart.file!"                  
          goto 301
c         EXIT                                                              
        ELSEIF (ios.lt.0) THEN                                              
          EXIT                                                              
        ELSEIF (ios.eq.0) THEN                                              
          WRITE(301+mirep,'(80a)') tmpStr                                          
        ENDIF                                                               
      ENDDO                                                                 
      CLOSE(201+mirep)                                                             
      CLOSE(301+mirep)                                                             
      return                                                                
      end                                                                   
                                                                            
                                                                            
      FUNCTION FileNmStr(movieFrameNumber)                                   
C ----------------------------------------------------------------------
C Tyson 011807 -                                                             
C  Function FileNmStr calculates the movie and restart file numbering   
C  scheme given the movieFrameNumber that is currently being reported   
C ----------------------------------------------------------------------
      IMPLICIT NONE                                                          
      INTEGER :: iremainder                                                  
      INTEGER,INTENT( IN ) :: movieFrameNumber                              
      CHARACTER(LEN=10) :: numbers                                           
      CHARACTER(LEN=9) :: tmpStr,FileNmStr                                   
      iremainder=0                                                           
      numbers(1:10)='0123456789'                                             
      tmpStr(1:9)='000000000'                                                
      if(movieFrameNumber.gt.99999999)then                                   
         iremainder=mod(movieFrameNumber,1000000000)/100000000               
         tmpStr(1:1)=                                                       
     &        numbers(iremainder+1                                          
     &        :iremainder+1)                                                
      endif                                                                 
      if(movieFrameNumber.gt.9999999)then                                   
         iremainder=mod(movieFrameNumber,100000000)/10000000                
         tmpStr(2:2)=                                                       
     &        numbers(iremainder+1                                          
     &        :iremainder+1)                                                
      endif                                                                  
      if(movieFrameNumber.gt.999999)then                                     
         iremainder=mod(movieFrameNumber,10000000)/1000000                   
         tmpStr(3:3)=                                                        
     &        numbers(iremainder+1                                           
     &        :iremainder+1)                                                 
      endif                                                                  
      if(movieFrameNumber.gt.99999)then                                      
         iremainder=mod(movieFrameNumber,1000000)/100000                     
         tmpStr(4:4)=                                                        
     &        numbers(iremainder+1                                           
     &        :iremainder+1)                                                 
      endif                                                                  
      if(movieFrameNumber.gt.9999)then                                       
         iremainder=mod(movieFrameNumber,100000)/10000                       
         tmpStr(5:5)=                                                        
     &        numbers(iremainder+1                                           
     &        :iremainder+1)                                                 
      endif                                                                  
      if(movieFrameNumber.gt.999)then                                        
         iremainder=mod(movieFrameNumber,10000)/1000                         
         tmpStr(6:6)=                                                        
     &        numbers(iremainder+1                                           
     &        :iremainder+1)                                                 
      endif                                                                  
      if(movieFrameNumber.gt.99)then                                         
         iremainder=mod(movieFrameNumber,1000)/100                           
         tmpStr(7:7)=                                                        
     &        numbers(iremainder+1                                           
     &        :iremainder+1)                                                 
      endif                                                                  
      if(movieFrameNumber.gt.9)then                                          
         iremainder=mod(movieFrameNumber,100)/10                             
         tmpStr(8:8)=                                                        
     &        numbers(iremainder+1                                           
     &        :iremainder+1)                                                 
      endif                                                                  
      if(movieFrameNumber.gt.0)then                                          
         iremainder=mod(movieFrameNumber,10)                                 
         tmpStr(9:9)=                                                        
     &        numbers(iremainder+1                                           
     &        :iremainder+1)                                                 
      endif                                                                  
      FileNmStr=tmpStr                                                       
      END FUNCTION FileNmStr                                                 
                                                                             
                                                                             
                                                                             

