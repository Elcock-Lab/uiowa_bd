

      subroutine calculate_treecode_elec(mythrd_loc,num_q_zs,
     &             miv_num_loc,miv_beq_loc,miv_enq_loc,
     &             crg_of_atm_cpz,perm_cpz,
     &             zpengtar,zforce,
     &             c_q_x_cpz,c_q_y_cpz,c_q_z_cpz,
     &             c_q_x_taz,c_q_y_taz,c_q_z_taz)

      use allocatable_arrays         
      implicit real(a-h,o-z)        

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
                                      
C first, make a private copy of all coordinates
C note that we make sure that all coords are from -xlen/2 to +xlen/2
C prior to passing them in the case of periodic boundaries (i_pbc=1)

      if(i_pbc.eq.1) then
        do jjj=1,num_q_as       
          iii=q_info(jjj)%ida
          c_q_x_cpz(jjj)=c_a_x(1,iii,1)-
     &        xlen*anint(c_a_x(1,iii,1)*xinv)
          c_q_y_cpz(jjj)=c_a_x(2,iii,1)-
     &        ylen*anint(c_a_x(2,iii,1)*yinv)
          c_q_z_cpz(jjj)=c_a_x(3,iii,1)-
     &        zlen*anint(c_a_x(3,iii,1)*zinv)
          c_q_x_taz(jjj)=c_q_x_cpz(jjj)
          c_q_y_taz(jjj)=c_q_y_cpz(jjj)
          c_q_z_taz(jjj)=c_q_z_cpz(jjj)
          crg_of_atm_cpz(jjj)=crg_of_atm(iii)
          perm_cpz(jjj)=jjj
        enddo
      else
        do jjj=1,num_q_as       
          iii=q_info(jjj)%ida
          c_q_x_cpz(jjj)=c_a_x(1,iii,1)
          c_q_y_cpz(jjj)=c_a_x(2,iii,1)
          c_q_z_cpz(jjj)=c_a_x(3,iii,1)
          c_q_x_taz(jjj)=c_a_x(1,iii,1)
          c_q_y_taz(jjj)=c_a_x(2,iii,1)
          c_q_z_taz(jjj)=c_a_x(3,iii,1)
          crg_of_atm_cpz(jjj)=crg_of_atm(iii)
          perm_cpz(jjj)=jjj
        enddo
      endif

C we should remember that we do not use periodic boundaries when we
C construct the tree. we *do* however use them when we calculate the
C interactions using the tree. I *think* that this is correct...

C note that each thread is calling the treecode independently - they're
C all passing all of the charges to the tree as sources, but each thread
C only asks for the forces/energies on charges ibeg_q_a to iend_q_a

C note also that I pass a lot of other stuff during the call to the treecode
C - a lot of these relate to periodic boundary conditions, salt
C concentration (in the case of a Debye-Huckel electrostatic model), and
C a grid-based Ewald electrostatics model. I checked my implementation
C of the latter worked so long as I stuck to assuming that only
C monopole-monopole interactions were important for periodic images. However,
C this was before I made a modified version of the treecode that used
C your Yukawa code - I therefore don't vouch for its correctness now. I
C vouch only for the correctness of the non-periodic cases and/or
C periodic-cases-that-use-only-minimum-image-interactions...

c           write(*,*)'calling treecode ',mythrd_loc
            call treecode3d_targ_openmp
     &                        (c_q_x_taz,
     &                         c_q_y_taz,
     &                         c_q_z_taz,
     &                         num_q_as,  
     &                         miv_num_loc,
     &                         miv_beq_loc,
     &                         miv_enq_loc,
     &                         zpengtar,
     &                         zforce,
     &                         c_q_x_cpz,
     &                         c_q_y_cpz,
     &                         c_q_z_cpz,
     &                         crg_of_atm_cpz,
     &                         perm_cpz,
     &                         num_q_as,
     &                         2,  ! iflag
     &                         treecode_order,
     &                         treecode_theta,
     &                         treecode_shrink,
     &                         treecode_maxatm,
     &                         0.0, ! dist_tol
     &                         timetree,
     &                         prinfo,  
     &                         num_q_as,
     &                         num_q_as,
     &                         mythrd_loc,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)
c           write(*,*)'done calling treecode ',mythrd_loc

      return
      end

