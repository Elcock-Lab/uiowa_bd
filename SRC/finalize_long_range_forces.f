
      subroutine finalize_long_range_forces(
     &             miv_num_loc,miv_beq_loc,miv_enq_loc,
     &             num_q_zs,num_tot_atmz,
     &             fff_ll_e_st,fff_ll_e_md,fff_ll_e_lg,
     &             zpengtar,zforce,perm_cpz,
     &             eee_e_st_loc,eee_e_md_loc,eee_e_lg_loc)

      use allocatable_arrays         
      implicit real(a-h,o-z)        

      integer, intent(in)     :: miv_num_loc
      integer, intent(in)     :: miv_beq_loc(miv_num_loc)
      integer, intent(in)     :: miv_enq_loc(miv_num_loc)
      integer, intent(in)     :: num_q_zs                        ! num_q_as
      integer, intent(in)     :: perm_cpz(1:num_q_zs)
      real,    intent(inout)  :: fff_ll_e_st(1:3,1:num_tot_atmz) ! fff_ll_e_st
      real,    intent(inout)  :: fff_ll_e_md(1:3,1:num_tot_atmz) ! fff_ll_e_md
      real,    intent(inout)  :: fff_ll_e_lg(1:3,1:num_tot_atmz) ! fff_ll_e_lg
      real,    intent(inout)  :: zpengtar(1:num_q_zs)            ! zpengtar
      real,    intent(inout)  :: zforce(1:3,1:num_q_zs)          ! zforce
      real,    intent(in)     :: eee_e_st_loc                    ! 
      real,    intent(in)     :: eee_e_md_loc                    ! 
      real,    intent(inout)  :: eee_e_lg_loc                    ! 

C keep the long-range force as the treecode force minus the
C short-range and medium-range electrostatic forces that were already
C calculated by the nonbonded lists.

      do nnn=1,miv_num_loc
        do jkl=miv_beq_loc(nnn),miv_enq_loc(nnn)

        iii=q_info(jkl)%ida

C note that we first need to calculate the electrostatic contribution
C from charges that are bonded to each other- these are calculated by 
C the treecode but are omitted from the regular short-range electrostatic 
C calculations. again, this might seem a bit weird, but it's usual in
C molecular mechanics simulations to not model electrostatic
C interactions between atoms that are bonded, e.g. in a hypothetical
C molecule such as:
C                             A-B-C-D-E
C                             1 2 3 4 5
C 
C the electrostatic interactions between A-B, A-C and A-D for example, would
C not be computed but those of A-E would be. (some MD force fields *do*
C calculate electrostatic interactions for A-D - so-called 1-4
C interactions, but we don't)

C note that 'iii' and 'jjj' are actual atom numbers here
 
C note also that we add the force to 'fff_ll_e_st','f_e_y_st','f_e_z_st'
C etc - we could have made a separate accumulator for this correction term 
C but it's cleaner just to add it to 'fff_ll_e_st' etc instead. note that we do 
C *not* intend to actually use these terms in the move step - we just need 
C to add it here so that we can *subtract* it from the treecode forces (and
C fff_ll_e_st isn't used after this point anyway)

C note also that even though we may be calculating the short-range
C interactions using our own direct code, the same interactions are
C calculated (possibly with minor errors) in the treecode.
C this means that those errors will remain within the treecode's
C contribution (to the long-range term). there's nothing to be done
C about that - and it's almost certainly inconsequential - but it's worth 
C being aware of it.

C finally, note that we should only do this for atoms that are
C not part of groups that don't 'see' other groups of atoms 
C (i.e. where 'i_res_no_nb groups' > 0) as these atoms will already have 
C their bonded terms cancelled out by upengtar...

        if(i_res_no_nb(iii).eq.0) then

C for all normal situations we will definitely reach this point...
C I have separate loops depending on whether we're using regular
C Coulombic interactions or Debye-Huckel screened interactions. this is
C another of those bifurcations that there's strictly no need for...

C note that here we do not make use of Newton's third law so we have to
C do the multiplying by 0.5 for the energy as we're adding it on twice

          if(implicit_ions) then
            do kkk=1,q_info(jkl)%num_bnd
              klm=q_info(jkl)%jda(kkk) ! this is a charge #
              jjj=q_info(klm)%ida      ! this is the atom #
              if(iii.eq.jjj) cycle ! needed for rings in DNA
              qqq=332.080*crg_of_atm(iii)*crg_of_atm(jjj)
              w1=c_a_x(1,iii,1)-c_a_x(1,jjj,1) 
              w2=c_a_x(2,iii,1)-c_a_x(2,jjj,1) 
              w3=c_a_x(3,iii,1)-c_a_x(3,jjj,1) 
              dist2=w1*w1+w2*w2+w3*w3
              dist=sqrt(dist2)     
              dist_inv=1.0/dist  
              temp=qqq*exp(-kappa*(dist-r_ion))
              u_v=temp*dielectric_inv*dist_inv*
     &            one_over_one_plus_kappa_rion
              f_v=u_v*(kappa+dist_inv)*dist_inv
              eee_e_lg_loc=eee_e_lg_loc-u_v*0.5
              fff_ll_e_st(1,iii)=fff_ll_e_st(1,iii)+f_v*w1    
              fff_ll_e_st(2,iii)=fff_ll_e_st(2,iii)+f_v*w2  
              fff_ll_e_st(3,iii)=fff_ll_e_st(3,iii)+f_v*w3   
            enddo

C here's the loop when we're using straight Coulombic interactions

          else 
            do kkk=1,q_info(jkl)%num_bnd
              klm=q_info(jkl)%jda(kkk) ! this is a charge #
              jjj=q_info(klm)%ida      ! this is the atom #
              if(iii.eq.jjj) cycle ! needed for rings in DNA
              qqq=332.080*crg_of_atm(iii)*crg_of_atm(jjj)
              w1=c_a_x(1,iii,1)-c_a_x(1,jjj,1) 
              w2=c_a_x(2,iii,1)-c_a_x(2,jjj,1) 
              w3=c_a_x(3,iii,1)-c_a_x(3,jjj,1) 
              dist2_inv=1.0/(w1*w1+w2*w2+w3*w3)
              dist_inv=sqrt(dist2_inv)       
              u_v=qqq*dielectric_inv*dist_inv
              f_v=u_v*dist2_inv
              eee_e_lg_loc=eee_e_lg_loc-u_v*0.5
              fff_ll_e_st(1,iii)=fff_ll_e_st(1,iii)+f_v*w1    
              fff_ll_e_st(2,iii)=fff_ll_e_st(2,iii)+f_v*w2  
              fff_ll_e_st(3,iii)=fff_ll_e_st(3,iii)+f_v*w3   
            enddo
          endif  

C close the 'if' statement over 'i_res_no_nb_grps'

        endif  

C if using ewald_elec_grid add each charge's interaction with its own images
C note that this has not been vetted for use with i_res_no_nb or ifree

        if(ewald_elec_grid) then
          u_p=ewald_elec_ei_with_iprime*
     &        crg_of_atm(iii)*crg_of_atm(iii)
        else
          u_p=0.0
        endif

C remember to remove the interactions from i_res_no_nb>0 atom pairs
C note that we don't do this for 'normal' atoms

        if(i_res_no_nb(iii).gt.0) then

          write(*,*)'AHE NEEDS TO PASS UFORCE TO SUBROUTINE'
          write(*,*)'SO THIS IS NOT YET FULLY IMPLEMENTED'
          stop ! TEMPORARY

          fff_ll_e_lg(1,iii)=                
     &      zforce(1,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      uforce(1,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      fff_ll_e_st(1,iii)-fff_ll_e_md(1,iii)
          fff_ll_e_lg(2,iii)=                
     &      zforce(2,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      uforce(2,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      fff_ll_e_st(2,iii)-fff_ll_e_md(2,iii)
          fff_ll_e_lg(3,iii)=                
     &      zforce(3,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      uforce(3,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      fff_ll_e_st(3,iii)-fff_ll_e_md(3,iii)
          eee_e_lg_loc=eee_e_lg_loc+
     &      zpengtar(jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      upengtar(jkl)*crg_of_atm(iii)*332.08*dielectric_inv+
     &      u_p
c       write(*,8301)iii,jkl,zpengtar(jkl),upengtar(jkl),eee_e_lg_loc,
c    &                       zpengtar(jkl)-upengtar(jkl)
8301    format('Check IRES ON  ',2i5,4f12.5) 

        else

          fff_ll_e_lg(1,iii)=                
     &      zforce(1,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      fff_ll_e_st(1,iii)-fff_ll_e_md(1,iii)
          fff_ll_e_lg(2,iii)=                
     &      zforce(2,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      fff_ll_e_st(2,iii)-fff_ll_e_md(2,iii)
          fff_ll_e_lg(3,iii)=                
     &      zforce(3,jkl)*crg_of_atm(iii)*332.08*dielectric_inv-
     &      fff_ll_e_st(3,iii)-fff_ll_e_md(3,iii)
          eee_e_lg_loc=eee_e_lg_loc+
     &      zpengtar(jkl)*crg_of_atm(iii)*332.08*dielectric_inv+
     &      u_p
c       write(*,8302)iii,jkl,zpengtar(jkl),upengtar(jkl),eee_e_lg_loc,
c    &                       zpengtar(jkl)-upengtar(jkl)
8302    format('Check IRES OFF ',2i5,4f12.5) 

c       write(*,5129)jkl,iii,    
c    &     zforce(1,jkl)*crg_of_atm(iii)*332.08*dielectric_inv,
c    &     zforce(2,jkl)*crg_of_atm(iii)*332.08*dielectric_inv,
c    &     zforce(3,jkl)*crg_of_atm(iii)*332.08*dielectric_inv
5129    format('check zforces ',2i6,3f12.5)
c       write(*,*)crg_of_atm(iii)
c       write(*,*)fff_ll_e_st(1,iii)
c       write(*,*)fff_ll_e_md(1,iii)
c       write(*,*)fff_ll_e_lg(1,iii)
c       write(*,*)zforce(1,jkl)
c       write(*,*)perm_cpz(iii)
c       write(*,5128)jkl,iii,
c    &     crg_of_atm(iii),fff_ll_e_st(1,iii),fff_ll_e_md(1,iii),
c    &                     fff_ll_e_lg(1,iii),
c    &     zforce(1,jkl)*crg_of_atm(iii)*332.08*dielectric_inv,
c    &     zforce(1,jkl),
c    &                perm_cpz(iii)
5128    format('check me out ',2i5,6f15.5,i5)

        endif
        
C increment the energy with the energy of this atom (note u_p is added)

        enddo
      enddo

c     write(*,*)'ENDING FINALIZE!!!',
c    &   eee_e_st_loc,eee_e_md_loc,eee_e_lg_loc

C finally, having added up the energies for all atoms we're going to 
C subtract the short-range and medium-range energies out here
C so that they're not added on multiple times

      eee_e_lg_loc=eee_e_lg_loc-eee_e_st_loc-eee_e_md_loc

      return
      end

