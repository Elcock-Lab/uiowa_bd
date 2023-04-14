      
      subroutine calculate_uncorrelated_randoms(mythrd_loc,mirep,
     &                      miv_num_loc,miv_typ_loc,miv_beg_loc,
     &                      miv_end_loc,miv_dyn_loc,miv_hyd_loc,
     &                      num_rnd_loc,ran_trm_loc)
                                      
      use allocatable_arrays         
      USE MKL_VSL_TYPE
      USE MKL_VSL

      implicit real(a-h,o-z)        

      real ran_trm_loc(1:num_rnd_loc)

      integer,              intent(in) :: miv_num_loc
      integer,              intent(in) :: miv_typ_loc(miv_num_loc)
      integer,              intent(in) :: miv_beg_loc(miv_num_loc)
      integer,              intent(in) :: miv_end_loc(miv_num_loc)
      integer,              intent(in) :: miv_dyn_loc(miv_num_loc)
      integer,              intent(in) :: miv_hyd_loc(miv_num_loc)

C move terms into 2D r_u_r array with num_hyd_atms*3 rows and num_hyd_stp cols
C at the end both the LD and BD r_u_r terms are in displacement form
C note that noHI terms are scaled by root.Dtrans ; HI are not since they 
C get multiplied by Smatrix later (this provides the root.D term)

C BD terms should be fine. LD terms need checking!!

      do nnn=1,miv_num_loc
        if(miv_typ_loc(nnn).eq.1) then
          if(miv_dyn_loc(nnn).eq.1) then 

C flex + BD + noHI

            if(miv_hyd_loc(nnn).eq.0) then
              ran_trm_loc=ran_trm_loc*d_hydro
              j=0
              do n=1,num_hyd_stp
                do i=miv_beg_loc(nnn)*3-2,miv_end_loc(nnn)*3
                  j=j+1
                  r_u_r(i,n,mirep)=ran_trm_loc(j)*d_t_a_root(i)
                enddo
              enddo

C flex + BD + HI

            elseif(miv_hyd_loc(nnn).ge.1) then
              kkk=miv_hyd_loc(nnn)
              ran_trm_loc=ran_trm_loc*d_hydro
              j=0
              do n=1,num_hyd_stp
                do i=miv_beg_loc(nnn)*3-2,miv_end_loc(nnn)*3
                  j=j+1
                  r_u_r(i,n,mirep)=ran_trm_loc(j)
                enddo
              enddo              

              ikk=atm_to_hyd_atm(miv_beg_loc(nnn))
              ik3=ikk*3
              ik3=ik3-miv_beg_loc(nnn)*3
              do n=1,num_hyd_stp
                do ii3=miv_beg_loc(nnn)*3-2,miv_end_loc(nnn)*3
                  h_u_r(ii3+ik3,n,kkk,mirep)=r_u_r(ii3,n,mirep)
                enddo
              enddo
            endif

          elseif(miv_dyn_loc(nnn).eq.2) then

C flex + LD +noHI
C (1) convert random term to undamped uncorrelated random force
C (2) convert to damped uncorrelated random force 
C (3) convert back to a displacement
C note that in the last step we use ri = fi * D * dt / kT

            if(miv_hyd_loc(nnn).eq.0) then
              ran_trm_loc=ran_trm_loc*f_hydro
              j=0
              do n=1,num_hyd_stp
                do i=miv_beg_loc(nnn)*3-2,miv_end_loc(nnn)*3
                  j=j+1
                  r_u_r(i,n,mirep)=ran_trm_loc(j)*d_t_a_root_inv(i)
                  f_r_r(i,mirep)=ld_fac1(i)*r_u_r(i,n,mirep)+
     &                           ld_fac2(i)*v_r_r(i,mirep)
                  v_r_r(i,mirep)=ld_fac3(i)*v_r_r(i,mirep)+
     &                           ld_fac4(i)*r_u_r(i,n,mirep)
                  r_u_r(i,n,mirep)=f_r_r(i,mirep)*d_t_a_3frm(i)*df_hydro
                enddo
              enddo

C flex + LD + HI
C as with the BD, note that the noHI term is effectively scaled by an
C additional sqrt.Dtrans term that is missing here...

            elseif(miv_hyd_loc(nnn).ge.1) then
              kkk=miv_hyd_loc(nnn)
              ran_trm_loc=ran_trm_loc*f_hydro
              j=0
              do n=1,num_hyd_stp
                do i=miv_beg_loc(nnn)*3-2,miv_end_loc(nnn)*3
                  j=j+1
                  r_u_r(i,n,mirep)=ran_trm_loc(j)*d_t_a_root_inv(i)
                  f_r_r(i,mirep)=ld_fac1(i)*r_u_r(i,n,mirep)+
     &                           ld_fac2(i)*v_r_r(i,mirep)
                  v_r_r(i,mirep)=ld_fac3(i)*v_r_r(i,mirep)+
     &                           ld_fac4(i)*r_u_r(i,n,mirep)
                  r_u_r(i,n,mirep)=f_r_r(i,mirep)*d_t_a_root(i)*df_hydro
                enddo
              enddo            

              ikk=atm_to_hyd_atm(miv_beg_loc(nnn))
              ik3=ikk*3
              ik3=ik3-miv_beg_loc(nnn)*3
              do n=1,num_hyd_stp
                do ii3=miv_beg_loc(nnn)*3-2,miv_end_loc(nnn)*3
                  h_u_r(ii3+ik3,n,kkk,mirep)=r_u_r(ii3,n,mirep)
                enddo
              enddo

            endif
          endif
        endif
      enddo

      return
      end

