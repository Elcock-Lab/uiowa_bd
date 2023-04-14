
C noe that there are separate blocks of code for each treatment of HI

      subroutine move_atoms(num_cnt,mirep,
     &                      miv_num_loc,miv_typ_loc,miv_beg_loc,
     &                      miv_end_loc,miv_dyn_loc,miv_hyd_loc)
                                      
      use allocatable_arrays         
      implicit real(a-h,o-z)        

      integer,              intent(in) :: num_cnt ! num_hyd_cnt
      integer,              intent(in) :: mirep
      integer,              intent(in) :: miv_num_loc
      integer,              intent(in) :: miv_typ_loc(miv_num_loc)
      integer,              intent(in) :: miv_beg_loc(miv_num_loc)
      integer,              intent(in) :: miv_end_loc(miv_num_loc)
      integer,              intent(in) :: miv_dyn_loc(miv_num_loc)
      integer,              intent(in) :: miv_hyd_loc(miv_num_loc)

C simple option for steepest descent here - assume we're doing BD
C without HI and this should work fine - note that the time_step sets
C the step size in the move...

C ----------------------------------------------------------------------
      if(steepest_descent) then
C ----------------------------------------------------------------------

        do nnn=1,miv_num_loc
          if(miv_typ_loc(nnn).eq.1) then
            if(miv_dyn_loc(nnn).eq.1) then
              if(miv_hyd_loc(nnn).eq.0) then
                do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                  ii3=iii*3-2
                  c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+
     &                         f_a_x(1,iii,mirep)*time_step
                  c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+
     &                         f_a_x(2,iii,mirep)*time_step
                  c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+
     &                         f_a_x(3,iii,mirep)*time_step
                enddo
              endif
            endif
          endif
        enddo

C ----------------------------------------------------------------------
      elseif(rpy_hydrodynamics.or.no_hydrodynamics) then
C ----------------------------------------------------------------------

C we assume here that miv_typ_loc(nnn)=1, i.e. flexible atoms...

        do nnn=1,miv_num_loc

C BD moves first

          if(miv_dyn_loc(nnn).eq.1) then

C flex + BD + noHI
C r_u_r was scaled by root.2Dt.Dtrans so no further scaling needed

            if(miv_hyd_loc(nnn).eq.0) then

              if(i_do_lincs) then
                do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                  c_a_x_old(1,iii,mirep)=c_a_x(1,iii,mirep)
                  c_a_x_old(2,iii,mirep)=c_a_x(2,iii,mirep)
                  c_a_x_old(3,iii,mirep)=c_a_x(3,iii,mirep)
                enddo
              endif
              do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                ii3=iii*3-2
                c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+
     &                       f_a_x(1,iii,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3,num_cnt,mirep)
                c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+
     &                       f_a_x(2,iii,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3+1,num_cnt,mirep)
                c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+
     &                       f_a_x(3,iii,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3+2,num_cnt,mirep)
              enddo

C flex + BD + HI
C h_c_r comes from Sij*r_u_r so no need for more scaling
C h_a_r     contains those elements of f_a_r that are for hydro atoms

            elseif(miv_hyd_loc(nnn).ge.1) then

C use kkk to store the identity of the current HI system

              kkk=miv_hyd_loc(nnn)

              do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                ii3=iii*3-2 
                ikk=atm_to_hyd_atm(iii)
                ik3=ikk*3-2
                si1_x=d_a_r(1,ik3  ,kkk,mirep)*h_a_r(1,kkk,mirep)
                si1_y=d_a_r(1,ik3+1,kkk,mirep)*h_a_r(1,kkk,mirep)
                si1_z=d_a_r(1,ik3+2,kkk,mirep)*h_a_r(1,kkk,mirep)
                do jk3=2,num_hyd_atms3(kkk)
                  si1_x=si1_x+d_a_r(jk3,ik3,kkk,mirep)*
     &                        h_a_r(jk3,kkk,mirep)
                  si1_y=si1_y+d_a_r(jk3,ik3+1,kkk,mirep)*
     &                        h_a_r(jk3,kkk,mirep)
                  si1_z=si1_z+d_a_r(jk3,ik3+2,kkk,mirep)*
     &                        h_a_r(jk3,kkk,mirep)
                enddo
                c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+s_hydro*si1_x+
     &                             h_c_r(ik3  ,num_cnt,kkk,mirep)
                c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+s_hydro*si1_y+
     &                             h_c_r(ik3+1,num_cnt,kkk,mirep)
                c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+s_hydro*si1_z+
     &                             h_c_r(ik3+2,num_cnt,kkk,mirep)

              enddo

            endif

C LD moves second

          elseif(miv_dyn_loc(nnn).eq.2) then

C flex + LD + noHI
C r_u_r needs to be checked for scaling - read Winter/Geyer

            if(miv_hyd_loc(nnn).eq.0) then
              if(i_do_lincs) then
                do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                  c_a_x_old(1,iii,mirep)=c_a_x(1,iii,mirep)
                  c_a_x_old(2,iii,mirep)=c_a_x(2,iii,mirep)
                  c_a_x_old(3,iii,mirep)=c_a_x(3,iii,mirep)
                enddo
              endif
              do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                ii3=iii*3-2
                c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+
     &                       f_e_r(ii3  ,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3,num_cnt,mirep)
                c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+
     &                       f_e_r(ii3+1,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3+1,num_cnt,mirep)
                c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+
     &                       f_e_r(ii3+2,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3+2,num_cnt,mirep)
              enddo

C flex + LD + HI
C r_c_r needs to be checked for scaling - read Winter/Geyer
C h_e_r     contains those elements of f_e_r that are for hydro atoms

            elseif(miv_hyd_loc(nnn).ge.1) then

C use kkk to store the identity of the current HI system

              kkk=miv_hyd_loc(nnn)

              do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                ii3=iii*3-2 
                ikk=atm_to_hyd_atm(iii)
                ik3=ikk*3-2
                si1_x=d_a_r(1,ik3  ,kkk,mirep)*h_e_r(1,kkk,mirep)
                si1_y=d_a_r(1,ik3+1,kkk,mirep)*h_e_r(1,kkk,mirep)
                si1_z=d_a_r(1,ik3+2,kkk,mirep)*h_e_r(1,kkk,mirep)
                do jk3=2,num_hyd_atms3(kkk)
                  si1_x=si1_x+d_a_r(jk3,ik3  ,kkk,mirep)*
     &                        h_e_r(jk3,kkk,mirep)
                  si1_y=si1_y+d_a_r(jk3,ik3+1,kkk,mirep)* 
     &                        h_e_r(jk3,kkk,mirep)
                  si1_z=si1_z+d_a_r(jk3,ik3+2,kkk,mirep)*
     &                        h_e_r(jk3,kkk,mirep)
                enddo
                c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+s_hydro*si1_x+
     &                       h_c_r(ik3  ,num_cnt,kkk,mirep)
                c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+s_hydro*si1_y+
     &                       h_c_r(ik3+1,num_cnt,kkk,mirep)
                c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+s_hydro*si1_z+
     &                       h_c_r(ik3+2,num_cnt,kkk,mirep)
              enddo

            endif ! close the if HI/noHI

          endif   ! close the if BD/LD

        enddo

C ----------------------------------------------------------------------
C here's where we do orientally averaged RPY instead
C ----------------------------------------------------------------------

C note that most of this is a rehash of the above block of code

      elseif(oarpy_hydrodynamics.or.nodiv_hydrodynamics) then

C we assume here that miv_typ_loc(nnn)=1, i.e. flexible atoms...

        do nnn=1,miv_num_loc

C BD moves first

          if(miv_dyn_loc(nnn).eq.1) then

C flex + BD + noHI
C r_u_r was scaled by root.2Dt.Dtrans so no further scaling needed

            if(miv_hyd_loc(nnn).eq.0) then

              if(i_do_lincs) then
                do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                  c_a_x_old(1,iii,mirep)=c_a_x(1,iii,mirep)
                  c_a_x_old(2,iii,mirep)=c_a_x(2,iii,mirep)
                  c_a_x_old(3,iii,mirep)=c_a_x(3,iii,mirep)
                enddo
              endif
              do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                ii3=iii*3-2
                c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+
     &                       f_a_x(1,iii,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3,num_cnt,mirep)
                c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+
     &                       f_a_x(2,iii,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3+1,num_cnt,mirep)
                c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+
     &                       f_a_x(3,iii,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3+2,num_cnt,mirep)
              enddo

C flex + BD + HI
C h_c_r comes from Sij*r_u_r so no need for more scaling
C h_a_r     contains those elements of f_a_r that are for hydro atoms

            elseif(miv_hyd_loc(nnn).ge.1) then

C use kkk to store the identity of the current HI system

              kkk=miv_hyd_loc(nnn)

              do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                ii3=iii*3-2 
                ikk=atm_to_hyd_atm(iii)
                ik3=ikk*3-2
                si1_x=d_a_r(1,ik3  ,kkk,mirep)*h_a_r(1,kkk,mirep)
                si1_y=d_a_r(1,ik3+1,kkk,mirep)*h_a_r(1,kkk,mirep)
                si1_z=d_a_r(1,ik3+2,kkk,mirep)*h_a_r(1,kkk,mirep)
                do jk3=2,num_hyd_atms3(kkk)
                  si1_x=si1_x+d_a_r(jk3,ik3,kkk,mirep)*
     &                        h_a_r(jk3,kkk,mirep)
                  si1_y=si1_y+d_a_r(jk3,ik3+1,kkk,mirep)*
     &                        h_a_r(jk3,kkk,mirep)
                  si1_z=si1_z+d_a_r(jk3,ik3+2,kkk,mirep)*
     &                        h_a_r(jk3,kkk,mirep)
                enddo
                c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+s_hydro*si1_x+
     &                             h_c_r(ik3  ,num_cnt,kkk,mirep)
                c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+s_hydro*si1_y+
     &                             h_c_r(ik3+1,num_cnt,kkk,mirep)
                c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+s_hydro*si1_z+
     &                             h_c_r(ik3+2,num_cnt,kkk,mirep)

C add on the div_x terms (need to multiply by timestep)
C note that we have the option of *not* adding this on using
C nodiv_hydrodynamics - this is only to demo the divergence term

                if(oarpy_hydrodynamics) then
                  c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+
     &                               div_x(iii)*time_step
                  c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+
     &                               div_y(iii)*time_step
                  c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+
     &                               div_z(iii)*time_step
                endif

              enddo

            endif

C LD moves second

          elseif(miv_dyn_loc(nnn).eq.2) then

C flex + LD + noHI
C r_u_r needs to be checked for scaling - read Winter/Geyer

            if(miv_hyd_loc(nnn).eq.0) then
              if(i_do_lincs) then
                do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                  c_a_x_old(1,iii,mirep)=c_a_x(1,iii,mirep)
                  c_a_x_old(2,iii,mirep)=c_a_x(2,iii,mirep)
                  c_a_x_old(3,iii,mirep)=c_a_x(3,iii,mirep)
                enddo
              endif
              do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                ii3=iii*3-2
                c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+
     &                       f_e_r(ii3  ,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3,num_cnt,mirep)
                c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+
     &                       f_e_r(ii3+1,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3+1,num_cnt,mirep)
                c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+
     &                       f_e_r(ii3+2,mirep)*d_t_a(iii)*s_hydro+
     &                       r_u_r(ii3+2,num_cnt,mirep)
              enddo

C flex + LD + HI
C r_c_r needs to be checked for scaling - read Winter/Geyer
C h_e_r     contains those elements of f_e_r that are for hydro atoms

            elseif(miv_hyd_loc(nnn).ge.1) then

C use kkk to store the identity of the current HI system

              kkk=miv_hyd_loc(nnn)

              do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
                ii3=iii*3-2 
                ikk=atm_to_hyd_atm(iii)
                ik3=ikk*3-2
                si1_x=d_a_r(1,ik3  ,kkk,mirep)*h_e_r(1,kkk,mirep)
                si1_y=d_a_r(1,ik3+1,kkk,mirep)*h_e_r(1,kkk,mirep)
                si1_z=d_a_r(1,ik3+2,kkk,mirep)*h_e_r(1,kkk,mirep)
                do jk3=2,num_hyd_atms3(kkk)
                  si1_x=si1_x+d_a_r(jk3,ik3  ,kkk,mirep)*
     &                        h_e_r(jk3,kkk,mirep)
                  si1_y=si1_y+d_a_r(jk3,ik3+1,kkk,mirep)* 
     &                        h_e_r(jk3,kkk,mirep)
                  si1_z=si1_z+d_a_r(jk3,ik3+2,kkk,mirep)*
     &                        h_e_r(jk3,kkk,mirep)
                enddo
                c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+s_hydro*si1_x+
     &                       h_c_r(ik3  ,num_cnt,kkk,mirep)
                c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+s_hydro*si1_y+
     &                       h_c_r(ik3+1,num_cnt,kkk,mirep)
                c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+s_hydro*si1_z+
     &                       h_c_r(ik3+2,num_cnt,kkk,mirep)

C add on the div_x terms (need to multiply by timestep)
C note that we have the option of *not* adding this on using
C nodiv_hydrodynamics - this is only to demo the divergence term

                if(oarpy_hydrodynamics) then
                  c_a_x(1,iii,mirep)=c_a_x(1,iii,mirep)+
     &                               div_x(iii)*time_step
                  c_a_x(2,iii,mirep)=c_a_x(2,iii,mirep)+
     &                               div_y(iii)*time_step
                  c_a_x(3,iii,mirep)=c_a_x(3,iii,mirep)+
     &                               div_z(iii)*time_step
                endif

              enddo

            endif ! close the if HI/noHI

          endif   ! close the if BD/LD

        enddo

C close the if statement regarding steepest descent or not

      endif

      return
      end

