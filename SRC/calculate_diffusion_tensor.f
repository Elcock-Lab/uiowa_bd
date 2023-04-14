
      subroutine calculate_diffusion_tensor(mithrd,mirep,
     &           miv_num_loc,miv_typ_loc,miv_beg_loc,
     &           miv_end_loc,miv_dyn_loc,miv_hyd_loc)

      use allocatable_arrays         
      implicit real(a-h,o-z)        

      integer                :: miv_num_loc
      integer                :: miv_typ_loc(miv_num_loc)
      integer                :: miv_beg_loc(miv_num_loc)
      integer                :: miv_end_loc(miv_num_loc)
      integer                :: miv_dyn_loc(miv_num_loc)
      integer                :: miv_hyd_loc(miv_num_loc)

C a few different options for calculating diffusion tensors here - for
C now, the only one that's properly checked is the one that is done 
C without any ewald corrections for periodic boundaries. it shouldn't be
C too difficult to reinstate all the other methods in the near future

      if(.not.rpy_hydrodynamics.and..not.oarpy_hydrodynamics.and.
     &   .not.nodiv_hydrodynamics) return

      if(rpy_hydrodynamics) then

        do mmm=1,miv_num_loc

          if(miv_hyd_loc(mmm).eq.0) cycle

C use kkk to store the identity of the current HI system

          kkk=miv_hyd_loc(mmm)

C iii & jjj are real atom numbers - used to provide coordinates
C ikk & jkk are hydro atom numbers - used for diffusion tensor

          do iii=miv_beg_loc(mmm),miv_end_loc(mmm)
            ikk=atm_to_hyd_atm(iii)
            ik3=ikk*3-2 ! note use of hyd atm # here

C first loop from j=1,i-1

            do jkk=1,ikk-1
              jk3=jkk*3-2 ! note use of hyd atm # here
              jjj=hyd_atm_to_atm(jkk,kkk)
              w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
              w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
              w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 

C 2022 don't apply minimum distance convention with HIs
C
CCC           if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
CCC             w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
CCC             w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
CCC             w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
CCC           endif                        ! PBC_NO PBC_YES

              dist2=w1**2+w2**2+w3**2
              rrijsq=1.0/dist2
              rij=sqrt(dist2)
              crap0=(r_t_a(iii)+r_t_a(jjj))/2.0
              if(rij.ge.(r_t_a(iii)+r_t_a(jjj))) then
                crap1=1.0+(2.0*crap0*crap0*rrijsq/3.0)
                crap2=1.0-(2.0*crap0*crap0*rrijsq)
              else
                crap1=rij*((8.0/3.0)-(0.75*rij/
     &         (crap0)))/(2.0*crap0)
                crap2=dist2/(8.0*crap0*crap0)
              endif
              crap3=kT_over_8pi_eta/rij
              d_a_r(jk3  ,ik3  ,kkk,mirep)=crap3*
     &                                    (crap1+crap2*w1*w1*rrijsq)
              d_a_r(jk3+1,ik3  ,kkk,mirep)=crap3*(crap2*w2*w1*rrijsq)
              d_a_r(jk3+2,ik3  ,kkk,mirep)=crap3*(crap2*w3*w1*rrijsq)
              d_a_r(jk3  ,ik3+1,kkk,mirep)=crap3*(crap2*w1*w2*rrijsq)
              d_a_r(jk3+1,ik3+1,kkk,mirep)=crap3*
     &                                    (crap1+crap2*w2*w2*rrijsq)
              d_a_r(jk3+2,ik3+1,kkk,mirep)=crap3*(crap2*w3*w2*rrijsq)
              d_a_r(jk3  ,ik3+2,kkk,mirep)=crap3*(crap2*w1*w3*rrijsq)
              d_a_r(jk3+1,ik3+2,kkk,mirep)=crap3*(crap2*w2*w3*rrijsq)
              d_a_r(jk3+2,ik3+2,kkk,mirep)=crap3*
     &                                    (crap1+crap2*w3*w3*rrijsq)
            enddo

C now loop (in effect) from j=i,i  

            d_a_r(ik3  ,ik3  ,kkk,mirep)=d_t_a(iii)                  
            d_a_r(ik3+1,ik3+1,kkk,mirep)=d_t_a(iii)               
            d_a_r(ik3+2,ik3+2,kkk,mirep)=d_t_a(iii)               

C now loop from j=i+1,n  

            do jkk=ikk+1,num_hyd_atms(kkk)
              jjj=hyd_atm_to_atm(jkk,kkk)
              jk3=jkk*3-2 ! note use of hyd atm # here
              w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
              w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
              w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 

C 2022 don't apply minimum distance convention with HIs
C
CCC           if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
CCC             w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
CCC             w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
CCC             w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
CCC           endif                        ! PBC_NO PBC_YES

              dist2=w1**2+w2**2+w3**2
              rrijsq=1.0/dist2
              rij=sqrt(dist2)
              crap0=(r_t_a(iii)+r_t_a(jjj))/2.0
              if(rij.ge.(r_t_a(iii)+r_t_a(jjj))) then
                crap1=1.0+(2.0*crap0*crap0*rrijsq/3.0)
                crap2=1.0-(2.0*crap0*crap0*rrijsq)
              else
                crap1=rij*((8.0/3.0)-(0.75*rij/
     &         (crap0)))/(2.0*crap0)
                crap2=dist2/(8.0*crap0*crap0)
              endif
              crap3=kT_over_8pi_eta/rij
              d_a_r(jk3  ,ik3  ,kkk,mirep)=crap3*
     &                                    (crap1+crap2*w1*w1*rrijsq)
              d_a_r(jk3+1,ik3  ,kkk,mirep)=crap3*(crap2*w2*w1*rrijsq)
              d_a_r(jk3+2,ik3  ,kkk,mirep)=crap3*(crap2*w3*w1*rrijsq)
              d_a_r(jk3  ,ik3+1,kkk,mirep)=crap3*(crap2*w1*w2*rrijsq)
              d_a_r(jk3+1,ik3+1,kkk,mirep)=crap3*
     &                                    (crap1+crap2*w2*w2*rrijsq)
              d_a_r(jk3+2,ik3+1,kkk,mirep)=crap3*(crap2*w3*w2*rrijsq)
              d_a_r(jk3  ,ik3+2,kkk,mirep)=crap3*(crap2*w1*w3*rrijsq)
              d_a_r(jk3+1,ik3+2,kkk,mirep)=crap3*(crap2*w2*w3*rrijsq)
              d_a_r(jk3+2,ik3+2,kkk,mirep)=crap3*
     &                                    (crap1+crap2*w3*w3*rrijsq)
            enddo

C end the do loop over the iii atoms

          enddo

C end the do loop over miv_num_loc

        enddo

      elseif(oarpy_hydrodynamics.or.nodiv_hydrodynamics) then

        do mmm=1,miv_num_loc

          if(miv_hyd_loc(mmm).eq.0) cycle

C use kkk to store the identity of the current HI system

          kkk=miv_hyd_loc(mmm)

          if(miv_typ_loc(mmm).eq.0) then 
            write(*,*)'calculate_diffusion_tensor not yet'
            write(*,*)'working for rigids'
            write(*,*)'quitting :('
            stop
          endif

C iii & jjj are real atom numbers - used to provide coordinates
C ikk & jkk are hydro atom numbers - used for diffusion tensor

          do iii=miv_beg_loc(mmm),miv_end_loc(mmm)
            ikk=atm_to_hyd_atm(iii)
            ik3=ikk*3-2 ! note use of hyd atm # here

C first loop from j=1,i-1

C PREAVERAGE - zero out the div_x etc terms here - then accumulate

            div_x(iii)=0.0
            div_y(iii)=0.0
            div_z(iii)=0.0
            do jkk=1,ikk-1
              jk3=jkk*3-2 ! note use of hyd atm # here
              jjj=hyd_atm_to_atm(jkk,kkk)
              w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
              w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
              w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 

C 2022 don't apply minimum distance convention with HIs
C
CCC           if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
CCC             w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
CCC             w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
CCC             w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
CCC           endif                        ! PBC_NO PBC_YES

              dist2=w1**2+w2**2+w3**2
              rrijsq=1.0/dist2
              rij=sqrt(dist2)
              crap0=(r_t_a(iii)+r_t_a(jjj))/2.0

C PREAVERAGE - much simpler expressions to use in this case

              rij3=rij*rij*rij
              if(rij.ge.(r_t_a(iii)+r_t_a(jjj))) then
                d_a_r(jk3  ,ik3  ,kkk,mirep)=kT_over_6pi_eta/rij
                d_a_r(jk3+1,ik3+1,kkk,mirep)=kT_over_6pi_eta/rij
                d_a_r(jk3+2,ik3+2,kkk,mirep)=kT_over_6pi_eta/rij
                div_x(iii)=div_x(iii)+kT_over_6pi_eta*w1/rij3
                div_y(iii)=div_y(iii)+kT_over_6pi_eta*w2/rij3
                div_z(iii)=div_z(iii)+kT_over_6pi_eta*w3/rij3
              else
                d_a_r(jk3  ,ik3  ,kkk,mirep)=kT_over_6pi_eta*   
     &               (1.0-0.25*rij/crap0)/crap0
                d_a_r(jk3+1,ik3+1,kkk,mirep)=kT_over_6pi_eta*   
     &               (1.0-0.25*rij/crap0)/crap0
                d_a_r(jk3+2,ik3+2,kkk,mirep)=kT_over_6pi_eta*   
     &               (1.0-0.25*rij/crap0)/crap0
                div_x(iii)=div_x(iii)+
     &                kT_over_6pi_eta*w1/rij/crap0/crap0/4.0
                div_y(iii)=div_y(iii)+
     &                kT_over_6pi_eta*w2/rij/crap0/crap0/4.0
                div_z(iii)=div_z(iii)+
     &                kT_over_6pi_eta*w3/rij/crap0/crap0/4.0
              endif
            enddo

C now loop (in effect) from j=i,i  

            d_a_r(ik3  ,ik3  ,kkk,mirep)=d_t_a(iii)                  
            d_a_r(ik3+1,ik3+1,kkk,mirep)=d_t_a(iii)               
            d_a_r(ik3+2,ik3+2,kkk,mirep)=d_t_a(iii)               

C now loop from j=i+1,n  

            do jkk=ikk+1,num_hyd_atms(kkk)
              jjj=hyd_atm_to_atm(jkk,kkk)
              jk3=jkk*3-2 ! note use of hyd atm # here
              w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
              w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
              w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 

C 2022 don't apply minimum distance convention with HIs
C
CCC           if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
CCC             w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
CCC             w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
CCC             w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
CCC           endif                        ! PBC_NO PBC_YES

              dist2=w1**2+w2**2+w3**2
              rrijsq=1.0/dist2
              rij=sqrt(dist2)
              crap0=(r_t_a(iii)+r_t_a(jjj))/2.0

C PREAVERAGE - much simpler expressions to use in this case

              rij3=rij*rij*rij
              if(rij.ge.(r_t_a(iii)+r_t_a(jjj))) then
                d_a_r(jk3  ,ik3  ,kkk,mirep)=kT_over_6pi_eta/rij
                d_a_r(jk3+1,ik3+1,kkk,mirep)=kT_over_6pi_eta/rij
                d_a_r(jk3+2,ik3+2,kkk,mirep)=kT_over_6pi_eta/rij
                div_x(iii)=div_x(iii)+kT_over_6pi_eta*w1/rij3
                div_y(iii)=div_y(iii)+kT_over_6pi_eta*w2/rij3
                div_z(iii)=div_z(iii)+kT_over_6pi_eta*w3/rij3
              else
                d_a_r(jk3  ,ik3  ,kkk,mirep)=kT_over_6pi_eta*   
     &               (1.0-0.25*rij/crap0)/crap0
                d_a_r(jk3+1,ik3+1,kkk,mirep)=kT_over_6pi_eta*   
     &               (1.0-0.25*rij/crap0)/crap0
                d_a_r(jk3+2,ik3+2,kkk,mirep)=kT_over_6pi_eta*   
     &               (1.0-0.25*rij/crap0)/crap0
                div_x(iii)=div_x(iii)+
     &                kT_over_6pi_eta*w1/rij/crap0/crap0/4.0
                div_y(iii)=div_y(iii)+
     &                kT_over_6pi_eta*w2/rij/crap0/crap0/4.0
                div_z(iii)=div_z(iii)+
     &                kT_over_6pi_eta*w3/rij/crap0/crap0/4.0
              endif
            enddo

C end the do loop over the iii atoms

          enddo

C end the do loop over miv_num_loc

        enddo

      endif

      return
      end
