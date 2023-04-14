
      subroutine add_up_forces(
     &             miv_num_loc,miv_beg_loc,miv_end_loc,
     &             fx_bnd,fx_sht,fx_med,fx_lng,fx_umb,fx_tot,mirep)

      use allocatable_arrays         
      implicit real(a-h,o-z)        

      integer   :: miv_num_loc
      integer   :: miv_beg_loc(miv_num_loc)
      integer   :: miv_end_loc(miv_num_loc)
      real      :: fx_bnd(1:3,1:num_tot_atms)
      real      :: fx_sht(1:3,1:num_tot_atms)
      real      :: fx_med(1:3,1:num_tot_atms)
      real      :: fx_lng(1:3,1:num_tot_atms)
      real      :: fx_umb(1:3,1:num_tot_atms)
      real      :: fx_tot(1:3,1:num_tot_atms)
                                      
      if(treecode_elec) then

        do nnn=1,miv_num_loc

          fx_tot(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))=
     &    fx_bnd(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))+
     &    fx_sht(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))+
     &    fx_med(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))+
     &    fx_lng(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))

        enddo

      else

        do nnn=1,miv_num_loc

          fx_tot(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))=
     &    fx_bnd(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))+
     &    fx_sht(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))+
     &    fx_med(1:3,miv_beg_loc(nnn):miv_end_loc(nnn))

        enddo
      endif

      return
      end

