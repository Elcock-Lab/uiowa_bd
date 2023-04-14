
      subroutine shrink_box

      use allocatable_arrays
      implicit real(a-h,o-z)

C 20-2-15 note that i just commented out a whole lot of stuff that has
C to do with resetting limits on grid sizes later in this code - i don't
C think that we need it...

      if(.not.sphere.and..not.cylinder.and.nstep.le.n_size) then
        xmax=xmax_final+
     &      (real(n_size-nstep)/real(n_size))*
     &      (r_size_fac-1.0)*xmax_final
        ymax=ymax_final+
     &      (real(n_size-nstep)/real(n_size))*
     &      (r_size_fac-1.0)*ymax_final
        zmax=zmax_final+
     &      (real(n_size-nstep)/real(n_size))*
     &      (r_size_fac-1.0)*zmax_final
        xmin=xmin_final+
     &      (real(n_size-nstep)/real(n_size))*
     &      (r_size_fac-1.0)*xmin_final
        ymin=ymin_final+
     &      (real(n_size-nstep)/real(n_size))*
     &      (r_size_fac-1.0)*ymin_final
        zmin=zmin_final+
     &      (real(n_size-nstep)/real(n_size))*
     &      (r_size_fac-1.0)*zmin_final
      else
        xmax=xmax_final
        ymax=ymax_final
        zmax=zmax_final
        xmin=xmin_final
        ymin=ymin_final
        zmin=zmin_final
      endif

      if(sphere) then
        if(nstep.le.n_size) then
          s_size=r_size+
     &          (real(n_size-nstep)/real(n_size))*
     &          (r_size_fac-1.0)*r_size
        else
          s_size=r_size
        endif
        s_size2=s_size**2
      endif

      if(cylinder) then
        if(nstep.le.n_size) then
          s_size=r_size+
     &          (real(n_size-nstep)/real(n_size))*
     &          (r_size_fac-1.0)*r_size
          m_size=l_size+
     &          (real(n_size-nstep)/real(n_size))*
     &          (r_size_fac-1.0)*l_size
        else
          s_size=r_size
          m_size=l_size
        endif
      endif

      if(afm) then
        afm_tip_cur=afm_tip_cur+1
        if(afm_tip_dir.eq.1) then
          afm_tip_x_crd=afm_tip_x_beg+afm_tip_x_dsp*real(afm_tip_cur)
          afm_tip_y_crd=afm_tip_y_beg+afm_tip_y_dsp*real(afm_tip_cur)
          afm_tip_z_crd=afm_tip_z_beg+afm_tip_z_dsp*real(afm_tip_cur)
        elseif(afm_tip_dir.eq.-1) then
          afm_tip_x_crd=afm_tip_x_end-afm_tip_x_dsp*real(afm_tip_cur)
          afm_tip_y_crd=afm_tip_y_end-afm_tip_y_dsp*real(afm_tip_cur)
          afm_tip_z_crd=afm_tip_z_end-afm_tip_z_dsp*real(afm_tip_cur)
        endif
        if(afm_tip_cur.eq.afm_tip_stp) then
          afm_tip_dir=afm_tip_dir*-1
          afm_tip_cur=0
        endif
      endif

      return
      end

