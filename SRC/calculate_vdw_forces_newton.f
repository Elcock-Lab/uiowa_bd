
C this code should work fine for both _st and _md forces

      subroutine calculate_vdw_forces_newton(n_ll_v_st,
     &             ida_ll_v_st,jda_ll_v_st,s11_ll_v_st,
     &             s22_ll_v_st,s33_ll_v_st,s44_ll_v_st,
     &             d11_ll_v_st,d22_ll_v_st,e00_ll_v_st,
     &             ene_ll_v_st,fff_ll_v_st,mirep)

      use allocatable_arrays
      implicit real(a-h,o-z)

      integer   :: ida_ll_v_st(1:n_ll_v_st)
      integer   :: jda_ll_v_st(1:n_ll_v_st)
      real      :: s11_ll_v_st(1:n_ll_v_st)
      real      :: s22_ll_v_st(1:n_ll_v_st)
      real      :: s33_ll_v_st(1:n_ll_v_st)
      real      :: s44_ll_v_st(1:n_ll_v_st)
      real      :: d11_ll_v_st(1:n_ll_v_st)
      real      :: d22_ll_v_st(1:n_ll_v_st)
      real      :: e00_ll_v_st(1:n_ll_v_st)
      real      :: fff_ll_v_st(1:3,1:num_tot_atms)

      do mmm=1,n_ll_v_st
        iii=ida_ll_v_st(mmm)
        jjj=jda_ll_v_st(mmm)
        s1a=s11_ll_v_st(mmm)
        s2a=s22_ll_v_st(mmm)
        s3a=s33_ll_v_st(mmm)
        s4a=s44_ll_v_st(mmm)
        d1a=d11_ll_v_st(mmm)
        d2a=d22_ll_v_st(mmm)
        e0a=e00_ll_v_st(mmm)
        w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
        w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
        w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 
        if(i_pbc.eq.1) then
          w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
          w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
          w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
        endif
        dist2=w1*w1+w2*w2+w3*w3
        dist2_inv=1.0/dist2
        if(dist2.gt.d2a.or.r_f_st.le.0.0) then
          f_v=s1a*(dist2_inv**7)-s2a*(dist2_inv**6)
          u_v=s3a*(dist2_inv**6)-s4a*(dist2_inv**5)
        else
          dist=sqrt(dist2)
          f_v=-r_f_st*(dist-d1a)/dist
          u_v= r_f_hf*(dist-d1a)**2+e0a
        endif
        ene_ll_v_st=ene_ll_v_st+u_v
        fff_ll_v_st(1,iii)=fff_ll_v_st(1,iii)+f_v*w1    
        fff_ll_v_st(2,iii)=fff_ll_v_st(2,iii)+f_v*w2  
        fff_ll_v_st(3,iii)=fff_ll_v_st(3,iii)+f_v*w3   
        fff_ll_v_st(1,jjj)=fff_ll_v_st(1,jjj)-f_v*w1    
        fff_ll_v_st(2,jjj)=fff_ll_v_st(2,jjj)-f_v*w2  
        fff_ll_v_st(3,jjj)=fff_ll_v_st(3,jjj)-f_v*w3   
      enddo

      return
      end

