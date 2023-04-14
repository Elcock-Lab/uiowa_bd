
      subroutine calculate_elec_forces_st_newton(n_ll_e_st,
     &             ida_ll_e_st,jda_ll_e_st,qqq_ll_e_st,
     &             ene_ll_e_st,fff_ll_e_st,ggg_ll_e_st,mirep)

      use allocatable_arrays
      implicit real(a-h,o-z)

      integer   :: ida_ll_e_st(1:n_ll_e_st)
      integer   :: jda_ll_e_st(1:n_ll_e_st)
      real      :: qqq_ll_e_st(1:n_ll_e_st)
      real      :: fff_ll_e_st(1:3,1:num_tot_atms) ! f_a_x_st
      real      :: ggg_ll_e_st(1:3,1:num_tot_atms) ! f_e_x_st

      if(cutoff_elec.and.treecode_elec) then

        if(implicit_ions) then
          do mmm=1,n_ll_e_st
            iii=ida_ll_e_st(mmm)
            jjj=jda_ll_e_st(mmm)
            qqq=qqq_ll_e_st(mmm)
            w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
            w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
            w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 
            if(i_pbc.eq.1) then
              w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
              w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
              w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
            endif
            dist2=w1**2+w2**2+w3**2
            dist=sqrt(dist2)       
            dist_inv=1.0/dist  
            temp=qqq*exp(-kappa*(dist-r_ion))
            u_v=temp*dielectric_inv*dist_inv*
     &          one_over_one_plus_kappa_rion
            f_v=u_v*(kappa+dist_inv)*dist_inv
            ene_ll_e_st=ene_ll_e_st+u_v
            fff_ll_e_st(1,iii)=fff_ll_e_st(1,iii)+f_v*w1    
            fff_ll_e_st(2,iii)=fff_ll_e_st(2,iii)+f_v*w2  
            fff_ll_e_st(3,iii)=fff_ll_e_st(3,iii)+f_v*w3   
            fff_ll_e_st(1,jjj)=fff_ll_e_st(1,jjj)-f_v*w1
            fff_ll_e_st(2,jjj)=fff_ll_e_st(2,jjj)-f_v*w2
            fff_ll_e_st(3,jjj)=fff_ll_e_st(3,jjj)-f_v*w3
            ggg_ll_e_st(1,iii)=ggg_ll_e_st(1,iii)+f_v*w1    
            ggg_ll_e_st(2,iii)=ggg_ll_e_st(2,iii)+f_v*w2  
            ggg_ll_e_st(3,iii)=ggg_ll_e_st(3,iii)+f_v*w3   
            ggg_ll_e_st(1,jjj)=ggg_ll_e_st(1,jjj)-f_v*w1    
            ggg_ll_e_st(2,jjj)=ggg_ll_e_st(2,jjj)-f_v*w2  
            ggg_ll_e_st(3,jjj)=ggg_ll_e_st(3,jjj)-f_v*w3   
          enddo
        else
          do mmm=1,n_ll_e_st
            iii=ida_ll_e_st(mmm)
            jjj=jda_ll_e_st(mmm)
            qqq=qqq_ll_e_st(mmm)
            w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
            w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
            w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 
            if(i_pbc.eq.1) then
              w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
              w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
              w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
            endif
            dist2=w1**2+w2**2+w3**2
            dist2_inv=1.0/(w1*w1+w2*w2+w3*w3)
            dist_inv=sqrt(dist2_inv)       
            u_v=qqq*dielectric_inv*dist_inv
            f_v=u_v*dist2_inv
            ene_ll_e_st=ene_ll_e_st+u_v
            fff_ll_e_st(1,iii)=fff_ll_e_st(1,iii)+f_v*w1    
            fff_ll_e_st(2,iii)=fff_ll_e_st(2,iii)+f_v*w2  
            fff_ll_e_st(3,iii)=fff_ll_e_st(3,iii)+f_v*w3   
            fff_ll_e_st(1,jjj)=fff_ll_e_st(1,jjj)-f_v*w1
            fff_ll_e_st(2,jjj)=fff_ll_e_st(2,jjj)-f_v*w2
            fff_ll_e_st(3,jjj)=fff_ll_e_st(3,jjj)-f_v*w3
            ggg_ll_e_st(1,iii)=ggg_ll_e_st(1,iii)+f_v*w1    
            ggg_ll_e_st(2,iii)=ggg_ll_e_st(2,iii)+f_v*w2  
            ggg_ll_e_st(3,iii)=ggg_ll_e_st(3,iii)+f_v*w3   
            ggg_ll_e_st(1,jjj)=ggg_ll_e_st(1,jjj)-f_v*w1    
            ggg_ll_e_st(2,jjj)=ggg_ll_e_st(2,jjj)-f_v*w2  
            ggg_ll_e_st(3,jjj)=ggg_ll_e_st(3,jjj)-f_v*w3   
          enddo
        endif

      elseif(cutoff_elec.and..not.treecode_elec) then

        if(implicit_ions) then
          do mmm=1,n_ll_e_st
            iii=ida_ll_e_st(mmm)
            jjj=jda_ll_e_st(mmm)
            qqq=qqq_ll_e_st(mmm)
            w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
            w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
            w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 
            if(i_pbc.eq.1) then
              w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
              w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
              w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
            endif
            dist2=w1**2+w2**2+w3**2
            dist=sqrt(dist2)       
            dist_inv=1.0/dist  
            temp=qqq*exp(-kappa*(dist-r_ion))
            u_v=temp*dielectric_inv*dist_inv*
     &          one_over_one_plus_kappa_rion
            f_v=u_v*(kappa+dist_inv)*dist_inv
            ene_ll_e_st=ene_ll_e_st+u_v
            fff_ll_e_st(1,iii)=fff_ll_e_st(1,iii)+f_v*w1    
            fff_ll_e_st(2,iii)=fff_ll_e_st(2,iii)+f_v*w2  
            fff_ll_e_st(3,iii)=fff_ll_e_st(3,iii)+f_v*w3   
            fff_ll_e_st(1,jjj)=fff_ll_e_st(1,jjj)-f_v*w1
            fff_ll_e_st(2,jjj)=fff_ll_e_st(2,jjj)-f_v*w2
            fff_ll_e_st(3,jjj)=fff_ll_e_st(3,jjj)-f_v*w3
          enddo
        else
          do mmm=1,n_ll_e_st
            iii=ida_ll_e_st(mmm)
            jjj=jda_ll_e_st(mmm)
            qqq=qqq_ll_e_st(mmm)
            w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep) 
            w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep) 
            w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep) 
            if(i_pbc.eq.1) then
              w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
              w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
              w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
            endif
            dist2=w1**2+w2**2+w3**2
            dist2_inv=1.0/(w1*w1+w2*w2+w3*w3)
            dist_inv=sqrt(dist2_inv)       
            u_v=qqq*dielectric_inv*dist_inv
            f_v=u_v*dist2_inv
            ene_ll_e_st=ene_ll_e_st+u_v
            fff_ll_e_st(1,iii)=fff_ll_e_st(1,iii)+f_v*w1    
            fff_ll_e_st(2,iii)=fff_ll_e_st(2,iii)+f_v*w2  
            fff_ll_e_st(3,iii)=fff_ll_e_st(3,iii)+f_v*w3   
            fff_ll_e_st(1,jjj)=fff_ll_e_st(1,jjj)-f_v*w1
            fff_ll_e_st(2,jjj)=fff_ll_e_st(2,jjj)-f_v*w2
            fff_ll_e_st(3,jjj)=fff_ll_e_st(3,jjj)-f_v*w3
          enddo
        endif

      endif

      return
      end

