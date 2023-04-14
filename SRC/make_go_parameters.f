
      subroutine make_go_parameters      
                                       
      use allocatable_arrays          
      implicit real(a-h,o-z)         
                                    
      character*80 superstring     
      character*5 atemp           
      integer jnk(10)            
      num_tot_go_pairs=0

C now read the list of go pairs once to find the number of contacts for each pair
C note that we're reading the m typ #'s,not m #'s                                      
C                                                                                           
      open(32,file=goparam_file,form='formatted',status='unknown')  
      read(32,789)atemp                                          

C put in a check to make sure that the go parameters have the right functional form

      if((atemp.eq.'12-10'.and.i_use_v_typ.ne.1).or.
     &   (atemp.eq.'12-06'.and.i_use_v_typ.ne.2).or.
     &   (atemp.eq.'08-04'.and.i_use_v_typ.ne.3)) then
        write(*,*)'go parameters do not match the desired LJ potential'
        stop
      endif
789   format(a5)           
799   format(1x)           
21    read(32,802,end=299)itp,iat,jtp,jat,t0,e1 
802   format(4i10,2f15.5)                      

C note that we use go_eps_low here so that only those go contacts with
C an epsilon above a threshold value get included in the calculation of
C Q values

      num_tot_go_pairs=num_tot_go_pairs+1
      if(e1.ge.go_eps_low) then
        num_des_go_pairs(itp,jtp)=num_des_go_pairs(itp,jtp)+1
      endif

      goto 21     
299   rewind(32)      

C now think about allocating memory for the go contacts of each m typ pair

C first allocate space for num_curr_conts - in order to do this we 
C need to know the number of molecules of each t_typ

      allocate(itp_tot_go_pair(1:num_tot_go_pairs))
      allocate(ida_tot_go_pair(1:num_tot_go_pairs))
      allocate(jtp_tot_go_pair(1:num_tot_go_pairs))
      allocate(jda_tot_go_pair(1:num_tot_go_pairs))
      allocate(go_t11(1:num_tot_go_pairs))
      allocate(go_t22(1:num_tot_go_pairs))
      allocate(go_t33(1:num_tot_go_pairs))
      allocate(go_t44(1:num_tot_go_pairs))
      allocate(go_dis(1:num_tot_go_pairs))
      allocate(go_di2(1:num_tot_go_pairs))
      allocate(go_elo(1:num_tot_go_pairs))
      allocate(go_d12(1:num_tot_go_pairs))
      allocate(go_mod(1:num_tot_go_pairs))
      allocate(go_eps(1:num_tot_go_pairs))

C now read through again and get all the appropriate parameters
C note that we assume that the go_mod = 1 if it's not given...

      num_tot_go_pairs=0                        
      itmp=0
      read(32,799)          
41    continue
      read(32,804,end=399)superstring           
804   format(a80)
      read(superstring,802)itp,iat,jtp,jat,t0,e1
      id=1
803   format(4i10,2f15.5,i10)          

      num_tot_go_pairs=num_tot_go_pairs+1
      t1=60.000*e1*t0**12                             
      t2=60.000*e1*t0**10                             
      t3= 5.000*e1*t0**12                             
      t4= 6.000*e1*t0**10                             
      itp_tot_go_pair(num_tot_go_pairs)=itp
      ida_tot_go_pair(num_tot_go_pairs)=iat
      jtp_tot_go_pair(num_tot_go_pairs)=jtp
      jda_tot_go_pair(num_tot_go_pairs)=jat
      go_t11(num_tot_go_pairs)=t1
      go_t22(num_tot_go_pairs)=t2
      go_t33(num_tot_go_pairs)=t3
      go_t44(num_tot_go_pairs)=t4
      go_dis(num_tot_go_pairs)=t0
      go_di2(num_tot_go_pairs)=t0**2
      go_mod(num_tot_go_pairs)=id
      go_eps(num_tot_go_pairs)=e1
      if(id.gt.go_mod_max) go_mod_max=id

C allow for a softer short-range interaction

      if(i_use_v_typ.eq.1) then
        if(t1.lt.1.000) then
          go_elo(num_tot_go_pairs)=0.0     
        else
          go_elo(num_tot_go_pairs)=t3/(t0**12)-t4/(t0**10)
        endif
      endif
      go_d12(num_tot_go_pairs)= t0**2*(1.200**2)
      goto 41                         
399   close(32)  

      return                         
      end                           
