
      subroutine make_charge_parameters                                                                         
                                                                                                      
      use allocatable_arrays                                                                          
      implicit real(a-h,o-z)                                                                          
                                                                                                      
      real,       allocatable :: re1(:)                                                        
      real,       allocatable :: re2(:)                                                        
      real,       allocatable :: re3(:)                                                        
      real,       allocatable :: re4(:)                                                        
      real,       allocatable :: rs0(:)                                                        
      real,       allocatable :: rs1(:)                                                        
      real,       allocatable :: rs2(:)                                                        
      real,       allocatable :: rs3(:)                                                        
      real,       allocatable :: rs4(:)                                                        
      integer,    allocatable :: myg(:)                                                        
      integer,    allocatable :: id_tmp(:)
      integer,    allocatable :: jd_tmp(:)
      integer,    allocatable :: kd_tmp(:,:) 
      character*3,allocatable :: fattyp(:)                                                               
      character*100 string                                                                  
      character*4 a                                                                      
      character*10 ajk1                                                                   
      data a /'ATOM'/                                                                    
      dimension coords(3)                                                                   
      character*6 name1                                                                     
      character*4 name2                                                                     
      character*3 name3                                                                     
      character*5 name4                                                                     
      character*2 name5                                                                     
      character*3 attyp                                                                    
      character*1 atn_typ                                                                  
      integer     ierror(100)                                                               

C open parameter file and read # of atom typs for allocating arrays                    
                                                                                            
      ntmp=0                                                                                
      open(30,file=parameter_file,form='formatted',status='old')                                   
14    read(30,12,end=13)
12    format(1x)                                                                   
      ntmp=ntmp+1                                                                           
      goto 14                                                                               
13    continue                                                                              
      allocate(re1(1:ntmp))                                                                 
      allocate(re2(1:ntmp))                                                                 
      allocate(re3(1:ntmp))                                                                 
      allocate(re4(1:ntmp))                                                                 
      allocate(rs0(1:ntmp))                                                                 
      allocate(rs1(1:ntmp))                                                                 
      allocate(rs2(1:ntmp))                                                                 
      allocate(rs3(1:ntmp))                                                                 
      allocate(rs4(1:ntmp))                                                                 
      allocate(myg(1:ntmp))                                                                 
      allocate(fattyp(1:ntmp))                                                                 
      rewind(30)                                                                            
      
C set # of atom typs actually found in the system to zero initially                      
                                                                                            
      ntpfnd=0                                                                         

C do a loop over all molecules                                                              
                                                                                            
      do 666 iii=1,num_t_ms                                                               

C temporary write out of types for checking purposes...

        if(i_debug) write(*,667)iii,i_t_tp(iii),i_f_tp(iii),i_r_tp(iii)
667     format('mol # and typs : ',4i8)

C first ask if we've seen this typ of molecule before                                      

        do mmm=iii-1,1,-1
          if(i_t_tp(iii).eq.i_t_tp(mmm)) goto 666
        enddo
 
C zero out the arrays                                                                       

        num_bonds(i_t_tp(iii))=0
        num_angls(i_t_tp(iii))=0
        num_dihes(i_t_tp(iii))=0
                      
C open up the charge (pdb) file                                                             
                     
        open(32,file=qef_file(i_t_tp(iii)),form='formatted',
     &       status='old')          

        qtot=0.00   
        natt=0     
                                                                
C come here to read a new atom

 1      continue  

        read(32,500,end=99)string    
        if(string(1:4).ne.'ATOM') goto 1

        read(string,501)coords(1),coords(2),coords(3),
     &                  atn_typ,q,r

501     format(30x,3f8.3,1x,a,3f10.3)

        read(string,502)name1,name2,name3,name4                
502     format(5x,a6,1x,a4,1x,a3,2x,a4)                                       

C zero the charge if i_use_no_elec

        if(i_use_no_elec) q=0.0                                                           

C note for this to work as intended we're going to have 3 'Q' atoms in
C each rigid molecule - these will be treated like 'Q' atoms in flexible
C molecules - we'll want to give all 3 atoms zero charge and only put a
C hydrodynamic radius on the 3rd atom...

        if(atn_typ.eq.'Q') then
          natt=natt+1  
          c_typ(i_t_tp(iii))%x(natt)=coords(1)                                        
          c_typ(i_t_tp(iii))%y(natt)=coords(2)                                        
          c_typ(i_t_tp(iii))%z(natt)=coords(3)                                        
          c_typ(i_t_tp(iii))%q(natt)=q
          c_typ(i_t_tp(iii))%r(natt)=r
          n_typ(i_t_tp(iii))%ch1(natt)=name1                                          
          n_typ(i_t_tp(iii))%ch2(natt)=name2                                          
          n_typ(i_t_tp(iii))%ch3(natt)=name3                                          
          n_typ(i_t_tp(iii))%ch4(natt)=name4                                          
          num_atms_pt(i_t_tp(iii))=num_atms_pt(i_t_tp(iii))+1

C here is where we can read in extra 'D' atoms for rigid molecules -
C note that these will be added to the list of 'Q' atoms for c_typ and
C n_typ arrays, but we'll keep a separate counter for them 

        elseif(atn_typ.eq.'D') then
          natt=natt+1  
          c_typ(i_t_tp(iii))%x(natt)=coords(1)                                        
          c_typ(i_t_tp(iii))%y(natt)=coords(2)                                        
          c_typ(i_t_tp(iii))%z(natt)=coords(3)                                        
          n_typ(i_t_tp(iii))%ch1(natt)=name1                                          
          n_typ(i_t_tp(iii))%ch2(natt)=name2                                          
          n_typ(i_t_tp(iii))%ch3(natt)=name3                                          
          n_typ(i_t_tp(iii))%ch4(natt)=name4                                          
          num_dyns_pt(i_t_tp(iii))=num_dyns_pt(i_t_tp(iii))+1
        endif

        qtot=qtot+q                                                                     
                                                                                            
C now, if this is a 'Q' atom of a flexible molecule, find vdw parameters     
                  
        if(ityp(iii).eq.1.and.atn_typ.eq.'Q') then
                                                                                            
C open up the parameter file                                                                
                                                                                            
          igotfound=0                                                                       
          rewind(30)                                                                        
919       read(30,*,end=934)attyp,radius,r_epsilon,igreasy                                   

C                                                                                           
C note that the van der Waals term that we're applying here                                 
C is a 12-6 term (at the moment)                                                            
C                                                                                           
C so: energy = 2.0*sqrt(epsi)*2.0*sqrt(epsj)*                                               
C                            [sigi**6  *  sigj**6/r**12 -                                                  
C              sigi**3  *  sigj**3/r**6]                                                    
C                                                                                           
C this looks like a somewhat crazy functional form,but it                                  
C makes a lot of sense when we come to calculate forces                                     
C                                                                                           
C note that assignment of parameters is based on the ement in col 14                      
C note also that we keep a list of all new typs found (ntpfnd)                        
C                                                                                           
          if(attyp.eq.string(14:16)) then                                                  
            notnew=0                                                                        
            do ijk=1,ntpfnd                                                            
              if(attyp.eq.fattyp(ijk)) then                                               
                c_typ(i_t_tp(iii))%a(num_atms_pt(i_f_tp(iii)))=ijk
                notnew=1                                                                    
              endif                                                                         
            enddo                                                                           
            if(notnew.eq.0) then                                                            
              ntpfnd=ntpfnd+1                                                     
              fattyp(ntpfnd)=attyp                                                   
              c_typ(i_t_tp(iii))%a(num_atms_pt(i_f_tp(iii)))=ntpfnd
C                                                                                           
C i_use_v_typ=1,2,3 for 12-10,12-06,08-04 respectively
C                                                                                           
              if(igreasy.ge.1) then                                                    [12-6 potential  ]
                if(i_use_v_typ.eq.1) then                                                [12-6 potential  ]
                  rs0(ntpfnd)=radius                                                 [12-6 potential  ]
                  rs1(ntpfnd)=sqrt(12.000)*radius**6                                 [12-6 potential  ]
                  rs2(ntpfnd)=sqrt(10.000)*radius**5                                             [12-6 potential  ]
                  rs3(ntpfnd)=radius**6                                             [12-6 potential  ]
                  rs4(ntpfnd)=radius**5                                             [12-6 potential  ]
                  re1(ntpfnd)=sqrt(5.000*r_epsilon)                                  [12-6 potential  ]
                  re2(ntpfnd)=sqrt(6.000*r_epsilon)                                   [12-6 potential  ]
                  myg(ntpfnd)=igreasy                                               [12-6 potential  ]
                elseif(i_use_v_typ.eq.2) then                                                                     [12-6 potential  ]
                  rs0(ntpfnd)=radius                                                 [12-6 potential  ]
                  rs1(ntpfnd)=sqrt(12.000)*radius**6                                 [12-6 potential  ]
                  rs2(ntpfnd)=sqrt(6.000)*radius**3                                             [12-6 potential  ]
                  rs3(ntpfnd)=radius**6                                             [12-6 potential  ]
                  rs4(ntpfnd)=radius**3                                             [12-6 potential  ]
                  re1(ntpfnd)=sqrt(1.000*r_epsilon)                                  [12-6 potential  ]
                  re2(ntpfnd)=sqrt(2.000*r_epsilon)                                   [12-6 potential  ]
                  myg(ntpfnd)=igreasy                                               [12-6 potential  ]
                elseif(i_use_v_typ.eq.3) then                                                                     [12-6 potential  ]
                  rs0(ntpfnd)=radius                                                 [12-6 potential  ]
                  rs1(ntpfnd)=sqrt(8.000)*radius**4                                 [12-6 potential  ]
                  rs2(ntpfnd)=sqrt(4.000)*radius**2                                             [12-6 potential  ]
                  rs3(ntpfnd)=radius**4                                             [12-6 potential  ]
                  rs4(ntpfnd)=radius**2                                             [12-6 potential  ]
                  re1(ntpfnd)=sqrt(1.000*r_epsilon)                                  [12-6 potential  ]
                  re2(ntpfnd)=sqrt(2.000*r_epsilon)                                   [12-6 potential  ]
                  myg(ntpfnd)=igreasy                                               [12-6 potential  ]
                endif                                                                    [12-6 potential  ]
              else                                                                     [12-6 potential  ]
C                                                                                           
C 27-05-07 note that unlike the 12-6 potentials,we're assuming that                        
C energy = eps*sig^12/dis^12 for repulsive-only interactions in order                       
C to be consistent with our PLoS paper. In the case of 12-6 potentials                      
C we just use the first component of the 12-6 potential unchanged...                        
C                                                                                           
                if(i_use_v_typ.eq.1) then                                                [12-6 potential  ]
                  rs0(ntpfnd)=radius                                                 [12-6 potential  ]
                  rs1(ntpfnd)=sqrt(12.000)*radius**6                                 [12-6 potential  ]
                  rs2(ntpfnd)=0.000                                                 [12-6 potential  ]
                  rs3(ntpfnd)=radius**6                                             [12-6 potential  ]
                  rs4(ntpfnd)=0.000                                                 [12-6 potential  ]
                  re1(ntpfnd)=sqrt(1.000*r_epsilon)                                  [12-6 potential  ]
                  re2(ntpfnd)=sqrt(1.000*r_epsilon)                                   [12-6 potential  ]
                  myg(ntpfnd)=0                                                     [12-6 potential  ]
                elseif(i_use_v_typ.eq.2) then                                         [12-6 potential  ]
                  rs0(ntpfnd)=radius                                                 [12-6 potential  ]
                  rs1(ntpfnd)=sqrt(12.000)*radius**6                                 [12-6 potential  ]
                  rs2(ntpfnd)=0.000                                                 [12-6 potential  ]
                  rs3(ntpfnd)=radius**6                                             [12-6 potential  ]
                  rs4(ntpfnd)=0.000                                                 [12-6 potential  ]
                  re1(ntpfnd)=sqrt(r_epsilon)                                  [12-6 potential  ]
C                                                                                           
C note that although the following looks like it should be zero it                          
C actually isn't a problem because of statements that follow later...                       
C                                                                                           
                  re2(ntpfnd)=sqrt(2.000*r_epsilon)                                   [12-6 potential  ]
                  myg(ntpfnd)=0                                                     [12-6 potential  ]
                elseif(i_use_v_typ.eq.3) then                                         [12-6 potential  ]
                  rs0(ntpfnd)=radius                                                 [12-6 potential  ]
                  rs1(ntpfnd)=sqrt(8.000)*radius**4                                 [12-6 potential  ]
                  rs2(ntpfnd)=0.000                                                 [12-6 potential  ]
                  rs3(ntpfnd)=radius**4                                             [12-6 potential  ]
                  rs4(ntpfnd)=0.000                                                 [12-6 potential  ]
                  re1(ntpfnd)=sqrt(r_epsilon)                                  [12-6 potential  ]
                  re2(ntpfnd)=sqrt(2.000*r_epsilon)                                   [12-6 potential  ]
                  myg(ntpfnd)=0                                                     [12-6 potential  ]
                endif                                                                    [12-6 potential  ]
              endif                                                                    [12-6 potential  ]
            endif                                                                           
            goto 939                                                                        
          endif                                                                             
          goto 919                                                                          
C                                                                                           
C if we couldn't find parameters for an atom we will quit                                              
C                                                                                           
934       rewind(30)                                                                        
          write(*,1001)string(14:16)                                                        
1001      format('could not find parameters for a ',a3)                                  
          stop                                                                              
C                                                                                           
C now carry on to the next atom                                                          
C                                                                                           
939       continue                                                                          
        endif                                                                               
C                                                                                           
C keep reading through the file                                                             
C                                                                                           
      goto 1                                                                                
 99   continue                                                                              
      close (32)                                                                         

C write out the total charge on the current molecule type

      write(*,199)i_t_tp(iii),qtot
199   format('Net charge on moltype ',i8,' = ',f15.5)

C close the loop over molecules                                                             
                                                                                            
666   continue                                                                              

      if(allocated(id_tmp)) deallocate(id_tmp)
      if(allocated(jd_tmp)) deallocate(jd_tmp)
      if(allocated(kd_tmp)) deallocate(kd_tmp)

C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                                   
C xx now figure out all typs of pairwise interactions xx                                    
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                                   
                                                                                            
      allocate(d11(1:ntpfnd,1:ntpfnd),stat=ierror(1))                                
      allocate(d22(1:ntpfnd,1:ntpfnd),stat=ierror(1))                                
      allocate(s11(1:ntpfnd,1:ntpfnd),stat=ierror(1))                                
      allocate(s22(1:ntpfnd,1:ntpfnd),stat=ierror(2))                                
      allocate(s33(1:ntpfnd,1:ntpfnd),stat=ierror(3))                                
      allocate(s44(1:ntpfnd,1:ntpfnd),stat=ierror(4))                                
      allocate(e11(1:ntpfnd,1:ntpfnd),stat=ierror(5))                                
      allocate(e22(1:ntpfnd,1:ntpfnd),stat=ierror(6))                                
      allocate(e00(1:ntpfnd,1:ntpfnd),stat=ierror(6))                                
      do ijk=1,ntpfnd                                                                  
        do jkl=1,ntpfnd                                                                
          d11(ijk,jkl)=sqrt(rs0(ijk)*rs0(jkl))                                                   
          d22(ijk,jkl)=rs0(ijk)*rs0(jkl)                                                   
          s11(ijk,jkl)=rs1(ijk)*rs1(jkl)                                                    
          s22(ijk,jkl)=rs2(ijk)*rs2(jkl)                                                    
          s33(ijk,jkl)=rs3(ijk)*rs3(jkl)                                                    
          s44(ijk,jkl)=rs4(ijk)*rs4(jkl)                                                    
          e11(ijk,jkl)=re1(ijk)*re1(jkl)                                                    
          e22(ijk,jkl)=re2(ijk)*re2(jkl)    

C 22-03-08 this sets the term that needs to be added to the energy for v interactions
C so that the 'short-range' harmonic term meshes smoothly with the 'lg-range' 12-10 potential

          if(i_use_v_typ.eq.1) then
            if(myg(ijk).ne.0.and.myg(jkl).ne.0) then
              e00(ijk,jkl)=-re1(ijk)*re1(jkl)/5.000
            else
              e00(ijk,jkl)=re1(ijk)*re1(jkl)
            endif
          endif
C                                                                                            
C the follow lines make hetero-greasy interactions non-attractive                           
C                                                                                           
          if(myg(ijk).gt.0.and.myg(jkl).gt.0.and.                                           
     &       myg(ijk).ne.myg(jkl)) s22(ijk,jkl)=0.000                                       
          if(myg(ijk).gt.0.and.myg(jkl).gt.0.and.                                           
     &       myg(ijk).ne.myg(jkl)) s44(ijk,jkl)=0.000                                       
C                                                                                            
C the follow lines make greasy-nongreasy interactions like non-non                          
C                                                                                           
          if(myg(ijk).eq.0.and.myg(jkl).gt.0) then                                          
            e00(ijk,jkl)=re1(ijk)*re1(ijk)
            e11(ijk,jkl)=re1(ijk)*re1(ijk)                                                  
            e22(ijk,jkl)=re2(ijk)*re2(ijk)                                                  
          endif                                                                             
          if(myg(ijk).gt.0.and.myg(jkl).eq.0) then                                          
            e00(ijk,jkl)=re1(jkl)*re1(jkl)
            e11(ijk,jkl)=re1(jkl)*re1(jkl)                                                  
            e22(ijk,jkl)=re2(jkl)*re2(jkl)                                                  
          endif                                                                             
C                                                                                            
C note that we multiply through by e11 & e22 here to avoid repetition                       
C so that they are not actually used any more in the simulation...                          
C note also since s33 and s44 will be zero for nongreasy interactions                       
C we don't need to worry about e22 not being zero (see earlier)                             
C                                                                                           
C 27-05-07 note that we need to use e11 & e22 differently depending                         
C on whether or not we're using 12-10 potentials...                                         
C                                                                                           
C 28-05-07 above issue no lger true - we've made our 12-6 potential                       
C more like our 12-10 potential                                                             
C                                                                                           
          s11(ijk,jkl)=s11(ijk,jkl)*e11(ijk,jkl)                                            
          s22(ijk,jkl)=s22(ijk,jkl)*e22(ijk,jkl)                                            
          s33(ijk,jkl)=s33(ijk,jkl)*e11(ijk,jkl)                                            
          s44(ijk,jkl)=s44(ijk,jkl)*e22(ijk,jkl)                                            
        enddo                                                                               
      enddo                                                                                 
C                                                                                            
C write out the parameters to output - worth checking these numbers                         
C note that although e11 & e22 are written,they're not used hence...                       
C                                    

C 32-07-12 all of this stuff is currently unused...
                                                       
      if(.not.i_limit_verbosity) then
      open(37,file='potential_check.file',form='formatted')                                      
      do ijk=1,ntpfnd                                                                  
        do jkl=1,ntpfnd                                                                
          write(*,555)fattyp(ijk),fattyp(jkl),                                           
     &                s11(ijk,jkl),s22(ijk,jkl),s33(ijk,jkl),                              
     &                s44(ijk,jkl),e11(ijk,jkl),e22(ijk,jkl),
     &                d11(ijk,jkl),e00(ijk,jkl)
555       format('for a typs ',a3,' & ',a3,1x,6e14.6,2f10.3)                                     
          do ipp=1,150
            ri=(real(ipp)/10.000)
            dist2=ri**2
            fconst1=                                                                           
     &       (s11(ijk,jkl)/(dist2**7) -                                              [12-6 potential  ]
     &        s22(ijk,jkl)/(dist2**6))                                               [12-6 potential  ]
            uv1=                                                                             
     &       (s33(ijk,jkl)/(dist2**6) -                                              [12-6 potential  ]
     &        s44(ijk,jkl)/(dist2**5))                                               [12-6 potential  ]
            dist=sqrt(dist2)                                                                             
            fconst2=-r_f_st*(dist-d11(ijk,jkl))/dist                    
            uv2=   r_f_hf*(dist-d11(ijk,jkl))**2+
     &                              e00(ijk,jkl)
            write(37,19)fattyp(ijk),fattyp(jkl),
     &                  ri,fconst1,uv1,fconst2,uv2,
     &                  d11(ijk,jkl),e00(ijk,jkl)
19          format('please check ',a3,1x,a3,7f12.5)
          enddo                                                                               
        enddo                                                                               
      enddo                       
      close(37)                                                               

      endif ! close if statement for i_limit_verbosity
                                                                                             
C write out the number of charges for flexible molecules                                     
                                                                                              
      if(.not.i_limit_verbosity) then
      open(37,file='charge_check.file',form='formatted')                                      

      do 777 iii=1,num_t_ms                                                               
  
C skip rigid molecules
 
        if(ityp(iii).eq.0) cycle

C ask if we've seen this typ of molecule before                                      
                                                                                            
        do mmm=iii-1,1,-1                                                                 
          if(i_f_tp(iii).eq.i_f_tp(mmm)) goto 777                               
        enddo                                                                             

        do jjj=1,num_atms_pt(i_f_tp(iii))                                           
          write(37,73)i_t_tp(iii),jjj,
     &      c_typ(i_t_tp(iii))%x(jjj),
     &      c_typ(i_t_tp(iii))%y(jjj),
     &      c_typ(i_t_tp(iii))%z(jjj),
     &      c_typ(i_t_tp(iii))%q(jjj)
        enddo                                                                             
777   enddo                                                                                 

      close(37)                                                                             
      endif

73    format('FLEX ',2i5,4f10.3)                                                            
74    format('RIGD ',2i5,4f10.3)                                                            
C                                                                                           
C some format statements                                                                    
C                                                                                           
 500  format(a100)                                                                             
 601  format('No of charges read for molecule ',i2,':',i5,/,                               
     &       'total charge                    :',f8.3)                                    
 602  format('No of surfas read for molecule ',i2,':',i5,/,                             
     &       'total r_scale factor                    :',f8.3)                              
 603  format('No of goas read for molecule ',i2,':',i5)                                  
C                                                                                           
C deallocate some arrays                                                                    
C                                                                                           
      deallocate(re1)     
      deallocate(re2)    
      deallocate(rs1)   
      deallocate(rs2)  
      deallocate(rs3) 
      deallocate(rs4)
      deallocate(myg)
      close(30)     
                              
C quit and go back                                                                          
                               
      return                    
      end                        
