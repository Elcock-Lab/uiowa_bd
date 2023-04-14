
      subroutine make_internal_parameters 
                                         
      use allocatable_arrays            
      implicit real(a-h,o-z)           
                                      
      character*128 string           
      character*4 a                 
      character*4 ajk1            
      data a /'ATOM'/             
      dimension coords(3)        
      character*6 name1         
      character*4 name2        
      character*3 name3       
      character*5 name4      
      character*2 name5     
      character*1 attyp    
      character*1 fattyp(20)         
      integer     ierr(100)         
C                                  
C do a loop over all molecules                                                                        
C                      
      do 666 iii=1,num_t_ms 
                             
C first ask if we've seen this typ of molecule before                                                

        do mmm=iii-1,1,-1     
          if(i_t_tp(iii).eq.i_t_tp(mmm)) goto 666   
        enddo                  

C if it's rigid simply read the generalized diffusion tensors

        if(ityp(iii).eq.0) then
          open(unit=46,file=bnd_file(i_t_tp(iii)),status='unknown')
          read(46,*) ! skip title line
          read(46,*)D_11(i_r_tp(iii)),a,b,c,d,e
          read(46,*)a,D_22(i_r_tp(iii)),a,b,c,d
          read(46,*)a,b,D_33(i_r_tp(iii)),a,b,c
          read(46,*)a,b,c,D_44(i_r_tp(iii)),a,b
          read(46,*)a,b,c,d,D_55(i_r_tp(iii)),a
          read(46,*)a,b,c,d,e,D_66(i_r_tp(iii))
          E_11(i_r_tp(iii))=sqrt(2.0*D_11(i_r_tp(iii))*time_step)
          E_22(i_r_tp(iii))=sqrt(2.0*D_22(i_r_tp(iii))*time_step)
          E_33(i_r_tp(iii))=sqrt(2.0*D_33(i_r_tp(iii))*time_step)
          E_44(i_r_tp(iii))=sqrt(2.0*D_44(i_r_tp(iii))*time_step)
          E_55(i_r_tp(iii))=sqrt(2.0*D_55(i_r_tp(iii))*time_step)
          E_66(i_r_tp(iii))=sqrt(2.0*D_66(i_r_tp(iii))*time_step)
          close(46)
          goto 666
        endif
                    
C if it's flexible, zero out the bond arrays in prep for next stage
                   
        num_bonds(i_t_tp(iii))=0  
        num_angls(i_t_tp(iii))=0 
        num_dihes(i_t_tp(iii))=0

666   continue                                                                                        
                                                                                                      
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                                     
C xx read the bonds,etc from the bnd_files xx                                     
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                                     
                                                                                                      
C note that i_done is just a counter that says whether we've identified                               
C the very first bond in the lg list of bonds for molecule 'ick'                                    
                                                                                                      
      do ick=1,num_t_ts                                                                       
        i_done_bonds(ick)=0                                                                         
        i_done_angls(ick)=0                                                                        
        i_done_dihes(ick)=0                                                                     
      enddo                                                                                         
      num_tot_bonds=0                                                                                 
      num_tot_angls=0                                                                                
      num_tot_dihes=0                                                                             
                                                                                                      
C do a loop over all molecules                                                                        
                                                                                                      
      do 777 iii=1,num_t_ms                                                                         

C skip this molecule if it's rigid

        if(ityp(iii).eq.0) goto 777                                                              

C skip this molecule if we've done it before

        do mmm=iii-1,1,-1                                                                           
          if(i_t_tp(iii).eq.i_t_tp(mmm)) goto 777                                         
        enddo                                                                                       

C otherwise, open the internal parameters file

        open(unit=46,file=bnd_file(i_t_tp(iii)),status='unknown')                                              
                                                                                                      
C 01-02-07 AHE : read the whole thing twice,the first time through                                   
C we're only interested in finding total number of bonds,angles,dihedrls                              
C so that we can efficiently allocate memory for them...                                              
                                                                                                      
        read(46,'(a4,3x,i8)')ajk1,ijk2                                                                    
        if(ajk1.eq.'bond') then                                                               
          do m=1,ijk2                                                                               
            read(46,1030)ic1,ic2,ic3,fc1,fc2                                                          

C 2019 - skip this bond if all atoms in it are static

            if(c_typ(ic1)%s(ic2).eq.1.and.
     &         c_typ(ic1)%s(ic3).eq.1) cycle

            num_tot_bonds=num_tot_bonds+1                                                                 
            if(i_done_bonds(ic1).eq.0) then                                                           
              id_beg_bonds(ic1)=num_tot_bonds-1                                                         
              i_done_bonds(ic1)=1                                                                     
            endif                                                                                     
          enddo                                                                                       
        endif                                                                                         
        read(46,'(a4,3x,i8)')ajk1,ijk2                                                                    
        if(ajk1.eq.'angl') then                                                               
          do m=1,ijk2                                                                               
            read(46,1031)ic1,ic2,ic3,ic4,fc1,fc2                                                      

C 2019 - skip this angl if all atoms in it are static

            if(c_typ(ic1)%s(ic2).eq.1.and.
     &         c_typ(ic1)%s(ic3).eq.1.and.
     &         c_typ(ic1)%s(ic4).eq.1) cycle

            num_tot_angls=num_tot_angls+1                                                               
            if(i_done_angls(ic1).eq.0) then                                                          
              id_beg_angls(ic1)=num_tot_angls-1                                                       
              i_done_angls(ic1)=1                                                                    
            endif                                                                                     
          enddo                                                                                       
        endif                                                                                         
        read(46,'(a4,3x,i8)')ajk1,ijk2                                                                    
        if(ajk1.eq.'dihe') then                                                               
          do m=1,ijk2                                                                               

C 2022 only read atom indices, 2 energies and 2 min-energy angles

            read(46,1032)ic1,ic2,ic3,ic4,ic5,fc1,fc2,fc3,fc4                                         
CCC         read(46,1032)ic1,ic2,ic3,ic4,ic5,fc1,fc2,fc3,fc4,                                        
CCC  &                   fc5,fc6,fc7,fc8,fc9,fc10                                                     

C 2019 - skip this dihe if all atoms in it are static

            if(c_typ(ic1)%s(ic2).eq.1.and.
     &         c_typ(ic1)%s(ic3).eq.1.and.
     &         c_typ(ic1)%s(ic4).eq.1.and.
     &         c_typ(ic1)%s(ic5).eq.1) cycle

            num_tot_dihes=num_tot_dihes+1                                                         
            if(i_done_dihes(ic1).eq.0) then                                                       
              id_beg_dihes(ic1)=num_tot_dihes-1                                                 
              i_done_dihes(ic1)=1                                                                 
            endif                                                                                     
          enddo                                                                                       
        endif                                                                                         
        close(46)                                                                                     
777   enddo                                                                                              
C                                                                                                     
C at this point we can allocate memory for the internal parameters                                    
C      
      write(*,*)'about to allocate memory for bonds etc'
      allocate(nbb1(1:num_tot_bonds),stat=ierr(1))                                                 
      allocate(nbb2(1:num_tot_bonds),stat=ierr(2))                                                 
      allocate(rkbb(1:num_tot_bonds),stat=ierr(3))                                                 
      allocate(r0bb(1:num_tot_bonds),stat=ierr(4))                                                 
      allocate(nba1(1:num_tot_angls),stat=ierr(5))                                                
      allocate(nba2(1:num_tot_angls),stat=ierr(6))                                                
      allocate(nba3(1:num_tot_angls),stat=ierr(7))                                                
      allocate(rkba(1:num_tot_angls),stat=ierr(8))                                                
      allocate(r0ba(1:num_tot_angls),stat=ierr(9))                                                
      allocate(nbd1(1:num_tot_dihes),stat=ierr(10))                                            
      allocate(nbd2(1:num_tot_dihes),stat=ierr(11))                                            
      allocate(nbd3(1:num_tot_dihes),stat=ierr(12))                                            
      allocate(nbd4(1:num_tot_dihes),stat=ierr(13))                                            
      allocate(r0bd(1:num_tot_dihes),stat=ierr(13))                                            
      allocate(tors1(4,1:num_tot_dihes),stat=ierr(14))                                         
      allocate(tors3(4,1:num_tot_dihes),stat=ierr(15))                                         
      write(*,*)'just allocated memory for bonds etc'
C                                                                                                     
C now do the whole thing all over again and read everything properly...                              
C                                                                                                     
      num_tot_bonds=0                                                                                 
      num_tot_angls=0                                                                                
      num_tot_dihes=0                                                                             
C                                                                                                     
C do a loop over all molecules                                                                        
C                                                                                                     
      do 888 iii=1,num_t_ms                                                                         
        if(ityp(iii).eq.0) goto 888                                                              
        do mmm=iii-1,1,-1                                                                           
          if(i_f_tp(iii).eq.i_f_tp(mmm)) goto 888                                         
        enddo

C 02-08-07 ic1 (i_f_tp) used to be read from each line of bnd_file
C          now,set ic1 = i_f_tp(iii) here,so ic1 is determined
C          independently of the bnd_file. this means that if we want to 
C          use the same molecule for a different simulation,in which it
C          might be the 4th f typ,instead of the 3rd,we don't need
C          to mess around and remake the bnd_file

        ic1=i_f_tp(iii)                                                                                       
        open(unit=46,file=bnd_file(i_t_tp(iii)),status='unknown')                                              
        read(46,'(a4,3x,i8)')ajk1,ijk2                                                                    
        if(ajk1.eq.'bond') then                                                               
          do m=1,ijk2                                                                               
            read(46,1030)itt,ic2,ic3,fc1,fc2                                                          

C 2019 - skip this bond if all atoms in it are static

            if(c_typ(ic1)%s(ic2).eq.1.and.
     &         c_typ(ic1)%s(ic3).eq.1) cycle

            num_tot_bonds=num_tot_bonds+1                                                                 
            num_bonds(ic1)=num_bonds(ic1)+1                                                           
            nbb1(num_bonds(ic1)+id_beg_bonds(ic1))=ic2
            nbb2(num_bonds(ic1)+id_beg_bonds(ic1))=ic3
            rkbb(num_bonds(ic1)+id_beg_bonds(ic1))=fc1
            r0bb(num_bonds(ic1)+id_beg_bonds(ic1))=fc2
1030        format(3i10,2f10.5)
          enddo              
        endif                
        read(46,'(a4,3x,i8)')ajk1,ijk2  
        if(ajk1.eq.'angl') then     
          do m=1,ijk2                    
            read(46,1031)itt,ic2,ic3,ic4,fc1,fc2 

C 2019 - skip this angl if all atoms in it are static

            if(c_typ(ic1)%s(ic2).eq.1.and.
     &         c_typ(ic1)%s(ic3).eq.1.and.
     &         c_typ(ic1)%s(ic4).eq.1) cycle

            num_tot_angls=num_tot_angls+1                                                               
            num_angls(ic1)=num_angls(ic1)+1                                                         
            nba1(num_angls(ic1)+id_beg_angls(ic1))=ic2                                              
            nba2(num_angls(ic1)+id_beg_angls(ic1))=ic3                                              
            nba3(num_angls(ic1)+id_beg_angls(ic1))=ic4                                              
            rkba(num_angls(ic1)+id_beg_angls(ic1))=fc1                                              
            r0ba(num_angls(ic1)+id_beg_angls(ic1))=fc2                                              
1031        format(4i10,2f10.5)                                                                       
          enddo                                                                                       
        endif                                                                                         
        read(46,'(a4,3x,i8)')ajk1,ijk2                                                                    
        if(ajk1.eq.'dihe') then                                                               
          do m=1,ijk2                                                                               

C 2022 only read atom indices, 2 energies and 2 min-energy angles

            read(46,1032)itt,ic2,ic3,ic4,ic5,fc1,fc2,fc3,fc4                                         
            fc5 =fc1
            fc6 =cos(fc3)
            fc7 =sin(fc3)
            fc8 =fc2 
            fc9 =cos(fc4)
            fc10=sin(fc4)
CCC         read(46,1032)itt,ic2,ic3,ic4,ic5,fc1,fc2,fc3,fc4,                                        
CCC  &                   fc5,fc6,fc7,fc8,fc9,fc10                                                     

C 2019 - skip this dihe if all atoms in it are static

            if(c_typ(ic1)%s(ic2).eq.1.and.
     &         c_typ(ic1)%s(ic3).eq.1.and.
     &         c_typ(ic1)%s(ic4).eq.1.and.
     &         c_typ(ic1)%s(ic5).eq.1) cycle

            num_tot_dihes=num_tot_dihes+1                                                         
            num_dihes(ic1)=num_dihes(ic1)+1                                                   
            nbd1(num_dihes(ic1)+id_beg_dihes(ic1)) =ic2                                       
            nbd2(num_dihes(ic1)+id_beg_dihes(ic1)) =ic3                                       
            nbd3(num_dihes(ic1)+id_beg_dihes(ic1)) =ic4                                       
            nbd4(num_dihes(ic1)+id_beg_dihes(ic1)) =ic5                                       
            r0bd(num_dihes(ic1)+id_beg_dihes(ic1)) =fc3                                       
            tors1(1,num_dihes(ic1)+id_beg_dihes(ic1))=fc5                                     
            tors1(3,num_dihes(ic1)+id_beg_dihes(ic1))=fc6                                     
            tors1(4,num_dihes(ic1)+id_beg_dihes(ic1))=fc7                                     
            tors3(1,num_dihes(ic1)+id_beg_dihes(ic1))=fc8                                     
            tors3(3,num_dihes(ic1)+id_beg_dihes(ic1))=fc9                                     
            tors3(4,num_dihes(ic1)+id_beg_dihes(ic1))=fc10                                    
1032      format(5i10,10f10.5)                                                                        
          enddo                                                                                       
        endif                                                                                         
        close(46)                                                                                     
888   enddo                                                                                         
1033  format(i5)                                                                                      
1034  format(2i5)                                                                                     
1035  format(10i5)                                                                                    

C now allocate and fill arrays that store the identities of the 
C atoms that are involved in a bonding interaction with each atom
C note that the following allocation (to 1000) seems to be too much 
C but in some cases may *not* be... it works fine for regular fine
C or coarse models but fails with 50 for cytoplasm models (hence 1000)
        
      do itp=1,num_f_ts                                                                               
        write(*,*)'now working on molecule type ',itp
        num_tmp=1000 ! assumes no more than 1000 bond exclusions per atom
        allocate(b_typ(itp)%id2(b_typ(itp)%num,1:num_tmp))
        do jjj=1,b_typ(itp)%num                                                                  
          b_typ(itp)%id1(jjj)=0
        enddo                                                                                         
        do kkk=1,num_bonds(itp)                                                                       
          jjj=nbb1(kkk+id_beg_bonds(itp))                                                             
          jj2=nbb2(kkk+id_beg_bonds(itp))                                                             
          b_typ(itp)%id1(jjj)=   
     &    b_typ(itp)%id1(jjj)+1
          b_typ(itp)%id2(jjj,
     &    b_typ(itp)%id1(jjj))=jj2
          b_typ(itp)%id1(jj2)=
     &    b_typ(itp)%id1(jj2)+1
          b_typ(itp)%id2(jj2,
     &    b_typ(itp)%id1(jj2))=jjj
        enddo                                                                                         
        do kkk=1,num_angls(itp)                                                                      
          jjj=nba1(kkk+id_beg_angls(itp))                                                            
          jj2=nba3(kkk+id_beg_angls(itp))                                                            
          b_typ(itp)%id1(jjj)=
     &    b_typ(itp)%id1(jjj)+1
          b_typ(itp)%id2(jjj,
     &    b_typ(itp)%id1(jjj))=jj2
          b_typ(itp)%id1(jj2)=
     &    b_typ(itp)%id1(jj2)+1
          b_typ(itp)%id2(jj2,
     &    b_typ(itp)%id1(jj2))=jjj
        enddo                                                                                         
        do kkk=1,num_dihes(itp)                                                                   
          jjj=nbd1(kkk+id_beg_dihes(itp))                                                         
          jj2=nbd4(kkk+id_beg_dihes(itp))                                                         
          b_typ(itp)%id1(jjj)=
     &    b_typ(itp)%id1(jjj)+1
          b_typ(itp)%id2(jjj,
     &    b_typ(itp)%id1(jjj))=jj2
          b_typ(itp)%id1(jj2)=
     &    b_typ(itp)%id1(jj2)+1
          b_typ(itp)%id2(jj2,
     &    b_typ(itp)%id1(jj2))=jjj
        enddo                                                                                         
      enddo                                                                                           
      if(.not.i_limit_verbosity) then
      do itp=1,num_f_ts                                                                               
        do jtp=1,b_typ(itp)%num                                                                         
          do ktp=1,b_typ(itp)%id1(jtp)                                                                    
            write(16,221)itp,jtp,b_typ(itp)%id1(jtp),
     &                           b_typ(itp)%id2(jtp,ktp)                              
          enddo                                                                                       
        enddo                                                                                       
      enddo                                                                                       
      endif

C 407-05-09
C convert bond exclusions into one massive 1D array
C num_tot_b_typ is the total # of types of exclusions - only one copy
C included for each molecule type

      num_tot_b_typ=0
      do itp=1,num_f_ts
        do jtp=1,b_typ(itp)%num                                                                         
          do ktp=1,b_typ(itp)%id1(jtp)                                                                    
            num_tot_b_typ=num_tot_b_typ+1
          enddo  
        enddo  
      enddo  
      allocate(itp_tot_b_typ(1:num_tot_b_typ))
      allocate(ida_tot_b_typ(1:num_tot_b_typ))
      allocate(jda_tot_b_typ(1:num_tot_b_typ))
      allocate(ind_tot_b_typ(1:num_flx_atms))
      allocate(jnd_tot_b_typ(1:num_flx_atms))

C for each unique atom prepare to store its first and last entries
C in the list that goes from 1:num_tot_b_typ
C if we have that information then we can look up the entries for
C each actual atom by using my_uniq_atm_num(1:num_flx_atms)

      allocate(ibeg_tmp_b_typ(1:num_uniqs))
      allocate(iend_tmp_b_typ(1:num_uniqs))
      ibeg_tmp_b_typ=-1000
      iend_tmp_b_typ=-1000

      write(*,*)'assigning 1D bond exclusion array ',num_tot_b_typ,
     &  num_uniqs

C note that we're just looking at each molecule type *once* here

      num_tot_b_typ=0
      num_uniq_atms_tmp=0
      do itp=1,num_f_ts
        do jtp=1,b_typ(itp)%num                                                                         
          num_uniq_atms_tmp=num_uniq_atms_tmp+1
          do ktp=1,b_typ(itp)%id1(jtp)                                                                    
            num_tot_b_typ=num_tot_b_typ+1

C to find the first and last entries associated with this unique atom we
C just use the first and last entries for ktp...

C AHE is not yest sure that this is perfect code but it has at least
C been fixed so that it works correctly with one-bond molecules

            if(ktp.eq.1) then
              ibeg_tmp_b_typ(num_uniq_atms_tmp)=num_tot_b_typ
            endif
            if(ktp.eq.b_typ(itp)%id1(jtp)) then
              iend_tmp_b_typ(num_uniq_atms_tmp)=num_tot_b_typ
            endif

            itp_tot_b_typ(num_tot_b_typ)=itp                     
            ida_tot_b_typ(num_tot_b_typ)=jtp                      
            jda_tot_b_typ(num_tot_b_typ)=b_typ(itp)%id2(jtp,ktp)
            if(.not.i_limit_verbosity)
     &      write(17,91)num_tot_b_typ,
     &                  itp_tot_b_typ(num_tot_b_typ),
     &                  ida_tot_b_typ(num_tot_b_typ),
     &                  jda_tot_b_typ(num_tot_b_typ)
91          format('bond exclusion pairs ',4i10)
          enddo  
        enddo  
      enddo  

C whereas now we're looking at each atom in the system

      write(*,*)'indexing 1D bond exclusion array for all atoms'
      do iii=1,num_flx_atms   
        kkk=iii-id_beg(mol_of_atm(iii))
        itp=i_f_tp(mol_of_atm(iii))
        ind_tot_b_typ(iii)=-1000 
        jnd_tot_b_typ(iii)=-1000 

C by defining ibeg_tmp_b_typ and iend_tmp_b_typ above we should be
C able to do a much more focused search loop here
C TEMP TEMP TEMP NOTE THAT IT DOESN"T WORK YET!!!
C UPDATE - IT MAY NOW WORK NOW THAT I've ADDED IN THE if(ibeg...)

        jkl=my_uniq_atm_num(iii)

C note that we only do this if this is an atom that has bond exclusions
C - it's quite possible that atoms won't (especially ions)

        if(ibeg_tmp_b_typ(jkl).ne.-1000) then
          do jjj=ibeg_tmp_b_typ(jkl),iend_tmp_b_typ(jkl)
ccc       do jjj=1,num_tot_b_typ           
            if(itp_tot_b_typ(jjj).lt.itp) cycle
            if(itp_tot_b_typ(jjj).gt.itp) goto 999
            if(ida_tot_b_typ(jjj).eq.kkk.and.
     &         ind_tot_b_typ(iii).eq.-1000)
     &         ind_tot_b_typ(iii)=jjj
            if(ida_tot_b_typ(jjj).eq.kkk)      
     &         jnd_tot_b_typ(iii)=jjj
          enddo
        endif
        if(ind_tot_b_typ(iii).eq.-1000) ind_tot_b_typ(iii)=0
        if(jnd_tot_b_typ(iii).eq.-1000) jnd_tot_b_typ(iii)=-1
999     continue
      enddo 
      write(*,*)'done with the indexing'
      if(.not.i_limit_verbosity) then
      open(unit=18,file='bond_exclusion.txt',status='unknown')
      do iii=1,num_flx_atms   
        write(18,81)iii,ind_tot_b_typ(iii),jnd_tot_b_typ(iii)
81      format('bond exclusion indices ',3i10)
      enddo
      close(18)
      endif
221   format(4i10)                                                                                
C                                                                                                     
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                              
C xx write out the bonds,angles,dihedrals to internal.parameters.final                              
C xx note that this section is largely superfluous now                                                
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                              
C
      if(.not.i_limit_verbosity) then                                                                                                     
      open(unit=47,file='internal.parameters.final',                                            [floop     beg]
     &     status='unknown')                                                                          
      ntotbond=    0                                                                                  
      ntotangle=   0                                                                                  
      ntotdihedral=0                                                                                  
      do itp=1,num_f_ts                                                                               
        ntotbond    =ntotbond    +num_bonds(itp)                                                      
        ntotangle   =ntotangle   +num_angls(itp)                                                     
        ntotdihedral=ntotdihedral+num_dihes(itp)                                                  
      enddo                                                                                           
      write(47,1129)ntotbond                                                                          
1129  format('bonds     ',i5)                                                                         
      do itp=1,num_f_ts                                                                               
        do jjj=1,num_bonds(itp)                                                                       
          write(47,1030)itp,nbb1(jjj+id_beg_bonds(itp)),                                             
     &                      nbb2(jjj+id_beg_bonds(itp)),                                             
     &                      rkbb(jjj+id_beg_bonds(itp)),                                             
     &                      r0bb(jjj+id_beg_bonds(itp))                                               
        enddo                                                                                         
      enddo                                                                                           
      write(47,1139)ntotangle                                                                         
1139  format('angles    ',i5)                                                                         
      do itp=1,num_f_ts                                                                               
        do jjj=1,num_angls(itp)                                                                      
          write(47,1031)itp,nba1(jjj+id_beg_angls(itp)),                                            
     &                      nba2(jjj+id_beg_angls(itp)),                                            
     &                      nba3(jjj+id_beg_angls(itp)),                                            
     &                      rkba(jjj+id_beg_angls(itp)),                                            
     &                      r0ba(jjj+id_beg_angls(itp))                                              
        enddo                                                                                         
      enddo                                                                                           
      write(47,1149)ntotdihedral                                                                      
1149  format('dihedrals ',i5)                                                                         
      do itp=1,num_f_ts                                                                               
        do jjj=1,num_dihes(itp)                                                                   
          write(47,1032)itp,nbd1(jjj+id_beg_dihes(itp)),                                         
     &                      nbd2(jjj+id_beg_dihes(itp)),                                         
     &                      nbd3(jjj+id_beg_dihes(itp)),                                         
     &                      nbd4(jjj+id_beg_dihes(itp))                                           
        enddo                                                                                         
      enddo                                                                                           
      close(47)
      endif ! i_limit_verbosity                                                                                  [floop     end]
C                                                                                                     
C some format statements                                                                              
C                                                                                                     
 500  format(a)                                                                                       
 501  format(30x,3f8.3,f15.8)                                                                         
 502  format(5x,a6,1x,a4,1x,a3,1x,a5,4x,3f8.3,18x,a2)                                                 
 601  format('No of charges read for molecule ',i2,':',i5,/,                                         
     &         'total charge                    :',f8.3)                                              
 602  format('No of surfas read for molecule ',i2,':',i5,/,                                       
     &         'total r_scale factor                    :',f8.3)                                        
 603  format('No of goas read for molecule ',i2,':',i5)                                            
C                                                                                                     
C quit and go back                                                                                    
C                                                                                                     
      return                                                                                          
      end                                                                                             
