
      subroutine write_movie_file(movieframe,crashed,mirep)
                                                     
      use allocatable_arrays                        
      implicit real(a-h,o-z)                       
                                                  
      dimension ax(3),ay(3),az(3)                
      character*10000 chain                     
      character*26 blah                        
C     character*36 filename                                                                           
      character*54 filename   ! new for 64-11-11                                                      
      character*5 poop 
      character*4 cack1
      character*3 atemp
      character(LEN=9) :: FileNmStr             
      blah='ABCDEFGHIJKLMNOPQRSTUVWXYZ'        

      filename='MOVIE/movie.XXX.'//FileNmStr(movieFrameNumber)//'.pdb'

      write(atemp,'(i3)')mirep
      do j=1,3
        if(atemp(j:j).eq.' ') atemp(j:j)='0'
      enddo
      filename(13:15)=atemp

      sclfac=1.00000000                       
      tmin=xmin
      if(ymin.lt.tmin) tmin=ymin
      if(zmin.lt.tmin) tmin=zmin
      tmin=abs(tmin) ! we'll add this on!
      tmin=0.0       ! no we won't!!!

      if(icrashed.eq.1) then    
        open(unit=45+mirep,file='MOVIE/crashed.pdb',status='unknown') 
      else 
301     continue ! should come back here if failed on open  
        write(*,900)movieframe                            
900     format("writing movie number = ",i10)                           
        open(unit=45+mirep,file=filename,status='unknown',err=301) 
      endif            

      write(45+mirep,40)tot_time,sclfac          
      write(45+mirep,11)xlen*sclfac,ylen*sclfac,zlen*sclfac
11    format('CRYST1',3f9.3,'  90.00  90.00  90.00 P 1           1')
      nawrite=0                                   
      do iii=1,num_t_ms                          

        write(45+mirep,88)iii                         
88      format('REMARK molecule number ',i8)   
        ncrap=0                               

C TEMP TEMP TEMP TEMP TEMP - just doing this to get *something* written                                                                                                     
C       if(ityp(iii).eq.1) then                                                                  
        if(ityp(iii).eq.1.or.ityp(iii).eq.0) then    

c         inumber=mod(i_t_tp(iii),26) ! do this for moltype in chainid
          inumber=mod(iii,26)                       
          if(inumber.eq.0) inumber=26              

          do jjj=1,num_atms_pt(i_t_tp(iii))       
            x1=c_a_x(1,id_beg(iii)+jjj,mirep)
            y1=c_a_x(2,id_beg(iii)+jjj,mirep)
            z1=c_a_x(3,id_beg(iii)+jjj,mirep)       

C note the new changes to the code here - all coordinates are shifted by
C half a box length when we write out - this means that their
C coords go from 0 to xlen etc so that they are aligned with pbc_box
C note that if wrap_molecules=2 then we figure out what change would be
C applied to the 1st atom of that molecule and then we apply this same
C change to all other atoms of the molecule

C note that if we're using i_debug then we write out coords as they are 
C read-in from the restart.file - no changes at all...

            if(.not.i_debug) then 
              if(wrap_molecules.eq.0) then
                x1=x1+xlen*0.5
                y1=y1+ylen*0.5
                z1=z1+zlen*0.5
              elseif(wrap_molecules.eq.1) then
                x1=x1-xlen*anint(x1*xinv)
                y1=y1-ylen*anint(y1*yinv)
                z1=z1-zlen*anint(z1*zinv)
                x1=x1+xlen*0.5
                y1=y1+ylen*0.5
                z1=z1+zlen*0.5
              elseif(wrap_molecules.eq.2) then
                if(jjj.eq.1) then
                  x1beg=x1
                  y1beg=y1
                  z1beg=z1
                  x1end=x1-xlen*anint(x1*xinv)
                  y1end=y1-ylen*anint(y1*yinv)
                  z1end=z1-zlen*anint(z1*zinv)
                  x1mov=x1end-x1beg
                  y1mov=y1end-y1beg
                  z1mov=z1end-z1beg
                endif
                x1=x1+x1mov+xlen*0.5
                y1=y1+y1mov+ylen*0.5
                z1=z1+z1mov+zlen*0.5
              endif
            endif
            nawrite=nawrite+1                
            ncrap=ncrap+1                   

C note that we write atoms that have not yet been released in a
C different format (using 'Z's as identifiers) and we put their
C coordinates in the same place as the last free atom - but note that we
C only do this if .not.i_debug...

C (for assigning bonds with conect statements it's easier if we still
C consider the unfree as atoms contributing to nawrite...)
C note that we only do this if 

            iresnum=mod(nawrite+10000,10000)       
            if(ifree(nawrite).eq.0) then  
              if(.not.i_debug) then
                x1=c_a_x(1,jfree(nawrite),mirep)   
                y1=c_a_x(2,jfree(nawrite),mirep)  
                z1=c_a_x(3,jfree(nawrite),mirep) 
                if(wrap_molecules.eq.0) then
                  x1=x1+xlen*0.5
                  y1=y1+ylen*0.5
                  z1=z1+zlen*0.5
                elseif(wrap_molecules.eq.1) then
                  x1=x1-xlen*anint(x1*xinv)
                  y1=y1-ylen*anint(y1*yinv)
                  z1=z1-zlen*anint(z1*zinv)
                  x1=x1+xlen*0.5
                  y1=y1+ylen*0.5
                  z1=z1+zlen*0.5
                elseif(wrap_molecules.eq.2) then
                  if(jjj.eq.1) then
                    x1beg=x1
                    y1beg=y1
                    z1beg=z1
                    x1end=x1-xlen*anint(x1*xinv)
                    y1end=y1-ylen*anint(y1*yinv)
                    z1end=z1-zlen*anint(z1*zinv)
                    x1mov=x1end-x1beg
                    y1mov=y1end-y1beg
                    z1mov=z1end-z1beg
                  endif
                  x1=x1+x1mov+xlen*0.5
                  y1=y1+y1mov+ylen*0.5
                  z1=z1+z1mov+zlen*0.5
                endif
              endif
              write(45+mirep,29) 
     &        nawrite,    
     &        n_typ(i_t_tp(iii))%ch2(jjj),
     &        n_typ(i_t_tp(iii))%ch3(jjj),
     &        'Z',                     
     &        n_typ(i_t_tp(iii))%ch4(jjj),
     &        x1+tmin,y1+tmin,z1+tmin,
     &        real(iii)/100.0,
     &        real(jjj)*50.000/real(num_atms_pt(i_t_tp(iii)))  
            else

C if doing protease then we need to check the atom names haven't changed

              cack1=n_typ(i_t_tp(iii))%ch2(jjj)
              if(i_do_protease) then
                if(name_override(nawrite).eq.-1) cack1(2:2)='O'
                if(name_override(nawrite).eq.1)  cack1(2:2)='N' 
              endif

              write(45+mirep,29)                           
     &        nawrite,                              
     &        cack1,                         
     &        n_typ(i_t_tp(iii))%ch3(jjj),
     &        blah(inumber:inumber),
     &        n_typ(i_t_tp(iii))%ch4(jjj),
     &        x1+tmin,y1+tmin,z1+tmin,
     &        real(iii)/100.0,
     &        real(jjj)*50.000/real(num_atms_pt(i_t_tp(iii)))
29    format('ATOM',i7,1x,a4,1x,a3,1x,a,a4,4x,3f8.3,2f6.2)  

            endif
          enddo                                            
        endif                                    
        write(45+mirep,99)                            
99      format('TER')                          
35    format('ATOM',i7,1x,a4,1x,a3,1x,a,a4,4x,3f12.5,'  1.00',f6.2)
30    format('ATOM',i7,1x,a4,1x,a3,1x,a,a4,4x,3f8.3,'  1.00',f6.2)
      enddo                                                     

C if afm write out the position of the tip

      if(afm) then
        nawrite=nawrite+1
        write(45+mirep,30)
     &        nawrite,
     &        ' X  ','XXX','?','   1',
     &        afm_tip_x_crd*sclfac,
     &        afm_tip_y_crd*sclfac,
     &        afm_tip_z_crd*sclfac,
     &        0.0
      endif

      nawriteorig=nawrite                       
      nawrite=0                                
      do iii=1,num_f_ms                       
        if(iii.gt.1) nawrite=nawrite+
     &               num_atms_pt(i_f_tp(iii-1))  
        do jjj=1,num_bonds(i_f_tp(iii))         

C note that we now skip atoms that have not yet been released

        if(c_f_m(iii)%free
     &    (nbb1(jjj+id_beg_bonds(i_f_tp(iii)))).eq.0) cycle
        if(c_f_m(iii)%free
     &    (nbb2(jjj+id_beg_bonds(i_f_tp(iii)))).eq.0) cycle

C new for 046-11-11 code - ask if the distance between the two atoms
C exceeds half of the box length. if so, then don't write a CONECT

        iat=nbb1(jjj+id_beg_bonds(i_f_tp(iii)))+nawrite
        jat=nbb2(jjj+id_beg_bonds(i_f_tp(iii)))+nawrite
        dist=sqrt((c_a_x(1,iat,mirep)-c_a_x(1,jat,mirep))**2+
     &            (c_a_x(2,iat,mirep)-c_a_x(2,jat,mirep))**2+
     &            (c_a_x(3,iat,mirep)-c_a_x(3,jat,mirep))**2)
        if(dist.gt.xlen*0.50) then
c         write(*,*)'Skipping CONECT for ',iat,' & ',jat
          cycle
        endif

        write(45+mirep,31)nbb1(jjj+id_beg_bonds(i_f_tp(iii)))+nawrite, 
     &              nbb2(jjj+id_beg_bonds(i_f_tp(iii)))+nawrite 
        enddo                                                  
      enddo                                                   
31    format('CONECT',2i5)                                   
      nawrite=nawriteorig                                   
      close(45+mirep)                                            
      if(movieframe.eq.0.and.mirep.eq.1) then                            
75003   open(unit=15,file='MOVIE/box.pdb',
     &       status='unknown',err=75003)
      nawrite=1                                          
      write(15,50)nawrite,nawrite,xmin,ymin,zmin
      nawrite=nawrite+1     
      write(15,50)nawrite,nawrite,xmin,ymin,zmax
      nawrite=nawrite+1     
      write(15,50)nawrite,nawrite,xmin,ymax,zmin
      nawrite=nawrite+1     
      write(15,50)nawrite,nawrite,xmin,ymax,zmax
      nawrite=nawrite+1     
      write(15,50)nawrite,nawrite,xmax,ymin,zmin
      nawrite=nawrite+1     
      write(15,50)nawrite,nawrite,xmax,ymin,zmax
      nawrite=nawrite+1     
      write(15,50)nawrite,nawrite,xmax,ymax,zmin
      nawrite=nawrite+1     
      write(15,50)nawrite,nawrite,xmax,ymax,zmax
      nawrite=nawrite+1     
      write(15,61)          
      write(15,62)          
      write(15,63)          
      write(15,64)          
      write(15,65)          
      write(15,66)          
      write(15,67)          
      write(15,68)                  
      write(15,69)                  
      write(15,70)                  
      write(15,71)                  
      write(15,72)                  
61    format('CONECT    1    2')    
62    format('CONECT    1    3')    
63    format('CONECT    1    5')    
64    format('CONECT    2    4')    
65    format('CONECT    2    6')    
66    format('CONECT    3    4')    
67    format('CONECT    3    7')    
68    format('CONECT    4    8')    
69    format('CONECT    5    6')    
70    format('CONECT    5    7')    
71    format('CONECT    6    8')    
72    format('CONECT    7    8')    
      close(15)                     
      endif                         
10    format('HEADER ',6f10.3)      
20    format(i5)                    
40    format('REMARK time= ',f15.5,' ps  sclfac =',f12.5)      
50    format('ATOM  ',i5,'  X   BOX',2x,i4,4x,3f8.3)  
      return                                          
      end                                             
