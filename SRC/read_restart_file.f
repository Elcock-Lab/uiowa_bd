
      subroutine read_restart_file  
                                      
      use allocatable_arrays         
      implicit real(a-h,o-z)        
      character*80             :: string                      
      character*16 restart_file_name
      character*3  atemp

      do ijk=1,num_reps

        restart_file_name='restart.file.XXX'
        write(atemp,'(i3)')ijk
        do j=1,3
          if(atemp(j:j).eq.' ') atemp(j:j)='0'
        enddo
        restart_file_name(14:16)=atemp

        write(*,*)'reading ',restart_file_name

      open(21,file=restart_file_name,status='unknown') 

      read(21,705)string 
705   format(a80)

C possible to specify what moviepdb # to start at on first line

      if(string(45:45).eq.' ') then
        read(string,710)tot_time_orig               
      else
        read(21,761)tot_time_orig,movieframenumber   
      endif
      write(*,*)'we will start pdbmovie at frame # ',movieframenumber

      tot_time=tot_time_orig                          
      tot_time=transfer(tot_time_orig,tot_time)   ! eh?               
715   continue                                                       

      read(21,770,end=725)kii,t1,t2,t3

      if(i_f_tp(kii).ge.1) then

        if(i_debug) write(*,901)kii,num_atms_pt(i_t_tp(kii))
901     format('reading flexible mol ',i8,
     &            ' ; it appears to contain ',i8,
     &            ' pseudoatoms')

        do lii=1,num_atms_pt(i_t_tp(kii))                    

          read(21,730)t1,t2,t3

C not sure that we ever use c_f_m but okay we'll keep it

          c_f_m(kii)%x(lii)=t1
          c_f_m(kii)%y(lii)=t2
          c_f_m(kii)%z(lii)=t3

          t_a_x(1,id_beg(kii)+lii,ijk)=c_f_m(kii)%x(lii)
          t_a_x(2,id_beg(kii)+lii,ijk)=c_f_m(kii)%y(lii)
          t_a_x(3,id_beg(kii)+lii,ijk)=c_f_m(kii)%z(lii)

        enddo                                                  

      elseif(i_r_tp(kii).ge.1) then

        if(i_debug) write(*,902)kii
902     format('reading rigid    mol ',i8)

        c_m_x(kii)=t1
        c_m_y(kii)=t2
        c_m_z(kii)=t3

        read(21,730)M11(kii),M12(kii),M13(kii)
        read(21,730)M21(kii),M22(kii),M23(kii)
        read(21,730)M31(kii),M32(kii),M33(kii)

C note that we use fact that 1st atom is at 1,0,0 in mol frame while 2nd
C atom is at 0,1,0 in mol frame and 3rd atom is at 0,0,0
C especially remember that M11 etc is the rotation matrix not atoms!

        cx1(kii)=M11(kii)+c_m_x(kii)
        cy1(kii)=M21(kii)+c_m_y(kii)
        cz1(kii)=M31(kii)+c_m_z(kii)
        cx2(kii)=M12(kii)+c_m_x(kii)
        cy2(kii)=M22(kii)+c_m_y(kii)
        cz2(kii)=M32(kii)+c_m_z(kii)
        cx3(kii)=c_m_x(kii)
        cy3(kii)=c_m_y(kii)
        cz3(kii)=c_m_z(kii)

C probably better to just put the coords the t_a_x array...

        t_a_x(1,id_beg(kii)+1,ijk)=cx1(kii)
        t_a_x(2,id_beg(kii)+1,ijk)=cy1(kii)
        t_a_x(3,id_beg(kii)+1,ijk)=cz1(kii)
        t_a_x(1,id_beg(kii)+2,ijk)=cx2(kii)
        t_a_x(2,id_beg(kii)+2,ijk)=cy2(kii)
        t_a_x(3,id_beg(kii)+2,ijk)=cz2(kii)
        t_a_x(1,id_beg(kii)+3,ijk)=cx3(kii)
        t_a_x(2,id_beg(kii)+3,ijk)=cy3(kii)
        t_a_x(3,id_beg(kii)+3,ijk)=cz3(kii)

      endif

      goto 715                                                    
725   continue                                                   

      close(21)                                                 

      write(*,*)
      write(*,*)'done reading ',restart_file_temp
      write(*,*)

C close the loop over replicas

      enddo

710   format(7x,f20.5)                                         
761   format(7x,f20.5,8x,i10)                                         
770   format(i8,8x,3f20.5)                                    
730   format(8x,8x,3f20.5)                                   

C quit and go back                                                                          
                               
      return                    
      end                        
