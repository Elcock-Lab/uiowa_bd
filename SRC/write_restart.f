
      subroutine write_restart(movieframe,mirep)
                                      
      use allocatable_arrays         
      implicit real(a-h,o-z)        
      character*16 restart_file_name
      character*3  atemp

C this block is for writing out back up restart.files

      if ((mod(movieframe-1,10)).eq.0.and.nstep.ne.0) then                 

        call WrtRstrtBk(movieframe-1,mirep)
73000   continue

        restart_file_name='restart.file.XXX'
        write(atemp,'(i3)')mirep
        do j=1,3
          if(atemp(j:j).eq.' ') atemp(j:j)='0'
        enddo
        restart_file_name(14:16)=atemp

        write(*,*)'opening up restart ',restart_file_name

        open(unit=101+mirep,file=restart_file_name,
     &       status='unknown',err=73000) 
        write(101+mirep,10)tot_time             
        do kkk=1,num_t_ms                
c         write(101+mirep,20)kkk,c_m_x(kkk),c_m_y(kkk),c_m_z(kkk) 
          write(101+mirep,20)kkk,0.0,0.0,0.0
          do lll=1,num_atms_pt(i_f_tp(kkk))               
731     format(2i8,3f20.5)
10      format('time = ',f20.5,' ps')
20      format(i8,'       0',3f20.5)
734       format(2i8,3f20.5,i8)     
            write(101+mirep,731)kkk,lll,     
     &        c_a_x(1,id_beg(kkk)+lll,mirep),  
     &        c_a_x(2,id_beg(kkk)+lll,mirep), 
     &        c_a_x(3,id_beg(kkk)+lll,mirep) 
          enddo                         
        enddo                          
        close(101+mirep)                     
      endif

C this block is for overwriting restart.file

      if((mod(movieframe,1)).eq.0) then                 

        restart_file_name='restart.file.XXX'
        write(atemp,'(i3)')mirep
        do j=1,3
          if(atemp(j:j).eq.' ') atemp(j:j)='0'
        enddo
        restart_file_name(14:16)=atemp

73001   continue
        open(unit=101+mirep,file=restart_file_name,
     &       status='unknown',err=73001)      
        write(101+mirep,10)tot_time                 
        do kkk=1,num_t_ms                   
c         write(101+mirep,20)kkk,c_m_x(kkk),c_m_y(kkk),c_m_z(kkk) 
          write(101+mirep,20)kkk,0.0,0.0,0.0
          do lll=1,num_atms_pt(i_f_tp(kkk))
              write(101+mirep,731)kkk,lll,
     &          c_a_x(1,id_beg(kkk)+lll,mirep),
     &          c_a_x(2,id_beg(kkk)+lll,mirep),
     &          c_a_x(3,id_beg(kkk)+lll,mirep)
          enddo
        enddo                

        close(101+mirep)           
      endif                                                    

      return
      end
