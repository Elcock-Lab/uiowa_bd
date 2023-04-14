                                                                        
      subroutine read_wall_functions(num_walls)
                                      
      use allocatable_arrays         
      implicit real(a-h,o-z)        

      open(unit=11,file=wall_file,status='unknown')

C first, count up # of uniq atom types (probably already done this)

      num_uniqs=0
      do iii=1,num_t_ts
        do jjj=1,num_atms_pt(iii)
          num_uniqs=num_uniqs+1
        enddo
      enddo

C allocate array for the #walls per uniq atom type

      allocate(wall_ptr(1:num_uniqs))
      do nnn=1,num_uniqs
        wall_ptr(nnn)%num=0
      enddo

C now read wall_file once to get #walls per uniq atom type

      allocate(iyz(1:100))
      nnn=0
      do iii=1,num_t_ts
        do jjj=1,num_atms_pt(iii)
          nnn=nnn+1
          read(11,*)iil,jjl,(iyz(kkk),kkk=1,num_walls)
          if(iii.ne.iil.or.jjj.ne.jjl) then
            write(*,*)'walls assigned to atom types in wrong order'
            write(*,*)'quitting :('
            stop
          endif
          do kkk=1,num_walls
            if(iyz(kkk).ne.0) then
              wall_ptr(nnn)%num=wall_ptr(nnn)%num+1
            endif
          enddo    
        enddo
      enddo               

C now that we know how many wall types each uniq atom sees we can
C allocate memory and fill it in
C note that when we read the entries all we're looking for is a '1' to
C indicate that particular wall is required, e.g.:
C     2 15 0 0 1 0 1 1 
C indicates that the 15th atom of the 2nd molecule type sees walls 3,5,6

      do nnn=1,num_uniqs
        allocate(wall_ptr(nnn)%iwl(1:wall_ptr(nnn)%num))
        wall_ptr(nnn)%num=0 ! reset to zero prior to recomputing
      enddo
      rewind(11)
      nnn=0
      do iii=1,num_t_ts
        do jjj=1,num_atms_pt(iii)
          nnn=nnn+1
          read(11,*)iil,jjl,(iyz(kkk),kkk=1,num_walls)
          do kkk=1,num_walls
            if(iyz(kkk).ne.0) then
              wall_ptr(nnn)%num=wall_ptr(nnn)%num+1
              wall_ptr(nnn)%iwl(wall_ptr(nnn)%num)=kkk
            endif
          enddo    
        enddo
      enddo               
      close(11)     

C quit and go back                                                                          
                               
      return                    
      end                        
