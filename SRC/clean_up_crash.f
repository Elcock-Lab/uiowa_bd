
      subroutine clean_up_crash(mirep)
                                      
      use allocatable_arrays         
      implicit real(a-h,o-z)        

      if(.not.i_look_for_crashes) then
        call write_movie_file(movieframenumber,icrashed)
      else

        c_a_x_fail=c_a_x 

        tot_time=tot_time-10.0*time_step
        do j=1,10
          do i=1,num_flx_atms
            c_a_x(1,i,mirep)=c_c_x(1,i,j,mirep) 
            c_a_x(2,i,mirep)=c_c_x(2,i,j,mirep) 
            c_a_x(3,i,mirep)=c_c_x(3,i,j,mirep) 
          enddo
          tot_time=tot_time+time_step
          call write_movie_file(movieframenumber,icrashed,mirep)
          movieframenumber=movieframenumber+1
        enddo

        c_a_x=c_a_x_fail
        call write_movie_file(movieframenumber,icrashed,mirep)

      endif

      return
      end
