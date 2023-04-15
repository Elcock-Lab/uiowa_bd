
      subroutine calculate_go_pairs

      use allocatable_arrays
      implicit real(a-h,o-z)

C note that this routine is currently called only by one thread in the
C entire run - we *could* call it from the master thread of each rep

      num_tot_intermol=0
      num_tot_intermod=0
      num_cur_Q=0.0 ! set all sums to zero

      do iii=1,num_g_ms
        do jjj=1,num_g_ms

C first analyze on a molecule-only basis

          if(num_req_Q(iii,jjj).ge.1) then
            do itp=1,num_threads
              kkk=mythrds_rep_curr(itp)
              num_cur_Q(iii,jjj,kkk)=
     &        num_cur_Q(iii,jjj,kkk)+
     &        num_cur_Q_loc(iii,jjj,itp)
            enddo
            do kkk=1,num_reps
              if(num_cur_Q(iii,jjj,kkk).gt.0) then
                qtmp=real(num_cur_Q(iii,jjj,kkk))/
     &               real(num_req_Q(iii,jjj))

C note that we only write out to the screen if jjj>=iii - there is no
C point writing out the same information twice - but we need to do the
C loop over jjj=1,num_g_ms in order to pick up i_do_NAM cases where the
C diffusing molecule has a higher moltype # than the target molecule

                if(jjj.ge.iii) then

                  num_tot_intermol=num_tot_intermol+1
                  write(*,82)tot_time,kkk,iii,jjj,
     &                       num_cur_Q(iii,jjj,kkk),
     &                       num_req_Q(iii,jjj),qtmp,
     &                       rep_go_scale(kkk)
                  ibin=int(qtmp*40.0)+1
                  if(ibin.gt.40) ibin=40
                  hist_of_Q(iii,jjj,ibin,kkk)=
     &            hist_of_Q(iii,jjj,ibin,kkk)+1
82      format('At time = ',f20.3,' replica# ',i6,
     &         ' Go mol  contacts ',4i10,' Q = ',f10.5,
     &         ' rep_go_scale = ',f8.5)

                endif

C check for termination (this code originally used earlier)

                if(iii.eq.mol_Q1.and.jjj.eq.mol_Q2) rat=qtmp

C 2023 v1.1 add support for unfolding simulations using fold_mode

                if(fold_mode.eq.1) then
                  if(rat.ge.Q_des) then
                    write(*,87)mol_Q1,mol_Q2,rat,tot_time
87                  format('for mols ',i5,' & ',i5,' Q = ',f10.5,
     &                     ' at time = ',f20.5,' FOLDED/ASSOCIATED')
                    call write_movie_file(movieframenumber,icrashed,kkk)
                    call write_restart(movieframenumber,kkk)
                    stop
                  endif
                elseif(fold_mode.eq.-1) then
                  if(rat.le.Q_des) then
                    write(*,88)mol_Q1,mol_Q2,rat,tot_time
88                  format('for mols ',i5,' & ',i5,' Q = ',f10.5,
     &                     ' at time = ',f20.5,' UNFOLDED/DISSOCIATED')
                    call write_movie_file(movieframenumber,icrashed,kkk)
                    call write_restart(movieframenumber,kkk)
                    stop
                  endif
                endif

              endif
            enddo
          endif
        enddo
      enddo

C write out the histogram of Q values

      open(unit=59,file='hist_of_Q.txt',status='unknown')

c     do mmm=1,num_reps
        do iii=1,num_g_ms
          do kkk=1,40
            rlo=real(kkk-1)*0.025
            rhi=real(kkk)*0.025
            rav=(rlo+rhi)/2.0
            write(59,911)iii,rav,
     &                (jjj,(hist_of_Q(iii,jjj,kkk,mmm),
     &                      jjj=iii,num_g_ms),mmm=1,num_reps)
          enddo
        enddo
c     enddo

911   format('for mol1 ',i4,' Q = ',f8.4,' with mol2 ',200(i4,i9,' * '))
      close(59)

      return
      end
