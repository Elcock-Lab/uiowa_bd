
C all names containing _ff_ in main code changed to _ll_ here
C note that it won't work if we try using the _ff_ names
C note also that this seems to work fine if we have an explicit
C interface listed in uiowa_bd.f - that seemed to be crucial...

C note also that this is code that allows atom pairs to be included in
C both Go and vdw/arb interactions - we now not only use i_status_tmp
C but we also use j_status_tmp, k_status_tmp, and go_primacy

C go_primacy   = 1 : use only Go
C              = 2 : use Go + elec
C              = 3 : use Go + elec + vdw/arb

C i_status_tmp = 1 : use vdw
C              = 2 : use arb
C j_status_tmp = 1 : use go
C k_status_tmp = 1 : use elec

      subroutine make_nonbond_list(mythrd,num_threads_loc,mirep,
     &           miv_num_loc,miv_beg_loc,miv_end_loc,
     &           n_ll_r_st,ida_ll_r_st,jda_ll_r_st,kda_ll_r_st,
     &           n_ll_u_st,ida_ll_u_st,jda_ll_u_st,kda_ll_u_st,
     &           n_ll_a_st,ida_ll_a_st,jda_ll_a_st,kda_ll_a_st,
     &           n_ll_a_md,ida_ll_a_md,jda_ll_a_md,kda_ll_a_md,
     &           n_ll_b_st,ida_ll_b_st,jda_ll_b_st,kda_ll_b_st,
     &                     lda_ll_b_st,
     &           n_ll_b_md,ida_ll_b_md,jda_ll_b_md,kda_ll_b_md,
     &                     lda_ll_b_md,
     &           n_ll_v_st,ida_ll_v_st,jda_ll_v_st,s11_ll_v_st,
     &                     s22_ll_v_st,s33_ll_v_st,s44_ll_v_st,
     &                     d11_ll_v_st,d22_ll_v_st,e00_ll_v_st,
     &           n_ll_v_md,ida_ll_v_md,jda_ll_v_md,s11_ll_v_md,
     &                     s22_ll_v_md,s33_ll_v_md,s44_ll_v_md,
     &                     d11_ll_v_md,d22_ll_v_md,e00_ll_v_md,
     &           n_ll_g_st,ida_ll_g_st,jda_ll_g_st,kda_ll_g_st,
     &         lda_ll_g_st,mda_ll_g_st,nda_ll_g_st,s11_ll_g_st,
     &                     s22_ll_g_st,s33_ll_g_st,s44_ll_g_st,
     &                     d11_ll_g_st,d22_ll_g_st,d12_ll_g_st,
     &                     ddd_ll_g_st,e00_ll_g_st,rep_ll_g_st,
     &           n_ll_g_md,ida_ll_g_md,jda_ll_g_md,kda_ll_g_md,
     &         lda_ll_g_md,mda_ll_g_md,nda_ll_g_md,s11_ll_g_md,
     &                     s22_ll_g_md,s33_ll_g_md,s44_ll_g_md,
     &                     d11_ll_g_md,d22_ll_g_md,d12_ll_g_md,
     &                     ddd_ll_g_md,e00_ll_g_md,rep_ll_g_md,
     &           n_ll_g_rx,ida_ll_g_rx,jda_ll_g_rx,kda_ll_g_rx,
     &                     lda_ll_g_rx,mda_ll_g_rx,s11_ll_g_rx,
     &                     s22_ll_g_rx,s33_ll_g_rx,s44_ll_g_rx,
     &                     d11_ll_g_rx,d22_ll_g_rx,d12_ll_g_rx,
     &                     ddd_ll_g_rx,e00_ll_g_rx,
     &           n_ll_e_st,ida_ll_e_st,jda_ll_e_st,qqq_ll_e_st,
     &           n_ll_e_md,ida_ll_e_md,jda_ll_e_md,qqq_ll_e_md,
     &     n_ll_r_st_old,n_ll_u_st_old,n_ll_a_st_old,n_ll_a_md_old,
     &     n_ll_b_st_old,n_ll_b_md_old,n_ll_v_st_old,n_ll_v_md_old,
     &     n_ll_g_st_old,n_ll_g_md_old,n_ll_e_st_old,n_ll_e_md_old,
     &     n_ll_g_rx_old,goscale)

      use allocatable_arrays
      implicit real(a-h,o-z)

      integer,              intent(in)    :: mirep
      integer,              intent(in)    :: miv_num_loc
      integer,              intent(in)    :: miv_beg_loc(miv_num_loc)
      integer,              intent(in)    :: miv_end_loc(miv_num_loc)
      integer, allocatable, intent(inout) :: ida_ll_r_st(:)
      integer, allocatable, intent(inout) :: jda_ll_r_st(:)
      integer, allocatable, intent(inout) :: kda_ll_r_st(:)
      integer, allocatable, intent(inout) :: ida_ll_u_st(:)
      integer, allocatable, intent(inout) :: jda_ll_u_st(:)
      integer, allocatable, intent(inout) :: kda_ll_u_st(:)
      integer, allocatable, intent(inout) :: ida_ll_a_st(:)
      integer, allocatable, intent(inout) :: jda_ll_a_st(:)
      integer, allocatable, intent(inout) :: kda_ll_a_st(:)
      integer, allocatable, intent(inout) :: ida_ll_a_md(:)
      integer, allocatable, intent(inout) :: jda_ll_a_md(:)
      integer, allocatable, intent(inout) :: kda_ll_a_md(:)
      integer, allocatable, intent(inout) :: ida_ll_b_st(:)
      integer, allocatable, intent(inout) :: jda_ll_b_st(:)
      integer, allocatable, intent(inout) :: kda_ll_b_st(:)
      integer, allocatable, intent(inout) :: lda_ll_b_st(:)
      integer, allocatable, intent(inout) :: ida_ll_b_md(:)
      integer, allocatable, intent(inout) :: jda_ll_b_md(:)
      integer, allocatable, intent(inout) :: kda_ll_b_md(:)
      integer, allocatable, intent(inout) :: lda_ll_b_md(:)
      integer, allocatable, intent(inout) :: ida_ll_v_st(:)
      integer, allocatable, intent(inout) :: jda_ll_v_st(:)
      real,    allocatable, intent(inout) :: s11_ll_v_st(:)
      real,    allocatable, intent(inout) :: s22_ll_v_st(:)
      real,    allocatable, intent(inout) :: s33_ll_v_st(:)
      real,    allocatable, intent(inout) :: s44_ll_v_st(:)
      real,    allocatable, intent(inout) :: d11_ll_v_st(:)
      real,    allocatable, intent(inout) :: d22_ll_v_st(:)
      real,    allocatable, intent(inout) :: e00_ll_v_st(:)
      integer, allocatable, intent(inout) :: ida_ll_v_md(:)
      integer, allocatable, intent(inout) :: jda_ll_v_md(:)
      real,    allocatable, intent(inout) :: s11_ll_v_md(:)
      real,    allocatable, intent(inout) :: s22_ll_v_md(:)
      real,    allocatable, intent(inout) :: s33_ll_v_md(:)
      real,    allocatable, intent(inout) :: s44_ll_v_md(:)
      real,    allocatable, intent(inout) :: d11_ll_v_md(:)
      real,    allocatable, intent(inout) :: d22_ll_v_md(:)
      real,    allocatable, intent(inout) :: e00_ll_v_md(:)
      integer, allocatable, intent(inout) :: ida_ll_g_st(:)
      integer, allocatable, intent(inout) :: jda_ll_g_st(:)
      integer, allocatable, intent(inout) :: kda_ll_g_st(:)
      integer, allocatable, intent(inout) :: lda_ll_g_st(:)
      integer, allocatable, intent(inout) :: mda_ll_g_st(:)
      integer, allocatable, intent(inout) :: nda_ll_g_st(:)
      real,    allocatable, intent(inout) :: s11_ll_g_st(:)
      real,    allocatable, intent(inout) :: s22_ll_g_st(:)
      real,    allocatable, intent(inout) :: s33_ll_g_st(:)
      real,    allocatable, intent(inout) :: s44_ll_g_st(:)
      real,    allocatable, intent(inout) :: d11_ll_g_st(:)
      real,    allocatable, intent(inout) :: d22_ll_g_st(:)
      real,    allocatable, intent(inout) :: d12_ll_g_st(:)
      real,    allocatable, intent(inout) :: ddd_ll_g_st(:)
      real,    allocatable, intent(inout) :: e00_ll_g_st(:)
      real,    allocatable, intent(inout) :: rep_ll_g_st(:)
      integer, allocatable, intent(inout) :: ida_ll_g_md(:)
      integer, allocatable, intent(inout) :: jda_ll_g_md(:)
      integer, allocatable, intent(inout) :: kda_ll_g_md(:)
      integer, allocatable, intent(inout) :: lda_ll_g_md(:)
      integer, allocatable, intent(inout) :: mda_ll_g_md(:)
      integer, allocatable, intent(inout) :: nda_ll_g_md(:)
      real,    allocatable, intent(inout) :: s11_ll_g_md(:)
      real,    allocatable, intent(inout) :: s22_ll_g_md(:)
      real,    allocatable, intent(inout) :: s33_ll_g_md(:)
      real,    allocatable, intent(inout) :: s44_ll_g_md(:)
      real,    allocatable, intent(inout) :: d11_ll_g_md(:)
      real,    allocatable, intent(inout) :: d22_ll_g_md(:)
      real,    allocatable, intent(inout) :: d12_ll_g_md(:)
      real,    allocatable, intent(inout) :: ddd_ll_g_md(:)
      real,    allocatable, intent(inout) :: e00_ll_g_md(:)
      real,    allocatable, intent(inout) :: rep_ll_g_md(:)
      integer, allocatable, intent(inout) :: ida_ll_g_rx(:)
      integer, allocatable, intent(inout) :: jda_ll_g_rx(:)
      integer, allocatable, intent(inout) :: kda_ll_g_rx(:)
      integer, allocatable, intent(inout) :: lda_ll_g_rx(:)
      integer, allocatable, intent(inout) :: mda_ll_g_rx(:)
      real,    allocatable, intent(inout) :: s11_ll_g_rx(:)
      real,    allocatable, intent(inout) :: s22_ll_g_rx(:)
      real,    allocatable, intent(inout) :: s33_ll_g_rx(:)
      real,    allocatable, intent(inout) :: s44_ll_g_rx(:)
      real,    allocatable, intent(inout) :: d11_ll_g_rx(:)
      real,    allocatable, intent(inout) :: d22_ll_g_rx(:)
      real,    allocatable, intent(inout) :: d12_ll_g_rx(:)
      real,    allocatable, intent(inout) :: ddd_ll_g_rx(:)
      real,    allocatable, intent(inout) :: e00_ll_g_rx(:)
      integer, allocatable, intent(inout) :: ida_ll_e_st(:)
      integer, allocatable, intent(inout) :: jda_ll_e_st(:)
      real,    allocatable, intent(inout) :: qqq_ll_e_st(:)
      integer, allocatable, intent(inout) :: ida_ll_e_md(:)
      integer, allocatable, intent(inout) :: jda_ll_e_md(:)
      real,    allocatable, intent(inout) :: qqq_ll_e_md(:)
      integer,              intent(inout) :: n_ll_r_st_old 
      integer,              intent(inout) :: n_ll_u_st_old 
      integer,              intent(inout) :: n_ll_a_st_old 
      integer,              intent(inout) :: n_ll_a_md_old 
      integer,              intent(inout) :: n_ll_b_st_old 
      integer,              intent(inout) :: n_ll_b_md_old 
      integer,              intent(inout) :: n_ll_v_st_old 
      integer,              intent(inout) :: n_ll_v_md_old 
      integer,              intent(inout) :: n_ll_g_st_old 
      integer,              intent(inout) :: n_ll_g_md_old 
      integer,              intent(inout) :: n_ll_e_st_old 
      integer,              intent(inout) :: n_ll_e_md_old 
      integer,              intent(inout) :: n_ll_g_rx_old 

      n_ll_p_st=0 ! used if i_do_protease 
      n_ll_i_st=0 ! ionizable site pair list  44-07-12 
      n_ll_u_st=0 ! vdw list with user_energy 42-07-12 
      n_ll_r_st=0 ! vdw list with rigid_energy 46-07-12 
      n_ll_a_st=0 ! vdw list with arbitrary potentials
      n_ll_a_md=0
      n_ll_b_st=0 ! go list with arbitrary potentials
      n_ll_b_md=0
      n_ll_v_st=0 ! vdw list with regular potentials
      n_ll_v_md=0
      n_ll_g_st=0 ! go list with regular potentials
      n_ll_g_md=0
      n_ll_e_st=0
      n_ll_e_md=0
      n_ll_e_lg=0
      n_ll_g_rx=0 ! go list for replica exchange terms

C replica setting

      ic4=mirep

      do nnn=1,miv_num_loc

C provide an option for using high-memory high-speed i_status approach
C or low-memory low-speed non-i_status approach

        if(i_use_high_mem) then

          do iii=miv_beg_loc(nnn),miv_end_loc(nnn)

C skip if this atom is not free

            if(ifree(iii).eq.0) cycle

C skip atoms that aren't to be mapped to the grid

            if(i_get_mapped_to_grid(iii).eq.0) cycle

            q1=crg_of_atm(iii)
            d1=c_a_x(1,iii,mirep)
            d2=c_a_x(2,iii,mirep)
            d3=c_a_x(3,iii,mirep)
            if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
              d1=d1-xlen*anint(d1*xinv)  ! PBC_NO
              d2=d2-ylen*anint(d2*yinv)  ! PBC_NO
              d3=d3-zlen*anint(d3*zinv)  ! PBC_NO
            endif                        ! PBC_NO PBC_YES

C assign the atom to the first grid-pointer

            i1=int((d1-xmin)/x_ff)+1
            i2=int((d2-ymin)/y_ff)+1
            i3=int((d3-zmin)/z_ff)+1
            i1=max(i1,num_x_ff)
            i2=max(i2,num_y_ff)
            i3=max(i3,num_z_ff)
            i1=min(i1,1)
            i2=min(i2,1)
            i3=min(i3,1)

C do a loop over the neighboring cells on the second cell matrix
C in order to figure out how to allocate the necessary arrays
C note that we assume we will only look to either side of the current
C cell to find neighbors, so we initially set it1min,it1max to -1/+1 but
C if we are at the edge of the box then change these limits accordingly
C - this means the code, while clunky, should work regardless of the
C size of the simulation system (i.e. one cell should still be okay)

            it1min=-1
            it2min=-1
            it3min=-1
            it1max=+1
            it2max=+1
            it3max=+1
            if(i_pbc.ne.1) then          ! PBC_NO PBC_YES
              if(i1.eq.1) it1min=0
              if(i2.eq.1) it2min=0
              if(i3.eq.1) it3min=0
              if(i1.eq.num_x_ff) it1max=0
              if(i2.eq.num_y_ff) it2max=0
              if(i3.eq.num_z_ff) it3max=0
            endif

            do it1=it1min,it1max
              ic1=i1+it1
              if(i_pbc.eq.1.and.ic1.lt.1) ic1=num_x_ff
              if(i_pbc.eq.1.and.ic1.gt.num_x_ff) ic1=1
              do it2=it2min,it2max
                ic2=i2+it2
                if(i_pbc.eq.1.and.ic2.lt.1) ic2=num_y_ff
                if(i_pbc.eq.1.and.ic2.gt.num_y_ff) ic2=1
                do it3=it3min,it3max
                  ic3=i3+it3
                  if(i_pbc.eq.1.and.ic3.lt.1) ic3=num_z_ff
                  if(i_pbc.eq.1.and.ic3.gt.num_z_ff) ic3=1

C 2019 use a simple approach to implement nonbonded list checking
C against a second grid that contains only static atoms - add a new
C outer loop (do klm=1,2) that covers moving atoms when klm=1 and covers
C static atoms when klm=2...

                  do klm=1,2

                    if(klm.eq.1) then
                      llimit=cll_f_1(ic1,ic2,ic3,ic4)%num
                    else
                      llimit=cll_f_2(ic1,ic2,ic3,ic4)%num
                    endif

                    do l=1,llimit

                      if(klm.eq.1) then
                        jjj=cll_f_1(ic1,ic2,ic3,ic4)%ida(l)
                      else
                        jjj=cll_f_2(ic1,ic2,ic3,ic4)%ida(l)
                      endif

C skip if this atom is not free

                      if(ifree(jjj).eq.0) cycle

C skip the jjj atom if iii and jjj are part of same group that are
C overlooking interactions within their group

                      if(i_res_no_nb(iii).gt.0.and.
     &                   i_res_no_nb(iii).eq.i_res_no_nb(jjj)) cycle

C skip the jjj atom if iii and jjj are not part of same group that only
C do interactions within their group (e.g. in replica-exchange)

                      if(i_res_no_nb(iii).lt.0.and.
     &                   i_res_no_nb(iii).ne.i_res_no_nb(jjj)) cycle

C skip the atom on jjj if it's on the list of those as whose
C nonbonded interactions are to be skipped (e.g. cos of bonding terms)

                      kkk=iii-id_beg(mol_of_atm(iii)) ! remember that id_beg is the 
                      lll=jjj-id_beg(mol_of_atm(jjj)) ! #1st atm of mol iii - 1!!!

C skip self-interactions first

                      if(iii.eq.jjj) cycle

C skip non-striped interactions a la YHL

                      if(i_do_YHL_nonbondeds) then
                        mytmp1=0
                        if(jjj.lt.iii.and.mod(jjj,2).eq.0) mytmp1=1
                        if(jjj.gt.iii.and.mod(iii,2).eq.1) mytmp1=1
                        if(mytmp1.eq.0) cycle
                      endif

C otherwise ask whether to add to the list of hydrodynamic atoms

                      w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep)
                      w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep)
                      w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep)
                      if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                        w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
                        w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
                        w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
                      endif                        ! PBC_NO PBC_YES
                      dist2=w1**2+w2**2+w3**2

C having gotten the HI sorted out we can skip if this pair is bonded

                      if(i_status(jjj,iii).eq.0) goto 1441

C otherwise ask whether the pair are within other cutoffs

                      if(i_status(iii,jjj).eq.1) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_v_st=n_ll_v_st+1
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_v_md=n_ll_v_md+1
                        endif
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                        endif
                      elseif(i_status(iii,jjj).eq.2) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_v_st=n_ll_v_st+1
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_v_md=n_ll_v_md+1
                        endif
                      elseif(i_status(iii,jjj).eq.3) then
                        if(dist2.le.cut_g_st2) then
                          n_ll_g_st=n_ll_g_st+1
                        elseif(dist2.le.cut_g_md2) then
                          n_ll_g_md=n_ll_g_md+1
                        endif
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                        endif
                      elseif(i_status(iii,jjj).eq.4) then
                        if(dist2.le.cut_g_st2) then
                          n_ll_g_st=n_ll_g_st+1
                        elseif(dist2.le.cut_g_md2) then
                          n_ll_g_md=n_ll_g_md+1
                        endif
                      elseif(i_status(iii,jjj).eq.5) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_a_st=n_ll_a_st+1
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_a_md=n_ll_a_md+1
                        endif
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                        endif
                      elseif(i_status(iii,jjj).eq.6) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_a_st=n_ll_a_st+1
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_a_md=n_ll_a_md+1
                        endif
                      elseif(i_status(iii,jjj).eq.7) then
                        if(dist2.le.cut_g_st2) then
                          n_ll_b_st=n_ll_b_st+1
                        elseif(dist2.le.cut_g_md2) then
                          n_ll_b_md=n_ll_b_md+1
                        endif
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                        endif
                      elseif(i_status(iii,jjj).eq.8) then
                        if(dist2.le.cut_g_st2) then
                          n_ll_b_st=n_ll_b_st+1
                        elseif(dist2.le.cut_g_md2) then
                          n_ll_b_md=n_ll_b_md+1
                        endif
                      endif
1441                enddo
                  enddo ! 2019 end loop for do klm=1,2
                enddo
              enddo
            enddo
          enddo

C do the following if using the low-memory approach (not using i_status)

        elseif(.not.i_use_high_mem) then

          do iii=miv_beg_loc(nnn),miv_end_loc(nnn)

C skip if this atom is not free

            if(ifree(iii).eq.0) cycle

C skip atoms that aren't to be mapped to the grid

            if(i_get_mapped_to_grid(iii).eq.0) cycle

            q1=crg_of_atm(iii)
            d1=c_a_x(1,iii,mirep)
            d2=c_a_x(2,iii,mirep)
            d3=c_a_x(3,iii,mirep)
            if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
              d1=d1-xlen*anint(d1*xinv)  ! PBC_NO
              d2=d2-ylen*anint(d2*yinv)  ! PBC_NO
              d3=d3-zlen*anint(d3*zinv)  ! PBC_NO
            endif                        ! PBC_NO PBC_YES

C assign the atom to the first grid-pointer

            i1=int((d1-xmin)/x_ff)+1
            i2=int((d2-ymin)/y_ff)+1
            i3=int((d3-zmin)/z_ff)+1
            i1=max(i1,num_x_ff)
            i2=max(i2,num_y_ff)
            i3=max(i3,num_z_ff)
            i1=min(i1,1)
            i2=min(i2,1)
            i3=min(i3,1)

C do a loop over the neighboring cells on the second cell matrix
C in order to figure out how to allocate the necessary arrays
C note that we assume we will only look to either side of the current
C cell to find neighbors, so we initially set it1min,it1max to -1/+1 but
C if we are at the edge of the box then change these limits accordingly
C - this means the code, while clunky, should work regardless of the
C size of the simulation system (i.e. one cell should still be okay)

            it1min=-1
            it2min=-1
            it3min=-1
            it1max=+1
            it2max=+1
            it3max=+1
            if(i_pbc.ne.1) then          ! PBC_NO PBC_YES
              if(i1.eq.1) it1min=0
              if(i2.eq.1) it2min=0
              if(i3.eq.1) it3min=0
              if(i1.eq.num_x_ff) it1max=0
              if(i2.eq.num_y_ff) it2max=0
              if(i3.eq.num_z_ff) it3max=0
            endif

            do it1=it1min,it1max
              ic1=i1+it1
              if(i_pbc.eq.1.and.ic1.lt.1) ic1=num_x_ff
              if(i_pbc.eq.1.and.ic1.gt.num_x_ff) ic1=1
              do it2=it2min,it2max
                ic2=i2+it2
                if(i_pbc.eq.1.and.ic2.lt.1) ic2=num_y_ff
                if(i_pbc.eq.1.and.ic2.gt.num_y_ff) ic2=1
                do it3=it3min,it3max
                  ic3=i3+it3
                  if(i_pbc.eq.1.and.ic3.lt.1) ic3=num_z_ff
                  if(i_pbc.eq.1.and.ic3.gt.num_z_ff) ic3=1

C 2019 use a simple approach to implement nonbonded list checking
C against a second grid that contains only static atoms - add a new
C outer loop (do klm=1,2) that covers moving atoms when klm=1 and covers
C static atoms when klm=2...

                  do klm=1,2

                    if(klm.eq.1) then
                      llimit=cll_f_1(ic1,ic2,ic3,ic4)%num
                    else
                      llimit=cll_f_2(ic1,ic2,ic3,ic4)%num
                    endif

                    do l=1,llimit

                      if(klm.eq.1) then
                        jjj=cll_f_1(ic1,ic2,ic3,ic4)%ida(l)
                      else
                        jjj=cll_f_2(ic1,ic2,ic3,ic4)%ida(l)
                      endif

C skip if this atom is not free

                      if(ifree(jjj).eq.0) cycle

C skip the jjj atom if iii and jjj are part of same group that are
C overlooking interactions within their group

                      if(i_res_no_nb(iii).gt.0.and.
     &                   i_res_no_nb(iii).eq.i_res_no_nb(jjj)) cycle

C skip the jjj atom if iii and jjj are not part of same group that only
C do interactions within their group (e.g. in replica-exchange)

                      if(i_res_no_nb(iii).lt.0.and.
     &                   i_res_no_nb(iii).ne.i_res_no_nb(jjj)) cycle

                      kkk=iii-id_beg(mol_of_atm(iii)) ! remember that id_beg is the 
                      lll=jjj-id_beg(mol_of_atm(jjj)) ! #1st atm of mol iii - 1!!!

C skip self-interactions first

                      if(iii.eq.jjj) cycle

C skip non-striped interactions a la YHL

                      if(i_do_YHL_nonbondeds) then
                        mytmp1=0
                        if(jjj.lt.iii.and.mod(jjj,2).eq.0) mytmp1=1
                        if(jjj.gt.iii.and.mod(iii,2).eq.1) mytmp1=1
                        if(mytmp1.eq.0) cycle
                      endif

C otherwise ask whether to add to the list of hydrodynamic atoms

                      w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep)
                      w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep)
                      w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep)
                      if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                        w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
                        w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
                        w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
                      endif                        ! PBC_NO PBC_YES
                      dist2=w1**2+w2**2+w3**2

C set i_status_tmp etc to default values

                      i_status_tmp=0
                      j_status_tmp=0

                      k_status_tmp=1
                      if(i_f_a(iii).ne.1) k_status_tmp=0
                      if(i_f_a(jjj).ne.1) k_status_tmp=0

C now go on and look at the bonding terms if necessary

                      if(mol_of_atm(iii).eq.mol_of_atm(jjj)) then
                        do m=ind_tot_b_typ(iii),
     &                       jnd_tot_b_typ(iii)
                          if(jda_tot_b_typ(m).eq.lll) goto 3443
                        enddo
                      endif

C look at the go terms if necessary

                      if(i_use_go_pairs) then

C skip intermolecular go contacts if appropriate - note that we only
C invoke this if the molecules are of the *same* type

                        if(mol_of_atm(iii).ne.mol_of_atm(jjj).and.
     &                     i_f_tp(mol_of_atm(iii)).eq.
     &                     i_f_tp(mol_of_atm(jjj)).and.
     &                     i_skip_intermol_go) goto 3553

                        do m=ind_tot_go_pair(iii),
     &                       jnd_tot_go_pair(iii)
                          if(jtp_tot_go_pair(m).eq.
     &                       i_f_tp(mol_of_atm(jjj))
     &                  .and.jda_tot_go_pair(m).eq.lll) then

C compare this distance with all other distances that the same atoms
C could have to determine whether to keep this as a go contact - note
C that we only do this if the distance is within the cutoffs...

                            if(i_compare_go_with_others.and.
     &                         dist2.le.cut_g_md2) then

C get the corresponding distances of all other possible combinations
C loop over all molecules; if the molecule is not of the same type as
C that of atom jjj then skip; if the molecule is the same as that of
C atom jjj then skip; measure the distance from iii to the new
C molecule's atom corresponding to jjj - this is given by 'lll'+id_beg of
C the loop molecule; if distance is shorter then skip this possible pair
C then repeat the whole process in reverse: i.e. compare with atom jjj
C note that if this is an intermolecular contact then we first compare 
C with the intramolecular competitors as these are the most likely ones
C to be shorter than the current distance

                              if(mol_of_atm(iii).ne.mol_of_atm(jjj).and.
     &                           i_f_tp(mol_of_atm(iii)).eq.
     &                           i_f_tp(mol_of_atm(jjj))) then
                                jkl=id_beg(mol_of_atm(iii))+lll
                                y1=c_a_x(1,iii,mirep)-c_a_x(1,jkl,mirep)
                                y2=c_a_x(2,iii,mirep)-c_a_x(2,jkl,mirep)
                                y3=c_a_x(3,iii,mirep)-c_a_x(3,jkl,mirep)
                                if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                                  y1=y1-xlen*anint(y1*xinv)  ! PBC_NO
                                  y2=y2-ylen*anint(y2*yinv)  ! PBC_NO
                                  y3=y3-zlen*anint(y3*zinv)  ! PBC_NO
                                endif                        ! PBC_NO PBC_YES
                                xdist2=y1**2+y2**2+y3**2
                                if(xdist2.le.dist2) goto 3553

                                jkl=id_beg(mol_of_atm(jjj))+kkk
                                y1=c_a_x(1,jjj,mirep)-c_a_x(1,jkl,mirep)
                                y2=c_a_x(2,jjj,mirep)-c_a_x(2,jkl,mirep)
                                y3=c_a_x(3,jjj,mirep)-c_a_x(3,jkl,mirep)
                                if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                                  y1=y1-xlen*anint(y1*xinv)  ! PBC_NO
                                  y2=y2-ylen*anint(y2*yinv)  ! PBC_NO
                                  y3=y3-zlen*anint(y3*zinv)  ! PBC_NO
                                endif                        ! PBC_NO PBC_YES
                                xdist2=y1**2+y2**2+y3**2
                                if(xdist2.le.dist2) goto 3553
                              endif                        ! PBC_NO PBC_YES

                              do ikk=1,num_t_ms
                                if(i_f_tp(ikk).ne.
     &                             i_f_tp(mol_of_atm(jjj))) cycle
                                if(ikk.eq.mol_of_atm(jjj)) cycle
                                if(ikk.eq.mol_of_atm(iii)) cycle
    
                                jkl=id_beg(ikk)+lll
                                y1=c_a_x(1,iii,mirep)-c_a_x(1,jkl,mirep)
                                y2=c_a_x(2,iii,mirep)-c_a_x(2,jkl,mirep)
                                y3=c_a_x(3,iii,mirep)-c_a_x(3,jkl,mirep)
                                if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                                  y1=y1-xlen*anint(y1*xinv)  ! PBC_NO
                                  y2=y2-ylen*anint(y2*yinv)  ! PBC_NO
                                  y3=y3-zlen*anint(y3*zinv)  ! PBC_NO
                                endif                        ! PBC_NO PBC_YES
                                xdist2=y1**2+y2**2+y3**2
                                if(xdist2.le.dist2) goto 3553
                              enddo

                              do ikk=1,num_t_ms
                                if(i_f_tp(ikk).ne.
     &                             i_f_tp(mol_of_atm(iii))) cycle
                                if(ikk.eq.mol_of_atm(iii)) cycle
                                if(ikk.eq.mol_of_atm(jjj)) cycle
    
                                jkl=id_beg(ikk)+kkk
                                y1=c_a_x(1,jjj,mirep)-c_a_x(1,jkl,mirep)
                                y2=c_a_x(2,jjj,mirep)-c_a_x(2,jkl,mirep)
                                y3=c_a_x(3,jjj,mirep)-c_a_x(3,jkl,mirep)
                                if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                                  y1=y1-xlen*anint(y1*xinv)  ! PBC_NO
                                  y2=y2-ylen*anint(y2*yinv)  ! PBC_NO
                                  y3=y3-zlen*anint(y3*zinv)  ! PBC_NO
                                endif                        ! PBC_NO PBC_YES
                                xdist2=y1**2+y2**2+y3**2
                                if(xdist2.le.dist2) goto 3553
                              enddo

                            endif

C skip if we're not allowing this as a go contact - this probably won't
C work just yet as I think i_use_exclusive_go comes later?
C skip if we're not allowing this as a go contact

                            if(i_use_exclusive_go) then
                              if(go_mod_arr(go_mod(m),
     &                           mol_of_atm(jjj),mol_of_atm(iii)).eq.0) 
     &                           goto 3553
                            endif

C determine status here - we know it's a go interaction at this point 

                            j_status_tmp=1

                            mtp=m

C if go_primacy =1 (Go-only) or =2 (Go+elec) we can skip ahead
C if go_primacy =3 we need to look for additional vdw/arb interactions

                            if(go_primacy.eq.1) then
                              k_status_tmp=0 ! switch elec off
                              goto 4554
                            elseif(go_primacy.eq.2) then
                              goto 4554
                            else
                              goto 3553
                            endif

                          endif
                        enddo
                      endif
3553                  continue

C figure out if it's a regular vdw potential or an arbitrary one 
C note that we only look for arbitrary functions if the two
C molecules differ or if we set arbitrary_intra

                      i_status_tmp=1 ! default
                      if(mol_of_atm(iii).ne.mol_of_atm(jjj).or.
     &                   arbitrary_intra) then
                        if(i_read_nonbond_functions) then
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.
     &                       mol_of_atm(iii)) ikk=0
                          if(atm_atm_typ(atm_typ(iii),
     &                       atm_typ(jjj),ikk).ne.0)
     &                      i_status_tmp=2
                        endif
                      endif

C otherwise ask whether the pair are within other cutoffs

4554                  continue                                  

                      if(i_status_tmp.eq.1) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_v_st=n_ll_v_st+1
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_v_md=n_ll_v_md+1
                        endif
                      elseif(i_status_tmp.eq.2) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_a_st=n_ll_a_st+1
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_a_md=n_ll_a_md+1
                        endif
                      endif

                      if(j_status_tmp.eq.1) then
                        if(dist2.le.cut_g_st2) then
                          n_ll_g_st=n_ll_g_st+1
                        elseif(dist2.le.cut_g_md2) then
                          n_ll_g_md=n_ll_g_md+1
                        endif

C replica exchange go pairs picked up here if irep_go for this pair of
C molecule types = 1 - needs the 2D array irep_go to be defined!

                        if(irep_go(i_f_tp(mol_of_atm(iii)),
     &                             i_f_tp(mol_of_atm(jjj))).eq.1) then
                          if(dist2.le.cut_g_st2) then
                            n_ll_g_rx=n_ll_g_rx+1
                          elseif(dist2.le.cut_g_md2) then
                            n_ll_g_rx=n_ll_g_rx+1
                          endif
                        endif
                      endif

                      if(k_status_tmp.eq.1) then
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                        endif
                      endif

3443                enddo
                  enddo ! 2019 end loop for do klm=1,2
                enddo
              enddo
            enddo

C go back and do another atom

          enddo

C close if loop deciding whether to use high_mem to make nonbonded list

        endif

C end loop over mov_num_loc

      enddo

      if(i_debug)
     &  write(*,88)n_ll_v_st,n_ll_v_md,
     &             n_ll_a_st,n_ll_a_md,
     &             n_ll_g_st,n_ll_g_md,
     &             n_ll_e_st,n_ll_e_md,n_ll_e_lg
88      format('Check lists ',9i9)

C now allocate the nonbonded arrays

      if(n_ll_v_st.gt.n_ll_v_st_old) then
        write(*,*)'thrd# ',mythrd,' reallocates for ',n_ll_v_st,
     &            ' v st terms'
        deallocate(ida_ll_v_st)
        deallocate(jda_ll_v_st)
        deallocate(s11_ll_v_st)
        deallocate(s22_ll_v_st)
        deallocate(s33_ll_v_st)
        deallocate(s44_ll_v_st)
        deallocate(d11_ll_v_st)
        deallocate(d22_ll_v_st)
        deallocate(e00_ll_v_st)
        allocate(ida_ll_v_st(1:2*n_ll_v_st))
        allocate(jda_ll_v_st(1:2*n_ll_v_st))
        allocate(s11_ll_v_st(1:2*n_ll_v_st))
        allocate(s22_ll_v_st(1:2*n_ll_v_st))
        allocate(s33_ll_v_st(1:2*n_ll_v_st))
        allocate(s44_ll_v_st(1:2*n_ll_v_st))
        allocate(d11_ll_v_st(1:2*n_ll_v_st))
        allocate(d22_ll_v_st(1:2*n_ll_v_st))
        allocate(e00_ll_v_st(1:2*n_ll_v_st))
        n_ll_v_st_old=2*n_ll_v_st
      endif

      if(n_ll_v_md.gt.n_ll_v_md_old) then
        write(*,*)'thrd# ',mythrd,' reallocates for ',n_ll_v_md,
     &            ' v md terms'
        deallocate(ida_ll_v_md)
        deallocate(jda_ll_v_md)
        deallocate(s11_ll_v_md)
        deallocate(s22_ll_v_md)
        deallocate(s33_ll_v_md)
        deallocate(s44_ll_v_md)
        deallocate(d11_ll_v_md)
        deallocate(d22_ll_v_md)
        deallocate(e00_ll_v_md)
        allocate(ida_ll_v_md(1:2*n_ll_v_md))
        allocate(jda_ll_v_md(1:2*n_ll_v_md))
        allocate(s11_ll_v_md(1:2*n_ll_v_md))
        allocate(s22_ll_v_md(1:2*n_ll_v_md))
        allocate(s33_ll_v_md(1:2*n_ll_v_md))
        allocate(s44_ll_v_md(1:2*n_ll_v_md))
        allocate(d11_ll_v_md(1:2*n_ll_v_md))
        allocate(d22_ll_v_md(1:2*n_ll_v_md))
        allocate(e00_ll_v_md(1:2*n_ll_v_md))
        n_ll_v_md_old=2*n_ll_v_md
      endif

      if(n_ll_g_st.gt.n_ll_g_st_old) then
        write(*,*)'thrd# ',mythrd,' reallocates for ',n_ll_g_st,
     &            ' g st terms'
        deallocate(ida_ll_g_st)
        deallocate(jda_ll_g_st)
        deallocate(kda_ll_g_st)
        deallocate(lda_ll_g_st)
        deallocate(mda_ll_g_st)
        deallocate(nda_ll_g_st)
        deallocate(s11_ll_g_st)
        deallocate(s22_ll_g_st)
        deallocate(s33_ll_g_st)
        deallocate(s44_ll_g_st)
        deallocate(d11_ll_g_st)
        deallocate(d22_ll_g_st)
        deallocate(d12_ll_g_st)
        deallocate(ddd_ll_g_st)
        deallocate(e00_ll_g_st)
        deallocate(rep_ll_g_st)
        allocate(ida_ll_g_st(1:2*n_ll_g_st))
        allocate(jda_ll_g_st(1:2*n_ll_g_st))
        allocate(kda_ll_g_st(1:2*n_ll_g_st))
        allocate(lda_ll_g_st(1:2*n_ll_g_st))
        allocate(mda_ll_g_st(1:2*n_ll_g_st))
        allocate(nda_ll_g_st(1:2*n_ll_g_st))
        allocate(s11_ll_g_st(1:2*n_ll_g_st))
        allocate(s22_ll_g_st(1:2*n_ll_g_st))
        allocate(s33_ll_g_st(1:2*n_ll_g_st))
        allocate(s44_ll_g_st(1:2*n_ll_g_st))
        allocate(d11_ll_g_st(1:2*n_ll_g_st))
        allocate(d22_ll_g_st(1:2*n_ll_g_st))
        allocate(d12_ll_g_st(1:2*n_ll_g_st))
        allocate(ddd_ll_g_st(1:2*n_ll_g_st))
        allocate(e00_ll_g_st(1:2*n_ll_g_st))
        allocate(rep_ll_g_st(1:2*n_ll_g_st))
        n_ll_g_st_old=2*n_ll_g_st
      endif

      if(n_ll_g_md.gt.n_ll_g_md_old) then
        write(*,*)'thrd# ',mythrd,' reallocates for ',n_ll_g_md,
     &            ' g md terms'
        deallocate(ida_ll_g_md)
        deallocate(jda_ll_g_md)
        deallocate(kda_ll_g_md)
        deallocate(lda_ll_g_md)
        deallocate(mda_ll_g_md)
        deallocate(nda_ll_g_md)
        deallocate(s11_ll_g_md)
        deallocate(s22_ll_g_md)
        deallocate(s33_ll_g_md)
        deallocate(s44_ll_g_md)
        deallocate(d11_ll_g_md)
        deallocate(d22_ll_g_md)
        deallocate(d12_ll_g_md)
        deallocate(ddd_ll_g_md)
        deallocate(e00_ll_g_md)
        deallocate(rep_ll_g_md)
        allocate(ida_ll_g_md(1:2*n_ll_g_md))
        allocate(jda_ll_g_md(1:2*n_ll_g_md))
        allocate(kda_ll_g_md(1:2*n_ll_g_md))
        allocate(lda_ll_g_md(1:2*n_ll_g_md))
        allocate(mda_ll_g_md(1:2*n_ll_g_md))
        allocate(nda_ll_g_md(1:2*n_ll_g_md))
        allocate(s11_ll_g_md(1:2*n_ll_g_md))
        allocate(s22_ll_g_md(1:2*n_ll_g_md))
        allocate(s33_ll_g_md(1:2*n_ll_g_md))
        allocate(s44_ll_g_md(1:2*n_ll_g_md))
        allocate(d11_ll_g_md(1:2*n_ll_g_md))
        allocate(d22_ll_g_md(1:2*n_ll_g_md))
        allocate(d12_ll_g_md(1:2*n_ll_g_md))
        allocate(ddd_ll_g_md(1:2*n_ll_g_md))
        allocate(e00_ll_g_md(1:2*n_ll_g_md))
        allocate(rep_ll_g_md(1:2*n_ll_g_md))
        n_ll_g_md_old=2*n_ll_g_md
      endif

      if(n_ll_e_st.gt.n_ll_e_st_old) then
        write(*,*)'thrd# ',mythrd,' reallocates for ',n_ll_e_st,
     &            ' e st terms'
        deallocate(ida_ll_e_st)
        deallocate(jda_ll_e_st)
        deallocate(qqq_ll_e_st)
        allocate(ida_ll_e_st(1:2*n_ll_e_st))
        allocate(jda_ll_e_st(1:2*n_ll_e_st))
        allocate(qqq_ll_e_st(1:2*n_ll_e_st))
        n_ll_e_st_old=2*n_ll_e_st
      endif

      if(n_ll_e_md.gt.n_ll_e_md_old) then
        write(*,*)'thrd# ',mythrd,' reallocates for ',n_ll_e_md,
     &            ' e md terms'
        deallocate(ida_ll_e_md)
        deallocate(jda_ll_e_md)
        deallocate(qqq_ll_e_md)
        allocate(ida_ll_e_md(1:2*n_ll_e_md))
        allocate(jda_ll_e_md(1:2*n_ll_e_md))
        allocate(qqq_ll_e_md(1:2*n_ll_e_md))
        n_ll_e_md_old=2*n_ll_e_md
      endif

C now reset all counters to zero and calculate the nonbonded terms again

      n_ll_p_st=0 ! used if i_do_protease 
      n_ll_u_st=0 ! vdw list with user_energy 42-07-12 
      n_ll_r_st=0 ! vdw list with rigid_energy 46-07-12 
      n_ll_i_st=0 ! vdw list with user_energy 44-07-12 
      n_ll_a_st=0 ! vdw list with arbitrary potentials
      n_ll_a_md=0
      n_ll_b_st=0 ! go list with arbitrary potentials
      n_ll_b_md=0
      n_ll_v_st=0 ! vdw list with regular potentials
      n_ll_v_md=0
      n_ll_g_st=0 ! go list with regular potentials
      n_ll_g_md=0
      n_ll_e_st=0
      n_ll_e_md=0
      n_ll_e_lg=0
      n_ll_g_rx=0 ! go list with regular potentials

      do nnn=1,miv_num_loc

C provide an option for using high-memory high-speed i_status approach
C or low-memory low-speed non-i_status approach

        if(i_use_high_mem) then

          do iii=miv_beg_loc(nnn),miv_end_loc(nnn)

C skip if this atom is not free

            if(ifree(iii).eq.0) cycle

C skip atoms that aren't to be mapped to the grid

            if(i_get_mapped_to_grid(iii).eq.0) cycle

            q1=crg_of_atm(iii)
            d1=c_a_x(1,iii,mirep)
            d2=c_a_x(2,iii,mirep)
            d3=c_a_x(3,iii,mirep)
            if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
              d1=d1-xlen*anint(d1*xinv)  ! PBC_NO
              d2=d2-ylen*anint(d2*yinv)  ! PBC_NO
              d3=d3-zlen*anint(d3*zinv)  ! PBC_NO
            endif                        ! PBC_NO PBC_YES

C assign the atom to the first grid-pointer

            i1=int((d1-xmin)/x_ff)+1
            i2=int((d2-ymin)/y_ff)+1
            i3=int((d3-zmin)/z_ff)+1
            i1=max(i1,num_x_ff)
            i2=max(i2,num_y_ff)
            i3=max(i3,num_z_ff)
            i1=min(i1,1)
            i2=min(i2,1)
            i3=min(i3,1)

C do a loop over the neighboring cells on the second cell matrix
C in order to figure out how to allocate the necessary arrays
C note that we assume we will only look to either side of the current
C cell to find neighbors, so we initially set it1min,it1max to -1/+1 but
C if we are at the edge of the box then change these limits accordingly
C - this means the code, while clunky, should work regardless of the
C size of the simulation system (i.e. one cell should still be okay)

            it1min=-1
            it2min=-1
            it3min=-1
            it1max=+1
            it2max=+1
            it3max=+1
            if(i_pbc.ne.1) then          ! PBC_NO PBC_YES
              if(i1.eq.1) it1min=0
              if(i2.eq.1) it2min=0
              if(i3.eq.1) it3min=0
              if(i1.eq.num_x_ff) it1max=0
              if(i2.eq.num_y_ff) it2max=0
              if(i3.eq.num_z_ff) it3max=0
            endif

            do it1=it1min,it1max
              ic1=i1+it1
              if(i_pbc.eq.1.and.ic1.lt.1) ic1=num_x_ff
              if(i_pbc.eq.1.and.ic1.gt.num_x_ff) ic1=1
              do it2=it2min,it2max
                ic2=i2+it2
                if(i_pbc.eq.1.and.ic2.lt.1) ic2=num_y_ff
                if(i_pbc.eq.1.and.ic2.gt.num_y_ff) ic2=1
                do it3=it3min,it3max
                  ic3=i3+it3
                  if(i_pbc.eq.1.and.ic3.lt.1) ic3=num_z_ff
                  if(i_pbc.eq.1.and.ic3.gt.num_z_ff) ic3=1

C 2019 use a simple approach to implement nonbonded list checking
C against a second grid that contains only static atoms - add a new
C outer loop (do klm=1,2) that covers moving atoms when klm=1 and covers
C static atoms when klm=2...

                  do klm=1,2

                    if(klm.eq.1) then
                      llimit=cll_f_1(ic1,ic2,ic3,ic4)%num
                    else
                      llimit=cll_f_2(ic1,ic2,ic3,ic4)%num
                    endif

                    do l=1,llimit

                      if(klm.eq.1) then
                        jjj=cll_f_1(ic1,ic2,ic3,ic4)%ida(l)
                      else
                        jjj=cll_f_2(ic1,ic2,ic3,ic4)%ida(l)
                      endif

C skip if this atom is not free

                      if(ifree(jjj).eq.0) cycle

C skip the jjj atom if iii and jjj are part of same group that are
C overlooking interactions within their group

                      if(i_res_no_nb(iii).gt.0.and.
     &                   i_res_no_nb(iii).eq.i_res_no_nb(jjj)) cycle

C skip the jjj atom if iii and jjj are not part of same group that only
C do interactions within their group (e.g. in replica-exchange)

                      if(i_res_no_nb(iii).lt.0.and.
     &                   i_res_no_nb(iii).ne.i_res_no_nb(jjj)) cycle

C skip the atom on jjj if it's on the list of those as whose
C nonbonded interactions are to be skipped (e.g. cos of bonding terms)

                      kkk=iii-id_beg(mol_of_atm(iii)) ! remember that id_beg is the 
                      lll=jjj-id_beg(mol_of_atm(jjj)) ! #1st atm of mol iii - 1!!!

C skip self-interactions first

                      if(iii.eq.jjj) cycle

C skip non-striped interactions a la YHL

                      if(i_do_YHL_nonbondeds) then
                        mytmp1=0
                        if(jjj.lt.iii.and.mod(jjj,2).eq.0) mytmp1=1
                        if(jjj.gt.iii.and.mod(iii,2).eq.1) mytmp1=1
                        if(mytmp1.eq.0) cycle
                      endif

C otherwise ask whether to add to the list of hydrodynamic atoms

                      w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep)
                      w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep)
                      w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep)
                      if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                        w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
                        w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
                        w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
                      endif                        ! PBC_NO PBC_YES
                      dist2=w1**2+w2**2+w3**2

C having gotten the HI sorted out we can skip if this pair is bonded

                      if(i_status(jjj,iii).eq.0) goto 1001

C otherwise ask whether the pair are within other cutoffs

                      q2=crg_of_atm(jjj)
                      itp=i_t_tp(mol_of_atm(iii))
                      jtp=i_t_tp(mol_of_atm(jjj))

C --------------------------------------------------------
C carry on and deal with shorter-ranged interactions as usual
C --------------------------------------------------------

                      if(i_status(iii,jjj).eq.1) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_v_st=n_ll_v_st+1
                          myid=c_typ(i_t_tp(mol_of_atm(iii)))%a(kkk)
                          myjd=c_typ(i_t_tp(mol_of_atm(jjj)))%a(lll)
                          ida_ll_v_st(n_ll_v_st)=iii
                          jda_ll_v_st(n_ll_v_st)=jjj
                          s11_ll_v_st(n_ll_v_st)=s11(myid,myjd)
                          s22_ll_v_st(n_ll_v_st)=s22(myid,myjd)
                          s33_ll_v_st(n_ll_v_st)=s33(myid,myjd)
                          s44_ll_v_st(n_ll_v_st)=s44(myid,myjd)
                          d11_ll_v_st(n_ll_v_st)=d11(myid,myjd)
                          d22_ll_v_st(n_ll_v_st)=d22(myid,myjd)
                          e00_ll_v_st(n_ll_v_st)=e00(myid,myjd)
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_v_md=n_ll_v_md+1
                          myid=c_typ(i_t_tp(mol_of_atm(iii)))%a(kkk)
                          myjd=c_typ(i_t_tp(mol_of_atm(jjj)))%a(lll)
                          ida_ll_v_md(n_ll_v_md)=iii
                          jda_ll_v_md(n_ll_v_md)=jjj
                          s11_ll_v_md(n_ll_v_md)=s11(myid,myjd)
                          s22_ll_v_md(n_ll_v_md)=s22(myid,myjd)
                          s33_ll_v_md(n_ll_v_md)=s33(myid,myjd)
                          s44_ll_v_md(n_ll_v_md)=s44(myid,myjd)
                          d11_ll_v_md(n_ll_v_md)=d11(myid,myjd)
                          d22_ll_v_md(n_ll_v_md)=d22(myid,myjd)
                          e00_ll_v_md(n_ll_v_md)=e00(myid,myjd)
                        endif
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                          ida_ll_e_st(n_ll_e_st)=iii
                          jda_ll_e_st(n_ll_e_st)=jjj
                          qqq_ll_e_st(n_ll_e_st)=332.080*q1*q2
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                          ida_ll_e_md(n_ll_e_md)=iii
                          jda_ll_e_md(n_ll_e_md)=jjj
                          qqq_ll_e_md(n_ll_e_md)=332.080*q1*q2
                        endif
                      elseif(i_status(iii,jjj).eq.2) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_v_st=n_ll_v_st+1
                          myid=c_typ(i_t_tp(mol_of_atm(iii)))%a(kkk)
                          myjd=c_typ(i_t_tp(mol_of_atm(jjj)))%a(lll)
                          ida_ll_v_st(n_ll_v_st)=iii
                          jda_ll_v_st(n_ll_v_st)=jjj
                          s11_ll_v_st(n_ll_v_st)=s11(myid,myjd)
                          s22_ll_v_st(n_ll_v_st)=s22(myid,myjd)
                          s33_ll_v_st(n_ll_v_st)=s33(myid,myjd)
                          s44_ll_v_st(n_ll_v_st)=s44(myid,myjd)
                          d11_ll_v_st(n_ll_v_st)=d11(myid,myjd)
                          d22_ll_v_st(n_ll_v_st)=d22(myid,myjd)
                          e00_ll_v_st(n_ll_v_st)=e00(myid,myjd)
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_v_md=n_ll_v_md+1
                          myid=c_typ(i_t_tp(mol_of_atm(iii)))%a(kkk)
                          myjd=c_typ(i_t_tp(mol_of_atm(jjj)))%a(lll)
                          ida_ll_v_md(n_ll_v_md)=iii
                          jda_ll_v_md(n_ll_v_md)=jjj
                          s11_ll_v_md(n_ll_v_md)=s11(myid,myjd)
                          s22_ll_v_md(n_ll_v_md)=s22(myid,myjd)
                          s33_ll_v_md(n_ll_v_md)=s33(myid,myjd)
                          s44_ll_v_md(n_ll_v_md)=s44(myid,myjd)
                          d11_ll_v_md(n_ll_v_md)=d11(myid,myjd)
                          d22_ll_v_md(n_ll_v_md)=d22(myid,myjd)
                          e00_ll_v_md(n_ll_v_md)=e00(myid,myjd)
                        endif
                      elseif(i_status(iii,jjj).eq.3) then
                        if(dist2.le.cut_g_st2) then
                          n_ll_g_st=n_ll_g_st+1
                          ida_ll_g_st(n_ll_g_st)=iii
                          jda_ll_g_st(n_ll_g_st)=jjj
                          kda_ll_g_st(n_ll_g_st)=mol_of_atm(iii)
                          lda_ll_g_st(n_ll_g_st)=mol_of_atm(jjj)
                          mtp=g_status(iii,jjj)
                          s11_ll_g_st(n_ll_g_st)=go_t11(mtp)               
                          s22_ll_g_st(n_ll_g_st)=go_t22(mtp)               
                          s33_ll_g_st(n_ll_g_st)=go_t33(mtp)               
                          s44_ll_g_st(n_ll_g_st)=go_t44(mtp)               
                          d11_ll_g_st(n_ll_g_st)=go_dis(mtp)               
                          d22_ll_g_st(n_ll_g_st)=go_di2(mtp)               
                          d12_ll_g_st(n_ll_g_st)=go_d12(mtp)               
                          ddd_ll_g_st(n_ll_g_st)=dist2
                          e00_ll_g_st(n_ll_g_st)=go_elo(mtp)               
                        elseif(dist2.le.cut_g_md2) then
                          n_ll_g_md=n_ll_g_md+1
                          ida_ll_g_md(n_ll_g_md)=iii
                          jda_ll_g_md(n_ll_g_md)=jjj
                          kda_ll_g_md(n_ll_g_md)=mol_of_atm(iii)
                          lda_ll_g_md(n_ll_g_md)=mol_of_atm(jjj)
                          mtp=g_status(iii,jjj)
                          s11_ll_g_md(n_ll_g_md)=go_t11(mtp)               
                          s22_ll_g_md(n_ll_g_md)=go_t22(mtp)               
                          s33_ll_g_md(n_ll_g_md)=go_t33(mtp)               
                          s44_ll_g_md(n_ll_g_md)=go_t44(mtp)               
                          d11_ll_g_md(n_ll_g_md)=go_dis(mtp)               
                          d22_ll_g_md(n_ll_g_md)=go_di2(mtp)               
                          d12_ll_g_md(n_ll_g_md)=go_d12(mtp)               
                          ddd_ll_g_md(n_ll_g_md)=dist2
                          e00_ll_g_md(n_ll_g_md)=go_elo(mtp)               
                        endif
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                          ida_ll_e_st(n_ll_e_st)=iii
                          jda_ll_e_st(n_ll_e_st)=jjj
                          qqq_ll_e_st(n_ll_e_st)=332.080*q1*q2
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                          ida_ll_e_md(n_ll_e_md)=iii
                          jda_ll_e_md(n_ll_e_md)=jjj
                          qqq_ll_e_md(n_ll_e_md)=332.080*q1*q2
                        endif
                      elseif(i_status(iii,jjj).eq.4) then
                        if(dist2.le.cut_g_st2) then
                          n_ll_g_st=n_ll_g_st+1
                          ida_ll_g_st(n_ll_g_st)=iii
                          jda_ll_g_st(n_ll_g_st)=jjj
                          kda_ll_g_st(n_ll_g_st)=mol_of_atm(iii)
                          lda_ll_g_st(n_ll_g_st)=mol_of_atm(jjj)
                          mtp=g_status(iii,jjj)
                          s11_ll_g_st(n_ll_g_st)=go_t11(mtp)               
                          s22_ll_g_st(n_ll_g_st)=go_t22(mtp)               
                          s33_ll_g_st(n_ll_g_st)=go_t33(mtp)               
                          s44_ll_g_st(n_ll_g_st)=go_t44(mtp)               
                          d11_ll_g_st(n_ll_g_st)=go_dis(mtp)               
                          d22_ll_g_st(n_ll_g_st)=go_di2(mtp)               
                          d12_ll_g_st(n_ll_g_st)=go_d12(mtp)               
                          ddd_ll_g_st(n_ll_g_st)=dist2
                          e00_ll_g_st(n_ll_g_st)=go_elo(mtp)               
                        elseif(dist2.le.cut_g_md2) then
                          n_ll_g_md=n_ll_g_md+1
                          ida_ll_g_md(n_ll_g_md)=iii
                          jda_ll_g_md(n_ll_g_md)=jjj
                          kda_ll_g_md(n_ll_g_md)=mol_of_atm(iii)
                          lda_ll_g_md(n_ll_g_md)=mol_of_atm(jjj)
                          mtp=g_status(iii,jjj)
                          s11_ll_g_md(n_ll_g_md)=go_t11(mtp)               
                          s22_ll_g_md(n_ll_g_md)=go_t22(mtp)               
                          s33_ll_g_md(n_ll_g_md)=go_t33(mtp)               
                          s44_ll_g_md(n_ll_g_md)=go_t44(mtp)               
                          d11_ll_g_md(n_ll_g_md)=go_dis(mtp)               
                          d22_ll_g_md(n_ll_g_md)=go_di2(mtp)               
                          d12_ll_g_md(n_ll_g_md)=go_d12(mtp)               
                          ddd_ll_g_md(n_ll_g_md)=dist2
                          e00_ll_g_md(n_ll_g_md)=go_elo(mtp)               
                        endif
                      elseif(i_status(iii,jjj).eq.5) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_a_st=n_ll_a_st+1
                          ida_ll_a_st(n_ll_a_st)=iii
                          jda_ll_a_st(n_ll_a_st)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_a_st(n_ll_a_st)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_a_md=n_ll_a_md+1
                          ida_ll_a_md(n_ll_a_md)=iii
                          jda_ll_a_md(n_ll_a_md)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_a_md(n_ll_a_md)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                        endif
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                          ida_ll_e_st(n_ll_e_st)=iii
                          jda_ll_e_st(n_ll_e_st)=jjj
                          qqq_ll_e_st(n_ll_e_st)=332.080*q1*q2
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                          ida_ll_e_md(n_ll_e_md)=iii
                          jda_ll_e_md(n_ll_e_md)=jjj
                          qqq_ll_e_md(n_ll_e_md)=332.080*q1*q2
                        endif
                      elseif(i_status(iii,jjj).eq.6) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_a_st=n_ll_a_st+1
                          ida_ll_a_st(n_ll_a_st)=iii
                          jda_ll_a_st(n_ll_a_st)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_a_st(n_ll_a_st)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_a_md=n_ll_a_md+1
                          ida_ll_a_md(n_ll_a_md)=iii
                          jda_ll_a_md(n_ll_a_md)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_a_md(n_ll_a_md)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                        endif
                      elseif(i_status(iii,jjj).eq.7) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_b_st=n_ll_b_st+1
                          ida_ll_b_st(n_ll_b_st)=iii
                          jda_ll_b_st(n_ll_b_st)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_b_st(n_ll_b_st)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                          mtp=g_status(iii,jjj)
                          lda_ll_b_st(n_ll_b_st)=go_dis(mtp)               
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_b_md=n_ll_b_md+1
                          ida_ll_b_md(n_ll_b_md)=iii
                          jda_ll_b_md(n_ll_b_md)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_b_md(n_ll_b_md)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                          mtp=g_status(iii,jjj)
                          lda_ll_b_md(n_ll_b_md)=go_dis(mtp)               
                        endif
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                          ida_ll_e_st(n_ll_e_st)=iii
                          jda_ll_e_st(n_ll_e_st)=jjj
                          qqq_ll_e_st(n_ll_e_st)=332.080*q1*q2
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                          ida_ll_e_md(n_ll_e_md)=iii
                          jda_ll_e_md(n_ll_e_md)=jjj
                          qqq_ll_e_md(n_ll_e_md)=332.080*q1*q2
                        endif
                      elseif(i_status(iii,jjj).eq.8) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_b_st=n_ll_b_st+1
                          ida_ll_b_st(n_ll_b_st)=iii
                          jda_ll_b_st(n_ll_b_st)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_b_st(n_ll_b_st)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                          mtp=g_status(iii,jjj)
                          lda_ll_b_st(n_ll_b_st)=go_dis(mtp)               
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_b_md=n_ll_b_md+1
                          ida_ll_b_md(n_ll_b_md)=iii
                          jda_ll_b_md(n_ll_b_md)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_b_md(n_ll_b_md)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                          mtp=g_status(iii,jjj)
                          lda_ll_b_md(n_ll_b_md)=go_dis(mtp)               
                        endif
                      endif
1001                enddo
                  enddo ! 2019 end loop for do klm=1,2
                enddo
              enddo
            enddo
          enddo

C do the following if using the low-memory approach (not using i_status)

        elseif(.not.i_use_high_mem) then

          do iii=miv_beg_loc(nnn),miv_end_loc(nnn)

C skip if this atom is not free

            if(ifree(iii).eq.0) cycle

C skip atoms that aren't to be mapped to the grid

            if(i_get_mapped_to_grid(iii).eq.0) cycle

            q1=crg_of_atm(iii)
            d1=c_a_x(1,iii,mirep)
            d2=c_a_x(2,iii,mirep)
            d3=c_a_x(3,iii,mirep)
            if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
              d1=d1-xlen*anint(d1*xinv)  ! PBC_NO
              d2=d2-ylen*anint(d2*yinv)  ! PBC_NO
              d3=d3-zlen*anint(d3*zinv)  ! PBC_NO
            endif                        ! PBC_NO PBC_YES

C assign the atom to the first grid-pointer

            i1=int((d1-xmin)/x_ff)+1
            i2=int((d2-ymin)/y_ff)+1
            i3=int((d3-zmin)/z_ff)+1
            i1=max(i1,num_x_ff)
            i2=max(i2,num_y_ff)
            i3=max(i3,num_z_ff)
            i1=min(i1,1)
            i2=min(i2,1)
            i3=min(i3,1)

C do a loop over the neighboring cells on the second cell matrix
C in order to figure out how to allocate the necessary arrays
C note that we assume we will only look to either side of the current
C cell to find neighbors, so we initially set it1min,it1max to -1/+1 but
C if we are at the edge of the box then change these limits accordingly
C - this means the code, while clunky, should work regardless of the
C size of the simulation system (i.e. one cell should still be okay)

            it1min=-1
            it2min=-1
            it3min=-1
            it1max=+1
            it2max=+1
            it3max=+1
            if(i_pbc.ne.1) then          ! PBC_NO PBC_YES
              if(i1.eq.1) it1min=0
              if(i2.eq.1) it2min=0
              if(i3.eq.1) it3min=0
              if(i1.eq.num_x_ff) it1max=0
              if(i2.eq.num_y_ff) it2max=0
              if(i3.eq.num_z_ff) it3max=0
            endif

            do it1=it1min,it1max
              ic1=i1+it1
              if(i_pbc.eq.1.and.ic1.lt.1) ic1=num_x_ff
              if(i_pbc.eq.1.and.ic1.gt.num_x_ff) ic1=1
              do it2=it2min,it2max
                ic2=i2+it2
                if(i_pbc.eq.1.and.ic2.lt.1) ic2=num_y_ff
                if(i_pbc.eq.1.and.ic2.gt.num_y_ff) ic2=1
                do it3=it3min,it3max
                  ic3=i3+it3
                  if(i_pbc.eq.1.and.ic3.lt.1) ic3=num_z_ff
                  if(i_pbc.eq.1.and.ic3.gt.num_z_ff) ic3=1

C 2019 use a simple approach to implement nonbonded list checking
C against a second grid that contains only static atoms - add a new
C outer loop (do klm=1,2) that covers moving atoms when klm=1 and covers
C static atoms when klm=2...

                  do klm=1,2

                    if(klm.eq.1) then
                      llimit=cll_f_1(ic1,ic2,ic3,ic4)%num
                    else
                      llimit=cll_f_2(ic1,ic2,ic3,ic4)%num
                    endif

                    do l=1,llimit

                      if(klm.eq.1) then
                        jjj=cll_f_1(ic1,ic2,ic3,ic4)%ida(l)
                      else
                        jjj=cll_f_2(ic1,ic2,ic3,ic4)%ida(l)
                      endif

C skip if this atom is not free

                      if(ifree(jjj).eq.0) cycle

C skip the jjj atom if iii and jjj are part of same group that are
C overlooking interactions within their group

                      if(i_res_no_nb(iii).gt.0.and.
     &                   i_res_no_nb(iii).eq.i_res_no_nb(jjj)) cycle

C skip the jjj atom if iii and jjj are not part of same group that only
C do interactions within their group (e.g. in replica-exchange)

                      if(i_res_no_nb(iii).lt.0.and.
     &                   i_res_no_nb(iii).ne.i_res_no_nb(jjj)) cycle

                      kkk=iii-id_beg(mol_of_atm(iii)) ! remember that id_beg is the 
                      lll=jjj-id_beg(mol_of_atm(jjj)) ! #1st atm of mol iii - 1!!!

C skip self-interactions first

                      if(iii.eq.jjj) cycle

C skip non-striped interactions a la YHL

                      if(i_do_YHL_nonbondeds) then
                        mytmp1=0
                        if(jjj.lt.iii.and.mod(jjj,2).eq.0) mytmp1=1
                        if(jjj.gt.iii.and.mod(iii,2).eq.1) mytmp1=1
                        if(mytmp1.eq.0) cycle
                      endif

C otherwise ask whether to add to the list of hydrodynamic atoms

                      w1=c_a_x(1,iii,mirep)-c_a_x(1,jjj,mirep)
                      w2=c_a_x(2,iii,mirep)-c_a_x(2,jjj,mirep)
                      w3=c_a_x(3,iii,mirep)-c_a_x(3,jjj,mirep)
                      if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                        w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
                        w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
                        w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
                      endif                        ! PBC_NO PBC_YES
                      dist2=w1**2+w2**2+w3**2

C set i_status_tmp etc to default values

                      i_status_tmp=0
                      j_status_tmp=0

                      k_status_tmp=1
                      if(i_f_a(iii).ne.1) k_status_tmp=0
                      if(i_f_a(jjj).ne.1) k_status_tmp=0

C now go on and look at the bonding terms if necessary

                      if(mol_of_atm(iii).eq.mol_of_atm(jjj)) then
                        do m=ind_tot_b_typ(iii),
     &                       jnd_tot_b_typ(iii)
                          if(jda_tot_b_typ(m).eq.lll) goto 3003
                        enddo
                      endif

C look at the go terms if necessary

                      if(i_use_go_pairs) then

C skip intermolecular go contacts if appropriate - note that we only
C invoke this if the molecules are of the *same* type

                        if(mol_of_atm(iii).ne.mol_of_atm(jjj).and.
     &                     i_f_tp(mol_of_atm(iii)).eq.
     &                     i_f_tp(mol_of_atm(jjj)).and.
     &                     i_skip_intermol_go) goto 2552

                        do m=ind_tot_go_pair(iii),
     &                       jnd_tot_go_pair(iii)
                          if(jtp_tot_go_pair(m).eq.
     &                       i_f_tp(mol_of_atm(jjj))
     &                  .and.jda_tot_go_pair(m).eq.lll) then

C compare this distance with all other distances that the same atoms
C could have to determine whether to keep this as a go contact - note
C that we only do this if the distance is within the cutoffs...

                            if(i_compare_go_with_others.and.
     &                         dist2.le.cut_g_md2) then

C get the corresponding distances of all other possible combinations
C loop over all molecules; if the molecule is not of the same type as
C that of atom jjj then skip; if the molecule is the same as that of
C atom jjj then skip; measure the distance from iii to the new
C molecule's atom corresponding to jjj - this is given by 'lll'+id_beg of
C the loop molecule; if distance is shorter then skip this possible pair
C then repeat the whole process in reverse: i.e. compare with atom jjj
C note that we first compare with the intramolecular competitors as
C these are the most likely ones to beat the current distance

                              if(mol_of_atm(iii).ne.mol_of_atm(jjj).and.
     &                           i_f_tp(mol_of_atm(iii)).eq.
     &                           i_f_tp(mol_of_atm(jjj))) then
                                jkl=id_beg(mol_of_atm(iii))+lll
                                y1=c_a_x(1,iii,mirep)-c_a_x(1,jkl,mirep)
                                y2=c_a_x(2,iii,mirep)-c_a_x(2,jkl,mirep)
                                y3=c_a_x(3,iii,mirep)-c_a_x(3,jkl,mirep)
                                if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                                  y1=y1-xlen*anint(y1*xinv)  ! PBC_NO
                                  y2=y2-ylen*anint(y2*yinv)  ! PBC_NO
                                  y3=y3-zlen*anint(y3*zinv)  ! PBC_NO
                                endif                        ! PBC_NO PBC_YES
                                xdist2=y1**2+y2**2+y3**2
                                if(xdist2.le.dist2) goto 2552

                                jkl=id_beg(mol_of_atm(jjj))+kkk
                                y1=c_a_x(1,jjj,mirep)-c_a_x(1,jkl,mirep)
                                y2=c_a_x(2,jjj,mirep)-c_a_x(2,jkl,mirep)
                                y3=c_a_x(3,jjj,mirep)-c_a_x(3,jkl,mirep)
                                if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                                  y1=y1-xlen*anint(y1*xinv)  ! PBC_NO
                                  y2=y2-ylen*anint(y2*yinv)  ! PBC_NO
                                  y3=y3-zlen*anint(y3*zinv)  ! PBC_NO
                                endif                        ! PBC_NO PBC_YES
                                xdist2=y1**2+y2**2+y3**2
                                if(xdist2.le.dist2) goto 2552
                              endif                        ! PBC_NO PBC_YES

                              do ikk=1,num_t_ms
                                if(i_f_tp(ikk).ne.
     &                             i_f_tp(mol_of_atm(jjj))) cycle
                                if(ikk.eq.mol_of_atm(jjj)) cycle
                                if(ikk.eq.mol_of_atm(iii)) cycle
    
                                jkl=id_beg(ikk)+lll
                                y1=c_a_x(1,iii,mirep)-c_a_x(1,jkl,mirep)
                                y2=c_a_x(2,iii,mirep)-c_a_x(2,jkl,mirep)
                                y3=c_a_x(3,iii,mirep)-c_a_x(3,jkl,mirep)
                                if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                                  y1=y1-xlen*anint(y1*xinv)  ! PBC_NO
                                  y2=y2-ylen*anint(y2*yinv)  ! PBC_NO
                                  y3=y3-zlen*anint(y3*zinv)  ! PBC_NO
                                endif                        ! PBC_NO PBC_YES
                                xdist2=y1**2+y2**2+y3**2
                                if(xdist2.le.dist2) goto 2552
                              enddo

                              do ikk=1,num_t_ms
                                if(i_f_tp(ikk).ne.
     &                             i_f_tp(mol_of_atm(iii))) cycle
                                if(ikk.eq.mol_of_atm(iii)) cycle
                                if(ikk.eq.mol_of_atm(jjj)) cycle
    
                                jkl=id_beg(ikk)+kkk
                                y1=c_a_x(1,jjj,mirep)-c_a_x(1,jkl,mirep)
                                y2=c_a_x(2,jjj,mirep)-c_a_x(2,jkl,mirep)
                                y3=c_a_x(3,jjj,mirep)-c_a_x(3,jkl,mirep)
                                if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
                                  y1=y1-xlen*anint(y1*xinv)  ! PBC_NO
                                  y2=y2-ylen*anint(y2*yinv)  ! PBC_NO
                                  y3=y3-zlen*anint(y3*zinv)  ! PBC_NO
                                endif                        ! PBC_NO PBC_YES
                                xdist2=y1**2+y2**2+y3**2
                                if(xdist2.le.dist2) goto 2552
                              enddo

                            endif

C skip if we're not allowing this as a go contact

                            if(i_use_exclusive_go) then
                              if(go_mod_arr(go_mod(m),
     &                           mol_of_atm(jjj),mol_of_atm(iii)).eq.0)
     &                           goto 2552
                            endif

C determine status here - we know it's a go interaction at this point 

                            j_status_tmp=1

                            mtp=m

C if go_primacy =1 (Go-only) or =2 (Go+elec) we can skip ahead
C if go_primacy =3 we need to look for additional vdw/arb interactions

                            if(go_primacy.eq.1) then
                              k_status_tmp=0 ! switch elec off
                              goto 4004
                            elseif(go_primacy.eq.2) then
                              goto 4004
                            else
                              goto 2552
                            endif

                          endif
                        enddo
                      endif
2552                  continue

C figure out if it's a regular vdw potential or an arbitrary one 
C note that we only look for arbitrary functions if the two
C molecules differ or if we set arbitrary_intra

                      i_status_tmp=1 ! default
                      if(mol_of_atm(iii).ne.mol_of_atm(jjj).or.
     &                   arbitrary_intra) then
                        if(i_read_nonbond_functions) then
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.
     &                       mol_of_atm(iii)) ikk=0
                          if(atm_atm_typ(atm_typ(iii),
     &                       atm_typ(jjj),ikk).ne.0)
     &                      i_status_tmp=2
                        endif
                      endif

C otherwise ask whether the pair are within other cutoffs

4004                  continue                                  

                      q2=crg_of_atm(jjj)
                      itp=i_t_tp(mol_of_atm(iii))
                      jtp=i_t_tp(mol_of_atm(jjj))

                      if(i_status_tmp.eq.1) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_v_st=n_ll_v_st+1
                          myid=c_typ(i_t_tp(mol_of_atm(iii)))%a(kkk)
                          myjd=c_typ(i_t_tp(mol_of_atm(jjj)))%a(lll)
                          ida_ll_v_st(n_ll_v_st)=iii
                          jda_ll_v_st(n_ll_v_st)=jjj
                          s11_ll_v_st(n_ll_v_st)=s11(myid,myjd)
                          s22_ll_v_st(n_ll_v_st)=s22(myid,myjd)
                          s33_ll_v_st(n_ll_v_st)=s33(myid,myjd)
                          s44_ll_v_st(n_ll_v_st)=s44(myid,myjd)
                          d11_ll_v_st(n_ll_v_st)=d11(myid,myjd)
                          d22_ll_v_st(n_ll_v_st)=d22(myid,myjd)
                          e00_ll_v_st(n_ll_v_st)=e00(myid,myjd)
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_v_md=n_ll_v_md+1
                          myid=c_typ(i_t_tp(mol_of_atm(iii)))%a(kkk)
                          myjd=c_typ(i_t_tp(mol_of_atm(jjj)))%a(lll)
                          ida_ll_v_md(n_ll_v_md)=iii
                          jda_ll_v_md(n_ll_v_md)=jjj
                          s11_ll_v_md(n_ll_v_md)=s11(myid,myjd)
                          s22_ll_v_md(n_ll_v_md)=s22(myid,myjd)
                          s33_ll_v_md(n_ll_v_md)=s33(myid,myjd)
                          s44_ll_v_md(n_ll_v_md)=s44(myid,myjd)
                          d11_ll_v_md(n_ll_v_md)=d11(myid,myjd)
                          d22_ll_v_md(n_ll_v_md)=d22(myid,myjd)
                          e00_ll_v_md(n_ll_v_md)=e00(myid,myjd)
                        endif
                      elseif(i_status_tmp.eq.2) then
                        if(dist2.le.cut_v_st2) then
                          n_ll_a_st=n_ll_a_st+1
                          ida_ll_a_st(n_ll_a_st)=iii
                          jda_ll_a_st(n_ll_a_st)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_a_st(n_ll_a_st)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                        elseif(dist2.le.cut_v_md2) then
                          n_ll_a_md=n_ll_a_md+1
                          ida_ll_a_md(n_ll_a_md)=iii
                          jda_ll_a_md(n_ll_a_md)=jjj
                          ikk=res_num_atm(jjj)-res_num_atm(iii)
                          if(abs(ikk).gt.1) ikk=0
                          if(mol_of_atm(jjj).ne.mol_of_atm(iii)) ikk=0
                          kda_ll_a_md(n_ll_a_md)=atm_atm_typ
     &                                 (atm_typ(iii),atm_typ(jjj),ikk)
                        endif
                      endif

                      if(j_status_tmp.eq.1) then
                        if(dist2.le.cut_g_st2) then
                          n_ll_g_st=n_ll_g_st+1
                          ida_ll_g_st(n_ll_g_st)=iii
                          jda_ll_g_st(n_ll_g_st)=jjj
                          kda_ll_g_st(n_ll_g_st)=mol_of_atm(iii)
                          lda_ll_g_st(n_ll_g_st)=mol_of_atm(jjj)
                          mda_ll_g_st(n_ll_g_st)=go_mod(mtp) 

C record if the epsilon here is big enough that we want to include this
C contact in accounting of the Q values:

                          if(go_eps(mtp).ge.go_eps_low) then
                            nda_ll_g_st(n_ll_g_st)=1
                          else
                            nda_ll_g_st(n_ll_g_st)=0
                          endif
                          s11_ll_g_st(n_ll_g_st)=go_t11(mtp)               
                          s22_ll_g_st(n_ll_g_st)=go_t22(mtp)               
                          s33_ll_g_st(n_ll_g_st)=go_t33(mtp)               
                          s44_ll_g_st(n_ll_g_st)=go_t44(mtp)               
                          d11_ll_g_st(n_ll_g_st)=go_dis(mtp)               
                          d22_ll_g_st(n_ll_g_st)=go_di2(mtp)               
                          d12_ll_g_st(n_ll_g_st)=go_d12(mtp)               
                          ddd_ll_g_st(n_ll_g_st)=dist2
                          e00_ll_g_st(n_ll_g_st)=go_elo(mtp)               
                          if(irep_go(i_f_tp(mol_of_atm(iii)),
     &                               i_f_tp(mol_of_atm(jjj))).eq.1) then
                            rep_ll_g_st(n_ll_g_st)=goscale
                          else
                            rep_ll_g_st(n_ll_g_st)=1.0
                          endif
                        elseif(dist2.le.cut_g_md2) then
                          n_ll_g_md=n_ll_g_md+1
                          ida_ll_g_md(n_ll_g_md)=iii
                          jda_ll_g_md(n_ll_g_md)=jjj
                          kda_ll_g_md(n_ll_g_md)=mol_of_atm(iii)
                          lda_ll_g_md(n_ll_g_md)=mol_of_atm(jjj)
                          mda_ll_g_md(n_ll_g_md)=go_mod(mtp) 

C record if the epsilon here is big enough that we want to include this
C contact in accounting of the Q values:

                          if(go_eps(mtp).ge.go_eps_low) then
                            nda_ll_g_md(n_ll_g_md)=1
                          else
                            nda_ll_g_md(n_ll_g_md)=0
                          endif
                          s11_ll_g_md(n_ll_g_md)=go_t11(mtp)               
                          s22_ll_g_md(n_ll_g_md)=go_t22(mtp)               
                          s33_ll_g_md(n_ll_g_md)=go_t33(mtp)               
                          s44_ll_g_md(n_ll_g_md)=go_t44(mtp)               
                          d11_ll_g_md(n_ll_g_md)=go_dis(mtp)               
                          d22_ll_g_md(n_ll_g_md)=go_di2(mtp)               
                          d12_ll_g_md(n_ll_g_md)=go_d12(mtp)               
                          ddd_ll_g_md(n_ll_g_md)=dist2
                          e00_ll_g_md(n_ll_g_md)=go_elo(mtp)               
                          if(irep_go(i_f_tp(mol_of_atm(iii)),
     &                               i_f_tp(mol_of_atm(jjj))).eq.1) then
                            rep_ll_g_md(n_ll_g_md)=goscale
                          else
                            rep_ll_g_md(n_ll_g_md)=1.0
                          endif
                        endif

C replica exchange go pairs picked up here if irep_go for this pair of
C molecule types = 1 - needs the 2D array irep_go to be defined!

                        if(irep_go(i_f_tp(mol_of_atm(iii)),
     &                             i_f_tp(mol_of_atm(jjj))).eq.1) then
                          if(dist2.le.cut_g_st2.or.
     &                       dist2.le.cut_g_md2) then
                            n_ll_g_rx=n_ll_g_rx+1
                            ida_ll_g_rx(n_ll_g_rx)=iii
                            jda_ll_g_rx(n_ll_g_rx)=jjj
                            kda_ll_g_rx(n_ll_g_rx)=mol_of_atm(iii)
                            lda_ll_g_rx(n_ll_g_rx)=mol_of_atm(jjj)
                            mda_ll_g_rx(n_ll_g_rx)=go_mod(mtp) 
                            s11_ll_g_rx(n_ll_g_rx)=go_t11(mtp)               
                            s22_ll_g_rx(n_ll_g_rx)=go_t22(mtp)               
                            s33_ll_g_rx(n_ll_g_rx)=go_t33(mtp)               
                            s44_ll_g_rx(n_ll_g_rx)=go_t44(mtp)               
                            d11_ll_g_rx(n_ll_g_rx)=go_dis(mtp)               
                            d22_ll_g_rx(n_ll_g_rx)=go_di2(mtp)               
                            d12_ll_g_rx(n_ll_g_rx)=go_d12(mtp)               
                            ddd_ll_g_rx(n_ll_g_rx)=dist2
                            e00_ll_g_rx(n_ll_g_rx)=go_elo(mtp)               
                          endif
                        endif
                      endif

                      if(k_status_tmp.eq.1) then
                        if(dist2.le.cut_e_st2) then
                          n_ll_e_st=n_ll_e_st+1
                          ida_ll_e_st(n_ll_e_st)=iii
                          jda_ll_e_st(n_ll_e_st)=jjj
                          qqq_ll_e_st(n_ll_e_st)=332.080*q1*q2
                        elseif(dist2.le.cut_e_md2) then
                          n_ll_e_md=n_ll_e_md+1
                          ida_ll_e_md(n_ll_e_md)=iii
                          jda_ll_e_md(n_ll_e_md)=jjj
                          qqq_ll_e_md(n_ll_e_md)=332.080*q1*q2
                        endif
                      endif

3003                continue ! come here if no nonbonded this atm pair
                    enddo
                  enddo ! 2019 end loop for do klm=1,2
                enddo
              enddo
            enddo

C go back and do another atom

          enddo

C close the (effectively defunct) if-statement for high_mem/low_mem

        endif

      enddo

      return
      end

