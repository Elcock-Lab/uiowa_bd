
      subroutine initialize_treecode_elec
                                                                
      use allocatable_arrays                                   
      implicit real(a-h,o-z)                                  

C for treecode (and other uses?) get the total number of charges in the
C system and assign information accordingly
C i_q_a stores, for each atom, its id in the list of charged atoms
C note that atoms only make it into this list if they have a charge and
C if they are free (i.e. are not subject to a 'grow' simulation)

      write(*,*)' initializing treecode_elec here '
      write(*,*)' note that we have completed make_charge_parameters '

      if(.not.allocated(q_info)) then
        write(*,*)' allocating charges (first time) '
      else
        write(*,*)' allocating charges (again) '
        deallocate(q_info)
      endif
      allocate(q_info(1:num_q_as))

      if(.not.allocated(p_info)) then
        write(*,*)' allocating protonatable sites (first time) '
      else
        write(*,*)' allocating protonatable (again) '
        deallocate(p_info)
      endif
      allocate(p_info(1:num_p_as))

      write(*,*)'There are ',num_q_as,' free-to-move charged atoms'
      write(*,*)'There are ',num_p_as,' free-to-move protonatable atoms'

      num_q_as=0
      num_p_as=0
      do iii=1,num_flx_atms

        if(i_q_a(iii).ne.0) then
          num_q_as=num_q_as+1
          q_info(num_q_as)%ida=iii ! store id of atom

C if found a charged atom, look for other bonded charges - first time
C through just find out how many there are

          num_bnd=0
          do jjj=ind_tot_b_typ(iii),
     &           jnd_tot_b_typ(iii)

C next line is tricky - need to make sure we add on id_beg here

            kkk=id_beg(mol_of_atm(iii))+jda_tot_b_typ(jjj)
            if(i_q_a(kkk).ne.0) num_bnd=num_bnd+1 
          enddo

C if there are bonded charges, then store their charge-ids - we can
C allocate memory for them and read them again now

          q_info(num_q_as)%num_bnd=num_bnd
          if(num_bnd.gt.0) then
            allocate(q_info(num_q_as)%jda(1:num_bnd))
            num_bnd=0
            do jjj=ind_tot_b_typ(iii),
     &             jnd_tot_b_typ(iii)

C next line is tricky: kkk is the atom number associated with the bond
C exclusions listed for iii - note that we have to add on the beginning
C atom number for the entire molecule for 1<kkk<num_flx_atms

              kkk=id_beg(mol_of_atm(iii))+jda_tot_b_typ(jjj)

              if(i_q_a(kkk).ne.0) then

C 2022 - before adding to the list we need to account for the fact that
C the bond exclusion list might have redundant entries - e.g. atoms 1
C and 3 might appear twice if they have an angle *and* a direct bond
C added between them, so before incrementing num_bnd, check that we
C haven't already stored this value of i_q_a(kkk)...

                iamnew=1
                do lll=1,num_bnd
                  if(q_info(num_q_as)%jda(lll).eq.i_q_a(kkk)) then
                    iamnew=0
                    exit
                  endif
                enddo
                if(iamnew.eq.1) then
                  num_bnd=num_bnd+1
                  q_info(num_q_as)%jda(num_bnd)=i_q_a(kkk)
                endif

              endif
            enddo
          endif

C 2022 - need to make sure %num_bnd is updated if some are skipped...

          q_info(num_q_as)%num_bnd=num_bnd

        endif

      enddo

C write out info to screen just to be sure

      if(i_debug) then
        do iii=1,num_flx_atms
          if(.not.i_do_pH.and.i_q_a(iii).ne.0)
     &    write(*,*)'charged & protonatable atoms are ',iii,i_q_a(iii)
        enddo
        do iii=1,num_p_as
          write(*,*)'protonatable sites are ',iii,p_info(iii)%ida,
     &                                            p_info(iii)%idq
        enddo
        do iii=1,num_q_as
          if(q_info(iii)%num_bnd.gt.0) then
            write(*,490)iii,q_info(iii)%ida,
     &               (q_info(iii)%jda(jjj),jjj=1,q_info(iii)%num_bnd)
490         format('charge # ',i8,' is atom # ',i8,' is bonded to ',
     &             'charge #s ',10i8)
          else
            write(*,491)iii,q_info(iii)%ida
491         format('charge # ',i8,' is atom # ',i8,' is bonded to ',
     &             'nothing ')
          endif
        enddo
      endif

C now allocate shared arrays that depend on num_q_as

      if(allocated(c_q_x_tar)) then
        deallocate(c_q_x_tar)
        deallocate(c_q_y_tar)
        deallocate(c_q_z_tar)
        deallocate(spengtar)
        deallocate(tpengtar)
        deallocate(upengtar)
        deallocate(vpengtar)
        deallocate(tforce)
        deallocate(uforce)
        deallocate(vforce)
      endif

      allocate(c_q_x_tar(1:num_q_as))  
      allocate(c_q_y_tar(1:num_q_as))  
      allocate(c_q_z_tar(1:num_q_as))  
      allocate(spengtar(1:num_q_as))   ! energy
      allocate(tpengtar(1:num_q_as))   ! energy
      allocate(upengtar(1:num_q_as))   ! energy
      allocate(vpengtar(1:num_q_as))   ! energy
      allocate(tforce(1:3,1:num_q_as)) ! force 
      allocate(uforce(1:3,1:num_q_as)) ! force 
      allocate(vforce(1:3,1:num_q_as)) ! force 

      write(*,*)'exiting initialize_treecode_elec'

      return
      end
      


