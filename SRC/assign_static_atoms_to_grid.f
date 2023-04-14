
      subroutine assign_static_atoms_to_grid(ndum1,mithrd,mirep,
     &           i1_lic,i2_lic,i3_lic)

      use allocatable_arrays
      implicit real(a-h,o-z)

      integer ndum1 ! how many atoms
      integer mirep ! myrep
      integer i1_lic(ndum1),i2_lic(ndum1),i3_lic(ndum1)

      i4=mirep
      cll_f_2(:,:,:,i4)%num=0

C do all current atom positions first
            
c     do iii=1,num_tot_atms
      do iii=num_mov_atms+1,num_tot_atms ! 2019 - only use static atoms

C skip atoms that aren't to be mapped to the grid

        if(i_get_mapped_to_grid(iii).eq.0) cycle

        w1=c_a_x(1,iii,mirep)
        w2=c_a_x(2,iii,mirep)
        w3=c_a_x(3,iii,mirep)

C if i_pbc=1 then wrap coordinates, otherwise check that the coordinates
C aren't out of the simulation box...

        if(i_pbc.eq.1) then          ! PBC_NO PBC_YES
          w1=w1-xlen*anint(w1*xinv)  ! PBC_NO
          w2=w2-ylen*anint(w2*yinv)  ! PBC_NO
          w3=w3-zlen*anint(w3*zinv)  ! PBC_NO
        else
          if(w1.lt.xmin.or.w1.gt.xmax.or.
     &       w2.lt.ymin.or.w2.gt.ymax.or.
     &       w3.lt.zmin.or.w3.gt.zmax) then
            write(*,911)iii,w1,w2,w3
911         format('atom# ',i8,' has coords ',3f10.5,
     &            ' these are out of the box & i_pbc=0; this is fatal')
            stop
          endif
        endif                        ! PBC_NO PBC_YES

C assign the atom to the first grid-pointer

        i1=int((w1-xmin)/x_ff)+1
        i2=int((w2-ymin)/y_ff)+1
        i3=int((w3-zmin)/z_ff)+1

C make sure that the limits are correctly applied:

        i1=max(i1,num_x_ff)
        i2=max(i2,num_y_ff)
        i3=max(i3,num_z_ff)
        i1=min(i1,1)
        i2=min(i2,1)
        i3=min(i3,1)

        i1_lic(iii)=i1
        i2_lic(iii)=i2
        i3_lic(iii)=i3
        cll_f_2(i1,i2,i3,i4)%num=cll_f_2(i1,i2,i3,i4)%num+1

      enddo

      do i3=1,num_z_ff
        do i2=1,num_y_ff
          do i1=1,num_x_ff
            deallocate(cll_f_2(i1,i2,i3,i4)%ida)
          enddo
        enddo
      enddo

      do i3=1,num_z_ff
        do i2=1,num_y_ff
          do i1=1,num_x_ff
            allocate(cll_f_2(i1,i2,i3,i4)%ida
     &              (cll_f_2(i1,i2,i3,i4)%num))
          enddo
        enddo
      enddo

      cll_f_2(:,:,:,i4)%num=0

c     do iii=1,num_tot_atms
      do iii=num_mov_atms+1,num_tot_atms ! 2019 - only use static atoms
        cll_f_2(i1_lic(iii),i2_lic(iii),i3_lic(iii),i4)%num=
     &  cll_f_2(i1_lic(iii),i2_lic(iii),i3_lic(iii),i4)%num+1
        cll_f_2(i1_lic(iii),i2_lic(iii),i3_lic(iii),i4)%ida
     & (cll_f_2(i1_lic(iii),i2_lic(iii),i3_lic(iii),i4)%num)=iii
      enddo                                           
 
      return
      end

