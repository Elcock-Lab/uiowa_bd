
      subroutine write_xtc
                                      
      use allocatable_arrays         
      implicit real(a-h,o-z)        

      do kkk=1,num_reps
      do iii=1,num_t_ms
        do jjj=1,num_atms_pt(i_t_tp(iii))
          x1=c_a_x(1,id_beg(iii)+jjj,kkk)
          y1=c_a_x(2,id_beg(iii)+jjj,kkk)
          z1=c_a_x(3,id_beg(iii)+jjj,kkk)
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
          ijk=id_beg(iii)+jjj
          x_xtc(ijk*3-2)=x1/10.0
          x_xtc(ijk*3-1)=y1/10.0
          x_xtc(ijk*3-0)=z1/10.0
        enddo
      enddo

      box(1)=xlen/10.0
      box(5)=ylen/10.0
      box(9)=zlen/10.0
      box(2)=0.0
      box(3)=0.0
      box(4)=0.0
      box(6)=0.0
      box(7)=0.0
      box(8)=0.0
      prec=1000.0 ! how does this affect file size?
      step=1
      if(kkk.eq.1) then
        num_xtc=num_xtc+1
        write(*,*)'writing xtc file frame# ',num_xtc
        ntcount=0
      endif

      if(nstep.eq.0) then
        call xdrfopen(xd2,xtc_file_out(kkk),"w",ret) ! overwrite
      else
        call xdrfopen(xd2,xtc_file_out(kkk),"a",ret) ! try append?
      endif

      call writextc(xd2,num_tot_atms,step,sngl(tot_time),
     &              box,x_xtc,prec,ret)
      call xdrfclose(xd2,ret)

      enddo

      return
      end
