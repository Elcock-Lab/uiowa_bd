
      subroutine read_position_restraints                        
                                                                
      use allocatable_arrays                                   
      implicit real(a-h,o-z)                                  
      character*180            :: string                                                                            

C i_restraint = 0 means no restraint
C             = 1 x,y,z-type restaint
C             = 2 radial restraint...
C             = 3 capsule-type restraint

C note the format here 
C if we want a conventional 1D.2D.3D restraint then we read a '1' followed by the atom number, 
C the 3D coords of the position and the 3 force constants
C if we want a radial restraint then we read a '2' followed by the atom number, the radius, and
C the radial force constant...

C 32-11-11 note that additional flag used here

      open(unit=44,file=pos_restraint_file,status='unknown')  
3     read(44,22,end=10)string
22    format(a180)
      if(string(1:1).eq.'1') then

C following section provided for backward compatibility - i.e. days when
C no flag was found in the 119'th column to indicate whether nonbonded
C interactions were to be skipped for this particular atom

        if(string(119:119).eq.' ') then
          read(string,23)i,xt,yt,zt,fxt,fyt,fzt
          i_restraint(i)=1
          cx_restraint(i)=xt
          cy_restraint(i)=yt
          cz_restraint(i)=zt
          fx_restraint(i)=fxt
          fy_restraint(i)=fyt
          fz_restraint(i)=fzt
          i_res_no_nb(i)=0 ! assume that we don't skip nonbondeds
        else
          read(string,23)i,xt,yt,zt,fxt,fyt,fzt,j
          i_restraint(i)=1
          cx_restraint(i)=xt
          cy_restraint(i)=yt
          cz_restraint(i)=zt
          fx_restraint(i)=fxt
          fy_restraint(i)=fyt
          fz_restraint(i)=fzt
          i_res_no_nb(i)=j  
          if(j.gt.num_i_res_no_nb_grps) num_i_res_no_nb_grps=j
        endif
      elseif(string(1:1).eq.'2') then
        if(string(49:49).eq.' ') then
          read(string,24)i,x1,y1,z1,xr,fr
          i_restraint(i)=2
          cx_restraint(i)=x1
          cy_restraint(i)=y1
          cz_restraint(i)=z1
          r_restraint(i)=xr ! radius
          fr_restraint(i)=fr  ! force constant
          i_res_no_nb(i)=0 ! assume that we don't skip nonbondeds
        else
          read(string,24)i,x1,y1,z1,xr,fr,j
          i_restraint(i)=2
          cx_restraint(i)=x1
          cy_restraint(i)=y1
          cz_restraint(i)=z1
          r_restraint(i)=xr ! radius
          fr_restraint(i)=fr  ! force constant
          i_res_no_nb(i)=j 
          if(j.gt.num_i_res_no_nb_grps) num_i_res_no_nb_grps=j
        endif
      elseif(string(1:1).eq.'3') then 
        if(string(49:49).eq.' ') then
          read(string,25)i,x1,y1,z1,x2,y2,z2,xr,fr
          i_restraint(i)=3
          cx_restraint(i)=x1
          cy_restraint(i)=y1
          cz_restraint(i)=z1
          dx_restraint(i)=x2
          dy_restraint(i)=y2
          dz_restraint(i)=z2
          r_restraint(i)=xr ! radius
          fr_restraint(i)=fr  ! force constant
          i_res_no_nb(i)=0 ! assume that we don't skip nonbondeds
        else
          read(string,25)i,x1,y1,z1,x2,y2,z2,xr,fr,j
          i_restraint(i)=3
          cx_restraint(i)=x1
          cy_restraint(i)=y1
          cz_restraint(i)=z1
          dx_restraint(i)=x2
          dy_restraint(i)=y2
          dz_restraint(i)=z2
          r_restraint(i)=xr ! radius
          fr_restraint(i)=fr  ! force constant
          i_res_no_nb(i)=j
          if(j.gt.num_i_res_no_nb_grps) num_i_res_no_nb_grps=j
        endif
      elseif(string(1:1).eq.'4') then 
        if(string(49:49).eq.' ') then
          read(string,25)i,x1,y1,z1,x2,y2,z2,xr,fr
          i_restraint(i)=4
          cx_restraint(i)=x1
          cy_restraint(i)=y1
          cz_restraint(i)=z1
          dx_restraint(i)=x2
          dy_restraint(i)=y2
          dz_restraint(i)=z2
          r_restraint(i)=xr ! radius
          fr_restraint(i)=fr  ! force constant
          i_res_no_nb(i)=0 ! assume that we don't skip nonbondeds
        else
          read(string,25)i,x1,y1,z1,x2,y2,z2,xr,fr,j
          i_restraint(i)=4
          cx_restraint(i)=x1
          cy_restraint(i)=y1
          cz_restraint(i)=z1
          dx_restraint(i)=x2
          dy_restraint(i)=y2
          dz_restraint(i)=z2
          r_restraint(i)=xr ! radius
          fr_restraint(i)=fr  ! force constant
          i_res_no_nb(i)=j
          if(j.gt.num_i_res_no_nb_grps) num_i_res_no_nb_grps=j
        endif
      elseif(string(1:1).eq.'5') then 
        if(string(49:49).eq.' ') then
          read(string,25)i,x1,y1,z1,x2,y2,z2,xr,fr
          i_restraint(i)=5
          cx_restraint(i)=x1
          cy_restraint(i)=y1
          cz_restraint(i)=z1
          dx_restraint(i)=x2
          dy_restraint(i)=y2
          dz_restraint(i)=z2
          r_restraint(i)=xr ! radius
          fr_restraint(i)=fr  ! force constant
          i_res_no_nb(i)=0 ! assume that we don't skip nonbondeds
        else
          read(string,25)i,x1,y1,z1,x2,y2,z2,xr,fr,j
          i_restraint(i)=5
          cx_restraint(i)=x1
          cy_restraint(i)=y1
          cz_restraint(i)=z1
          dx_restraint(i)=x2
          dy_restraint(i)=y2
          dz_restraint(i)=z2
          r_restraint(i)=xr ! radius
          fr_restraint(i)=fr  ! force constant
          i_res_no_nb(i)=j
          if(j.gt.num_i_res_no_nb_grps) num_i_res_no_nb_grps=j
        endif
      endif
23    format(1x,i8,3f20.5,3f15.5,i5)
24    format(1x,i8,3f20.5,2f15.5,i5)
25    format(1x,i8,6f20.5,2f15.5,i5)
      goto 3
10    close(44)
      return
      end
