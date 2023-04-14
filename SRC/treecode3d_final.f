
C AHE - this is the same code as used in other directories with the
C exception that tforce has been redimensioned (1:3,1:num_q_as) instead
C of the other way around...

C
C   Primary author:  Hans Johnston  (hansjohn@umich.edu) 
C   Department of Mathematics
C   University of Michigan, Ann Arbor
C
C   Copyright (c) 2003. The Regents of the University of Michigan.
C   All Rights Reserved.
C
C   This program is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   This program is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software Foundation,
C   Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE treecode3d_targ_procedures
      IMPLICIT NONE

C r8 is 8-byte (double precision) real

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(4)

C global variables for taylor expansions

      INTEGER :: torder,torderlim
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: cf,cf1,cf2,cf3
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:,:) :: a,b

C global variables to track tree levels 
 
      INTEGER :: minlevel,maxlevel
      
C global variables used when computing potential/force

      INTEGER :: orderoffset
      REAL(KIND=r8),DIMENSION(3) :: tarpos
      REAL(KIND=r8) :: thetasq,gdist_tolsq
      
!$omp threadprivate(torder,torderlim,cf,cf1,cf2,cf3,a,b,
!$omp$              minlevel,maxlevel,
!$omp$              orderoffset,tarpos,thetasq,gdist_tolsq)

CC global variables for postition and charge storage
CC
CC      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xcopy,ycopy,zcopy,qcopy

C node pointer and node type declarations

      TYPE tnode_pointer
           TYPE(tnode), POINTER :: p_to_tnode
      END TYPE tnode_pointer
      TYPE tnode
           INTEGER          :: numpar,ibeg,iend
           REAL(KIND=r8)    :: x_min,y_min,z_min
           REAL(KIND=r8)    :: x_max,y_max,z_max
           REAL(KIND=r8)    :: x_mid,y_mid,z_mid
           REAL(KIND=r8)    :: radius,sqradius,aspect
           INTEGER          :: level,num_children,exist_ms
           REAL(KIND=r8),DIMENSION(:,:,:),POINTER :: ms
           TYPE(tnode_pointer), DIMENSION(8) :: child
      END TYPE tnode

      CONTAINS
!!!!!!!!!!!!!!!
      SUBROUTINE SETUP(x,y,z,q,numpars,order,theta,iflag,dist_tol,
     &                 xyzminmax,arrdim)
      IMPLICIT NONE
C
C SETUP allocates and initializes arrays needed for the Taylor expansion.
C Also, global variables are set and the Cartesian coordinates of
C the smallest box containing the particles is determined. The particle
C postions and charges are copied so that they can be restored upon exit.
C
      INTEGER,INTENT(IN) :: numpars,order,iflag,arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q
      REAL(KIND=r8),INTENT(INOUT),DIMENSION(6) :: xyzminmax
      REAL(KIND=r8),INTENT(IN) :: theta,dist_tol

C local variables

      INTEGER :: err,i
      REAL(KIND=r8) :: t1

C global integers and reals:  TORDER, TORDERLIM and THETASQ

      torder=order
      IF (iflag .EQ. 1) THEN
          orderoffset=0
      ELSE
          orderoffset=1
      END IF
      torderlim=torder+orderoffset
      thetasq=theta*theta
      gdist_tolsq=dist_tol*dist_tol

C allocate global Taylor expansion variables

      ALLOCATE(cf(0:torder), STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocationg Taylor variables! '
         STOP
      END IF

      ALLOCATE(cf1(torderlim),cf2(torderlim),cf3(torderlim),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating Taylor variables! '
         STOP
      END IF

      ALLOCATE(a(-2:torderlim,-2:torderlim,-2:torderlim),
     &         b(-2:torderlim,-2:torderlim,-2:torderlim),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating Taylor variables! '
         STOP
      END IF

C initialize arrays for Taylor sums and coeffs

      a=0.0_r8
      b=0.0_r8

      DO i=0,torder
         cf(i) = REAL(i,KIND=r8)+1.0_r8
      END DO

      DO i=1,torderlim
         t1=1.0_r8/REAL(i,KIND=r8)
         cf1(i)=t1
         cf2(i)=1.0_r8-0.5_r8*t1
         cf3(i)=1.0_r8-t1
      END DO

C find bounds of Cartesian box enclosing the particles

      xyzminmax(1)=MINVAL(x(1:numpars))
      xyzminmax(2)=MAXVAL(x(1:numpars))
      xyzminmax(3)=MINVAL(y(1:numpars))
      xyzminmax(4)=MAXVAL(y(1:numpars))
      xyzminmax(5)=MINVAL(z(1:numpars))
      xyzminmax(6)=MAXVAL(z(1:numpars))


CC Allocate arrays and copy variables. Also create and initialize orderarr.

C      ALLOCATE(xcopy(numpars),ycopy(numpars),zcopy(numpars), 
C     &         qcopy(numpars),STAT=err)
C      IF (err .NE. 0) THEN
C         WRITE(6,*) 'Error allocating copy variables! '
C         STOP
C      END IF    
C
C      xcopy=x(1:numpars)
C      ycopy=y(1:numpars)
C      zcopy=z(1:numpars)
C      qcopy=q(1:numpars)

      RETURN
      END SUBROUTINE SETUP
!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE CREATE_TREE(p,ibeg,iend,x,y,z,q,perm,shrink,
     &                                 maxparnode,xyzmm,level,arrdim)
      IMPLICIT NONE
C
C CREATE_TREE recursively create the tree structure. Node P is
C input, which contains particles indexed from IBEG to IEND. After
C the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
C Real array XYZMM contains the min and max values of the coordinates
C of the particle in P, thus defining the box.   
C
      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: ibeg,iend,shrink,level,maxparnode,arrdim
      INTEGER,DIMENSION(arrdim),INTENT(INOUT) :: perm      
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z,q
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzmm

C local variables

      REAL(KIND=r8) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax,t1,t2,t3
      INTEGER, DIMENSION(8,2) :: ind
      REAL(KIND=r8), DIMENSION(6,8) :: xyzmms
      INTEGER :: i,j,limin,limax,err,loclev,numposchild
      REAL(KIND=r8), DIMENSION(6) ::  lxyzmm
     
C allocate pointer 

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF

C set node fields: number of particles, exist_ms
C and xyz bounds 

      p%numpar=iend-ibeg+1
      p%exist_ms=0

      IF (shrink .EQ. 1) THEN   
         p%x_min=MINVAL(x(ibeg:iend))
         p%x_max=MAXVAL(x(ibeg:iend))
         p%y_min=MINVAL(y(ibeg:iend))
         p%y_max=MAXVAL(y(ibeg:iend))
         p%z_min=MINVAL(z(ibeg:iend))
         p%z_max=MAXVAL(z(ibeg:iend))
      ELSE
         p%x_min=xyzmm(1)
         p%x_max=xyzmm(2)
         p%y_min=xyzmm(3)
         p%y_max=xyzmm(4)
         p%z_min=xyzmm(5)
         p%z_max=xyzmm(6)        
      END IF

C compute aspect ratio

      xl=p%x_max-p%x_min
      yl=p%y_max-p%y_min
      zl=p%z_max-p%z_min

      lmax=MAX(xl,yl,zl)
      t1=lmax
      t2=MIN(xl,yl,zl)
      IF (t2 .NE. 0.0_r8) THEN
         p%aspect=t1/t2
      ELSE
         p%aspect=0.0_r8
      END IF

C midpoint coordinates , RADIUS and SQRADIUS 

      p%x_mid=(p%x_max+p%x_min)/2.0_r8
      p%y_mid=(p%y_max+p%y_min)/2.0_r8
      p%z_mid=(p%z_max+p%z_min)/2.0_r8
      t1=p%x_max-p%x_mid
      t2=p%y_max-p%y_mid
      t3=p%z_max-p%z_mid
      p%sqradius=t1*t1+t2*t2+t3*t3
      p%radius=SQRT(p%sqradius)

C set particle limits, tree level of node, and nullify children pointers

      p%ibeg=ibeg
      p%iend=iend
      p%level=level
      IF (maxlevel .LT. level) THEN
         maxlevel=level
      END IF
      p%num_children=0
      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO
      
C compute moments - NOT HERE, but WHEN NEEDED 

!      IF (p%exist_ms .EQ. 0) THEN
!         ALLOCATE(p%ms(0:torder,0:torder,0:torder),STAT=err)
!          IF (err .NE. 0) THEN
!              WRITE(6,*) 'Error allocating node moments! '
!              STOP
!          END IF
!          CALL COMP_MS(p,x,y,z,q,arrdim)
!          p%exist_ms=1
!       END IF   

      IF (p%numpar .GT. maxparnode) THEN
C
C set IND array to 0 and then call PARTITION routine.  IND array holds indices
C of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
C
         xyzmms(1,1)=p%x_min
         xyzmms(2,1)=p%x_max
         xyzmms(3,1)=p%y_min
         xyzmms(4,1)=p%y_max
         xyzmms(5,1)=p%z_min
         xyzmms(6,1)=p%z_max
         ind(1,1)=ibeg
         ind(1,2)=iend
         x_mid=p%x_mid
         y_mid=p%y_mid
         z_mid=p%z_mid

         CALL PARTITION_8(x,y,z,q,perm,xyzmms,xl,yl,zl,lmax,
     &                    numposchild,x_mid,y_mid,z_mid,ind,arrdim)
C
C create children if indicated and store info in parent
C
         loclev=level+1
         DO i=1,numposchild
            IF (ind(i,1) .LE. ind(i,2)) THEN
               p%num_children=p%num_children+1
               lxyzmm=xyzmms(:,i)
               CALL CREATE_TREE(p%child(p%num_children)%p_to_tnode,
     &                          ind(i,1),ind(i,2),x,y,z,q,perm,shrink,
     &                          maxparnode,lxyzmm,loclev,arrdim)
            END IF
         END DO
      ELSE
         IF (level .LT. minlevel) THEN
            minlevel=level
         END IF
      END IF   

      END SUBROUTINE CREATE_TREE      
!!!!!!!!!!!!!!!
      SUBROUTINE PARTITION_8(x,y,z,q,perm,xyzmms,xl,yl,zl,lmax,
     &                       numposchild,x_mid,y_mid,z_mid,ind,arrdim)
      IMPLICIT NONE
C
C PARTITION_8 determines the particle indices of the eight sub boxes
C containing the particles after the box defined by particles I_BEG
C to I_END is divided by its midpoints in each coordinate direction.
C The determination of the indices is accomplished by the subroutine
C PARTITION. A box is divided in a coordinate direction as long as the
C resulting aspect ratio is not too large. This avoids the creation of
C "narrow" boxes in which Talyor expansions may become inefficient.
C On exit the INTEGER array IND (dimension 8 x 2) contains
C the indice limits of each new box (node) and NUMPOSCHILD the number 
C of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
C that box J is empty.
C
      INTEGER, INTENT(IN) :: arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z,q
      INTEGER,DIMENSION(arrdim),INTENT(INOUT) :: perm   
      INTEGER, DIMENSION(8,2),INTENT(INOUT) :: ind
      REAL(KIND=r8),DIMENSION(6,8),INTENT(INOUT) :: xyzmms
      REAL(KIND=r8), INTENT(IN) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax
      INTEGER,INTENT(INOUT) :: numposchild

C local variables

      INTEGER :: temp_ind,i
      REAL(KIND=r8) :: critlen

      numposchild=1
      critlen=lmax/sqrt(2.0_r8)

      IF (xl .GE. critlen) THEN
         CALL PARTITION(x,y,z,q,perm,ind(1,1),ind(1,2),
     &                  x_mid,temp_ind,arrdim)

         ind(2,1)=temp_ind+1
         ind(2,2)=ind(1,2)
         ind(1,2)=temp_ind
         xyzmms(:,2)=xyzmms(:,1)
         xyzmms(2,1)=x_mid
         xyzmms(1,2)=x_mid
         numposchild=2*numposchild
      END IF 
 
      IF (yl .GE. critlen) THEN
         DO i=1,numposchild
            CALL PARTITION(y,x,z,q,perm,ind(i,1),ind(i,2),
     &                     y_mid,temp_ind,arrdim)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(4,i)=y_mid
            xyzmms(3,numposchild+i)=y_mid
         END DO
         numposchild=2*numposchild
      END IF

      IF (zl .GE. critlen) THEN
         DO i=1,numposchild
            CALL PARTITION(z,x,y,q,perm,ind(i,1),ind(i,2),
     &                     z_mid,temp_ind,arrdim)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(6,i)=z_mid
            xyzmms(5,numposchild+i)=z_mid
         END DO
         numposchild=2*numposchild
      END IF

      RETURN 
      END SUBROUTINE PARTITION_8
!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE TREE_COMPFP(p,xtar,ytar,ztar,numtars,
     &                       pnum_beg,pnum_end,tpengtar,
     &                       x,y,z,q,tforce,numpars,farrdim,arrdim,
     &                       i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

      IMPLICIT NONE
C
C TREE_COMPF is the driver routine which calls COMPF_TREE for each target
C particle, setting the global variables TARPOS before the call. 
C P is the root node of the tree. 
C

      INTEGER,INTENT(IN) :: numpars,farrdim,arrdim,numtars,
     &                      pnum_beg,pnum_end,i_pbc
      REAL(KIND=r8),INTENT(IN) :: xlen,ylen,zlen,xinv,yinv,zinv

      LOGICAL,INTENT(IN) ::       ewald_elec_grid
      LOGICAL,INTENT(IN) ::       implicit_ions     
      INTEGER,INTENT(IN) ::       ewald_elec_grid_ndata,
     &                            ewald_elec_grid_nx
      REAL(KIND=r8),INTENT(IN) :: dielectric_inv,
     &                            kappa,
     &                            r_ion,
     &                            one_over_one_plus_kappa_rion,
     &                            ewald_elec_grid_rx,
     &                            ewald_elec_grid_rx_inv
      REAL(KIND=r8),DIMENSION(ewald_elec_grid_ndata),
     &              INTENT(IN) :: ewald_elec_grid_ff

      TYPE(tnode),POINTER :: p  
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q
      REAL(KIND=r8),DIMENSION(numtars),INTENT(IN) :: xtar,ytar,ztar 
      REAL(KIND=r8),DIMENSION(3,farrdim),INTENT(INOUT) :: tforce
      REAL(KIND=r8),DIMENSION(numtars),INTENT(INOUT) :: tpengtar

C local variables

      INTEGER :: i,j
      REAL(KIND=r8),DIMENSION(3) :: flocal,f
      REAL(KIND=r8) :: penglocal,peng

      IF (p%num_children .EQ. 0) THEN
         DO i = pnum_beg,pnum_end
            tarpos(1)=xtar(i)
            tarpos(2)=ytar(i)
            tarpos(3)=ztar(i)
            CALL COMPFP_TREE(p,peng,f,x,y,z,q,arrdim,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

            tforce(1:3,i) = -f(1:3) 
!            fx(i)=-f(1)
!            fy(i)=-f(2)
!            fz(i)=-f(3)
            tpengtar(i)=peng
         END DO
      ELSE
         DO i = pnum_beg,pnum_end
            tarpos(1)=xtar(i)
            tarpos(2)=ytar(i)
            tarpos(3)=ztar(i)
            penglocal=0.0_r8
            flocal=0.0_r8
            DO j=1,p%num_children
               CALL COMPFP_TREE(p%child(j)%p_to_tnode,peng,
     &                         f,x,y,z,q,arrdim,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

               flocal=flocal+f
               penglocal=penglocal+peng     
            END DO
            tforce(1:3,i) = -flocal(1:3)            
!            fx(i)=-flocal(1)
!            fy(i)=-flocal(2)
!            fz(i)=-flocal(3)
            tpengtar(i)=penglocal
         END DO
      END IF

      RETURN
      END SUBROUTINE TREE_COMPFP
!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE COMPFP_TREE(p,peng,f,x,y,z,q,arrdim,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: arrdim,i_pbc
      REAL(KIND=r8),INTENT(IN) :: xlen,ylen,zlen,xinv,yinv,zinv

      LOGICAL,INTENT(IN) ::       ewald_elec_grid
      LOGICAL,INTENT(IN) ::       implicit_ions     
      INTEGER,INTENT(IN) ::       ewald_elec_grid_ndata,
     &                            ewald_elec_grid_nx
      REAL(KIND=r8),INTENT(IN) :: dielectric_inv,
     &                            kappa,
     &                            r_ion,
     &                            one_over_one_plus_kappa_rion,
     &                            ewald_elec_grid_rx,
     &                            ewald_elec_grid_rx_inv
      REAL(KIND=r8),DIMENSION(ewald_elec_grid_ndata),
     &              INTENT(IN) :: ewald_elec_grid_ff

      TYPE(tnode),POINTER :: p      
      REAL(KIND=r8),INTENT(INOUT) :: peng
      REAL(KIND=r8),DIMENSION(3),INTENT(INOUT) :: f
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q

C local variables

      REAL(KIND=r8),DIMENSION(3) :: flocal
      REAL(KIND=r8) :: tx,ty,tz,distsq,t1,t2,penglocal,cfj,cfk
      INTEGER :: i,j,k,err

      INTEGER :: ijk,jkl,klm,lmn,lmn16a,lmn16b
      REAL(KIND=r8) :: dist,dist2,dist_inv,u_v,f_v
      REAL(KIND=r8) :: e1,e2,e3
      REAL(KIND=r8) :: g1,g2,g3
      REAL(KIND=r8) :: o1,o2,o3
      REAL(KIND=r8) :: w1,w2,w3
      REAL(KIND=r8) :: u_g,f_g1,f_g2,f_g3
      REAL(KIND=r8) :: u000,u001,u010,u011,u100,u101,u110,u111
      REAL(KIND=r8) :: dx,dy,dz,dxdy,dxdz,dydz,dxdydz
      REAL(KIND=r8) :: dy_dydz,dz_dydz,dx_dxdy,dz_dxdz
  
C determine DISTSQ for MAC test

      tx=tarpos(1)-p%x_mid
      ty=tarpos(2)-p%y_mid
      tz=tarpos(3)-p%z_mid
c     write(*,701)tarpos(1),tarpos(2),tarpos(3)
c     write(*,702)p%x_mid,p%y_mid,p%z_mid
701   format('tarpos ',3f12.5)
702   format('p%xyz  ',3f12.5)

C AHE adds correction for minimum image here

      if(i_pbc.eq.1) then
        tx=tx-xlen*anint(tx*xinv)
        ty=ty-ylen*anint(ty*yinv)
        tz=tz-zlen*anint(tz*zinv)
      endif

      distsq=tx*tx+ty*ty+tz*tz

C intialize potential energy and force 

      peng=0.0_r8
      f=0.0_r8

C If MAC is accepted and there is more than 1 particle in the 
C box use the expansion for the approximation.

C AHE wonders whether to *force* the routine to go through this section even
C when there is only one particle in the box we just comment out the
C second line of this if statement... (but it doesn't work)

      IF ((p%sqradius .LT. distsq*thetasq) .AND.
     &    (p%sqradius .NE. 0.0_r8)) THEN
         IF (p%exist_ms .EQ. 0) THEN
             ALLOCATE(p%ms(0:torder,0:torder,0:torder),STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error allocating node moments! '
                STOP
             END IF
             CALL COMP_MS(p,x,y,z,q,arrdim)
             p%exist_ms=1
         END IF   
         penglocal=0.0_r8   
         flocal=0.0_r8
         CALL COMP_TCOEFF(p,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

C add in separate blocks for Debye-Huckel electrostatics here
C note that in comparison with our 'old' approach we no 
C longer multiply by f_v or u_v - these were previously used as 
C some kind of conversion between pure coulomb and DH - not needed any
C more - therefore the code is much cleaner

C also note that we subtract the force terms rather than adding them as
C in the old code - presumably Hans switched signs on us again...

C finally, note that the energy terms are HALVED - multliplied by 0.5 -
C this is because atom i and j of an interacting pair are *both* called

c     write(*,901)tarpos(1),tarpos(2),tarpos(3)
901   format('using multipole for atom at pos ',3f12.5)

         if(implicit_ions) then

           dist2=tx**2+ty**2+tz**2
           dist=sqrt(dist2)
           dist_inv=1.0_r8/dist

           DO k=0,torder
              cfk=cf(k)
              DO j=0,torder-k
                 cfj=cf(j)
                 DO i=0,torder-k-j

                    flocal(1)=flocal(1)-cf(i)*a(i+1,j,k)
     &                       *p%ms(i,j,k)
                    flocal(2)=flocal(2)-cfj*a(i,j+1,k)
     &                       *p%ms(i,j,k)
                    flocal(3)=flocal(3)-cfk*a(i,j,k+1)
     &                       *p%ms(i,j,k)
                    peng=peng+a(i,j,k)*p%ms(i,j,k)*0.5

c                   write(*,88)i,j,k,cf(i),cfj,cfk,
c    &  a(i+1,j,k),a(i,j+1,k),a(i,j,k+1),p%ms(i,j,k),
c    &  a(i,j,k)*p%ms(i,j,k)*0.5
88      format('help me out ',3i8,12f10.5)
          
                 END DO
              END DO
           END DO

         else                    

           DO k=0,torder
              cfk=cf(k)
              DO j=0,torder-k
                 cfj=cf(j)
                 DO i=0,torder-k-j
  
                    flocal(1)=flocal(1)-cf(i)*a(i+1,j,k)
     &                       *p%ms(i,j,k)
                    flocal(2)=flocal(2)-cfj*a(i,j+1,k)
     &                       *p%ms(i,j,k)
                    flocal(3)=flocal(3)-cfk*a(i,j,k+1)
     &                       *p%ms(i,j,k)
                    peng=peng+a(i,j,k)*p%ms(i,j,k)*0.5

c                   write(*,88)i,j,k,cf(i),cfj,cfk,
c    &  a(i+1,j,k),a(i,j+1,k),a(i,j,k+1),p%ms(i,j,k),
c    &  a(i,j,k)*p%ms(i,j,k)*0.5
          
                 END DO
              END DO
           END DO
         endif                   

         f=f+flocal

      ELSE

C If MAC fails check to see if the are children. If not, perform direct 
C calculation.  If there are children, call routine recursively for each.
C
         IF (p%num_children .EQ. 0) THEN
            CALL COMPFP_DIRECT(penglocal,flocal,p%ibeg,p%iend,
     &                         x,y,z,q,arrdim,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

            peng=penglocal
            f=flocal
         ELSE
            DO i=1,p%num_children
               CALL COMPFP_TREE(p%child(i)%p_to_tnode,penglocal,
     &                         flocal,x,y,z,q,arrdim,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

               peng=peng+penglocal
               f=f+flocal
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPFP_TREE
!!!!!!!!
      SUBROUTINE COMP_TCOEFF(p,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

      IMPLICIT NONE
C
C COMP_TCOEFF computes the Taylor coefficients of the potential
C using a recurrence formula.  The center of the expansion is the
C midpoint of the node P.  TARPOS and TORDERLIM are globally defined.  
C
      INTEGER,INTENT(IN) :: i_pbc
      REAL(KIND=r8),INTENT(IN) :: xlen,ylen,zlen,xinv,yinv,zinv

      LOGICAL,INTENT(IN) ::       ewald_elec_grid
      LOGICAL,INTENT(IN) ::       implicit_ions     
      INTEGER,INTENT(IN) ::       ewald_elec_grid_ndata,
     &                            ewald_elec_grid_nx
      REAL(KIND=r8),INTENT(IN) :: dielectric_inv,
     &                            kappa,
     &                            r_ion,
     &                            one_over_one_plus_kappa_rion,
     &                            ewald_elec_grid_rx,
     &                            ewald_elec_grid_rx_inv
      REAL(KIND=r8),DIMENSION(ewald_elec_grid_ndata),
     &              INTENT(IN) :: ewald_elec_grid_ff

      TYPE(tnode),POINTER :: p      

C local varaibles

      REAL(KIND=r8) :: dx,dy,dz,ddx,ddy,ddz,dist,fac
      REAL(KIND=r8) :: kappax,kappay,kappaz
      INTEGER :: i,j,k,tp1

C setup variables

      dx=tarpos(1)-p%x_mid
      dy=tarpos(2)-p%y_mid
      dz=tarpos(3)-p%z_mid

C AHE adds correction for minimum image here

      if(i_pbc.eq.1) then
        dx=dx-xlen*anint(dx*xinv)
        dy=dy-ylen*anint(dy*yinv)
        dz=dz-zlen*anint(dz*zinv)
      endif

      ddx=2.0_r8*dx
      ddy=2.0_r8*dy
      ddz=2.0_r8*dz

      kappax=kappa*dx
      kappay=kappa*dy
      kappaz=kappa*dz

      dist=dx*dx+dy*dy+dz*dz
      fac=1.0_r8/dist
      dist=SQRT(dist)

C 0th coeff or function val

      b(0,0,0)=EXP(-kappa*dist)
      a(0,0,0)=b(0,0,0)/dist

C 2 indices are 0

      b(1,0,0)=kappax*a(0,0,0)
      b(0,1,0)=kappay*a(0,0,0)
      b(0,0,1)=kappaz*a(0,0,0)

      a(1,0,0)=fac*dx*(a(0,0,0)+kappa*b(0,0,0))
      a(0,1,0)=fac*dy*(a(0,0,0)+kappa*b(0,0,0))
      a(0,0,1)=fac*dz*(a(0,0,0)+kappa*b(0,0,0))

      DO i=2,torderlim
         b(i,0,0)=cf1(i)*kappa*(dx*a(i-1,0,0)-a(i-2,0,0))
         b(0,i,0)=cf1(i)*kappa*(dy*a(0,i-1,0)-a(0,i-2,0))
         b(0,0,i)=cf1(i)*kappa*(dz*a(0,0,i-1)-a(0,0,i-2))

         a(i,0,0)=fac*(ddx*cf2(i)*a(i-1,0,0)-cf3(i)*a(i-2,0,0)+
     &            cf1(i)*kappa*(dx*b(i-1,0,0)-b(i-2,0,0)))
         a(0,i,0)=fac*(ddy*cf2(i)*a(0,i-1,0)-cf3(i)*a(0,i-2,0)+
     &            cf1(i)*kappa*(dy*b(0,i-1,0)-b(0,i-2,0)))
         a(0,0,i)=fac*(ddz*cf2(i)*a(0,0,i-1)-cf3(i)*a(0,0,i-2)+
     &            cf1(i)*kappa*(dz*b(0,0,i-1)-b(0,0,i-2)))
      END DO

C 1 index 0, 1 index 1, other >=1

      b(1,1,0)=kappax*a(0,1,0)
      b(1,0,1)=kappax*a(0,0,1)
      b(0,1,1)=kappay*a(0,0,1)

      a(1,1,0)=fac*(dx*a(0,1,0)+ddy*a(1,0,0)+kappax*b(0,1,0))
      a(1,0,1)=fac*(dx*a(0,0,1)+ddz*a(1,0,0)+kappax*b(0,0,1))
      a(0,1,1)=fac*(dy*a(0,0,1)+ddz*a(0,1,0)+kappay*b(0,0,1))

      DO i=2,torderlim-1
         b(1,0,i)=kappax*a(0,0,i)
         b(0,1,i)=kappay*a(0,0,i)
         b(0,i,1)=kappaz*a(0,i,0)
         b(1,i,0)=kappax*a(0,i,0)
         b(i,1,0)=kappay*a(i,0,0)
         b(i,0,1)=kappaz*a(i,0,0)

         a(1,0,i)=fac*(dx*a(0,0,i)+ddz*a(1,0,i-1)-a(1,0,i-2)+
     &            kappax*b(0,0,i))
         a(0,1,i)=fac*(dy*a(0,0,i)+ddz*a(0,1,i-1)-a(0,1,i-2)+
     &            kappay*b(0,0,i))
         a(0,i,1)=fac*(dz*a(0,i,0)+ddy*a(0,i-1,1)-a(0,i-2,1)+
     &            kappaz*b(0,i,0))
         a(1,i,0)=fac*(dx*a(0,i,0)+ddy*a(1,i-1,0)-a(1,i-2,0)+
     &            kappax*b(0,i,0))
         a(i,1,0)=fac*(dy*a(i,0,0)+ddx*a(i-1,1,0)-a(i-2,1,0)+
     &            kappay*b(i,0,0))
         a(i,0,1)=fac*(dz*a(i,0,0)+ddx*a(i-1,0,1)-a(i-2,0,1)+
     &            kappaz*b(i,0,0))
      END DO

C 1 index 0, others >= 2

      DO i=2,torderlim-2
         DO j=2,torderlim-i
            b(i,j,0)=cf1(i)*kappa*(dx*a(i-1,j,0)-a(i-2,j,0))
            b(i,0,j)=cf1(i)*kappa*(dx*a(i-1,0,j)-a(i-2,0,j))
            b(0,i,j)=cf1(i)*kappa*(dy*a(0,i-1,j)-a(0,i-2,j))

            a(i,j,0)=fac*(ddx*cf2(i)*a(i-1,j,0)+ddy*a(i,j-1,0)
     &               -cf3(i)*a(i-2,j,0)-a(i,j-2,0)+
     &               cf1(i)*kappa*(dx*b(i-1,j,0)-b(i-2,j,0)))
            a(i,0,j)=fac*(ddx*cf2(i)*a(i-1,0,j)+ddz*a(i,0,j-1)
     &               -cf3(i)*a(i-2,0,j)-a(i,0,j-2)+
     &               cf1(i)*kappa*(dx*b(i-1,0,j)-b(i-2,0,j)))
            a(0,i,j)=fac*(ddy*cf2(i)*a(0,i-1,j)+ddz*a(0,i,j-1)
     &               -cf3(i)*a(0,i-2,j)-a(0,i,j-2)+
     &               cf1(i)*kappa*(dy*b(0,i-1,j)-b(0,i-2,j)))
         END DO
      END DO

C 2 indices 1, other >= 1
C b(1,1,1) is correct, but a little tricky!
C      b(1,1,1)=5.0*dz*fac*b(1,1,0)

      b(1,1,1)=kappax*a(0,1,1)
      a(1,1,1)=fac*(dx*a(0,1,1)+ddy*a(1,0,1)+ddz*a(1,1,0)+
     &         kappax*b(0,1,1))

      DO i=2,torderlim-2
         b(1,1,i)=kappax*a(0,1,i)
         b(1,i,1)=kappax*a(0,i,1)
         b(i,1,1)=kappay*a(i,0,1)

         a(1,1,i)=fac*(dx*a(0,1,i)+ddy*a(1,0,i)+ddz*a(1,1,i-1)
     &           -a(1,1,i-2)+kappax*b(0,1,i))
         a(1,i,1)=fac*(dx*a(0,i,1)+ddy*a(1,i-1,1)+ddz*a(1,i,0)
     &           -a(1,i-2,1)+kappax*b(0,i,1))
         a(i,1,1)=fac*(dy*a(i,0,1)+ddx*a(i-1,1,1)+ddz*a(i,1,0)
     &           -a(i-2,1,1)+kappay*b(i,0,1))
      END DO

C 1 index 1, others >=2

      DO i=2,torderlim-3
         DO j=2,torderlim-i
            b(1,i,j)=kappax*a(0,i,j)
            b(i,1,j)=kappay*a(i,0,j)
            b(i,j,1)=kappaz*a(i,j,0)

            a(1,i,j)=fac*(dx*a(0,i,j)+ddy*a(1,i-1,j)+ddz*a(1,i,j-1)
     &              -a(1,i-2,j)-a(1,i,j-2)+kappax*b(0,i,j))
            a(i,1,j)=fac*(dy*a(i,0,j)+ddx*a(i-1,1,j)+ddz*a(i,1,j-1)
     &              -a(i-2,1,j)-a(i,1,j-2)+kappay*b(i,0,j))
            a(i,j,1)=fac*(dz*a(i,j,0)+ddx*a(i-1,j,1)+ddy*a(i,j-1,1)
     &              -a(i-2,j,1)-a(i,j-2,1)+kappaz*b(i,j,0))

         END DO
      END DO

C all indices >=2

      DO k=2,torderlim-4
         DO j=2,torderlim-2-k
            DO i=2,torderlim-k-j
               b(i,j,k)=cf1(i)*kappa*(dx*a(i-1,j,k)-a(i-2,j,k))

               a(i,j,k)=fac*(ddx*cf2(i)*a(i-1,j,k)+ddy*a(i,j-1,k)
     &                 +ddz*a(i,j,k-1)-cf3(i)*a(i-2,j,k)
     &                 -a(i,j-2,k)-a(i,j,k-2)+
     &                 cf1(i)*kappa*(dx*b(i-1,j,k)-b(i-2,j,k)))
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE COMP_TCOEFF    
!!!!!!!!!!!!!!!
      SUBROUTINE COMP_MS(p,x,y,z,q,arrdim)
      IMPLICIT NONE
C
C COMP_MS computes the moments for node P needed in the Taylor approximation
C
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p 
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q

C local variables

      INTEGER :: i,k1,k2,k3
      REAL(KIND=r8) :: dx,dy,dz,tx,ty,tz
     
      p%ms=0.0_r8
      DO i=p%ibeg,p%iend
         dx=x(i)-p%x_mid
         dy=y(i)-p%y_mid
         dz=z(i)-p%z_mid
         tz=1.0_r8
         DO k3=0,torder
            ty=1.0_r8
            DO k2=0,torder-k3
               tx=1.0_r8 
               DO k1=0,torder-k3-k2
                  p%ms(k1,k2,k3)=p%ms(k1,k2,k3)+q(i)*tx*ty*tz
                  tx=tx*dx
              END DO
              ty=ty*dy
            END DO
            tz=tz*dz
         END DO
      END DO
         
      RETURN
      END SUBROUTINE COMP_MS
!!!!!!!!!!!!!!!!
      SUBROUTINE COMPFP_DIRECT(peng,f,ibeg,iend,x,y,z,q,arrdim,
     &        i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

      IMPLICIT NONE
C
C COMPF_DIRECT directly computes the force on the current target
C particle determined by the global variable TARPOS. 
C
      INTEGER,INTENT(IN) :: ibeg,iend,arrdim,i_pbc
      REAL(KIND=r8),INTENT(IN) :: xlen,ylen,zlen,xinv,yinv,zinv

      LOGICAL,INTENT(IN) ::       ewald_elec_grid
      LOGICAL,INTENT(IN) ::       implicit_ions     
      INTEGER,INTENT(IN) ::       ewald_elec_grid_ndata,
     &                            ewald_elec_grid_nx
      REAL(KIND=r8),INTENT(IN) :: dielectric_inv,
     &                            kappa,
     &                            r_ion,
     &                            one_over_one_plus_kappa_rion,
     &                            ewald_elec_grid_rx,
     &                            ewald_elec_grid_rx_inv
      REAL(KIND=r8),DIMENSION(ewald_elec_grid_ndata),
     &              INTENT(IN) :: ewald_elec_grid_ff

      REAL(KIND=r8),DIMENSION(3),INTENT(OUT) :: f
      REAL(KIND=r8),INTENT(OUT) :: peng
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q

C local variables

      INTEGER :: i
      REAL(KIND=r8) :: t1,t2,tx,ty,tz

      INTEGER :: ijk,jkl,klm,lmn,lmn16a,lmn16b
      REAL(KIND=r8) :: dist,dist2,dist_inv,u_v,f_v
      REAL(KIND=r8) :: e1,e2,e3
      REAL(KIND=r8) :: g1,g2,g3
      REAL(KIND=r8) :: o1,o2,o3
      REAL(KIND=r8) :: w1,w2,w3
      REAL(KIND=r8) :: u_g,f_g1,f_g2,f_g3
      REAL(KIND=r8) :: u000,u001,u010,u011,u100,u101,u110,u111
      REAL(KIND=r8) :: dx,dy,dz,dxdy,dxdz,dydz,dxdydz
      REAL(KIND=r8) :: dy_dydz,dz_dydz,dx_dxdy,dz_dxdz
  
      peng=0.0_r8
      f=0.0_r8

C AHE adds separate loops for with/without pbc here
C note that Hans calculates the negative of the force here

C do the following if no periodic boundaries
C note that there's two options here: with/without Debye-Huckel
C the Debye-Huckel code is the same as that in the BD code

      if(i_pbc.ne.1) then

        if(implicit_ions) then
          DO i=ibeg,iend
             tx=x(i)-tarpos(1)
             ty=y(i)-tarpos(2)
             tz=z(i)-tarpos(3)
             dist2=tx*tx+ty*ty+tz*tz
             IF (dist2 .GT. gdist_tolsq) THEN
                dist=sqrt(dist2)
                dist_inv=1.0_r8/dist
                u_v=exp(-kappa*(dist-r_ion))*dist_inv*
     &              one_over_one_plus_kappa_rion
                f_v=u_v*(kappa+dist_inv)*dist_inv
                peng=peng+q(i)*u_v*0.5
                f(1)=f(1)+q(i)*tx*f_v
                f(2)=f(2)+q(i)*ty*f_v
                f(3)=f(3)+q(i)*tz*f_v
             END IF 
          END DO   
        else

C do the following if not using Debye-Huckel

          DO i=ibeg,iend
             tx=x(i)-tarpos(1)
             ty=y(i)-tarpos(2)
             tz=z(i)-tarpos(3)
             t2=tx*tx+ty*ty+tz*tz
             IF (t2 .GT. gdist_tolsq) THEN
                t2=1.0_r8/SQRT(t2)
                t1=t2*t2*t2
                peng=peng+q(i)*t2*0.5
                f(1)=f(1)+q(i)*tx*t1
                f(2)=f(2)+q(i)*ty*t1
                f(3)=f(3)+q(i)*tz*t1
             END IF 
          END DO   
        endif 

C do the following if periodic boundaries and minimum image
C note that there's two options here: with/without Debye-Huckel
C the Debye-Huckel code is the same as that in the BD code

      elseif(i_pbc.eq.1.and..not.ewald_elec_grid) then

        if(implicit_ions) then
          DO i=ibeg,iend
             tx=x(i)-tarpos(1)
             ty=y(i)-tarpos(2)
             tz=z(i)-tarpos(3)
             tx=tx-xlen*anint(tx*xinv)
             ty=ty-ylen*anint(ty*yinv)
             tz=tz-zlen*anint(tz*zinv)
             dist2=tx*tx+ty*ty+tz*tz
             IF (dist2 .GT. gdist_tolsq) THEN
                dist=sqrt(dist2)
                dist_inv=1.0/dist
                u_v=exp(-kappa*(dist-r_ion))*dist_inv*
     &              one_over_one_plus_kappa_rion
                f_v=u_v*(kappa+dist_inv)*dist_inv
                peng=peng+q(i)*u_v*0.5
                f(1)=f(1)+q(i)*tx*f_v
                f(2)=f(2)+q(i)*ty*f_v
                f(3)=f(3)+q(i)*tz*f_v
             END IF 
          END DO   
        else

C do the following if not using Debye-Huckel

          DO i=ibeg,iend
             tx=x(i)-tarpos(1)
             ty=y(i)-tarpos(2)
             tz=z(i)-tarpos(3)
             tx=tx-xlen*anint(tx*xinv)
             ty=ty-ylen*anint(ty*yinv)
             tz=tz-zlen*anint(tz*zinv)
             t2=tx*tx+ty*ty+tz*tz
             IF (t2 .GT. gdist_tolsq) THEN
                t2=1.0_r8/SQRT(t2)
                t1=t2*t2*t2
                peng=peng+q(i)*t2*0.5
                f(1)=f(1)+q(i)*tx*t1
                f(2)=f(2)+q(i)*ty*t1
                f(3)=f(3)+q(i)*tz*t1
             END IF 
          END DO   
        endif 

      endif

      RETURN
      END SUBROUTINE COMPFP_DIRECT
!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEANUP(p)
      IMPLICIT NONE
C
C CLEANUP deallocates allocated global variables and then
C calls recursive routine REMOVE_NODE to delete the tree.
C
      TYPE(tnode),POINTER :: p      

C local variables
  
      INTEGER :: err

      DEALLOCATE(cf,cf1,cf2,cf3,a,b, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating Taylor variables! '
         STOP
      END IF      

CC      DEALLOCATE(xcopy,ycopy,zcopy,qcopy,STAT=err)
CC      IF (err .NE. 0) THEN
CC         WRITE(6,*) 'Error deallocating copy variables! '
CC         STOP
CC      END IF  

      CALL REMOVE_NODE(p)
      DEALLOCATE(p, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating root node! '
         STOP
      END IF 
      NULLIFY(p)         

      RETURN
      END SUBROUTINE CLEANUP
!!!!!!!!!!!
      RECURSIVE SUBROUTINE REMOVE_NODE(p)
      IMPLICIT NONE
C
C REMOVE_NODE recursively removes each node from the
C tree and deallocates its memory for MS array if it
C exits.
C
      TYPE(tnode),POINTER :: p 

C local variables

      INTEGER :: i,err

      IF (p%exist_ms .EQ. 1) THEN
         DEALLOCATE(p%ms,STAT=err)
         IF (err .NE. 0) THEN
            WRITE(6,*) 'Error deallocating node MS! '
            STOP
         END IF               
      END IF

      IF (p%num_children .GT. 0) THEN
          DO i=1,p%num_children
            CALL REMOVE_NODE(p%child(i)%p_to_tnode)
            DEALLOCATE(p%child(i)%p_to_tnode,STAT=err)
            IF (err .NE. 0) THEN
               WRITE(6,*) 'Error deallocating node child! '
               STOP
            END IF                           
          END DO
      END IF 

      RETURN                
      END SUBROUTINE REMOVE_NODE      

      END MODULE treecode3d_targ_procedures
!!!!!!!!!!!!!!!
      SUBROUTINE TREECODE3D_TARG_OpenMP(xtar,ytar,ztar,numtars,
     &                           miv_num_loc,miv_beq_loc,miv_enq_loc,
     &                           tpengtar,tforce,x,y,z,q,perm,numpars,
     &                           iflag,order,theta,shrink,
     &                           maxparnode,dist_tol,
     &                           timetree,prinfo,farrdim,arrdim,myid,
     &                           i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

      USE treecode3d_targ_procedures
      IMPLICIT NONE
C
C TREECODE3D approximates the electrostatic force on particles using a hierarchical
C oct-tree algorithm with body-cell forces approximated via Taylor expansions. 
C
C Input:
C
C XTAR,YTAR,ZTAR : (REAL one-dimensional arrays of length NUMTARS) 
C                  Cartesian coordinates of the target particles 
C
C NUMTARS : (INTEGER)  number of target particles
C
C DIST_TOL : (REAL) Potential energy/force contribution is excluded
C            if any particle lies within DIST_TOL of a target particle.
C            (NOTE: this parameter is only used in the case of a 
C                   direct computation)
C
C X,Y,Z : (REAL one-dimensional arrays of length ARRDIM) Cartesian coordinates
C         of the particles 
C
C Q   :  (REAL one-dimensional array of length ARRDIM) Corresponding charge 
C         of each particle
C
C NUMPARS : (INTEGER)  number of particles
C
C IFLAG   : (INTEGER)  flag indicating what is to be computer:
C                  1 - potential energy, 2 - potential energy and force
C
C ORDER :  (INTEGER) order of Taylor expansion used in the body-cell approximation
C
C THETA : (REAL) opening angle used in the BMAX MAC 
C
C SHRINK : (INTEGER) switch used in the oct-tree construction. If equal to 1
C          then upon subdivision of a node the node bounds are taken
C          the smallest Cartesian box containing the particles in the node.
C
C MAXPARNODE : (INTEGER) maximum number of particles allowed in a leaf.
C
C PRINFO : (INTEGER) swith for outputting information in this routine.  If set
C          to 1 then information is writen to stdout.
C
C FARRDIM : (INTEGER) declared length of arrays FX, FY, and FZ. 
C
C ARRDIM  : (INTEGER) declared length of arrays X, Y, Z, and Q.
C
C Output:
C
C TPENGTAR    : (REAL one-dimensional arrays of length NUMTARS) 
C               potential energy ith target particle.
C
C FX,FY,FZ : (REAL one-dimensional arrays of length FARRDIM) force on ith particle
C            is given by (FX(i),FY(i),FZ(i)).
C
C TREETIME : (REAL) time (in seconds) for the tree algorithm computation.
C
C
      integer,intent(in) :: miv_num_loc
      integer,intent(in) :: miv_beq_loc(miv_num_loc)
      integer,intent(in) :: miv_enq_loc(miv_num_loc)
      integer            :: nnn,pnum_beg,pnum_end

      INTEGER,INTENT(IN) :: numpars,iflag,order,farrdim,arrdim,
     &                      shrink,prinfo,maxparnode,numtars,
     &                      myid,i_pbc
      INTEGER,DIMENSION(arrdim),INTENT(INOUT) :: perm
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z,q
      REAL(KIND=r8),DIMENSION(numtars),INTENT(IN)   :: xtar,ytar,ztar
      REAL(KIND=r8),DIMENSION(numtars),INTENT(INOUT) :: tpengtar
      REAL(KIND=r8),DIMENSION(3,farrdim),INTENT(INOUT) :: tforce
      REAL(KIND=r8),INTENT(IN) :: theta,dist_tol
      REAL(KIND=r8),INTENT(IN) :: xlen,ylen,zlen,xinv,yinv,zinv

      LOGICAL,INTENT(IN) ::       ewald_elec_grid
      LOGICAL,INTENT(IN) ::       implicit_ions     
      INTEGER,INTENT(IN) ::       ewald_elec_grid_ndata,
     &                            ewald_elec_grid_nx
      REAL(KIND=r8),INTENT(IN) :: dielectric_inv,
     &                            kappa,
     &                            r_ion,
     &                            one_over_one_plus_kappa_rion,
     &                            ewald_elec_grid_rx,
     &                            ewald_elec_grid_rx_inv
      REAL(KIND=r8),DIMENSION(ewald_elec_grid_ndata),
     &              INTENT(IN) :: ewald_elec_grid_ff

      REAL(KIND=r8),INTENT(OUT) :: timetree

C local variables

      TYPE(tnode),POINTER :: troot
      INTEGER :: i,level,err
      REAL(KIND=r8), DIMENSION(6) :: xyzminmax

C variables needed for f90 DATE_AND_TIME intrinsic

      INTEGER,DIMENSION(8) :: time1,time2 
      CHARACTER (LEN=8)  :: datec
      CHARACTER (LEN=10) :: timec
      CHARACTER (LEN=5)  :: zonec
      REAL(KIND=r8)      :: totaltime


C Call SETUP to allocate arrays for Taylor expansions
C and setup global variables. Also, copy variables into global copy arrays. 

      CALL SETUP(x,y,z,q,numpars,order,theta,iflag,dist_tol,xyzminmax,
     &           arrdim)

C nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(troot)  

      CALL DATE_AND_TIME(datec,timec,zonec,time1)

C set global variables to track tree levels during construction

      level=0
      minlevel=50000
      maxlevel=0

      IF (prinfo .EQ. 1) THEN
         WRITE(6,*) ' '
         WRITE(6,*) 'Creating tree '
      END IF

      CALL CREATE_TREE(troot,1,numpars,x,y,z,q,perm,shrink,
     &                 maxparnode,xyzminmax,level,arrdim)

      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,totaltime)

C print tree information to stdout if PRINFO=1

      IF (prinfo .EQ. 1) THEN

         WRITE(6,*) ' '
         WRITE(6,*) 'Tree created '
         WRITE(6,*) ' '
         WRITE(6,*) 'Tree creation time (secs) : ',totaltime
         WRITE(6,*) ' '
         WRITE(6,*) 'Run synopsis : '
         WRITE(6,*) ' '
         WRITE(6,*) 'numpar :',troot%numpar
         WRITE(6,*) 'x_mid  :',troot%x_mid
         WRITE(6,*) 'y_mid  :',troot%y_mid
         WRITE(6,*) 'z_mid  :',troot%z_mid
         WRITE(6,*) 'radius :',troot%radius   
         WRITE(6,*) 'order           :',order
         WRITE(6,*) 'theta           :',theta
         WRITE(6,*) 'shrink          :',shrink
         WRITE(6,*) 'maxparnode      :',maxparnode
         WRITE(6,*) 'iflag           :',iflag
         WRITE(6,*) 'tree minlevel   :',minlevel
         WRITE(6,*) 'tree maxlevel   :',maxlevel
         WRITE(6,*) ' '
         WRITE(6,*) 'numtars         :',numtars
         WRITE(6,*) 'dist_tol        :',dist_tol
         WRITE(6,*) ' '
      END IF
 
      CALL DATE_AND_TIME(datec,timec,zonec,time1)

C AHE removed option of only calculating potential energy (TREE_COMPP)

      IF (iflag .NE. 1) THEN
         IF (prinfo .EQ. 1) THEN
            WRITE(6,*) ' '
            WRITE(6,*) 'Computing potential energy and force with tree '
         END IF

         do nnn=1,miv_num_loc

           pnum_beg=miv_beq_loc(nnn)
           pnum_end=miv_enq_loc(nnn)

         CALL TREE_COMPFP(troot,xtar,ytar,ztar,numtars,
     &                    pnum_beg,pnum_end,
     &                    tpengtar,
     &                    x,y,z,q,tforce,numpars,farrdim,arrdim,
     &                    i_pbc,xlen,ylen,zlen,xinv,yinv,zinv,
     &        ewald_elec_grid,
     &        dielectric_inv, 
     &        implicit_ions,
     &        kappa,
     &        r_ion,
     &        one_over_one_plus_kappa_rion,
     &        ewald_elec_grid_ndata,
     &        ewald_elec_grid_rx,
     &        ewald_elec_grid_rx_inv,
     &        ewald_elec_grid_nx,
     &        ewald_elec_grid_ff)

          enddo

      END IF

      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,totaltime)
      timetree=totaltime

      IF (prinfo .EQ. 1) THEN      
          WRITE(6,*) ' '
          WRITE(6,*) 'thread ',myid,'Treecode time (secs) : ',timetree
      END IF

CC restore x, y, z and q to their original values
CC
CC      x(1:numpars)=xcopy
CC      y(1:numpars)=ycopy
CC      z(1:numpars)=zcopy
CC      q(1:numpars)=qcopy

C Call CLEANUP to deallocate global variables and tree structure.

      IF (prinfo .EQ. 1) THEN
         WRITE(6,*) ' '
         WRITE(6,*) 'Deallocating tree structure!'
         WRITE(6,*) ' '
      END IF

      CALL CLEANUP(troot)

      END SUBROUTINE TREECODE3D_TARG_OpenMP
!!!!!!!!!!!!!
      SUBROUTINE PARTITION(a,b,c,q,perm,ibeg,iend,val,midind,arrdim)
      IMPLICIT NONE
C
C PARTITION determines the index MIDIND, after partitioning
C in place the  arrays A,B,C and Q,  such that 
C A(IBEG:MIDIND) <= VAL and  A(MIDIND+1:IEND) > VAL. 
C If on entry IBEG > IEND or  A(IBEG:IEND) > VAL then MIDIND
C is returned as IBEG-1. 
C 
      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(4)
      INTEGER, INTENT(IN) :: arrdim,ibeg,iend
      INTEGER,DIMENSION(arrdim),INTENT(INOUT) :: perm
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: a,b,c,q
      INTEGER, INTENT(INOUT) :: midind   
      REAL(KIND=r8) val

C local variables

      REAL(KIND=r8) ta,tb,tc,tq
      INTEGER lower,upper,tperm

      IF (ibeg .LT. iend) THEN

C temporarily store IBEG entries and set A(IBEG)=VAL for 
C the partitoning algorithm.  

         ta=a(ibeg)
         tb=b(ibeg)
         tc=c(ibeg)
         tq=q(ibeg)
         tperm = perm(ibeg)
         a(ibeg)=val 
         upper=ibeg
         lower=iend

         DO WHILE (upper .NE. lower)
            DO WHILE ((upper .LT. lower) .AND. (val .LT. a(lower)))
                  lower=lower-1
            END DO
            IF (upper .NE. lower) THEN
               a(upper)=a(lower)
               b(upper)=b(lower)
               c(upper)=c(lower)
               q(upper)=q(lower)
               perm(upper) = perm(lower)
            END IF
            DO WHILE ((upper .LT. lower) .AND. (val .GE. a(upper)))
                  upper=upper+1
            END DO
            IF (upper .NE. lower) THEN
               a(lower)=a(upper)
               b(lower)=b(upper)
               c(lower)=c(upper)
               q(lower)=q(upper)
               perm(lower) = perm(upper)
            END IF
         END DO
         midind=upper

C replace TA in position UPPER and change MIDIND if TA > VAL 

         IF (ta .GT. val) THEN
            midind=upper-1
         END IF
         a(upper)=ta
         b(upper)=tb
         c(upper)=tc
         q(upper)=tq
         perm(upper) = tperm

      ELSEIF (ibeg .EQ. iend) THEN
         IF (a(ibeg) .LE. val) THEN
            midind=ibeg
         ELSE
            midind=ibeg-1
         END IF
      ELSE
         midind=ibeg-1
      END IF

      RETURN
      END SUBROUTINE PARTITION
!!!!!!!!!!!!!!!!!!
      SUBROUTINE TTIME(timebeg,timeend,totaltime)
      IMPLICIT NONE
C
C TTIME computes the time difference in seconds between
C the timestamps TIMEBEG and TIMEEND returned by the 
C f90 intrinsic DATE_AND_TIME
C
      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(4)
      INTEGER,DIMENSION(8),INTENT(INOUT) :: timebeg,timeend
      REAL(KIND=r8),INTENT(OUT) :: totaltime

C TIMEEND is modifed by borrowing in case each of its fields
C are not .GE. to the corresponding field in TIMEBEG (up to
C and including days) 

      IF (timeend(8) .LT. timebeg(8)) THEN
          timeend(8)=timeend(8)+1000
          timeend(7)=timeend(7)-1
      END IF
      IF (timeend(7) .LT. timebeg(7)) THEN
          timeend(7)=timeend(7)+60
          timeend(6)=timeend(6)-1
      END IF
      IF (timeend(6) .LT. timebeg(6)) THEN
          timeend(6)=timeend(6)+60
          timeend(5)=timeend(5)-1
      END IF
      IF (timeend(5) .LT. timebeg(5)) THEN
          timeend(5)=timeend(5)+24
          timeend(3)=timeend(3)-1
      END IF

      totaltime=  REAL(timeend(8)-timebeg(8),KIND=r8) +
     &       1000.0_r8*( REAL(timeend(7)-timebeg(7),KIND=r8) + 
     &         60.0_r8*( REAL(timeend(6)-timebeg(6),KIND=r8) +
     &         60.0_r8*( REAL(timeend(5)-timebeg(5),KIND=r8) +
     &         24.0_r8*( REAL(timeend(3)-timebeg(3),KIND=r8)))))
      totaltime=totaltime/1000.0_r8

     
      RETURN
      END SUBROUTINE TTIME



