
      subroutine calculate_bonded_forces
     &          (natm,nmol,ibeg,iend,
     &           miv_num_loc,miv_beg_loc,miv_end_loc,
     &           nbnd,ibnd,dbnd,ebnd,
     &           nang,iang,dang,eang,
     &           ndih,idih,ddih,edih,
     &           epos,ewll,eumb,
     &           i_die_now,master,
     &           n_ll_rods,ida_ll_rods,jda_ll_rods,
     &           s11_ll_rods,s22_ll_rods,s33_ll_rods,s44_ll_rods,
     &           d11_ll_rods,d22_ll_rods,e00_ll_rods,
     &           ufx1,ufy1,ufz1,
     &           ufx2,ufy2,ufz2,
     &           fx,nosm,posm,mirep,mithrd)

      use allocatable_arrays

      implicit real(a-h,o-z)                                                                   

      integer miv_num_loc
      integer miv_beg_loc(miv_num_loc)
      integer miv_end_loc(miv_num_loc)
      integer nbnd
      integer nmol
      integer ibeg
      integer iend
      integer ibnd(nbnd)
      integer iang(nang)
      integer idih(ndih)
      integer nosm  ! =0 if no walls; =6 if walls
      integer mirep
      integer mithrd
      logical i_die_now
      logical master
      real    dbnd(nbnd)
      real    dang(nang)
      real    ddih(ndih)
      real    ebnd(nmol)
      real    eang(nmol)
      real    edih(nmol)
      real    epos(nmol)
      real    ewll(nmol)
      real    fx(3,natm)
      real    posm(nosm) ! osmotic pressure
      integer ida_ll_rods(1:n_ll_rods)
      integer jda_ll_rods(1:n_ll_rods)
      real    s11_ll_rods(1:n_ll_rods)
      real    s22_ll_rods(1:n_ll_rods)
      real    s33_ll_rods(1:n_ll_rods)
      real    s44_ll_rods(1:n_ll_rods)
      real    d11_ll_rods(1:n_ll_rods)
      real    d22_ll_rods(1:n_ll_rods)
      real    e00_ll_rods(1:n_ll_rods)

C note that here we have updated all reference to the coords to include
C the extra dimension for myrep; no change is required to the forces as
C these are effectively private to each thread at least for bonded terms

C [1] do the bonds

      do ijk=1,nbnd

        iii=ibnd(ijk)
        m1=nbb1_new(iii)
        m2=nbb2_new(iii)

        xij1=c_a_x(1,m1,mirep)-c_a_x(1,m2,mirep)
        xij2=c_a_x(2,m1,mirep)-c_a_x(2,m2,mirep)
        xij3=c_a_x(3,m1,mirep)-c_a_x(3,m2,mirep)

        if(i_pbc.eq.1.and.periodic_bonds_okay) then
          xij1=xij1-xlen*anint(xij1*xinv)
          xij2=xij2-ylen*anint(xij2*yinv)
          xij3=xij3-zlen*anint(xij3*zinv)
        endif

        rij=sqrt(xij1**2+xij2**2+xij3**2)

C if bond_dev_quit>0 and bond length deviates by more than the 
C threshold then we will send a message to quit later
C if i_look_for_clashes then we will write out 10 pdb files
C that are for consecutive steps leading to the clash 

        if(abs(rij-r0bb_new(iii)).gt.bond_dev_quit.and.
     &         bond_dev_quit.gt.0.0) then

!$omp atomic
          i_need_to_die_now=i_need_to_die_now+1

          write(*,185)tot_time,mirep,m1,m2
          write(*,184)r0bb_new(iii),rij
          write(*,183)
          write(*,186)m1,mirep,
     &        c_a_x(1,m1,mirep),c_a_x(2,m1,mirep),c_a_x(3,m1,mirep)
          write(*,186)m2,mirep,
     &        c_a_x(1,m2,mirep),c_a_x(2,m2,mirep),c_a_x(3,m2,mirep)
185       format(' at time = ',f20.5,' for replica # ',i6,
     &           ' the bond for atoms ',i8,' & ',i8,' has exploded')
184       format(' it should be ',f10.3,' but is ',f10.3,
     &           ' so i will quit ')
183       format(' if you want to carry',
     &           ' on instead, increase bond_dev_quit ',
     &           ' or set it to <0 for no monitoring')
186       format('coords ',i8,' replica # ',i8,3f15.5)

        endif

C do the following if we're using harmonic function for the bond

        fbond=rkbb_new(iii)*(rij-r0bb_new(iii))/rij
        ebnd(iiib_new(iii))=
     &  ebnd(iiib_new(iii))+
     &  0.500*rkbb_new(iii)*(rij-r0bb_new(iii))**2*dbnd(ijk)

        fx(1,m1)=fx(1,m1)-fbond*xij1
        fx(2,m1)=fx(2,m1)-fbond*xij2
        fx(3,m1)=fx(3,m1)-fbond*xij3
        fx(1,m2)=fx(1,m2)+fbond*xij1
        fx(2,m2)=fx(2,m2)+fbond*xij2
        fx(3,m2)=fx(3,m2)+fbond*xij3

      enddo

C [2] do the angles

      do ijk=1,nang

        iii=iang(ijk)
        m1=nba1_new(iii)
        m2=nba2_new(iii)
        m3=nba3_new(iii)

        xia=c_a_x(1,m1,mirep)
        yia=c_a_x(2,m1,mirep)
        zia=c_a_x(3,m1,mirep)
        xib=c_a_x(1,m2,mirep)
        yib=c_a_x(2,m2,mirep)
        zib=c_a_x(3,m2,mirep)
        xic=c_a_x(1,m3,mirep)
        yic=c_a_x(2,m3,mirep)
        zic=c_a_x(3,m3,mirep)
        force=rkba_new(iii)
        xab=xia-xib
        yab=yia-yib
        zab=zia-zib
        xcb=xic-xib
        ycb=yic-yib
        zcb=zic-zib
        if(i_pbc.eq.1.and.periodic_bonds_okay) then
          xab=xab-xlen*anint(xab*xinv)  ! PBC_NO
          yab=yab-ylen*anint(yab*yinv)  ! PBC_NO
          zab=zab-zlen*anint(zab*zinv)  ! PBC_NO
          xcb=xcb-xlen*anint(xcb*xinv)  ! PBC_NO
          ycb=ycb-ylen*anint(ycb*yinv)  ! PBC_NO
          zcb=zcb-zlen*anint(zcb*zinv)  ! PBC_NO
        endif
        rab2=xab*xab+yab*yab+zab*zab
        rcb2=xcb*xcb+ycb*ycb+zcb*zcb
        if(rab2.eq.0.000 .or. rcb2.eq.0.000) cycle
        xp1 = ycb*zab - zcb*yab
        yp1 = zcb*xab - xcb*zab
        zp1 = xcb*yab - ycb*xab
        rp = sqrt(xp1*xp1 + yp1*yp1 + zp1*zp1)

C place change limits on rp so that it is less likely to blow up
C - this is possibly unnecessary but this is an attempt to make
C the limits more similar to the dihedral limits

        rp = max(rp,0.0001)
        dot = xab*xcb + yab*ycb + zab*zcb
        cosine = dot / sqrt(rab2*rcb2)
        cosine = min(1.000,max(-1.000,cosine))
        angle = acos(cosine)

C do the following if we're using harmonic function for the angle

        dt = angle - r0ba_new(iii)
        dt2 = dt * dt
        eang(iiia_new(iii))=
     &  eang(iiia_new(iii))+0.500*force*dt2*dang(ijk)
        deddt = force * dt

        terma=-deddt/(rab2*rp)
        termc= deddt/(rcb2*rp)
        dedxia=terma*(yab*zp1-zab*yp1)
        dedyia=terma*(zab*xp1-xab*zp1)
        dedzia=terma*(xab*yp1-yab*xp1)
        dedxic=termc*(ycb*zp1-zcb*yp1)
        dedyic=termc*(zcb*xp1-xcb*zp1)
        dedzic=termc*(xcb*yp1-ycb*xp1)
        dedxib=-dedxia-dedxic
        dedyib=-dedyia-dedyic
        dedzib=-dedzia-dedzic
        fx(1,m1)=fx(1,m1)-dedxia
        fx(2,m1)=fx(2,m1)-dedyia
        fx(3,m1)=fx(3,m1)-dedzia
        fx(1,m2)=fx(1,m2)-dedxib
        fx(2,m2)=fx(2,m2)-dedyib
        fx(3,m2)=fx(3,m2)-dedzib
        fx(1,m3)=fx(1,m3)-dedxic
        fx(2,m3)=fx(2,m3)-dedyic
        fx(3,m3)=fx(3,m3)-dedzic
      enddo

C [3] do the dihedrals

      do ijk=1,ndih

        iii=idih(ijk)
        m1=nbd1_new(iii)
        m2=nbd2_new(iii)
        m3=nbd3_new(iii)
        m4=nbd4_new(iii)

        xia=c_a_x(1,m1,mirep)
        yia=c_a_x(2,m1,mirep)
        zia=c_a_x(3,m1,mirep)
        xib=c_a_x(1,m2,mirep)
        yib=c_a_x(2,m2,mirep)
        zib=c_a_x(3,m2,mirep)
        xic=c_a_x(1,m3,mirep)
        yic=c_a_x(2,m3,mirep)
        zic=c_a_x(3,m3,mirep)
        xid=c_a_x(1,m4,mirep)
        yid=c_a_x(2,m4,mirep)
        zid=c_a_x(3,m4,mirep)
        xba=xib-xia
        yba=yib-yia
        zba=zib-zia
        xcb=xic-xib
        ycb=yic-yib
        zcb=zic-zib
        xdc=xid-xic
        ydc=yid-yic
        zdc=zid-zic
        if(i_pbc.eq.1.and.periodic_bonds_okay) then
          xba=xba-xlen*anint(xba*xinv)  ! PBC_NO
          yba=yba-ylen*anint(yba*yinv)  ! PBC_NO
          zba=zba-zlen*anint(zba*zinv)  ! PBC_NO
          xcb=xcb-xlen*anint(xcb*xinv)  ! PBC_NO
          ycb=ycb-ylen*anint(ycb*yinv)  ! PBC_NO
          zcb=zcb-zlen*anint(zcb*zinv)  ! PBC_NO
          xdc=xdc-xlen*anint(xdc*xinv)  ! PBC_NO
          ydc=ydc-ylen*anint(ydc*yinv)  ! PBC_NO
          zdc=zdc-zlen*anint(zdc*zinv)  ! PBC_NO
        endif
        xt=yba*zcb-ycb*zba
        yt=zba*xcb-zcb*xba
        zt=xba*ycb-xcb*yba
        xu=ycb*zdc-ydc*zcb
        yu=zcb*xdc-zdc*xcb
        zu=xcb*ydc-xdc*ycb
        xtu=yt*zu-yu*zt
        ytu=zt*xu-zu*xt
        ztu=xt*yu-xu*yt
        rt2=xt*xt+yt*yt+zt*zt
        ru2=xu*xu+yu*yu+zu*zu

C skip this dihedral if angles are too close to being linear

        if(rt2.lt.0.001.or.
     &     ru2.lt.0.001) then
          cycle
        endif

        rtru=sqrt(rt2*ru2)
        if(rtru.eq.0.000) cycle ! redundant cos of above lines
        rcb=sqrt(xcb*xcb+ycb*ycb+zcb*zcb)
        cosine1=(xt*xu+yt*yu+zt*zu)/rtru
        sine1=(xcb*xtu+ycb*ytu+zcb*ztu)/(rcb*rtru)

C do the following if we're using cosine function for the dihedral

        v1=tors1_new(1,iii)
        c1=tors1_new(3,iii)
        s1=tors1_new(4,iii)
        v3=tors3_new(1,iii)
        c3=tors3_new(3,iii)
        s3=tors3_new(4,iii)
        cosine2=cosine1*cosine1-sine1*sine1
        sine2=2.000*cosine1*sine1
        cosine3=cosine1*cosine2-sine1*sine2
        sine3=cosine1*sine2+sine1*cosine2

C AHE notes that the energy here is not yet completely figured out - the
C forces are, I think, all correct but the energy for the period=3 term
C is definitely not - hopefully will fix this soon

C I'm pretty sure that the following two lines *only* need changing -
C the forces must be correct because we correctly sample around
C dihedrals - but we must need to correct for the fact that we are now
C reading in the equilibrium phi, not the phase angle... the
C relationships between these are:

C phase1 = phi-min1 + pi
C phase3 = phi-min3 + pi
C      = 3*phi-min1 + pi

C hypothesis is that the following are calculating:

C 1 + cos (NPhi - phi-minN)

C when they should be calculating:

C 1 + cos (NPhi - phaseN)

C if so, then let's first calculate:

C cos (NPhi - Phi-minN)
C then cos(above-pi) = cos(above)cos(pi) - sin(above)sin(pi)

C ORIGINAL CODE
C       phi1=1.000+(cosine1*c1+sine1*s1)
C       phi3=1.000+(cosine3*c3+sine3*s3)

C cos(3Phi - 3Phi0) = cos(3Phi)cos(3Phi0) + 
C    &                sin(3Phi)sin(3Phi0)

        cos_tmp1=(cosine1*c1+sine1*s1)
        sin_tmp1=(sine1*c1-cosine1*s1)
        phi1=1.000+(cos_tmp1*cos(pi)-sin_tmp1*sin(pi))
        cos_tmp3=(cosine3*c3+sine3*s3)
        sin_tmp3=(sine3*c3-cosine3*s3)
        phi3=1.000+(cos_tmp3*cos(pi)-sin_tmp3*sin(pi))
 
        dphi1=(cosine1*s1-sine1*c1)
        dphi3=3.000*(cosine3*s3-sine3*c3)
        dedphi=(v1*dphi1+v3*dphi3)
        edih(iiid_new(iii))=
     &  edih(iiid_new(iii))+
     &             ((v1*phi1)+(v3*phi3))*ddih(ijk)

        xxca=xic-xia
        yyca=yic-yia
        zzca=zic-zia
        xdb=xid-xib
        ydb=yid-yib
        zdb=zid-zib
        if(i_pbc.eq.1) then
          xdb=xdb-xlen*anint(xdb*xinv)  ! PBC_NO
          ydb=ydb-ylen*anint(ydb*yinv)  ! PBC_NO
          zdb=zdb-zlen*anint(zdb*zinv)  ! PBC_NO
          xxca=xxca-xlen*anint(xxca*xinv)  ! PBC_NO
          yyca=yyca-ylen*anint(yyca*yinv)  ! PBC_NO
          zzca=zzca-zlen*anint(zzca*zinv)  ! PBC_NO
        endif
        dedxt=dedphi*(yt*zcb-ycb*zt)/(rt2*rcb)
        dedyt=dedphi*(zt*xcb-zcb*xt)/(rt2*rcb)
        dedzt=dedphi*(xt*ycb-xcb*yt)/(rt2*rcb)
        dedxu=-dedphi*(yu*zcb-ycb*zu)/(ru2*rcb)
        dedyu=-dedphi*(zu*xcb-zcb*xu)/(ru2*rcb)
        dedzu=-dedphi*(xu*ycb-xcb*yu)/(ru2*rcb)
        dedxia=zcb*dedyt-ycb*dedzt
        dedyia=xcb*dedzt-zcb*dedxt
        dedzia=ycb*dedxt-xcb*dedyt
        dedxib=yyca*dedzt-zzca*dedyt+zdc*dedyu-ydc*dedzu
        dedyib=zzca*dedxt-xxca*dedzt+xdc*dedzu-zdc*dedxu
        dedzib=xxca*dedyt-yyca*dedxt+ydc*dedxu-xdc*dedyu
        dedxic=zba*dedyt-yba*dedzt+ydb*dedzu-zdb*dedyu
        dedyic=xba*dedzt-zba*dedxt+zdb*dedxu-xdb*dedzu
        dedzic=yba*dedxt-xba*dedyt+xdb*dedyu-ydb*dedxu
        dedxid=zcb*dedyu-ycb*dedzu
        dedyid=xcb*dedzu-zcb*dedxu
        dedzid=ycb*dedxu-xcb*dedyu
        fx(1,m1)=fx(1,m1)+dedxia             
        fx(2,m1)=fx(2,m1)+dedyia             
        fx(3,m1)=fx(3,m1)+dedzia             
        fx(1,m2)=fx(1,m2)+dedxib             
        fx(2,m2)=fx(2,m2)+dedyib             
        fx(3,m2)=fx(3,m2)+dedzib             
        fx(1,m3)=fx(1,m3)+dedxic             
        fx(2,m3)=fx(2,m3)+dedyic             
        fx(3,m3)=fx(3,m3)+dedzic             
        fx(1,m4)=fx(1,m4)+dedxid             
        fx(2,m4)=fx(2,m4)+dedyid             
        fx(3,m4)=fx(3,m4)+dedzid             
      enddo

C [4] do the position restraints

      do nnn=1,miv_num_loc
        do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
          if(i_restraint(iii).eq.1) then
            xij1=c_a_x(1,iii,mirep)-cx_restraint(iii)
            xij2=c_a_x(2,iii,mirep)-cy_restraint(iii)
            xij3=c_a_x(3,iii,mirep)-cz_restraint(iii)
            if(i_pbc.eq.1) then
              xij1=xij1-xlen*anint(xij1*xinv)  ! PBC_NO
              xij2=xij2-ylen*anint(xij2*yinv)  ! PBC_NO
              xij3=xij3-zlen*anint(xij3*zinv)  ! PBC_NO
            endif
            epos(mol_of_atm(iii))=epos(mol_of_atm(iii))+
     &        0.5*fx_restraint(iii)*xij1**2+
     &        0.5*fy_restraint(iii)*xij2**2+
     &        0.5*fz_restraint(iii)*xij3**2
            fx(1,iii)=fx(1,iii)-fx_restraint(iii)*xij1
            fx(2,iii)=fx(2,iii)-fy_restraint(iii)*xij2
            fx(3,iii)=fx(3,iii)-fz_restraint(iii)*xij3

C spherial can capsule forces NOW implemented

          elseif(i_restraint(iii).eq.2) then
            rad2=(c_a_x(1,iii,mirep)-cx_restraint(iii))**2+
     &           (c_a_x(2,iii,mirep)-cy_restraint(iii))**2+
     &           (c_a_x(3,iii,mirep)-cz_restraint(iii))**2
            rad=1.0/sqrt(rad2)
            rad=rad*r_restraint(iii)
            yij1=(c_a_x(1,iii,mirep)-cx_restraint(iii))*rad
            yij1=yij1+cx_restraint(iii)
            yij2=(c_a_x(2,iii,mirep)-cy_restraint(iii))*rad
            yij2=yij2+cy_restraint(iii)
            yij3=(c_a_x(3,iii,mirep)-cz_restraint(iii))*rad
            yij3=yij3+cz_restraint(iii)
            xij1=c_a_x(1,iii,mirep)-yij1
            xij2=c_a_x(2,iii,mirep)-yij2
            xij3=c_a_x(3,iii,mirep)-yij3
            fx(1,iii)=fx(1,iii)-fr_restraint(iii)*xij1
            fx(2,iii)=fx(2,iii)-fr_restraint(iii)*xij2
            fx(3,iii)=fx(3,iii)-fr_restraint(iii)*xij3

C capsule forces here - harmonic restraint to a desired radius

          elseif(i_restraint(iii).eq.3) then

C if the x coord is low then we're in the left-hand hemisphere

            if(c_a_x(1,iii,mirep).lt.cx_restraint(iii)) then
              rad2=(c_a_x(1,iii,mirep)-cx_restraint(iii))**2+
     &             (c_a_x(2,iii,mirep)-cy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-cz_restraint(iii))**2
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=(c_a_x(1,iii,mirep)-cx_restraint(iii))*rad
              yij1=yij1+cx_restraint(iii)
              yij2=(c_a_x(2,iii,mirep)-cy_restraint(iii))*rad
              yij2=yij2+cy_restraint(iii)
              yij3=(c_a_x(3,iii,mirep)-cz_restraint(iii))*rad
              yij3=yij3+cz_restraint(iii)
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
c             rij=sqrt(xij1**2+xij2**2+xij3**2)

C if the x coord is high then we're in the right-hand hemisphere

            elseif(c_a_x(1,iii,mirep).gt.dx_restraint(iii)) then
              rad2=(c_a_x(1,iii,mirep)-dx_restraint(iii))**2+
     &             (c_a_x(2,iii,mirep)-dy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-dz_restraint(iii))**2
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=(c_a_x(1,iii,mirep)-dx_restraint(iii))*rad
              yij1=yij1+dx_restraint(iii)
              yij2=(c_a_x(2,iii,mirep)-dy_restraint(iii))*rad
              yij2=yij2+dy_restraint(iii)
              yij3=(c_a_x(3,iii,mirep)-dz_restraint(iii))*rad
              yij3=yij3+dz_restraint(iii)
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
c             rij=sqrt(xij1**2+xij2**2+xij3**2)

C if the x coord is in-between the we're in the cylindrical part

            else
              rad2=(c_a_x(2,iii,mirep)-cy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-cz_restraint(iii))**2
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=c_a_x(1,iii,mirep)
              yij2=c_a_x(2,iii,mirep)*rad
              yij3=c_a_x(3,iii,mirep)*rad
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
c             rij=sqrt(xij1**2+xij2**2+xij3**2)
c             write(*,819)iii,c_a_x(1,iii,mirep),c_a_x(2,iii,mirep),
c    &                        c_a_x(3,iii,mirep),yij1,yij2,yij3,
c    &                        rij,sqrt(rad2),
c    &        cx_restraint(iii),cy_restraint(iii),cz_restraint(iii)
819           format('CHECK 1 ',i8,11f15.5)
            endif
            fx(1,iii)=fx(1,iii)-fr_restraint(iii)*xij1
            fx(2,iii)=fx(2,iii)-fr_restraint(iii)*xij2
            fx(3,iii)=fx(3,iii)-fr_restraint(iii)*xij3

C capsule forces - one-sided potentials here - wall is on outside

          elseif(i_restraint(iii).eq.4) then
            if(c_a_x(1,iii,mirep).lt.cx_restraint(iii)) then
              rad2=(c_a_x(1,iii,mirep)-cx_restraint(iii))**2+
     &             (c_a_x(2,iii,mirep)-cy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-cz_restraint(iii))**2
              if(rad2.lt.r_restraint(iii)**2) cycle
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=(c_a_x(1,iii,mirep)-cx_restraint(iii))*rad
              yij1=yij1+cx_restraint(iii)
              yij2=(c_a_x(2,iii,mirep)-cy_restraint(iii))*rad
              yij2=yij2+cy_restraint(iii)
              yij3=(c_a_x(3,iii,mirep)-cz_restraint(iii))*rad
              yij3=yij3+cz_restraint(iii)
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
c             rij=sqrt(xij1**2+xij2**2+xij3**2)

C if the x coord is high then we're in the right-hand hemisphere

            elseif(c_a_x(1,iii,mirep).gt.dx_restraint(iii)) then
              rad2=(c_a_x(1,iii,mirep)-dx_restraint(iii))**2+
     &             (c_a_x(2,iii,mirep)-dy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-dz_restraint(iii))**2
              if(rad2.lt.r_restraint(iii)**2) cycle
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=(c_a_x(1,iii,mirep)-dx_restraint(iii))*rad
              yij1=yij1+dx_restraint(iii)
              yij2=(c_a_x(2,iii,mirep)-dy_restraint(iii))*rad
              yij2=yij2+dy_restraint(iii)
              yij3=(c_a_x(3,iii,mirep)-dz_restraint(iii))*rad
              yij3=yij3+dz_restraint(iii)
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
c             rij=sqrt(xij1**2+xij2**2+xij3**2)

C if the x coord is in-between the we're in the cylindrical part

            else
              rad2=(c_a_x(2,iii,mirep)-cy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-cz_restraint(iii))**2
              if(rad2.lt.r_restraint(iii)**2) cycle
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=c_a_x(1,iii,mirep)
              yij2=c_a_x(2,iii,mirep)*rad
              yij3=c_a_x(3,iii,mirep)*rad
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
            endif
            fx(1,iii)=fx(1,iii)-fr_restraint(iii)*xij1
            fx(2,iii)=fx(2,iii)-fr_restraint(iii)*xij2
            fx(3,iii)=fx(3,iii)-fr_restraint(iii)*xij3

C capsule forces - one-sided potentials here - wall is on inside

          elseif(i_restraint(iii).eq.5) then
            if(c_a_x(1,iii,mirep).lt.cx_restraint(iii)) then
              rad2=(c_a_x(1,iii,mirep)-cx_restraint(iii))**2+
     &             (c_a_x(2,iii,mirep)-cy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-cz_restraint(iii))**2
              if(rad2.gt.r_restraint(iii)**2) cycle
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=(c_a_x(1,iii,mirep)-cx_restraint(iii))*rad
              yij1=yij1+cx_restraint(iii)
              yij2=(c_a_x(2,iii,mirep)-cy_restraint(iii))*rad
              yij2=yij2+cy_restraint(iii)
              yij3=(c_a_x(3,iii,mirep)-cz_restraint(iii))*rad
              yij3=yij3+cz_restraint(iii)
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
c             rij=sqrt(xij1**2+xij2**2+xij3**2)

C if the x coord is high then we're in the right-hand hemisphere

            elseif(c_a_x(1,iii,mirep).gt.dx_restraint(iii)) then
              rad2=(c_a_x(1,iii,mirep)-dx_restraint(iii))**2+
     &             (c_a_x(2,iii,mirep)-dy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-dz_restraint(iii))**2
              if(rad2.gt.r_restraint(iii)**2) cycle
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=(c_a_x(1,iii,mirep)-dx_restraint(iii))*rad
              yij1=yij1+dx_restraint(iii)
              yij2=(c_a_x(2,iii,mirep)-dy_restraint(iii))*rad
              yij2=yij2+dy_restraint(iii)
              yij3=(c_a_x(3,iii,mirep)-dz_restraint(iii))*rad
              yij3=yij3+dz_restraint(iii)
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
c             rij=sqrt(xij1**2+xij2**2+xij3**2)

C if the x coord is in-between the we're in the cylindrical part

            else
              rad2=(c_a_x(2,iii,mirep)-cy_restraint(iii))**2+
     &             (c_a_x(3,iii,mirep)-cz_restraint(iii))**2
              if(rad2.gt.r_restraint(iii)**2) cycle
              rad=1.0/sqrt(rad2)
              rad=rad*r_restraint(iii)
              yij1=c_a_x(1,iii,mirep)
              yij2=c_a_x(2,iii,mirep)*rad
              yij3=c_a_x(3,iii,mirep)*rad
              xij1=c_a_x(1,iii,mirep)-yij1
              xij2=c_a_x(2,iii,mirep)-yij2
              xij3=c_a_x(3,iii,mirep)-yij3
            endif
            fx(1,iii)=fx(1,iii)-fr_restraint(iii)*xij1
            fx(2,iii)=fx(2,iii)-fr_restraint(iii)*xij2
            fx(3,iii)=fx(3,iii)-fr_restraint(iii)*xij3
          endif
        enddo
      enddo

C [5] do wall terms

C note that +1 means wall is on the positive x-side of the atoms
C           -1 means wall is on the negative x-side of the atoms
C           +2 means wall is on the positive y-side of the atoms
C           -2 means wall is on the negative y-side of the atoms
C           +3 means wall is on the positive z-side of the atoms
C           -3 means wall is on the negative z-side of the atoms

C note that we could also add in other kinds of confinement forces, e.g.
C by putting in 0...

C we also accumulate the *absolute* value of the force on each wall
C for purposes of calculating the osmotic pressure

      if(walls) then
        do nnn=1,miv_num_loc
          do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
            ijk=my_uniq_atm_num(iii)
            do jkl=1,wall_ptr(ijk)%num
              n=wall_ptr(ijk)%iwl(jkl) ! get wall #
              if(crd_wall(n).eq.1) then
                if(c_a_x(1,iii,mirep).gt.pos_wall(n)) then
                  xij1=c_a_x(1,iii,mirep)-pos_wall(n)
                  fx(1,iii)=fx(1,iii)-fct_wall(n)*xij1
                  posm(1)=posm(1)+fct_wall(n)*xij1
                endif
              elseif(crd_wall(n).eq.2) then
                if(c_a_x(2,iii,mirep).gt.pos_wall(n)) then
                  xij1=c_a_x(2,iii,mirep)-pos_wall(n)
                  fx(2,iii)=fx(2,iii)-fct_wall(n)*xij1
                  posm(2)=posm(2)+fct_wall(n)*xij1
                endif
              elseif(crd_wall(n).eq.3) then 
                if(c_a_x(3,iii,mirep).gt.pos_wall(n)) then
                  xij1=c_a_x(3,iii,mirep)-pos_wall(n)
                  fx(3,iii)=fx(3,iii)-fct_wall(n)*xij1
                  posm(3)=posm(3)+fct_wall(n)*xij1
                endif
              elseif(crd_wall(n).eq.-1) then
                if(c_a_x(1,iii,mirep).lt.pos_wall(n)) then
                  xij1=c_a_x(1,iii,mirep)-pos_wall(n)
                  fx(1,iii)=fx(1,iii)-fct_wall(n)*xij1
                  posm(4)=posm(4)-fct_wall(n)*xij1
                endif
              elseif(crd_wall(n).eq.-2) then
                if(c_a_x(2,iii,mirep).lt.pos_wall(n)) then
                  xij1=c_a_x(2,iii,mirep)-pos_wall(n)
                  fx(2,iii)=fx(2,iii)-fct_wall(n)*xij1
                  posm(5)=posm(5)-fct_wall(n)*xij1
                endif
              elseif(crd_wall(n).eq.-3) then
                if(c_a_x(3,iii,mirep).lt.pos_wall(n)) then
                  xij1=c_a_x(3,iii,mirep)-pos_wall(n)
                  fx(3,iii)=fx(3,iii)-fct_wall(n)*xij1
                  posm(6)=posm(6)-fct_wall(n)*xij1
                endif
              endif
              ewll(mol_of_atm(iii))=ewll(mol_of_atm(iii))+
     &          0.5*fct_wall(n)*xij1**2
            enddo
          enddo
        enddo

      endif

C [8] do confinement terms

      if(sphere) then
        s_size2=s_size**2
        do nnn=1,miv_num_loc
          do iii=miv_beg_loc(nnn),miv_end_loc(nnn)
            rad2=c_a_x(1,iii,mirep)**2+
     &           c_a_x(2,iii,mirep)**2+
     &           c_a_x(3,iii,mirep)**2
            if(rad2.le.s_size2) cycle
            rad=1.0/sqrt(rad2)

C get coords of sphere surface along the center->atom vector
C note that we also accumulate force for osmotic pressure calc

            yij1=c_a_x(1,iii,mirep)*s_size*rad
            yij2=c_a_x(2,iii,mirep)*s_size*rad
            yij3=c_a_x(3,iii,mirep)*s_size*rad
            xij1=c_a_x(1,iii,mirep)-yij1
            xij2=c_a_x(2,iii,mirep)-yij2
            xij3=c_a_x(3,iii,mirep)-yij3
            fx(1,iii)=fx(1,iii)-f_size*xij1
            fx(2,iii)=fx(2,iii)-f_size*xij2
            fx(3,iii)=fx(3,iii)-f_size*xij3
            fmag=sqrt(f_size*f_size*(xij1**2+xij2**2+xij3**2))
            posm(7)=posm(7)+fmag
          enddo
        enddo
      endif

      if(cylinder) then
        s_size2=s_size**2
        m_size2=m_size**2
        m_sizeh=m_size*0.5
        do nnn=1,miv_num_loc
          do iii=miv_beg_loc(nnn),miv_end_loc(nnn)

C c_a_x(1,iii) is the x coord of atom iii
C c_a_x(2,iii) is the y coord of atom iii
C c_a_x(3,iii) is the z coord of atom iii
C
C yij1         is the x coord of the nearest point on the capsule
C              surface to atom iii
C yij2         is the y coord of the nearest point on the capsule
C              surface to atom iii
C yij3         is the z coord of the nearest point on the capsule
C              surface to atom iii

C if the x coord is low then we're in the left-hand hemisphere

            if(c_a_x(1,iii,mirep).lt.-m_sizeh) then
              rad2=(c_a_x(1,iii,mirep)+m_sizeh)**2+
     &              c_a_x(2,iii,mirep)**2+c_a_x(3,iii,mirep)**2
              if(rad2.le.s_size2) cycle
              rad=1.0/sqrt(rad2)
              yij1=(c_a_x(1,iii,mirep)+m_sizeh)*s_size*rad
              yij1=yij1-m_sizeh
              yij2=c_a_x(2,iii,mirep)*s_size*rad
              yij3=c_a_x(3,iii,mirep)*s_size*rad

C if the x coord is high then we're in the right-hand hemisphere

            elseif(c_a_x(1,iii,mirep).gt.m_sizeh) then
              rad2=(c_a_x(1,iii,mirep)-m_sizeh)**2+
     &              c_a_x(2,iii,mirep)**2+c_a_x(3,iii,mirep)**2
              if(rad2.le.s_size2) cycle
              rad=1.0/sqrt(rad2)
              yij1=(c_a_x(1,iii,mirep)-m_sizeh)*s_size*rad
              yij1=yij1+m_sizeh
              yij2=c_a_x(2,iii,mirep)*s_size*rad
              yij3=c_a_x(3,iii,mirep)*s_size*rad

C if the x coord is in-between the we're in the cylindrical part

            else
              rad2=c_a_x(2,iii,mirep)**2+c_a_x(3,iii,mirep)**2
              if(rad2.le.s_size2) cycle
              rad=1.0/sqrt(rad2)
              yij1=c_a_x(1,iii,mirep)
              yij2=c_a_x(2,iii,mirep)*s_size*rad
              yij3=c_a_x(3,iii,mirep)*s_size*rad
            endif
            xij1=c_a_x(1,iii,mirep)-yij1
            xij2=c_a_x(2,iii,mirep)-yij2
            xij3=c_a_x(3,iii,mirep)-yij3
            fx(1,iii)=fx(1,iii)-f_size*xij1
            fx(2,iii)=fx(2,iii)-f_size*xij2
            fx(3,iii)=fx(3,iii)-f_size*xij3
            fmag=sqrt(f_size*f_size*(xij1**2+xij2**2+xij3**2))
            posm(7)=posm(7)+fmag
c           write(*,999)iii,c_a_x(1,iii,mirep),
c    &                  c_a_x(2,iii,mirep),c_a_x(3,iii,mirep),
c    &                    fx(1,iii),fx(2,iii),fx(3,iii)
999         format('CYL FORCE ',i8,6f12.3)
          enddo
        enddo
      endif

      end


