!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
! The following Parameters control static dimensioning for SGSIM:
!
!   MAXX      maximum nodes in X
!   MAXY      maximum nodes in Y
!   MAXZ      maximum nodes in Z
!   MAXCTX    maximum X points in covariance table (ODD number)
!   MAXCTY    maximum Y points in covariance table (ODD number)
!   MAXCTZ    maximum Z points in covariance table (ODD number)
!   MAXSBX    maximum super block nodes in X direction
!   MAXSBY    maximum super block nodes in Y direction
!   MAXSBZ    maximum super block nodes in Z direction
!   MAXDAT    maximum number of data
!   MAXNOD    maximum number of previously simulated nodes to use
!   MAXSAM    maximum number of data for one simulation
!   MAXNST    maximum number of nested structures for variogram
!
!-----------------------------------------------------------------------
!
! User Adjustable Parameters:
!
      parameter(MAXX   =   700, MAXY   = 620, MAXZ   = 1,
     +          MAXCTX =    71, MAXCTY =  71, MAXCTZ = 21,
     +          MAXSBX =    21, MAXSBY =  21, MAXSBZ = 11,
     +          MAXDAT = 10000, MAXNOD =  100,MAXSAM = 100,
     +          MAXNST =     5  )
!
! Fixed Parameters:
!
      parameter(MAXKR1=MAXNOD+MAXSAM+1,MAXKR2=MAXKR1*MAXKR1,
     +          MAXROT=MAXNST+1,MAXCXY=MAXCTX*MAXCTY,
     +          MAXXYZ=MAXCTX*MAXCTY*MAXCTZ,MAXSB=MAXSBX*MAXSBY*MAXSBZ,
     +          MXYZ=MAXX*MAXY*MAXZ,UNEST=-99.0,EPSLON=1.0e-20,
     +          VERSION=2.000)
!
! ACORN parameters:
!
      parameter(KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
!
! Variable Declaration:
!
      real      x(MAXDAT),y(MAXDAT),z(MAXDAT),vr(MAXDAT),wt(MAXDAT),
     +          vrtr(MAXDAT),vrgtr(MAXDAT),close(MAXDAT),sec(MAXDAT),
     +          sim(MXYZ),lvm(MXYZ),tmp(MXYZ),order(MXYZ),
     +          c0(1),cc(MAXNST),aa(MAXNST),ang1(MAXNST),ang2(MAXNST),
     +          ang3(MAXNST),anis1(MAXNST),anis2(MAXNST),
     +          covtab(MAXCTX,MAXCTY,MAXCTZ),
     +          cnodex(MAXNOD),cnodey(MAXNOD),cnodez(MAXNOD),
     +          cnodev(MAXNOD),vra(MAXKR1),vrea(MAXKR1),ltpar
      real*8    rotmat(MAXROT,3,3),r(MAXKR1),rr(MAXKR1),s(MAXKR1),
     +          a(MAXKR2)
      integer   ltail,utail,sstrat,mults,nisb(MAXSB),nst(1),
     +          it(MAXNST),icnode(MAXNOD)
      integer*2 ixnode(MAXXYZ),iynode(MAXXYZ),iznode(MAXXYZ),
     +          ixsbtosr(8*MAXSB),iysbtosr(8*MAXSB),izsbtosr(8*MAXSB)
!
! Common Blocks:
!
      common /iaco/   ixv(MAXOP1)

      common /generl/ nx,ny,nz,xsiz,ysiz,zsiz,xmn,ymn,zmn,nxy,nxyz,
     +                nd,x,y,z,vr,wt,itrans,ntr,vrtr,vrgtr,sec,idbg,
     +                lin,lout,ldbg,llvm,icollvm,nvaril,ktype,colocorr,
     +                ltail,ltpar,utail,utpar,zmin,zmax,varred

      common /simula/ nsim,sim,lvm

      common /search/ radius,radsqd,sang1,sang2,sang3,sanis1,sanis2,
     +                noct,nclose,close,ndmin,ndmax,sstrat,mults,nmult

      common /cova3d/ nst,it,cmax,c0,cc,aa,ang1,ang2,ang3,
     +                anis1,anis2,isrot

      common /clooku/ nctx,ncty,nctz,covtab,nlooku,ncnode,icnode,
     +                cnodex,cnodey,cnodez,cnodev,nodmax,ixnode,iynode,
     +                iznode

      common /krigev/ rotmat,r,rr,s,a,vra,cbb
