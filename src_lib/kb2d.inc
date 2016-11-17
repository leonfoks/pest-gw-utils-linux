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
!                       2-D Kriging Program
!                       *******************
!
! The following Parameters control static dimensioning within okb2d:
!
!   MAXX      maximum nodes in X
!   MAXY      maximum nodes in Y
!   MAXDAT    maximum number of data points
!   MAXSAM    maximum number of data points to use in one kriging system
!   MAXDIS    maximum number of discretization points per block
!   MAXNST    maximum number of nested structures
!
!-----------------------------------------------------------------------
!
! User Adjustable:
!
      parameter(MAXDAT =10000,  &
                MAXSAM =  500,  &
                MAXDIS =   64,  &
                MAXNST =    5)
!
! Fixed
!
      parameter(MAXKD=MAXSAM+1,MAXKRG=MAXKD*MAXKD,UNEST=-999.,  &
                EPSLON=0.0000001,VERSION=2.000)
!
! Variable Declaration:
!
      real      x(MAXDAT),y(MAXDAT),vr(MAXDAT),aa(MAXNST),cc(MAXNST),  &
                ang(MAXNST),anis(MAXNST)
      integer   it(MAXNST)
!
! Common Blocks:
!
      common /datagv/ nd,x,y,vr

      common /parmfk/ nx,ny,xmn,ymn,xsiz,ysiz,nxdis,nydis,ndmin,    &
                      ndmax,radius,nst,c0,it,aa,cc,ang,anis,ktype,  &
                      skmean,idbg,lout,ldbg
 