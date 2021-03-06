!     Last change:  J    22 Apr 2002    4:44 pm

      subroutine readparm(ncall,k_ktype,r_radius,n_nx,x_xmn,x_xsiz,
     +n_ny,y_ymn,y_ysiz,n_ndmin,n_ndmax,n_ncnode,n_nst,c_c0,i_it,c_cc,
     +a_ang,a_aa,a_anis,ncond,n_nd,eastdat,northdat,valdat,ss,iseed)

      implicit integer(i-n), real(a-h, o-z)                    !jd

      include  'sgsim.inc'
      real      var(50)
      real*8    p,acorni,cp,oldcp,w
      logical   testfl

      integer i_it(n_nst)
      real valdat(n_nd),eastdat(n_nd),northdat(n_nd)
      real c_cc(n_nst),a_ang(n_nst),a_aa(n_nst),a_anis(n_nst)


! -- Some checks are made on array dimensions.

      call comsize(n_nx,MAXX,'MAXX')
      call comsize(n_ny,MAXY,'MAXY')
      call comsize(n_nd,MAXDAT,'MAXDAT')


      ixl=0
      iyl=0
      izl=0
      ivrl=0
      iwt=0
      isecvr=0

      tmin=-1.0e30
      tmax= 1.0e30

      itrans=0
      ismooth=0
      isvr=0
      iswt=0
      zmin=-1.030
      zmax=1.0e30
      ltail=1
      ltpar=0.0
      utail=1
      utpar=0.0

      ldbg=55                                        !debug
      idbg=0                                         !debug

      nsim=1

      nx=n_nx
      xmn=x_xmn
      xsiz=x_xsiz
      ny=n_ny
      ymn=y_ymn
      ysiz=y_ysiz
      nz=1
      zmn=0.0
      zsiz=sqrt(xsiz*ysiz)
      nxy  = nx*ny
      nxyz = nx*ny*nz

      if(ncall.eq.1) then
!        ixv(1)=3248593
!        ixv(1)=6984956
         ixv(1)=iseed
        do i=1,1000
           p = acorni(idum)
        end do
      end if

      ndmin= n_ndmin
      ndmax= n_ndmax
      nodmax=n_ncnode
      sstrat=0
      if(sstrat.eq.1) ndmax = 0

      mults=0
      nmult=0

      noct=0
      radius=r_radius
      radius1=r_radius
      radius2=r_radius
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      sang1=0.0
      sang2=0.0
      sang3=0.0

      ktype=k_ktype
      colocorr = 0.0
      icollvm=0

      nst(1)=n_nst
      c0(1)=c_c0
      sill = c0(1)

      do i=1,nst(1)
         it(i)=i_it(i)
         cc(i)=c_cc(i)
         ang1(i)=a_ang(i)
         ang2(i)=0.0
         ang3(i)=0.0
         aa(i)=a_aa(i)
         anis1(i)=a_anis(i)
         anis2(i)=0.0
         sill     = sill + cc(i)
         if(it(i).eq.4) then
              write(*,*) ' A power model is NOT allowed '
              write(*,*) ' Choose a different model and re start '
              stop
         endif
      end do

!
! Warn the user if the sill is different than 1.0:
!
!jd      if(sill.gt.(1.0+EPSLON).or.sill.lt.(1.0-EPSLON)) then
!jd            write(*,*) 'WARNING the sill of your variogram is not 1.0!'
!jd            write(*,*) '        the sill = ',sill
!jd            write(*,*)
!jd      end if
      nd = 0
      av = 0.0
      ss = 0.0
      if(ncond.eq.0)then
            ndmin  = 0
            ndmax  = 0
            sstrat = 1
      else
            twt = 0.0
            nd  = 0
            nt  = 0
            do idat=1,ncond
              nd = nd + 1
              vr(nd)=valdat(idat)
              x(nd)=eastdat(idat)
              y(nd)=northdat(idat)
              z(nd) = zmn
              wt(nd) = 1.0
              sec(nd) = UNEST
              twt = twt + wt(nd)
              av  = av  + vr(nd)*wt(nd)
              ss  = ss  + vr(nd)*vr(nd)*wt(nd)
            end do
!
! Compute the averages and variances as an error check for the user:
!
            av = av / max(twt,EPSLON)
            ss =(ss / max(twt,EPSLON)) - av * av
!jd            write(ldbg,111) nd,nt,av,ss
!jd 111  format(/,' Data for SGSIM: Number of acceptable data  = ',i8,/,
!jd     +         '                 Number trimmed             = ',i8,/,
!jd     +         '                 Weighted Average           = ',f12.4,/,
!jd     +         '                 Weighted Variance          = ',f12.4,/)
            ss=sqrt(ss)                                                  !jd
      endif

      return

      end



      subroutine sgsim(maxtempdim,temparray,avval,astdev)
!-----------------------------------------------------------------------
!
!           Conditional Simulation of a 3-D Rectangular Grid
!           ************************************************
!
! This subroutine generates 3-D realizations of a Gaussian process with
! a given autocovariance model, and conditional to input Gaussian data.
! The conditional simulation is achieved by sequential simulation of all
! the nodes visited by a random path.
!
!
!
! PROGRAM NOTES:
!
!  1. The three dimensional anisotropy parameters, i.e., of the search
!     ellipse and variogram ranges are described in section 2.3 of the
!     manual.   The variogram parameters are described in the same place
!
!  2. The original data and previously simulated grid nodes can be
!     searched separately.  There can be a different maximum number of
!     each and a minimum number of original data can be specified
!     to restrict simulation beyond the limits of the data.  The
!     closeness of previously simulated grid nodes is measured according
!     to the variogram structural distance.
!
!
!
! INPUT VARIABLES:
!
!   nd               Number of data (no missing values)
!   x,y,z(nd)        coordinates of the data
!   vr(nd)           gaussian data (normal scores)
!
!   nx,ny,nz         Number of blocks in X,Y, and Z
!   xmn,ymn,zmn      Coordinate at the center of the first Block
!   xsiz,ysiz,zsiz   spacing of the grid nodes (block size)
!
!   nsim             number of simulations
!   ktype            =1, ordinary kriging; =0, simple kriging
!   sim              the current realization
!   idbg             integer debugging level (0=none,2=normal,4=serious)
!   ldbg             unit number for the debugging output
!   lout             unit number for the output
!
!   radius           Maximum search radius
!   sang1            Azimuth angle of the principal search direction
!   sang2            Dip angle of the principal search direction
!   sang3            Third rotation angle of the search ellipse
!   sanis1           Anisotropy for the dip angle
!   sanis2           Anisotropy for the plunge angle
!   ndmin            Minimum number of data required before sim
!   ndmax            Maximum number of samples for simulation
!   noct             Maximum number per octant if an octant search is
!                      desired (if <= 0, then no octant search)
!
!   nodmax           Maximum number of previously simulated grid nodes
!                      to consider in the simulation.  The structural
!                      variogram distance is used to identify close ones
!
!   c0               Nugget constant (isotropic).
!   cc(nst)          Multiplicative factor of each nested structure.
!   aa(nst)          Parameter "a" of each nested structure.
!   it(nst)          Type of nested structures (1=sph,2=exp,3=gau,4=pow)
!   ang1(nst)        Azimuth angle for the principal direction
!   ang2(nst)        Dip angle for the principal direction
!   ang3(nst)        Third rotation angle to rotate the two minor
!                      directions around the principal direction
!   anis1(nst)       Anisotropy (radius in minor direction at 90
!                      degrees from "ang1" divided by the principal
!                      radius in direction "ang1")
!   anis2(nst)       Anisotropy (radius in minor direction at 90 degrees
!                      vertical from "ang1" divided by the principal
!                      radius in direction "ang1")
!
!
! OUTPUT VARIABLES:  Simulated Values are written to "lout"
!
!
!
! EXTERNAL REFERENCES:
!
!   super            Sets up the super block search of original data
!   search           Search for nearby data values
!   ctable           Builds a covariance table and "spiral" search
!   srchnd           Search for nearby simulated grid nodes
!   sqdist           computes anisotropic squared distance
!   sortem           sorts multiple arrays in ascending order (separate)
!   cova3            Calculates the covariance given a variogram model
!   krige            Sets up and solves either the SK or OK system
!   ksol             Linear system solver using Gaussian elimination
!
!
!
! Concepts taken from F. Alabert and E. Isaaks
!
!-----------------------------------------------------------------------
      include  'sgsim.inc'
      real      randnu(1),var(10)
      real*8    p,acorni,cp,oldcp,w
      logical   testind
      integer max_nx,max_ny
      real temparray(maxtempdim)
!
! Set up the rotation/anisotropy matrices that are needed for the
! variogram and search.
!
!jd      write(*,*) 'Setting up rotation matrices for variogram and search'
      do is=1,nst(1)
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),
     +                  is,MAXROT,rotmat)
      end do
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
!
! Set up the super block search:
!
      if(sstrat.eq.0) then
!jd            write(*,*) 'Setting up super block search strategy'
            nsec = 1
            call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
     +                   vr,wt,nsec,sec,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,
     +                   nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +                   nzsup,zmnsup,zsizsup)
            call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
     +                   isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
     +                   iysbtosr,izsbtosr)
      end if
!
! Set up the covariance table and the spiral search:
!
      call ctable
!
! MAIN LOOP OVER ALL THE SIMULAUTIONS:
!
      do isim=1,nsim
!
! Read in the secondary data distribution for this realization:
!
            if(isim.gt.1.and.ktype.eq.4) then
!jd                  write(*,*)
!jd                  write(*,*) ' Reading next secondary model'
                  index = 0
                  do iz=1,nz
                        do iy=1,ny
                              do ix=1,nx
                                 index = index + 1
                                 read(llvm,*,end=977)(var(j),j=1,nvaril)
                                 lvm(index) = var(icollvm)
                                 sim(index) = real(index)
                              end do
                        end do
                  end do
!jd                  write(*,*) ' Building CDF from  secondary model'
                  call sortem(1,nxyz,lvm,1,sim,c,d,e,f,g,h)
                  oldcp = 0.0
                  cp    = 0.0
                  do i=1,nxyz
                        cp =  cp + dble(1.0/real(nxyz))
                        w  = (cp + oldcp)/2.0
                        call gauinv(w,lvm(i),ierr)
                        lvm(i) = lvm(i) * varred
                        oldcp  =  cp
                  end do
!jd                  write(*,*) ' Restoring order of secondary model'
                  call sortem(1,nxyz,sim,1,lvm,c,d,e,f,g,h)
 977              continue
!jd                  write(*,*)
            end if
!
! Work out a random path for this realization:
!
            do ind=1,nxyz
                  sim(ind)   = real(acorni(idum))
                  order(ind) = ind
            end do
!
! The multiple grid search works with multiples of 4 (yes, that is
! somewhat arbitrary):
!
            if(mults.eq.1) then
                  do imult=1,nmult
                        nnz = max(1,nz/(imult*4))
                        nny = max(1,ny/(imult*4))
                        nnx = max(1,nx/(imult*4))
                        jz  = 1
                        jy  = 1
                        jx  = 1
                        do iz=1,nnz
                           if(nnz.gt.1) jz = iz*imult*4
                           do iy=1,nny
                              if(nny.gt.1) jy = iy*imult*4
                              do ix=1,nnx
                                 if(nnx.gt.1) jx = ix*imult*4
                                 index = jx + (jy-1)*nx + (jz-1)*nxy
                                 sim(index) = sim(index) + imult
                              end do
                           end do
                        end do
                  end do
            end if
            call sortem(1,nxyz,sim,1,order,c,d,e,f,g,h)
!
! Initialize the simulation:
!
            do ind=1,nxyz
                  sim(ind) = UNEST
            end do
!jd            write(*,*)
!jd            write(*,*) 'Working on realization number ',isim
!
! Assign the data to the closest grid node:
!
            TINY = 0.0001
            do id=1,nd
                  call getindx(nx,xmn,xsiz,x(id),ix,testind)
                  call getindx(ny,ymn,ysiz,y(id),iy,testind)
                  call getindx(nz,zmn,zsiz,z(id),iz,testind)
                  ind = ix + (iy-1)*nx + (iz-1)*nxy
                  xx  = xmn + real(ix-1)*xsiz
                  yy  = ymn + real(iy-1)*ysiz
                  zz  = zmn + real(iz-1)*zsiz
                  test = abs(xx-x(id)) + abs(yy-y(id)) + abs(zz-z(id))
!
! Assign this data to the node (unless there is a closer data):
!
                  if(sstrat.eq.1) then
                        if(sim(ind).ge.0.0) then
                              id2 = int(sim(ind)+0.5)
                              test2 = abs(xx-x(id2)) + abs(yy-y(id2))
     +                                               + abs(zz-z(id2))
                              if(test.le.test2) sim(ind) = real(id)
!jd                              write(ldbg,102) id,id2
                        else
                              sim(ind) = real(id)
                        end if
                  end if
!
! Assign a flag so that this node does not get simulated:
!
                  if(sstrat.eq.0.and.test.le.TINY) sim(ind)=10.0*UNEST
            end do
 102        format(' WARNING data values ',2i5,' are both assigned to ',
     +           /,'         the same node - taking the closest')
!
! Now, enter data values into the simulated grid:
!
            do ind=1,nxyz
                  id = int(sim(ind)+0.5)
                  if(id.gt.0) sim(ind) = vr(id)
            end do
            irepo = max(1,min((nxyz/10),10000))
!
! MAIN LOOP OVER ALL THE NODES:
!
            do in=1,nxyz
!jd                  if((int(in/irepo)*irepo).eq.in) write(*,103) in
!jd 103              format('   currently on node ',i9)
!
! Figure out the location of this point and make sure it has
! not been assigned a value already:
!
                  index = int(order(in)+0.5)
                  if(sim(index).gt.(UNEST+EPSLON).or.
     +               sim(index).lt.(UNEST*2.0)) go to 5
                  iz = int((index-1)/nxy) + 1
                  iy = int((index-(iz-1)*nxy-1)/nx) + 1
                  ix = index - (iz-1)*nxy - (iy-1)*nx
                  xx = xmn + real(ix-1)*xsiz
                  yy = ymn + real(iy-1)*ysiz
                  zz = zmn + real(iz-1)*zsiz
!
! Now, we'll simulate the point ix,iy,iz.  First, get the close data
! and make sure that there are enough to actually simulate a value,
! we'll only keep the closest "ndmax" data, and look for previously
! simulated grid nodes:
!
                  if(sstrat.eq.0) then
                        call srchsupr(xx,yy,zz,radsqd,isrot,MAXROT,
     +                          rotmat,nsbtosr,ixsbtosr,iysbtosr,
     +                          izsbtosr,noct,nd,x,y,z,wt,nisb,nxsup,
     +                          xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +                          nzsup,zmnsup,zsizsup,nclose,close,
     +                          infoct)
                        if(nclose.lt.ndmin) go to 5
                        if(nclose.gt.ndmax) nclose = ndmax
                  endif
                  call srchnd(ix,iy,iz)
!
! Calculate the conditional mean and standard deviation.  This will be
! done with kriging if there are data, otherwise, the global mean and
! standard deviation will be used:
!
                  if(ktype.eq.2) then
                        gmean = lvm(index)
                  else
!jd                        gmean = 0.0
                        gmean = avval                !jd
                  end if
                  if((nclose+ncnode).lt.1) then
                        cmean  = gmean
!jd                        cstdev = 1.0
                        cstdev=astdev                !jd
                  else
!
! Perform the kriging.  Note that if there are fewer than four data
! then simple kriging is prefered so that the variance of the
! realization does not become artificially inflated:
!
                        lktype = ktype
                        if(ktype.eq.1.and.(nclose+ncnode).lt.4)lktype=0
                        call krige(ix,iy,iz,xx,yy,zz,lktype,
     +                             gmean,cmean,cstdev,astdev)   !jd
                  endif
!
! Draw a random number and assign a value to this node:
!
                  p = acorni(idum)
                  call gauinv(p,xp,ierr)
                  sim(index) = xp * cstdev + cmean
                  if(idbg.ge.3) write(ldbg,141) p,sim(index)
 141              format(' random number ',f6.4,' realization ',f7.4)
!
! Quick check for far out results:
!
!jd                  if(abs(cmean).gt.5.0.or.abs(cstdev).gt.5.0.or.
!jd     +               abs(sim(index)).gt.6.0) then
!jd                  write(ldbg,104) ix,iy,iz,cmean,cstdev,sim(index)
!jd  104             format('WARNING: grid node location: ',3i5,/,
!jd     +                   '         conditional mean:   ',f12.5,/,
!jd     +                   '         conditional stdev:  ',f12.5,/,
!jd     +                   '         simulated value:    ',f12.5)
!jd                  endif
!
! END MAIN LOOP OVER NODES:
!
 5                continue
            end do
!
! Do we need to reassign the data to the grid nodes?
!
            if(sstrat.eq.0) then
                  do id=1,nd
                        call getindx(nx,xmn,xsiz,x(id),ix,testind)
                        call getindx(ny,ymn,ysiz,y(id),iy,testind)
                        call getindx(nz,zmn,zsiz,z(id),iz,testind)
                        xx  = xmn + real(ix-1)*xsiz
                        yy  = ymn + real(iy-1)*ysiz
                        zz  = zmn + real(iz-1)*zsiz
                        ind = ix + (iy-1)*nx + (iz-1)*nxy
                        test=abs(xx-x(id))+abs(yy-y(id))+abs(zz-z(id))
                        if(test.le.TINY) sim(ind) = vr(id)
                  end do
            end if
!
! Back transform each value and write results:
!
            ne = 0
            av = 0.0
            ss = 0.0
            do ind=1,nxyz
                  simval = sim(ind)
!jd                  if(simval.gt.-9.0.and.simval.lt.9.0) then
                        ne = ne + 1
                        av = av + simval
                        ss = ss + simval*simval
!jd                  end if
                  if(itrans.eq.1.and.simval.gt.(UNEST+EPSLON)) then
                        simval = backtr(simval,ntr,vrtr,vrgtr,zmin,
     +                                  zmax,ltail,ltpar,utail,utpar)
                        if(simval.lt.zmin) simval = zmin
                        if(simval.gt.zmax) simval = zmax
                  end if
!jd                  write(lout,'(f12.4)') simval
!jd                  write(lout,*) simval                !jd
            end do
            av = av / max(real(ne),1.0)
            ss =(ss / max(real(ne),1.0)) - av * av
!jd            write(ldbg,111) isim,ne,av,ss
!jd            write(*,   111) isim,ne,av,ss             !field stats  !jd
 111        format(/,' Realization ',i3,': number   = ',i8,/,
     +               '                  mean     = ',f12.4,
     +               ' (close to 0.0?)',/,
     +               '                  variance = ',f12.4,
     +               ' (close to gammabar(V,V)? approx. 1.0)',/)
!
! END MAIN LOOP OVER SIMULATIONS:
!
      end do
!
! Return to the main program:
!

      do ind=1,nxyz
        temparray(ind)=sim(ind)
      end do

      return
      end



      subroutine ctable
!-----------------------------------------------------------------------
!
!               Establish the Covariance Look up Table
!               **************************************
!
! The idea is to establish a 3-D network that contains the covariance
! value for a range of grid node offsets that should be at as large
! as twice the search radius in each direction.  The reason it has to
! be twice as large as the search radius is because we want to use it
! to compute the data covariance matrix as well as the data-point
! covariance matrix.
!
! Secondly, we want to establish a search for nearby nodes that
! in order of closeness as defined by the variogram.
!
!
!
! INPUT VARIABLES:
!
!   xsiz,ysiz,zsiz  Definition of the grid being considered
!   MAXCTX,Y,Z      Number of blocks in covariance table
!
!   covariance table parameters
!
!
!
! OUTPUT VARIABLES:  covtab()         Covariance table
!
! EXTERNAL REFERENCES:
!
!   sqdist          Computes 3-D anisotropic squared distance
!   sortem          Sorts multiple arrays in ascending order
!   cova3           Computes the covariance according to a 3-D model
!
!
!
!-----------------------------------------------------------------------
      parameter(TINY=1.0e-10)
      include  'sgsim.inc'
      real*8    hsqd,sqdist
      logical   first
!
! Size of the look-up table:
!
      nctx = min(((MAXCTX-1)/2),(nx-1))
      ncty = min(((MAXCTY-1)/2),(ny-1))
      nctz = min(((MAXCTZ-1)/2),(nz-1))
!
! Debugging output:
!
!jd      write(ldbg,*)
!jd      write(ldbg,*) 'Covariance Look up table and search for previously'
!jd      write(ldbg,*) 'simulated grid nodes.  The maximum range in each '
!jd      write(ldbg,*) 'coordinate direction for covariance look up is:'
!jd      write(ldbg,*) '          X direction: ',nctx*xsiz
!jd      write(ldbg,*) '          Y direction: ',ncty*ysiz
!jd      write(ldbg,*) '          Z direction: ',nctz*zsiz
!jd      write(ldbg,*) 'Node Values are not searched beyond this distance!'
!jd      write(ldbg,*)
!
! NOTE: If dynamically allocating memory, and if there is no shortage
!       it would a good idea to go at least as far as the radius and
!       twice that far if you wanted to be sure that all covariances
!       in the left hand covariance matrix are within the table look-up.
!
! Initialize the covariance subroutine and cbb at the same time:
!
      call cova3(0.0,0.0,0.0,0.0,0.0,0.0,1,nst,MAXNST,c0,it,cc,aa,
     +           1,MAXROT,rotmat,cmax,cbb)
!
! Now, set up the table and keep track of the node offsets that are
! within the search radius:
!
      nlooku = 0
      do i=-nctx,nctx
      xx = i * xsiz
      ic = nctx + 1 + i
      do j=-ncty,ncty
      yy = j * ysiz
      jc = ncty + 1 + j
      do k=-nctz,nctz
      zz = k * zsiz
      kc = nctz + 1 + k
            call cova3(0.0,0.0,0.0,xx,yy,zz,1,nst,MAXNST,c0,it,cc,aa,
     +                 1,MAXROT,rotmat,cmax,covtab(ic,jc,kc))
            hsqd = sqdist(0.0,0.0,0.0,xx,yy,zz,isrot,
     +                          MAXROT,rotmat)
            if(real(hsqd).le.radsqd) then
                  nlooku         = nlooku + 1
!
! We want to search by closest variogram distance (and use the
! anisotropic Euclidean distance to break ties:
!
                  tmp(nlooku)   = - (covtab(ic,jc,kc)-TINY*real(hsqd))
                  order(nlooku) = real((kc-1)*MAXCXY+(jc-1)*MAXCTX+ic)
            endif
      end do
      end do
      end do
!
! Finished setting up the look-up table, now order the nodes such
! that the closest ones, according to variogram distance, are searched
! first. Note: the "loc" array is used because I didn't want to make
! special allowance for 2 byte integers in the sorting subroutine:
!
      call sortem(1,nlooku,tmp,1,order,c,d,e,f,g,h)
      do il=1,nlooku
            loc = int(order(il))
            iz  = int((loc-1)/MAXCXY) + 1
            iy  = int((loc-(iz-1)*MAXCXY-1)/MAXCTX) + 1
            ix  = loc-(iz-1)*MAXCXY - (iy-1)*MAXCTX
            iznode(il) = int(iz)
            iynode(il) = int(iy)
            ixnode(il) = int(ix)
      end do
      if(nodmax.gt.MAXNOD) then
!jd            write(ldbg,*)
!jd            write(ldbg,*) 'The maximum number of close nodes = ',nodmax
!jd            write(ldbg,*) 'this was reset from your specification due '
!jd            write(ldbg,*) 'to storage limitations.'
            nodmax = MAXNOD
      endif
!
! Debugging output if requested:
!
      if(idbg.lt.2) return
      write(ldbg,*)
      write(ldbg,*) 'There are ',nlooku,' nearby nodes that will be '
      write(ldbg,*) 'checked until enough close data are found.'
      write(ldbg,*)
      if(idbg.lt.14) return
      do i=1,nlooku
            xx = (ixnode(i) - nctx - 1) * xsiz
            yy = (iynode(i) - ncty - 1) * ysiz
            zz = (iznode(i) - nctz - 1) * zsiz
            write(ldbg,100) i,xx,yy,zz
      end do
 100  format('Point ',i3,' at ',3f12.4)
!
! All finished:
!
      return
      end



      subroutine srchnd(ix,iy,iz)
!-----------------------------------------------------------------------
!
!               Search for nearby Simulated Grid nodes
!               **************************************
!
! The idea is to spiral away from the node being simulated and note all
! the nearby nodes that have been simulated.
!
!
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   sim             the realization so far
!   nodmax          the maximum number of nodes that we want
!   nlooku          the number of nodes in the look up table
!   i[x,y,z]node    the relative indices of those nodes.
!   [x,y,z]mn       the origin of the global grid netwrok
!   [x,y,z]siz      the spacing of the grid nodes.
!
!
!
! OUTPUT VARIABLES:
!
!   ncnode          the number of close nodes
!   icnode()        the number in the look up table
!   cnode[x,y,z]()  the location of the nodes
!   cnodev()        the values at the nodes
!
!
!
!-----------------------------------------------------------------------
      include  'sgsim.inc'
      integer   ninoct(8)
!
! Consider all the nearby nodes until enough have been found:
!
      ncnode = 0
      if(noct.gt.0) then
            do i=1,8
                  ninoct(i) = 0
            end do
      end if
      do 2 il=2,nlooku
            if(ncnode.eq.nodmax) return
            i = ix + (int(ixnode(il))-nctx-1)
            j = iy + (int(iynode(il))-ncty-1)
            k = iz + (int(iznode(il))-nctz-1)
            if(i.lt. 1.or.j.lt. 1.or.k.lt. 1) go to 2
            if(i.gt.nx.or.j.gt.ny.or.k.gt.nz) go to 2
            ind = i + (j-1)*nx + (k-1)*nxy
            if(sim(ind).gt.UNEST) then
!
! Check the number of data already taken from this octant:
!
                  if(noct.gt.0) then
                        idx = ix - i
                        idy = iy - j
                        idz = iz - k
                        if(idz.gt.0) then
                              iq = 4
                              if(idx.le.0 .and. idy.gt.0) iq = 1
                              if(idx.gt.0 .and. idy.ge.0) iq = 2
                              if(idx.lt.0 .and. idy.le.0) iq = 3
                        else
                              iq = 8
                              if(idx.le.0 .and. idy.gt.0) iq = 5
                              if(idx.gt.0 .and. idy.ge.0) iq = 6
                              if(idx.lt.0 .and. idy.le.0) iq = 7
                        end if
                        ninoct(iq) = ninoct(iq) + 1
                        if(ninoct(iq).gt.noct) go to 2
                  end if
                  ncnode = ncnode + 1
                  icnode(ncnode) = il
                  cnodex(ncnode) = xmn + real(i-1)*xsiz
                  cnodey(ncnode) = ymn + real(j-1)*ysiz
                  cnodez(ncnode) = zmn + real(k-1)*zsiz
                  cnodev(ncnode) = sim(ind)
            endif
 2    continue
!
! Return to calling program:
!
      return
      end



      subroutine krige(ix,iy,iz,xx,yy,zz,lktype,gmean,cmean,cstdev,
     +astdev)                                                        !jd
!-----------------------------------------------------------------------
!
!            Builds and Solves the SK or OK Kriging System
!            *********************************************
!
! INPUT VARIABLES:
!
!   ix,iy,iz        index of the point currently being simulated
!   xx,yy,zz        location of the point currently being simulated
!
!
!
! OUTPUT VARIABLES:
!
!   cmean           kriged estimate
!   cstdev          kriged standard deviation
!
!
!
! EXTERNAL REFERENCES: ksol   Gaussian elimination system solution
!
!
!
!-----------------------------------------------------------------------
      include 'sgsim.inc'
      logical first
!
! Size of the kriging system:
!
      first = .false.
      na    = nclose + ncnode
      if(lktype.eq.0) neq = na
      if(lktype.eq.1) neq = na + 1
      if(lktype.eq.2) neq = na
      if(lktype.eq.3) neq = na + 2
      if(lktype.eq.4) neq = na + 1
!
! Set up kriging matrices:
!
      in=0
      do j=1,na
!
! Sort out the actual location of point "j"
!
            if(j.le.nclose) then
                  index  = int(close(j))
                  x1     = x(index)
                  y1     = y(index)
                  z1     = z(index)
                  vra(j) = vr(index)
                  vrea(j)= sec(index)
                  if(lktype.eq.2) vra(j) = vra(j) - vrea(j)
            else
!
! It is a previously simulated node (keep index for table look-up):
!
                  index  = j-nclose
                  x1     = cnodex(index)
                  y1     = cnodey(index)
                  z1     = cnodez(index)
                  vra(j) = cnodev(index)
                  ind    = icnode(index)
                  ix1    = ix + (int(ixnode(ind))-nctx-1)
                  iy1    = iy + (int(iynode(ind))-ncty-1)
                  iz1    = iz + (int(iznode(ind))-nctz-1)
                  index  = ix1 + (iy1-1)*nx + (iz1-1)*nxy
                  vrea(j)= lvm(index)
                  if(lktype.eq.2) vra(j) = vra(j) - vrea(j)
            endif
            do i=1,j
!
! Sort out the actual location of point "i"
!
                  if(i.le.nclose) then
                        index  = int(close(i))
                        x2     = x(index)
                        y2     = y(index)
                        z2     = z(index)
                  else
!
! It is a previously simulated node (keep index for table look-up):
!
                        index  = i-nclose
                        x2     = cnodex(index)
                        y2     = cnodey(index)
                        z2     = cnodez(index)
                        ind    = icnode(index)
                        ix2    = ix + (int(ixnode(ind))-nctx-1)
                        iy2    = iy + (int(iynode(ind))-ncty-1)
                        iz2    = iz + (int(iznode(ind))-nctz-1)
                  endif
!
! Now, get the covariance value:
!
                  in = in + 1
!
! Decide whether or not to use the covariance look-up table:
!
                  if(j.le.nclose.or.i.le.nclose) then
                        call cova3(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,c0,it,
     +                             cc,aa,1,MAXROT,rotmat,cmax,cov)
                        a(in) = dble(cov)
                  else
!
! Try to use the covariance look-up (if the distance is in range):
!
                        ii = nctx + 1 + (ix1 - ix2)
                        jj = ncty + 1 + (iy1 - iy2)
                        kk = nctz + 1 + (iz1 - iz2)
                        if(ii.lt.1.or.ii.gt.MAXCTX.or.
     +                     jj.lt.1.or.jj.gt.MAXCTY.or.
     +                     kk.lt.1.or.kk.gt.MAXCTZ) then
                              call cova3(x1,y1,z1,x2,y2,z2,1,nst,MAXNST,
     +                             c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                        else
                              cov = covtab(ii,jj,kk)
                        endif
                        a(in) = dble(cov)
                  endif
            end do
!
! Get the RHS value (possibly with covariance look-up table):
!
            if(j.le.nclose) then
                  call cova3(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,cc,aa,
     +                       1,MAXROT,rotmat,cmax,cov)
                  r(j) = dble(cov)
            else
!
! Try to use the covariance look-up (if the distance is in range):
!
                  ii = nctx + 1 + (ix - ix1)
                  jj = ncty + 1 + (iy - iy1)
                  kk = nctz + 1 + (iz - iz1)
                  if(ii.lt.1.or.ii.gt.MAXCTX.or.
     +               jj.lt.1.or.jj.gt.MAXCTY.or.
     +               kk.lt.1.or.kk.gt.MAXCTZ) then
                        call cova3(xx,yy,zz,x1,y1,z1,1,nst,MAXNST,c0,it,
     +                             cc,aa,1,MAXROT,rotmat,cmax,cov)
                  else
                        cov = covtab(ii,jj,kk)
                  endif
                  r(j) = dble(cov)
            endif
            rr(j) = r(j)
      end do
!
! Addition of OK constraint:
!
      if(lktype.eq.1.or.lktype.eq.3) then
            do i=1,na
                  in    = in + 1
                  a(in) = 1.0
            end do
            in       = in + 1
            a(in)    = 0.0
            r(na+1)  = 1.0
            rr(na+1) = 1.0
      endif
!
! Addition of the External Drift Constraint:
!
      if(lktype.eq.3) then
            edmin =  999999.
            edmax = -999999.
            do i=1,na
                  in    = in + 1
                  a(in) = vrea(i)
                  if(a(in).lt.edmin) edmin = a(in)
                  if(a(in).gt.edmax) edmax = a(in)
            end do
            in       = in + 1
            a(in)    = 0.0
            in       = in + 1
            a(in)    = 0.0
            ind      = ix + (iy-1)*nx + (iz-1)*nxy
            r(na+2)  = dble(lvm(ind))
            rr(na+2) = r(na+2)
            if((edmax-edmin).lt.EPSLON) neq = neq - 1
      endif
!
! Addition of Collocated Cosimulation Constraint:
!
      if(lktype.eq.4) then
            sfmin =  1.0e21
            sfmax = -1.0e21
            do i=1,na
                  in    = in + 1
                  a(in) = dble(colocorr)*r(i)
                  if(a(in).lt.sfmin) sfmin = a(in)
                  if(a(in).gt.sfmax) sfmax = a(in)
            end do
            in    = in + 1
            a(in) = 1.0
            ii    = na + 1
            r(ii) = dble(colocorr)
            rr(ii)= r(ii)
!           if((sfmax-sfmin).lt.EPSLON) neq = neq - 1
      end if
!
! Write out the kriging Matrix if Seriously Debugging:
!
      if(idbg.ge.3) then
            write(ldbg,100) ix,iy,iz
            is = 1
            do i=1,neq
                  ie = is + i - 1
                  write(ldbg,101) i,r(i),(a(j),j=is,ie)
                  is = is + i
            end do
 100        format(/,'Kriging Matrices for Node: ',3i4,' RHS first')
 101        format('    r(',i2,') =',f7.4,'  a= ',99f7.4)
      endif
!
! Solve the Kriging System:
!
      if(neq.eq.1.and.lktype.ne.3) then
            s(1)  = r(1) / a(1)
            ising = 0
      else
            call ksol(1,neq,1,a,r,s,ising)
      endif
!
! Write a warning if the matrix is singular:
!
      if(ising.ne.0) then
!jd            if(idbg.ge.1) then
!jd                  write(ldbg,*) 'WARNING SGSIM: singular matrix'
!jd                  write(ldbg,*) '               for node',ix,iy,iz
!jd            endif
            cmean  = gmean
            cstdev = astdev             !jd
            return
      endif
!
! Compute the estimate and kriging variance.  Recall that kriging type
!     0 = Simple Kriging:
!     1 = Ordinary Kriging:
!     2 = Locally Varying Mean:
!     3 = External Drift:
!     4 = Collocated Cosimulation:
!
      cmean  = 0.0
      cstdev = cbb
      sumwts = 0.0
      do i=1,na
            cmean  = cmean  + real(s(i))*vra(i)
            cstdev = cstdev - real(s(i)*rr(i))
            sumwts = sumwts + real(s(i))
      end do

      if(lktype.eq.0) cmean = cmean + (1.0 - sumwts) * gmean   !jd

      if(lktype.eq.1) cstdev = cstdev - real(s(na+1))

      if(lktype.eq.2) cmean  = cmean + gmean

      if(lktype.eq.4) then
            ind    = ix + (iy-1)*nx + (iz-1)*nxy
            cmean  = cmean  + real(s(na+1))*lvm(ind)
            cstdev = cstdev - real(s(na+1) *rr(na+1))
      end if
!
! Error message if negative variance:
!
      if(cstdev.lt.0.0) then
!jd            write(ldbg,*) 'ERROR: Negative Variance: ',cstdev
            cstdev = 0.0
      endif
      cstdev = sqrt(max(cstdev,0.0))
!
! Write out the kriging Weights if Seriously Debugging:
!
      if(idbg.ge.3) then
            do i=1,na
                  write(ldbg,140) i,vra(i),s(i)
            end do
 140        format(' Data ',i4,' value ',f8.4,' weight ',f8.4)
            if(lktype.eq.4) write(ldbg,141) lvm(ind),s(na+1)
 141        format(' Sec Data  value ',f8.4,' weight ',f8.4)
            write(ldbg,142) gmean,cmean,cstdev
 142        format(' Global mean ',f8.4,' conditional ',f8.4,
     +             ' std dev ',f8.4)
      end if
!
! Finished Here:
!
      return
      end


      subroutine comsize(n,m,achar)

        integer n,m
        character*(*) achar

        if(n.gt.m)then
          write(6,10) trim(achar),n
10        format(/,' *** Increase parameter ',a,' to at least ',i4,
     +    ' and re-compile program ***',/)
          stop
        end if

      end

      double precision function acorni(idum)
c-----------------------------------------------------------------------
c
c Fortran implementation of ACORN random number generator of order less
c than or equal to 12 (higher orders can be obtained by increasing the
c parameter value MAXORD).
c
c
c NOTES: 1. The variable idum is a dummy variable. The common block
c           IACO is used to transfer data into the function.
c
c        2. Before the first call to ACORN the common block IACO must
c           be initialised by the user, as follows. The values of
c           variables in the common block must not subsequently be
c           changed by the user.
c
c             KORDEI - order of generator required ( must be =< MAXORD)
c
c             MAXINT - modulus for generator, must be chosen small
c                      enough that 2*MAXINT does not overflow
c
c             ixv(1) - seed for random number generator
c                      require 0 < ixv(1) < MAXINT
c
c             (ixv(I+1),I=1,KORDEI)
c                    - KORDEI initial values for generator
c                      require 0 =< ixv(I+1) < MAXINT
c
c        3. After initialisation, each call to ACORN generates a single
c           random number between 0 and 1.
c
c        4. An example of suitable values for parameters is
c
c             KORDEI   = 10
c             MAXINT   = 2**30
c             ixv(1)   = an odd integer in the (approximate) range 
c                        (0.001 * MAXINT) to (0.999 * MAXINT)
c             ixv(I+1) = 0, I=1,KORDEI
c
c
c
c Author: R.S.Wikramaratna,                           Date: October 1990
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      common/iaco/ ixv(MAXOP1)
      do i=1,KORDEI
            ixv(i+1)=(ixv(i+1)+ixv(i))
            if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
      end do
      acorni=dble(ixv(KORDEI+1))/MAXINT
      return
      end


      subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c              Sets up an Anisotropic Rotation Matrix
c              **************************************
c
c Sets up the matrix to transform cartesian coordinates to coordinates
c accounting for angles and anisotropy (see manual for a detailed
c definition):
c
c
c INPUT PARAMETERS:
c
c   ang1             Azimuth angle for principal direction
c   ang2             Dip angle for principal direction
c   ang3             Third rotation angle
c   anis1            First anisotropy ratio
c   anis2            Second anisotropy ratio
c   ind              matrix indicator to initialize
c   MAXROT           maximum number of rotation matrices dimensioned
c   rotmat           rotation matrices
c
c
c NO EXTERNAL REFERENCES
c
c
c-----------------------------------------------------------------------
      parameter(DEG2RAD=3.141592654/180.0,EPSLON=1.e-20)
      real*8    rotmat(MAXROT,3,3),afac1,afac2,sina,sinb,sint,
     +          cosa,cosb,cost
c
c Converts the input angles to three angles which make more
c  mathematical sense:
c
c         alpha   angle between the major axis of anisotropy and the
c                 E-W axis. Note: Counter clockwise is positive.
c         beta    angle between major axis and the horizontal plane.
c                 (The dip of the ellipsoid measured positive down)
c         theta   Angle of rotation of minor axis about the major axis
c                 of the ellipsoid.
c
      if(ang1.ge.0.0.and.ang1.lt.270.0) then
            alpha = (90.0   - ang1) * DEG2RAD
      else
            alpha = (450.0  - ang1) * DEG2RAD
      endif
      beta  = -1.0 * ang2 * DEG2RAD
      theta =        ang3 * DEG2RAD
c
c Get the required sines and cosines:
c
      sina  = dble(sin(alpha))
      sinb  = dble(sin(beta))
      sint  = dble(sin(theta))
      cosa  = dble(cos(alpha))
      cosb  = dble(cos(beta))
      cost  = dble(cos(theta))
c
c Construct the rotation matrix in the required memory:
c
      afac1 = 1.0 / dble(max(anis1,EPSLON))
      afac2 = 1.0 / dble(max(anis2,EPSLON))
      rotmat(ind,1,1) =       (cosb * cosa)
      rotmat(ind,1,2) =       (cosb * sina)
      rotmat(ind,1,3) =       (-sinb)
      rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(ind,2,3) = afac1*( sint * cosb)
      rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(ind,3,3) = afac2*(cost * cosb)
c
c Return to calling program:
c
      return
      end


      subroutine setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
     +                   vr,tmp,nsec,sec1,sec2,sec3,MAXSBX,MAXSBY,
     +                   MAXSBZ,nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,
     +                   ysizsup,nzsup,zmnsup,zsizsup)
c-----------------------------------------------------------------------
c
c           Establish Super Block Search Limits and Sort Data
c           *************************************************
c
c This subroutine sets up a 3-D "super block" model and orders the data
c by super block number.  The limits of the super block is set to the
c minimum and maximum limits of the grid; data outside are assigned to
c the nearest edge block.
c
c The idea is to establish a 3-D block network that contains all the
c relevant data.  The data are then sorted by their index location in
c the search network, i.e., the index location is given after knowing
c the block index in each coordinate direction (ix,iy,iz):
c          ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix
c An array, the same size as the number of super blocks, is constructed
c that contains the cumulative number of data in the model.  With this
c array it is easy to quickly check what data are located near any given
c location.
c
c
c
c INPUT VARIABLES:
c
c   nx,xmn,xsiz      Definition of the X grid being considered
c   ny,ymn,ysiz      Definition of the Y grid being considered
c   nz,zmn,zsiz      Definition of the Z grid being considered
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   vr(nd)           Variable at each location.
c   tmp(nd)          Temporary storage to keep track of the super block
c                      index associated to each data (uses the same
c                      storage already allocated for the simulation)
c   nsec             Number of secondary variables to carry with vr
c   sec1(nd)         First secondary variable (if nsec >= 1)
c   sec2(nd)         Second secondary variable (if nsec >= 2)
c   sec3(nd)         Third secondary variable (if nsec = 3)
c   MAXSB[X,Y,Z]     Maximum size of super block network
c
c
c
c OUTPUT VARIABLES:
c
c   nisb()                Array with cumulative number of data in each
c                           super block.
c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
c   nysup,ymnsup,ysizsup  Definition of the Y super block grid
c   nzsup,zmnsup,zsizsup  Definition of the Z super block grid
c
c
c
c EXTERNAL REFERENCES:
c
c   sortem           Sorting routine to sort the data
c
c
c
c-----------------------------------------------------------------------
      real    x(*),y(*),z(*),vr(*),tmp(*),sec1(*),sec2(*),sec3(*)
      integer nisb(*)
      logical inflag
c
c Establish the number and size of the super blocks:
c
      nxsup   = min(nx,MAXSBX)
      nysup   = min(ny,MAXSBY)
      nzsup   = min(nz,MAXSBZ)
      xsizsup = real(nx)*xsiz/real(nxsup)
      ysizsup = real(ny)*ysiz/real(nysup)
      zsizsup = real(nz)*zsiz/real(nzsup)
      xmnsup  = (xmn-0.5*xsiz)+0.5*xsizsup
      ymnsup  = (ymn-0.5*ysiz)+0.5*ysizsup
      zmnsup  = (zmn-0.5*zsiz)+0.5*zsizsup
c
c Initialize the extra super block array to zeros:
c
      do i=1,nxsup*nysup*nzsup
            nisb(i) = 0
      end do
c
c Loop over all the data assigning the data to a super block and
c accumulating how many data are in each super block:
c
      do i=1,nd
            call getindx(nxsup,xmnsup,xsizsup,x(i),ix,inflag)
            call getindx(nysup,ymnsup,ysizsup,y(i),iy,inflag)
            call getindx(nzsup,zmnsup,zsizsup,z(i),iz,inflag)
            ii = ix + (iy-1)*nxsup + (iz-1)*nxsup*nysup
            tmp(i)   = ii
            nisb(ii) = nisb(ii) + 1
      end do
c
c Sort the data by ascending super block number:
c
      nsort = 4 + nsec
      call sortem(1,nd,tmp,nsort,x,y,z,vr,sec1,sec2,sec3)
c
c Set up array nisb with the starting address of the block data:
c
      do i=1,(nxsup*nysup*nzsup-1)
            nisb(i+1) = nisb(i) + nisb(i+1)
      end do
c
c Finished:
c
      return
      end

      subroutine picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
     +                   irot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
     +                   iysbtosr,izsbtosr)
c-----------------------------------------------------------------------
c
c             Establish Which Super Blocks to Search
c             **************************************
c
c This subroutine establishes which super blocks must be searched given
c that a point being estimated/simulated falls within a super block
c centered at 0,0,0.
c
c
c
c INPUT VARIABLES:
c
c   nxsup,xsizsup    Definition of the X super block grid
c   nysup,ysizsup    Definition of the Y super block grid
c   nzsup,zsizsup    Definition of the Z super block grid
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   radsqd           squared search radius
c
c
c
c OUTPUT VARIABLES:
c
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c
c
c
c-----------------------------------------------------------------------
      real*8  rotmat(MAXROT,3,3),hsqd,sqdist,shortest
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
c
c MAIN Loop over all possible super blocks:
c
      nsbtosr = 0
      do i=-(nxsup-1),(nxsup-1)
      do j=-(nysup-1),(nysup-1)
      do k=-(nzsup-1),(nzsup-1)
            xo = real(i)*xsizsup
            yo = real(j)*ysizsup
            zo = real(k)*zsizsup
c
c Find the closest distance between the corners of the super blocks:
c
            shortest = 1.0e21
            do i1=-1,1
            do j1=-1,1
            do k1=-1,1
                  do i2=-1,1
                  do j2=-1,1
                  do k2=-1,1
                        if(i1.ne.0.and.j1.ne.0.and.k1.ne.0.and.
     +                     i2.ne.0.and.j2.ne.0.and.k2.ne.0) then
                              xdis = real(i1-i2)*0.5*xsizsup + xo
                              ydis = real(j1-j2)*0.5*ysizsup + yo
                              zdis = real(k1-k2)*0.5*zsizsup + zo
                              hsqd = sqdist(0.0,0.0,0.0,xdis,ydis,zdis,
     +                                      irot,MAXROT,rotmat)
                              if(hsqd.lt.shortest) shortest = hsqd
                        end if
                  end do
                  end do
                  end do
            end do
            end do
            end do
c
c Keep this super block if it is close enoutgh:
c
            if(real(shortest).le.radsqd) then
                  nsbtosr = nsbtosr + 1
                  ixsbtosr(nsbtosr) = i
                  iysbtosr(nsbtosr) = j
                  izsbtosr(nsbtosr) = k
            end if
      end do
      end do
      end do
c
c Finished:
c
      return
      end


      subroutine sortem(ib,ie,a,iperm,b,c,d,e,f,g,h)
c-----------------------------------------------------------------------
c
c                      Quickersort Subroutine
c                      **********************
c
c This is a subroutine for sorting a real array in ascending order. This
c is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
c in collected algorithms of the ACM.
c
c The method used is that of continually splitting the array into parts
c such that all elements of one part are less than all elements of the
c other, with a third part in the middle consisting of one element.  An
c element with value t is chosen arbitrarily (here we choose the middle
c element). i and j give the lower and upper limits of the segment being
c split.  After the split a value q will have been found such that 
c a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
c performs operations on the two segments (i,q-1) and (q+1,j) as follows
c The smaller segment is split and the position of the larger segment is
c stored in the lt and ut arrays.  If the segment to be split contains
c two or fewer elements, it is sorted and another segment is obtained
c from the lt and ut arrays.  When no more segments remain, the array
c is completely sorted.
c
c
c INPUT PARAMETERS:
c
c   ib,ie        start and end index of the array to be sorteda
c   a            array, a portion of which has to be sorted.
c   iperm        0 no other array is permuted.
c                1 array b is permuted according to array a
c                2 arrays b,c are permuted.
c                3 arrays b,c,d are permuted.
c                4 arrays b,c,d,e are permuted.
c                5 arrays b,c,d,e,f are permuted.
c                6 arrays b,c,d,e,f,g are permuted.
c                7 arrays b,c,d,e,f,g,h are permuted.
c               >7 no other array is permuted.
c
c   b,c,d,e,f,g,h  arrays to be permuted according to array a.
c
c OUTPUT PARAMETERS:
c
c    a      = the array, a portion of which has been sorted.
c
c    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
c
c NO EXTERNAL ROUTINES REQUIRED:
c
c-----------------------------------------------------------------------
      dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
c
c The dimensions for lt and ut have to be at least log (base 2) n
c
      integer   lt(64),ut(64),i,j,k,m,p,q
c
c Initialize:
c
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
c
c If this segment has more than two elements  we split it
c
 10   if (j-i-1) 100,90,15
c
c p is the position of an arbitrary element in the segment we choose the
c middle element. Under certain circumstances it may be advantageous
c to choose p at random.
c
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
c
c Start at the beginning of the segment, search for k such that a(k)>t
c
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
c
c Such an element has now been found now search for a q such that a(q)<t
c starting at the end of the segment.
c
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
c
c a(q) has now been found. we interchange a(q) and a(k)
c
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
c
c Update q and search for another pair to interchange:
c
      q = q-1
      go to 20
 50   q = k-1
 60   continue
c
c The upwards search has now met the downwards search:
c
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
c
c The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
c store the position of the largest segment in lt and ut
c
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
c
c Update m and split the new smaller segment
c
 80   m = m+1
      go to 10
c
c We arrive here if the segment has  two elements we test to see if
c the segment is properly ordered if not, we perform an interchange
c
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
c
c If lt and ut contain more segments to be sorted repeat process:
c
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      return
      end

      subroutine gauinv(p,xp,ierr)
c-----------------------------------------------------------------------
c
c Computes the inverse of the standard normal cumulative distribution
c function with a numerical approximation from : Statistical Computing,
c by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95.
c
c
c
c INPUT/OUTPUT:
c
c   p    = double precision cumulative probability value: dble(psingle)
c   xp   = G^-1 (p) in single precision
c   ierr = 1 - then error situation (p out of range), 0 - OK
c
c
c-----------------------------------------------------------------------
      real*8 p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim,p
      save   p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim
c
c Coefficients of approximation:
c
      data lim/1.0e-10/
      data p0/-0.322232431088/,p1/-1.0/,p2/-0.342242088547/,
     +     p3/-0.0204231210245/,p4/-0.0000453642210148/
      data q0/0.0993484626060/,q1/0.588581570495/,q2/0.531103462366/,
     +     q3/0.103537752850/,q4/0.0038560700634/
c
c Check for an error situation:
c
      ierr = 1
      if(p.lt.lim) then
            xp = -1.0e10
            return
      end if
      if(p.gt.(1.0-lim)) then
            xp =  1.0e10
            return
      end if
      ierr = 0      
c
c Get k for an error situation:
c
      pp   = p
      if(p.gt.0.5) pp = 1 - pp
      xp   = 0.0
      if(p.eq.0.5) return
c
c Approximate the function:
c
      y  = dsqrt(dlog(1.0/(pp*pp)))
      xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) /
     +               ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
      if(real(p).eq.real(pp)) xp = -xp
c
c Return with G^-1(p):
c
      return
      end

      subroutine getindx(n,min,siz,loc,index,inflag)
c-----------------------------------------------------------------------
c
c     Gets the coordinate index location of a point within a grid
c     ***********************************************************
c
c
c n       number of "nodes" or "cells" in this coordinate direction
c min     origin at the center of the first cell
c siz     size of the cells
c loc     location of the point being considered
c index   output index within [1,n]
c inflag  true if the location is actually in the grid (false otherwise
c         e.g., if the location is outside then index will be set to
c         nearest boundary
c
c
c
c-----------------------------------------------------------------------
      integer   n,index
      real      min,siz,loc
      logical   inflag
c
c Compute the index of "loc":
c
      index = int( (loc-min)/siz + 1.5 )
c
c Check to see if in or out:
c
      if(index.lt.1) then
            index  = 1
            inflag = .false.
      else if(index.gt.n) then
            index  = n
            inflag = .false.
      else
            inflag = .true.
      end if
c
c Return to calling program:
c
      return
      end

      subroutine srchsupr(xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat,
     +                    nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd,
     +                    x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,
     +                    nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,
     +                    nclose,close,infoct)
c-----------------------------------------------------------------------
c
c              Search Within Super Block Search Limits
c              ***************************************
c
c
c This subroutine searches through all the data that have been tagged in
c the super block subroutine.  The close data are passed back in the
c index array "close".  An octant search is allowed.
c
c
c
c INPUT VARIABLES:
c
c   xloc,yloc,zloc   location of point being estimated/simulated
c   radsqd           squared search radius
c   irot             index of the rotation matrix for searching
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c   nsbtosr          Number of super blocks to search
c   ixsbtosr         X offsets for super blocks to search
c   iysbtosr         Y offsets for super blocks to search
c   izsbtosr         Z offsets for super blocks to search
c   noct             If >0 then data will be partitioned into octants
c   nd               Number of data
c   x(nd)            X coordinates of the data
c   y(nd)            Y coordinates of the data
c   z(nd)            Z coordinates of the data
c   tmp(nd)          Temporary storage to keep track of the squared
c                      distance associated with each data
c   nisb()                Array with cumulative number of data in each
c                           super block.
c   nxsup,xmnsup,xsizsup  Definition of the X super block grid
c   nysup,ymnsup,ysizsup  Definition of the X super block grid
c   nzsup,zmnsup,zsizsup  Definition of the X super block grid
c
c
c
c OUTPUT VARIABLES:
c
c   nclose           Number of close data
c   close()          Index of close data
c   infoct           Number of informed octants (only computes if
c                      performing an octant search)
c
c
c
c EXTERNAL REFERENCES:
c
c   sqdist           Computes anisotropic squared distance
c   sortem           Sorts multiple arrays in ascending order
c
c
c
c-----------------------------------------------------------------------
      real    x(*),y(*),z(*),tmp(*),close(*)
      real*8  rotmat(MAXROT,3,3),hsqd,sqdist
      integer nisb(*),inoct(8)
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
      logical inflag
c
c Determine the super block location of point being estimated:
c
      call getindx(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
c
c Loop over all the possible Super Blocks:
c
      nclose = 0
      do 1 isup=1,nsbtosr
c
c Is this super block within the grid system:
c
            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)
            if(ixsup.le.0.or.ixsup.gt.nxsup.or.
     +         iysup.le.0.or.iysup.gt.nysup.or.
     +         izsup.le.0.or.izsup.gt.nzsup) go to 1
c
c Figure out how many samples in this super block:
c
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif
c
c Loop over all the data in this super block:
c
            do 2 ii=1,nums
                  i = i + 1
c
c Check squared distance:
c
                  hsqd = sqdist(xloc,yloc,zloc,x(i),y(i),z(i),irot,
     +                          MAXROT,rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
c
c Accept this sample:
c
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = real(hsqd)
 2          continue
 1    continue
c
c Sort the nearby samples by distance to point being estimated:
c
      call sortem(1,nclose,tmp,1,close,c,d,e,f,g,h)
c
c If we aren't doing an octant search then just return:
c
      if(noct.le.0) return
c
c PARTITION THE DATA INTO OCTANTS:
c
      do i=1,8
            inoct(i) = 0
      end do
c
c Now pick up the closest samples in each octant:
c
      nt = 8*noct
      na = 0
      do j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
c
c Keep this sample if the maximum has not been exceeded:
c
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
c
c End of data selection. Compute number of informed octants and return:
c
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
c
c Finished:
c
      return
      end

      real function backtr(vrgs,nt,vr,vrg,zmin,zmax,ltail,ltpar,
     +                     utail,utpar)
c-----------------------------------------------------------------------
c
c           Back Transform Univariate Data from Normal Scores
c           *************************************************
c
c This subroutine backtransforms a standard normal deviate from a
c specified back transform table and option for the tails of the
c distribution.  Call once with "first" set to true then set to false
c unless one of the options for the tail changes.
c
c
c
c INPUT VARIABLES:
c
c   vrgs             normal score value to be back transformed
c   nt               number of values in the back transform tbale
c   vr(nt)           original data values that were transformed
c   vrg(nt)          the corresponding transformed values
c   zmin,zmax        limits possibly used for linear or power model
c   ltail            option to handle values less than vrg(1):
c   ltpar            parameter required for option ltail
c   utail            option to handle values greater than vrg(nt):
c   utpar            parameter required for option utail
c
c
c
c-----------------------------------------------------------------------
      parameter(EPSLON=1.0e-20)
      dimension vr(nt),vrg(nt)
      real      ltpar,utpar,lambda
      integer   ltail,utail
c
c Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
c
      if(vrgs.le.vrg(1)) then
            backtr = vr(1)
            cdflo  = gcum(vrg(1))
            cdfbt  = gcum(vrgs)
            if(ltail.eq.1) then
                  backtr = powint(0.0,cdflo,zmin,vr(1),cdfbt,1.0)
            else if(ltail.eq.2) then
                  cpow   = 1.0 / ltpar
                  backtr = powint(0.0,cdflo,zmin,vr(1),cdfbt,cpow)
            endif
c
c Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
c
      else if(vrgs.ge.vrg(nt)) then
            backtr = vr(nt)
            cdfhi  = gcum(vrg(nt))
            cdfbt  = gcum(vrgs)
            if(utail.eq.1) then
                  backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,1.0)
            else if(utail.eq.2) then
                  cpow   = 1.0 / utpar
                  backtr = powint(cdfhi,1.0,vr(nt),zmax,cdfbt,cpow)
            else if(utail.eq.4) then
                  lambda = (vr(nt)**utpar)*(1.0-gcum(vrg(nt)))
                  backtr = (lambda/(1.0-gcum(vrgs)))**(1.0/utpar)
            endif
      else
c
c Value within the transformation table:
c
            call locate(vrg,nt,1,nt,vrgs,j)
            j = max(min((nt-1),j),1)
            backtr = powint(vrg(j),vrg(j+1),vr(j),vr(j+1),vrgs,1.0)
      endif
      return
      end



      subroutine cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c                      5. hole effect model
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested 
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c-----------------------------------------------------------------------
      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-10)
      integer   nst(*),it(*)
      real      c0(*),cc(*),aa(*)
      real*8    rotmat(MAXROT,3,3),hsqd,sqdist
c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c
      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
            h = real(dsqrt(hsqd))
c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then
!jd                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
                  cova = cova + cc(ist)*exp(-1.0*h/aa(ist))
c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
!jd                  cova = cova + cc(ist)*exp(-(3.0*h/aa(ist))
!jd     +                                      *(3.0*h/aa(ist)))
                  cova = cova + cc(ist)*exp(-(1.0*h/aa(ist))
     +                                      *(1.0*h/aa(ist)))

c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*(h**aa(ist))
c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            endif
      end do
c
c Finished:
c
      return
      end


      real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c    Squared Anisotropic Distance Calculation Given Matrix Indicator
c    ***************************************************************
c
c This routine calculates the anisotropic distance between two points
c  given the coordinates of each point and a definition of the
c  anisotropy.
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   ind              The rotation matrix to use
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c
c OUTPUT VARIABLES:
c
c   sqdist           The squared distance accounting for the anisotropy
c                      and the rotation of coordinates (if any).
c
c
c NO EXTERNAL REFERENCES
c
c
c-----------------------------------------------------------------------
      real*8 rotmat(MAXROT,3,3),cont,dx,dy,dz
c
c Compute component distance vectors and the squared distance:
c
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx
     +             + rotmat(ind,i,2) * dy
     +             + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
      end do
      return
      end


      subroutine ksol(nright,neq,nsb,a,r,s,ising)
c-----------------------------------------------------------------------
c
c                Solution of a System of Linear Equations
c                ****************************************
c
c
c
c INPUT VARIABLES:
c
c   nright,nsb       number of columns in right hand side matrix.
c                      for KB2D: nright=1, nsb=1
c   neq              number of equations
c   a()              upper triangular left hand side matrix (stored 
c                      columnwise)
c   r()              right hand side matrix (stored columnwise)
c                      for kb2d, one column per variable
c
c
c
c OUTPUT VARIABLES:
c
c   s()              solution array, same dimension as  r  above.
c   ising            singularity indicator
c                      0,  no singularity problem
c                     -1,  neq .le. 1
c                      k,  a null pivot appeared at the kth iteration
c
c
c
c PROGRAM NOTES:
c
c   1. Requires the upper triangular left hand side matrix.
c   2. Pivots are on the diagonal.
c   3. Does not search for max. element for pivot.
c   4. Several right hand side matrices possible.
c   5. USE for ok and sk only, NOT for UK.
c
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
c
c If there is only one equation then set ising and return:
c
      if(neq.le.1) then
            ising = -1
            return
      endif
c
c Initialize:
c
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
c
c Start triangulation:
c
      do k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            do iv=1,nright
                  nm1=nm*(iv-1)
                  ii=kk+nn*(iv-1)
                  piv=1./a(ii)
                  lp=0
                  do i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
                        end do
                        do llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
                        end do
                  end do
            end do
      end do
c
c Error checking - singular matrix:
c
      ijm=ij-nn*(nright-1)
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
c
c Finished triangulation, start solving back:
c
      do iv=1,nright
            nm1=nm*(iv-1)
            ij=ijm+nn*(iv-1)
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
      end do
c
c Finished solving back, return:
c
      return
      end


      real function gcum(x)
c-----------------------------------------------------------------------
c
c Evaluate the standard normal cdf given a normal deviate x.  gcum is
c the area under a unit normal curve to the left of x.  The results are
c accurate only to about 5 decimal places.
c
c
c-----------------------------------------------------------------------
      z = x
      if(z.lt.0.) z = -z
      t    = 1./(1.+ 0.2316419*z)
      gcum = t*(0.31938153   + t*(-0.356563782 + t*(1.781477937 +
     +       t*(-1.821255978 + t*1.330274429))))
      e2   = 0.
c
c  6 standard deviations out gets treated as infinity:
c
      if(z.le.6.) e2 = exp(-z*z/2.)*0.3989422803
      gcum = 1.0- e2 * gcum
      if(x.ge.0.) return
      gcum = 1.0 - gcum
      return
      end


      real function powint(xlow,xhigh,ylow,yhigh,xval,pow)
c-----------------------------------------------------------------------
c
c Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
c                 for a value of x and a power pow.
c
c-----------------------------------------------------------------------
      parameter(EPSLON=1.0e-20)

      if((xhigh-xlow).lt.EPSLON) then
            powint = (yhigh+ylow)/2.0
      else
            powint = ylow + (yhigh-ylow)* 
     +               (((xval-xlow)/(xhigh-xlow))**pow)
      end if

      return
      end

      subroutine locate(xx,n,is,ie,x,j)
c-----------------------------------------------------------------------
c
c Given an array "xx" of length "n", and given a value "x", this routine
c returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
c must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
c returned to indicate that x is out of range.
c
c Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
c-----------------------------------------------------------------------
      dimension xx(n)
c
c Initialize lower and upper methods:
c
      if(is.le.0) is = 1
      jl = is-1
      ju = ie
      if(xx(n).le.x) then
            j = ie
            return
      end if
c
c If we are not done then compute a midpoint:
c
 10   if(ju-jl.gt.1) then
            jm = (ju+jl)/2
c
c Replace the lower or upper limit with the midpoint:
c
            if((xx(ie).gt.xx(is)).eqv.(x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
      endif
c
c Return with the array index:
c
      j = jl
      return
      end


