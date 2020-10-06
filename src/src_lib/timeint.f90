subroutine time_interp(ifail,nbore,ndays,nsecs,value,intday,intsec, &
rnear,rconst,valinterp,extrap,direction)

! -- Subroutine time_interp interpolates an array of times and values to a
! -- user-supplied time.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition arises
!       nbore:     number of times and corresponding values to be interpolated
!       ndays:     elapsed days corresponding to each value
!       nsecs:     elpased seconds corresponding to each value
!       value:     array of time-based values to be interpolated
!       intday:    the day to be interpolated to, expressed as elapsed days
!       intsec:    the time to be intepolated to, expressed as elapsed seconds
!       rnear:     maximum permitted days to nearest sample
!       rconst:    maximum days to nearest sample if interpolation cannot take
!                  place
!       valinterp: interpolated value
!       extrap:    'yes' if use linear extrapolation to (within rconst) if
!                  interpolation cannot take place
!       direction: 'lo' if extrapolation from two previous points,
!                  'hi' if extrapolation from two following points,
!                  'med' if interpolation if possible, otherwise extrapolation
!                  (note: 'med' is the default)

! -- Revision history:-
!       June-November, 1995: version 1.
!       Mid Year, 1997: incorporated "extrap" argument for SMPCAL
!	October, 1997: incorporated "direction" input for SMPCAL

	use defn
	use inter

	integer, intent(out)                    :: ifail
	integer, intent(in)                     :: nbore
	integer, intent(in), dimension(nbore)   :: ndays,nsecs
	double precision, intent(in), dimension(nbore)      :: value
	integer, intent(in)                     :: intday,intsec
	real, intent(in)                        :: rnear,rconst
	double precision, intent(out)           :: valinterp
	character (len=*), intent(in), optional :: extrap
	character (len=*), intent(in), optional :: direction

	integer                                 :: i,ie,id
	double precision                        :: secfac,diff,diff1,dentime
	character (len=3)                       :: atemp


	ie=0
	if(present(extrap)) then
	  atemp=extrap
	  call casetrans(atemp,'lo')
	  if(atemp.eq.'yes')then
	    ie=1
	  else if(atemp.eq.'no') then
	    ie=0
	  else
	    call sub_error('TIME_INTERP')
	  end if
	end if

	id=0
	if(present(direction))then
	  atemp=direction
	  call casetrans(atemp,'lo')
	  if(atemp.eq.'lo')then
	    id=-1
	  else if(atemp.eq.'hi')then
	    id=1
	  else if(atemp.eq.'med')then
	    id=0
	  else
	    call sub_error('TIME_INTERP')
	  end if
	end if

	if((id.ne.0).and.(ie.eq.0))then
	  call sub_error('TIME_INTERP')
	end if

	ifail=0
	secfac=1.0d0/86400.0d0
	if(nbore.eq.1) then
	  diff=dble(intday-ndays(1))+dble(intsec-nsecs(1))*secfac
	  if(abs(diff).le.rconst)then
	    valinterp=value(1)
	  else
	    if(diff.gt.0)then
	      valinterp=-9.1e37
	    else
	      valinterp=-8.1e37
	    end if
	  end if
	  return
	end if

	do i=1,nbore-1
	  if((ndays(i).gt.ndays(i+1)).or. &
	    ((ndays(i).eq.ndays(i+1)).and.(nsecs(i).ge.nsecs(i+1))))then
	    ifail=1
	    return
	  end if
	end do

	do i=1,nbore
	  diff=dble(ndays(i)-intday)+dble(nsecs(i)-intsec)*secfac
	  if(diff.ge.0)then
	    if(i.eq.1)then
	      if(diff.le.rconst)then
	        if(ie.eq.1)then
	          if((value(1).lt.-1.0e38).or.(value(2).lt.-1.0e38))then
	            valinterp=value(1)
	          else
	            dentime=dble(ndays(i+1)-ndays(i))+ &
                            dble(nsecs(i+1)-nsecs(i))*secfac 
	            if(dentime.le.0) then
	              ifail=1
	              return
	            else
	              valinterp=value(i)-(value(i+1)-value(i))*diff/dentime
	            end if
	          end if
	        else
		  valinterp=value(1)
	        end if
	      else
		valinterp=-8.1e37
	      end if
	      return
	    end if

	    if(id.eq.-1)then
	      if(i.eq.2)then
	        diff1=dble(intday-ndays(1))+dble(intsec-nsecs(1))*secfac
	        if(diff1.gt.rnear)then               !note - not rconst
	          valinterp=-7.1e37
	        else
	          valinterp=value(1)
	        end if
	        if(value(1).lt.-1.0e38) valinterp=-1.1e38
	      else
	        dentime=dble(ndays(i-1)-ndays(i-2))+ &
                        dble(nsecs(i-1)-nsecs(i-2))*secfac
	        if(dentime.lt.0.0d0)then
	          ifail=1
	          return
	        else
	          diff1=dble(intday-ndays(i-1))+  &
                  dble(intsec-nsecs(i-1))*secfac
	          if(diff1.gt.rnear)then
	            valinterp=-7.1e37
	          else
	            if(value(i-1).lt.-1.0e38)then
	              valinterp=-1.1e38
	            else if(value(i-2).lt.-1.0e38) then
	              valinterp=value(i-1)
	            else
	              valinterp=value(i-1)+ &
                      (value(i-1)-value(i-2))/dentime*diff1
	            end if
	          end if
	        end if
	      end if
	      return
	    else if(id.eq.1)then
	      if(i.eq.nbore)then
	        if(diff.gt.rnear)then
	          valinterp=-7.1e37
	        else
	          valinterp=value(i)
	        end if
	        if(value(i).lt.-1.0e38)valinterp=-1.1e38
	      else
	        dentime=dble(ndays(i+1)-ndays(i))+      &
	                dble(nsecs(i+1)-nsecs(i))*secfac
	        if(dentime.le.0)then
	          ifail=1
	          return
	        else
	          if(diff.gt.rnear)then
	            valinterp=-7.1e37
	          else
	            if(value(i).lt.-1.0e38)then
	              valinterp=-1.1e38
	            else if(value(i+1).lt.-1.0e38)then
	              valinterp=value(i)
	            else
	              valinterp=value(i)-  &
                      (value(i+1)-value(i))/dentime*diff
	            end if
	          end if
	        end if
	      end if
	      return
	    else 

	      dentime=dble(ndays(i)-ndays(i-1))+ &
	      dble(nsecs(i)-nsecs(i-1))*secfac
	      if(dentime.le.0)then
	        ifail=1
	        return
	      else 
	        diff1=dentime-diff
	        if((diff1.gt.rnear).and.(diff.gt.rnear))then
		  valinterp=-7.1e37
	        else
		  valinterp=value(i-1)+(value(i)-value(i-1))/dentime*diff1
	        end if
	        if(value(i).lt.-1.0e38)then
		  if(diff1.le.rconst)then
		    valinterp=value(i-1)
		  else
		    valinterp=-1.1e38
		  end if
	        else if(value(i-1).lt.-1.0e38)then
		  if(diff.le.rconst)then
		    valinterp=value(i)
		  else
		    valinterp=-1.1e38
		  end if
	        end if
	      end if
	    return
	    end if
	  end if
	end do

	diff1=dble(intday-ndays(nbore))+dble(intsec-nsecs(nbore))*secfac
	if(diff1.le.rconst)then
	  if(ie.eq.1) then
	    if((value(nbore).lt.-1.0e38).or.(value(nbore-1).lt.-1.0e38)) then
	      valinterp=value(nbore)
	    else
	      dentime=dble(ndays(nbore)-ndays(nbore-1))    &
                     +dble(nsecs(nbore)-nsecs(nbore-1))*secfac
	      if(dentime.le.0) then
	        ifail=1
	        return
	      else
	        valinterp=value(nbore)+(value(nbore)-value(nbore-1))*    &
                diff1/dentime
	      end if
	    end if
	  else
	    valinterp=value(nbore)
	  end if
	else
	  valinterp=-9.1e37
	end if

	return
	end
 