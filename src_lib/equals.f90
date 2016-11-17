
!*****************************************************************************
! subroutines comprising the generic subroutine EQUALS ------->
!*****************************************************************************

! -- Subroutine equals compares two numbers of the same kind.

! -- Arguments are as follows:-
!       r1:   the first number
!       r2:   the second number

! -- It returns .TRUE. or .FALSE.


logical function equals_int(r1,r2)

        integer, intent(in)   :: r1
        integer, intent(in)   :: r2

        equals_int=(r1.eq.r2)

end function equals_int



logical function equals_real(r1,r2)

          real, intent(in)      :: r1
          real, intent(in)      :: r2

          real                  :: rtemp

          rtemp=abs(3.0*spacing(r1))
          if(abs(r1-r2).lt.rtemp)then
            equals_real=.true.
          else
            equals_real=.false.
          end if

end function equals_real


logical function equals_dbl(r1,r2)

          real (kind (1.0d0)), intent(in)      :: r1
          real (kind (1.0d0)), intent(in)      :: r2

          real (kind (1.0d0))                  :: rtemp

          rtemp=abs(3.0*spacing(r1))
          if(abs(r1-r2).lt.rtemp)then
            equals_dbl=.true.
          else
            equals_dbl=.false.
          end if

end function equals_dbl


