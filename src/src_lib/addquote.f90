!     Last change:  JD   28 Dec 2000    9:08 pm

subroutine addquote(afile,aqfile)

! -- Subroutine ADDQUOTE adds quotes to a filename if it has a space in it.

! -- Arguments are as follows:-
!        afile:       the name of the file
!        aqfile:      the name of the file with quotes added

        character (len=*), intent(in)   :: afile
        character (len=*), intent(out)  :: aqfile
        integer nbb

        if(index(trim(afile),' ').eq.0)then
          aqfile=afile
        else
          aqfile(1:1)='"'
          aqfile(2:)=trim(afile)
          nbb=len_trim(aqfile)+1
          aqfile(nbb:nbb)='"'
        end if

        return
end subroutine addquote


