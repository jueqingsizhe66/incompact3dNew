  subroutine primefactors(num, factors, nfact)

    implicit none

    integer, intent(IN) :: num
    integer, intent(OUT), dimension(*) :: factors
    integer, intent(INOUT) :: nfact

    integer :: i, n

    i = 2
    nfact = 1
    n = num
    do
       if (mod(n,i) == 0) then
          factors(nfact) = i
          nfact = nfact + 1
          n = n / i
       else
          i = i + 1
       end if
       if (n == 1) then
          nfact = nfact - 1
          exit
       end if
    end do

    return

  end subroutine primefactors
