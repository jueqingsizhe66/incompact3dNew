  subroutine findfactor(num, factors, nfact)

    implicit none

    integer, intent(IN) :: num
    integer, intent(OUT), dimension(*) :: factors
    integer, intent(OUT) :: nfact
    integer :: i, m

    ! find the factors <= sqrt(num)
    m = int(sqrt(real(num)))
    nfact = 1
    do i=1,m
       if (num/i*i == num) then
          factors(nfact) = i
          nfact = nfact + 1
       end if
    end do
    nfact = nfact - 1

    ! derive those > sqrt(num)
    if (factors(nfact)**2/=num) then
       do i=nfact+1, 2*nfact
          factors(i) = num / factors(2*nfact-i+1)
       end do
       nfact = nfact * 2
    else
       do i=nfact+1, 2*nfact-1
          factors(i) = num / factors(2*nfact-i)
       end do
       nfact = nfact * 2 - 1
    endif

    return

  end subroutine findfactor
