c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c  IN THIS FILE THERE ARE A BUNCH OF TRIDIAGONAL MATRIX INVERSION
C  ROUTINES.
C  THESE CAN BE PERIODIC OR NOT
C  THESE ROUTINES RELAY ON LU MATRIX DECOMPOSITIONS
C  THE ROUTINES INCREASE THE EFFICIENCY BY DOING
C  INNER LLOPS VERY LONG
C  THIS IS THE REASON WHY IN THE RHS ONLY ONE INDEX IS
C  USED
C  THE USER INTERESTED TO VERY EFFICIENT CODES SHOULD CHANGE THESE
C  ROUTINES WHY THE MAJOR PART OF THE CPU TIME IS SPENT HERE
C  MAY BE MORE EFFICIENTS ROUTINE ARE AVAILABLE IN LITERATURE
C
c************************************************************************
      subroutine brtrj(a,b,c,m1,m2,m3)
c     a() sottodiagonale
c     b() diagonale
c     c() sopradiagonale
c     r() vettore termine noti e soluzione
c     m1,m2,m3   dimensioni della matrice
      include 'paramdi.f'
      common/rhsc/rhs(mm1*mm2*mm3)
      dimension a(mm2),b(mm2),c(mm2)
      dimension gam(mm2)
c***
c     if(m2.gt.mg) go to 8888
c***
      bet=b(1)
      subet=1./bet
      do k=1,m3
       do i=1,m1
        iadd=i+(k-1)*m1*m2
        rhs(iadd)=rhs(iadd)*subet
       end do
      end do
c***
      do j=2,m2
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        subet=1./bet
       do k=1,m3
        do i=1,m1
         iaddc=i+(j-1)*m1+(k-1)*m1*m2
         iaddm=i+(j-2)*m1+(k-1)*m1*m2
         rhs(iaddc)=(rhs(iaddc)-a(j)*rhs(iaddm))*subet
        end do
       end do
      end do
      do k=1,m3
       do j=m2-1,1,-1
        do i=1,m1
         iaddc=i+(j-1)*m1+(k-1)*m1*m2
         iaddp=i+(j  )*m1+(k-1)*m1*m2
         rhs(iaddc)=rhs(iaddc)-gam(j+1)*rhs(iaddp)
        end do
       end do
      end do
c***
      return
8888  print*,'mg too small'
      stop
      end
c************************************************************************
      subroutine brtrk(a,b,c,m1,m2,m3)
c     a() sottodiagonale
c     b() diagonale
c     c() sopradiagonale
c     r() vettore termine noti e soluzione
c     m1,m2,m3   dimensioni della matrice
      include 'paramdi.f'
      common/rhsc/rhs(mm1*mm2*mm3)
      dimension a(mm3),b(mm3),c(mm3)
      dimension gam(mm3)
c***
c     if(m2.gt.mg) go to 8888
c***
      bet=b(1)
      subet=1./bet
      do j=1,m2
       do i=1,m1
        iadd=i+(j-1)*m1
        rhs(iadd)=rhs(iadd)*subet
       end do
      end do
c***
      do k=2,m3
        gam(k)=c(k-1)/bet
        bet=b(k)-a(k)*gam(k)
        subet=1./bet
       do j=1,m2
        do i=1,m1
         iaddc=i+(j-1)*m1+(k-1)*m1*m2
         iaddm=i+(j-1)*m1+(k-2)*m1*m2
         rhs(iaddc)=(rhs(iaddc)-a(k)*rhs(iaddm))*subet
        end do
       end do
      end do
      do j=1,m2
       do k=m3-1,1,-1
        do i=1,m1
         iaddc=i+(j-1)*m1+(k-1)*m1*m2
         iaddp=i+(j-1)*m1+(k  )*m1*m2
         rhs(iaddc)=rhs(iaddc)-gam(k+1)*rhs(iaddp)
        end do
       end do
      end do
c***
      return
8888  print*,'mg too small'
      stop
      end
c***********************************************************************
      subroutine tripvmy(n1i,n1f,n2i,n2f,n3i,n3f,nn1,nn2)
      include'param.f'
c
c     vectorized for right hand side and coefficients
c
      dimension fn(m2),p(m2)
      ia = n1i + 1
      ii = n1i + n1f
      do k=n3i,n3f
      do 20 j=n2i,n2f
        q(j,n1i) = -api(j)/aci(j)
        s(j,n1i) = -ami(j)/aci(j)
        iaddf = n1f+(j-1)*nn1+(k-1)*nn1*nn2
        iaddi = n1i+(j-1)*nn1+(k-1)*nn1*nn2
        fn(j) = rhs(iaddf)
        rhs(iaddi) = rhs(iaddi)/aci(j)
   20 continue
c
c     forward elimination sweep
c
      do 10 i=ia,n1f
        do 21 j=n2i,n2f
          iadd = i+(j-1)*nn1+(k-1)*nn1*nn2
          iaddm = i-1+(j-1)*nn1+(k-1)*nn1*nn2
          p(j) =1./( aci(j) + ami(j)*q(j,i-1))
          q(j,i) = - api(j)*p(j)
          s(j,i) = - ami(j)*s(j,i-1)*p(j)
          rhs(iadd) = ( rhs(iadd) - ami(j)*rhs(iaddm))*p(j)
   21   continue
   10 continue
c
c     backward pass
c
      do 22 j=n2i,n2f
        s(j,n1f) = 1.
        fei(j,n1f) = 0.
   22 continue
      do 11 l=ia,n1f
        i = ii - l
        do 23 j=n2i,n2f
          iadd = i+(j-1)*nn1+(k-1)*nn1*nn2
          s(j,i) = s(j,i) + q(j,i)*s(j,i+1)
          fei(j,i) = rhs(iadd) + q(j,i)*fei(j,i+1)
   23   continue
   11 continue
      do 24 j=n2i,n2f
        iaddf = n1f+(j-1)*nn1+(k-1)*nn1*nn2
        rhs(iaddf)=(fn(j)-api(j)*fei(j,n1i) -
     %      ami(j)*fei(j,n1f-1))/(api(j)*s(j,n1i) +
     %      ami(j)*s(j,n1f-1)+aci(j))
   24 continue
c
c     backward elimination pass
c
      do 12 l=ia,n1f
        i = ii -l
        do 25 j=n2i,n2f
          iadd = i+(j-1)*nn1+(k-1)*nn1*nn2
          iaddf = n1f+(j-1)*nn1+(k-1)*nn1*nn2
          rhs(iadd) = rhs(iaddf)*s(j,i) + fei(j,i)
   25   continue
   12 continue
      end do
      return
      end

