c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c  
c  IN THIS FILE THERE ARE A BUNCH OF TRIDIAGONAL MATRIX INVERSION
C  ROUTINES.
C  THESE CAN BE PERIODIC OR NOT
C  THESE ROUTINES RELAY ON LU MATRIX DECOMPOSITIONS
C  SOME OF THE COEFFICIENTS ARE EVALUATED ONLY ONCE WHEN
C  IT POSSIBLE TO SAVE COMPUTATIONAL TIME
C  THE ROUTINES INCREASE THE EFFICIENCY BY DOING
C  INNER LLOPS VERY LONG
C  THIS IS THE REASON WHY SOMETIME IN THE RHS ONLY ONE INDEX IS
C  USED
C  THE USER INTERESTED TO VERY EFFICIENT CODES SHOULD CHANGE THESE
C  ROUTINES WHY THE MAJOR PART OF THE CPU TIME IS SPENT HERE
C  MAY BE MORE EFFICIENTS ROUTINE ARE AVAILABLE IN LITERATURE
C
c************************************************************************
c                                                                       *
c  ****************************** subrout ctpv1ij  ***********************
c                                                                       *
c************************************************************************
      subroutine ctpv1ij(ji,jf,ni,nf,nst)
      include 'param.f'
c
c     ROUTINE FOR THE INVERSION OF PERIODIC TRIDIAGONAL MATRICES
c     VECTORIZED FOR RIGHT HAND SIDE AND COEFFICIENTS
c
      ja = ji + 1
      jj = ji + jf
      do 20 k=ni,nf
      qqf(k,ji,nst) = -apif(k,ji,nst)/acif(k,ji,nst)
      ssf(k,ji,nst) = -amif(k,ji,nst)/acif(k,ji,nst)
   20 continue
c
c     forward elimination sweep
c
      do 10 j=ja,jf
      do 21 k=ni,nf
      pf(k) =1./( acif(k,j,nst) + amif(k,j,nst)*qqf(k,j-1,nst))
      qqf(k,j,nst) = - apif(k,j,nst)*pf(k)
      ssf(k,j,nst) = - amif(k,j,nst)*ssf(k,j-1,nst)*pf(k)
   21 continue
   10 continue
c
c     backward pass
c
      do 22 k=ni,nf
      ssf(k,jf,nst) = 1.
   22 continue
      do 11 i=ja,jf
      j = jj - i
      do 23 k=ni,nf
      ssf(k,j,nst) = ssf(k,j,nst) + qqf(k,j,nst)*ssf(k,j+1,nst)
   23 continue
   11 continue
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout trpv1ij  ***********************
c                                                                       *
c************************************************************************
      subroutine trpv1ij(ji,jf,ni,nf,nst)
      include 'param.f'
c
c     ROUTINE FOR THE INVERSION OF PERIODIC TRIDIAGONAL MATRICES
c     VECTORIZED FOR RIGHT HAND SIDE AND COEFFICIENTS
c
      ja = ji + 1
      jj = ji + jf
      do 20 k=ni,nf
      fnn(k) = fi(k,jf)
      fi(k,ji) = fi(k,ji)/acif(k,ji,nst)
   20 continue
c
c     forward elimination sweep
c
      do 10 j=ja,jf
      do 21 k=ni,nf
      pf(k) =1./( acif(k,j,nst) + amif(k,j,nst)*qqf(k,j-1,nst))
      fi(k,j) = ( fi(k,j) - amif(k,j,nst)*fi(k,j-1))*pf(k)
   21 continue
   10 continue
c
c     backward pass
c
      do 22 k=ni,nf
      fei(k,jf) = 0.
   22 continue
      do 11 i=ja,jf
      j = jj - i
      do 23 k=ni,nf
      fei(k,j) = fi(k,j) + qqf(k,j,nst)*fei(k,j+1)
   23 continue
   11 continue
      do 24 k=ni,nf
      fi(k,jf)=(fnn(k) - apif(k,jf,nst)*fei(k,ji)
     &       - amif(k,jf,nst)*fei(k,jf-1))
     &       /(apif(k,jf,nst)*ssf(k,ji,nst)
     &       + amif(k,jf,nst)*ssf(k,jf-1,nst)  +acif(k,jf,nst))
   24 continue
c
c     backward elimination pass
c
      do 12 i=ja,jf
      j = jj -i
      do 25 k=ni,nf
      fi(k,j) = fi(k,jf)*ssf(k,j,nst) + fei(k,j)
   25 continue
   12 continue
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout ctpv2ij  ***********************
c                                                                       *
c************************************************************************
      subroutine ctpv2ij(ji,jf,ni,nf,nst)
      include 'param.f'
c
c     ROUTINE FOR THE INVERSION OF PERIODIC TRIDIAGONAL MATRICES
c     VECTORIZED FOR RIGHT HAND SIDE AND COEFFICIENTS
c
      ja = ji + 1
      jj = ji + jf
      do 20 k=ni,nf
      qqv(k,ji,nst) = -apiv(k,ji,nst)/aciv(k,ji,nst)
      ssv(k,ji,nst) = -amiv(k,ji,nst)/aciv(k,ji,nst)
   20 continue
c
c
c     forward elimination sweep
c
      do 10 j=ja,jf
      do 21 k=ni,nf
      pv(k) =1./( aciv(k,j,nst) + amiv(k,j,nst)*qqv(k,j-1,nst))
      qqv(k,j,nst) = - apiv(k,j,nst)*pv(k)
      ssv(k,j,nst) = - amiv(k,j,nst)*ssv(k,j-1,nst)*pv(k)
   21 continue
   10 continue
c
c     backward pass
c
      do 22 k=ni,nf
      ssv(k,jf,nst) = 1.
   22 continue
      do 11 i=ja,jf
      j = jj - i
      do 23 k=ni,nf
      ssv(k,j,nst) = ssv(k,j,nst) + qqv(k,j,nst)*ssv(k,j+1,nst)
   23 continue
   11 continue
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout trpv2ij  ***********************
c                                                                       *
c************************************************************************
      subroutine trpv2ij(ji,jf,ni,nf,nst)
      include 'param.f'
c
c     ROUTINE FOR THE INVERSION OF PERIODIC TRIDIAGONAL MATRICES
c     VECTORIZED FOR RIGHT HAND SIDE AND COEFFICIENTS
c
      ja = ji + 1
      jj = ji + jf
      do 20 k=ni,nf
      fnn(k) = fi(k,jf)
      fi(k,ji) = fi(k,ji)/aciv(k,ji,nst)
   20 continue
c
c     forward elimination sweep
c
      do 10 j=ja,jf
      do 21 k=ni,nf
      pv(k) =1./( aciv(k,j,nst) + amiv(k,j,nst)*qqv(k,j-1,nst))
      fi(k,j) = ( fi(k,j) - amiv(k,j,nst)*fi(k,j-1))*pv(k)
   21 continue
   10 continue
c
c     backward pass
c
      do 22 k=ni,nf
      fei(k,jf) = 0.
   22 continue
      do 11 i=ja,jf
      j = jj - i
      do 23 k=ni,nf
      fei(k,j) = fi(k,j) + qqv(k,j,nst)*fei(k,j+1)
   23 continue
   11 continue
      do 24 k=ni,nf
      fi(k,jf)=(fnn(k) - apiv(k,jf,nst)*fei(k,ji)
     &       - amiv(k,jf,nst)*fei(k,jf-1))
     &       /(apiv(k,jf,nst)*ssv(k,ji,nst)
     &       + amiv(k,jf,nst)*ssv(k,jf-1,nst)  +aciv(k,jf,nst))
   24 continue
c
c     backward elimination pass
c
      do 12 i=ja,jf
      j = jj -i
      do 25 k=ni,nf
      fi(k,j) = fi(k,jf)*ssv(k,j,nst) + fei(k,j)
   25 continue
   12 continue
      return
      end
c************************************************************************
      subroutine btrjik(a,b,c,mm1,mm2,mm3)
      include'param.f'
c     a() sottodiagonale
c     b() diagonale
c     c() sopradiagonale
c     r() vettore termine noti e soluzione
c     m1,m2,m3   dimensioni della matrice
      dimension a(m2),b(m2),c(m2)
      dimension ggm(m2)
c***
c***
      mm1mm3=mm1*mm3
      bgt=b(1)
      subet=1./bgt
      do iadd=1,mm1mm3
        rhs(iadd)=rhs(iadd)*subet
      end do
c***
      do j=2,mm2
        ggm(j)=c(j-1)/bgt
        bgt=b(j)-a(j)*ggm(j)
        subet=1./bgt
       do ki=1,mm1mm3
         iaddc=ki+(j-1)*mm1mm3
         iaddm=ki+(j-2)*mm1mm3
         rhs(iaddc)=(rhs(iaddc)-a(j)*rhs(iaddm))*subet
       end do
      end do
      do j=mm2-1,1,-1
       do ki=1,mm1mm3
         iaddc=ki+(j-1)*mm1mm3
         iaddp=ki+(j  )*mm1mm3
         rhs(iaddc)=rhs(iaddc)-ggm(j+1)*rhs(iaddp)
       end do
      end do
c***
      return
      stop
      end
c***********************************************************************
      subroutine tripijk(n1i,n1f,n2i,n2f,n3i,n3f,nn1,nn2)
      include'param.f'
c
c     vectorized for right hand side and coefficients
c
      dimension fn(m2),p(m2),qi(m2,m1),si(m2,m1)
      ia = n1i + 1
      ii = n1i + n1f
      do k=n3i,n3f
      do 20 j=n2i,n2f
        qi(j,n1i) = -api(j)/aci(j)
        si(j,n1i) = -ami(j)/aci(j)
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
          p(j) =1./( aci(j) + ami(j)*qi(j,i-1))
          qi(j,i) = - api(j)*p(j)
          si(j,i) = - ami(j)*si(j,i-1)*p(j)
          rhs(iadd) = ( rhs(iadd) - ami(j)*rhs(iaddm))*p(j)
   21   continue
   10 continue
c
c     backward pass
c
      do 22 j=n2i,n2f
        si(j,n1f) = 1.
        fei(j,n1f) = 0.
   22 continue
      do 11 l=ia,n1f
        i = ii - l
        do 23 j=n2i,n2f
          iadd = i+(j-1)*nn1+(k-1)*nn1*nn2
          si(j,i) = si(j,i) + qi(j,i)*si(j,i+1)
          fei(j,i) = rhs(iadd) + qi(j,i)*fei(j,i+1)
   23   continue
   11 continue
      do 24 j=n2i,n2f
        iaddf = n1f+(j-1)*nn1+(k-1)*nn1*nn2
        rhs(iaddf)=(fn(j)-api(j)*fei(j,n1i) -
     %      ami(j)*fei(j,n1f-1))/(api(j)*si(j,n1i) +
     %      ami(j)*si(j,n1f-1)+aci(j))
   24 continue
c
c     backward elimination pass
c
      do 12 l=ia,n1f
        i = ii -l
        do 25 j=n2i,n2f
          iadd = i+(j-1)*nn1+(k-1)*nn1*nn2
          iaddf = n1f+(j-1)*nn1+(k-1)*nn1*nn2
          rhs(iadd) = rhs(iaddf)*si(j,i) + fei(j,i)
   25   continue
   12 continue
      end do
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout tripkji ********************* *
c                                                                       *
c************************************************************************
      subroutine tripkji(k1,k2,ni,nf,i1,i2,nn1,nn2)
      include 'param.f'
c
c     soluzione di una tridiagonale periodica 
c     con il metodo di Sherman-Morrison
c
c     amk() sottodiagonale
c     ack() diagonale
c     apk() sopradiagonale
c     rhs() vettore termine noti e soluzione
c     
      parameter (mg=200)
      dimension gmm(mg),r2(mg)
c***
c      rhs per correzione periodicita'
c***
      r2(1)=-amk(1)
      r2(k2-1)=-apk(k2-1)
      do k=2,k2-2
       r2(k)=0.
      end do
c
      betk=ack(1)
       subet=1./betk
       mm1mm2=nn1*nn2
      do ijadd=1,mm1mm2
        rhs(ijadd)=rhs(ijadd)*subet
      end do
      r2(1)=r2(1)/betk
c***
      do k=2,k2-1
       gmm(k)=apk(k-1)/betk
       betk=ack(k)-amk(k)*gmm(k)
       subet=1./betk
        do ijadd=1,mm1mm2
         kaddc= ijadd+(k-1)*mm1mm2
         kaddm= ijadd+(k-2)*mm1mm2
         rhs(kaddc)=(rhs(kaddc)-amk(k)*rhs(kaddm))*subet
        end do
      end do
      do k=k2-2,1,-1
        do ijadd=1,mm1mm2
         kaddc=ijadd+(k-1)*mm1mm2
         kaddp=ijadd+(k  )*mm1mm2
         rhs(kaddc)=rhs(kaddc)-gmm(k+1)*rhs(kaddp)
        end do
      end do
c
c     correzione periodicita'
      betk=ack(1)
      do k=2,k2-1
        gmm(k)=apk(k-1)/betk
        betk=ack(k)-amk(k)*gmm(k)
        r2(k)=(r2(k)-amk(k)*r2(k-1))/betk
      end do
      do k=k2-2,1,-1
        r2(k)=r2(k)-gmm(k+1)*r2(k+1)
      end do
c
c      
      betk=1./(ack(k2)+amk(k2)*r2(k2-1)+apk(k2)*r2(1))
        do ijadd=1,mm1mm2
         kaddk2=ijadd+(k2-1)*mm1mm2
         kaddk2m=ijadd+(k2-2)*mm1mm2
         rhs(kaddk2)=(rhs(kaddk2)-apk(k2)*rhs(ijadd)
     1              -amk(k2)*rhs(kaddk2m))*betk
        end do
      do k=1,k2-1
        do ijadd=1,mm1mm2
         kaddk2=ijadd+(k2-1)*mm1mm2
         kaddc =ijadd+(k-1)*mm1mm2
         rhs(kaddc)=rhs(kaddc)+r2(k)*rhs(kaddk2)
         end do
      end do
c***
      return
      end
