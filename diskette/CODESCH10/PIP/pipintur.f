c***********************************************************************
c                                                                      *
c     this routine gives the initial condition for the pipe            *
c     starting with the Poiseuille profile plus a
c     random disturbance with amplitude reduced near the wall
c                                                                      * 
c***********************************************************************
      subroutine intur
      include 'param.f'
      common/pertq3/duinf(m1)
      write(6,*)  '  '
      write(32,*) '  '
      write(6,*) ' I N I T I A L   C O N D I T I O N '
      write(32,*) ' I N I T I A L   C O N D I T I O N '
      write(6,*)  '  '
      write(32,*) '  '
c
c   Poiseuille profile
c
      call invela
c*******************    axial velocity          ************************
      qout=0.
      qinf=0.
            q3mai=0.
      do  j=1,n2m
      darea=g2rm(j)/dx2*rm(j)/dx1
      do  i=1,n1m
         q3mai=max(uinfth(i,j),q3mai)
      qout=qout+uinfth(i,j)*darea
      qinf=qinf+uinfth(i,j)*darea
      do  k=1,n3
      q3(i,j,k)=uinfth(i,j)
      enddo
      enddo
      enddo
      write(42,*)'init. cond q3',q3mai,qinf,qout
      write(6,*)'init. cond q3',q3mai,qinf,qout
c     q1 comp.                                                          
c*******************    azimuthal velocity    ************************
c
            q1mai=0.
            q3mai=0.
            q2mai=0.
      do  j=1,n2m
      do  i=1,n1m
      im=imv(i)
      q1(i,j,1)=0.
      do  k=2,n3m
      km=kmv(k)
      q1(i,j,k)=0.
      q1mai=max(abs(q1(i,j,k)),q1mai)
      enddo
      enddo
      enddo
c     q2 comp.                                                          
c*******************      axial  velocity    random disturbance  *
      do  i=1,n1m
      duinf(i)=0.
      do  j=1,n2m
      ydw=1.-rm(j)
      vperq=vper
ccccccccccccccccccccccccccccccccccccc
      if(ydw.lt.0.025) vperq=vper/5.
ccccccccccccccccccccccccccccccccccccc
      do  k=1,n3m
      q3per=vperq*(-1.+2.*rand())
c     q3per=vperq*(-1.+2.*ranf())
      duinf(i)=max(q3per,duinf(i))
      q3(i,j,k)=q3(i,j,k)+q3per  
      q3mai=max(abs(q3(i,j,k)),q3mai)
      enddo
      enddo
      enddo
c
c    radial velocity from continuity equation
c
      do  k=1,n3m
      kp=kpv(k)
      do  i=1,n1m
      q2(i,1,k)=0.
      q2(i,n2,k)=0.
      do  j=1,n2m
      jp=j+1
      q2(i,jp,k)=q2(i,j,k)-rm(j)*g2rm(j)*dx3/dx2*
     1           (q3(i,j,kp)-q3(i,j,k))
      q2mai=max(abs(q2(i,jp,k)),q2mai)
      enddo
      enddo
      enddo
      write(42,*)'init. cond q1,q2,q3',q1mai,q2mai,q3mai
      write(6,*)'init. cond q1,q2,q3',q1mai,q2mai,q3mai
      return                                                            
      end                                                               

c***********************************************************************
c                                                                      *
c  ***  inizialize with laminar Poiseuille q3 profile 
c                                                                      *
c***********************************************************************
      subroutine invela
      include 'param.f'
      dimension va(m2)
      common/pertq3/duinf(m1)
c                                                                       
      pi=2.*asin(1.)                                                    
c                                                                       
       uvpma=0.
      open(13,file='velinf.out')
      open(16,file='vorinf.out')
      open(17,file='thper.out')
      do  j=1,n2m
       xi=rm(j)
       y=r0-xi
       uinf(j)=r0-xi**2
      enddo
      do  j=1,n2m
      do  i=1,n1m                                                    
       uinfth(i,j)=uinf(j) 
       enddo
       enddo
       vamax=0.
      do  j=1,n2m
       amo=amo+(1.-uinf(j))*uinf(j)/dx2*g2rm(j)
      enddo
       amo=1./amo
       write(6,*)'amo=',amo
      do  j=2,n2m
       xi=rc(j)  
       va(j)=-(uinf(j)-uinf(j-1))*dx2/g2rc(j)
       vamax=max(abs(va(j)),vamax)
       write(13,201) xi,uinf(j)
  201  format(4e12.5)
       write(16,201) xi,va(j)
       enddo
       va(1)=va(2)
       va(n2)=va(n2m)
       close(13)
       close(16)
       close(17)
      return                                                            
      end                                                               

