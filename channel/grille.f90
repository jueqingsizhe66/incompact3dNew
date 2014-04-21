! **************************************************************
!
   program grille
!
! **************************************************************
!
     implicit none
!
     integer,parameter ::nx=129,ny=65,nz=32
     integer,parameter :: my=ny+1
!
     real,dimension(nx) :: xx(nx)
     real,dimension(ny) :: yy(ny)
     real,dimension(nz) :: zz(nz)
     integer :: i,j,k,nclx,ncly,nclz,istret
     real :: xlx,yly,zlz,dx,dy,dz,x,y,z
     real,dimension(ny) :: ppy,pp2y,pp3y,pp4y,fux1,d1,d
     real,dimension(my) :: yp,yeta
     real :: yinf,beta,den,xnum,alpha,xcx,den1,den2,den3,den4,xnum1,cst
     real :: pi,fourpi
!
     xlx=20.
     yly=12.
     zlz=6.
     istret=3
     beta=0.80!2.531023707461935
     nclx=2
     ncly=1
     nclz=0
!
   if (nclx==0) dx=xlx/nx 
   if (nclx==1 .or. nclx==2) dx=xlx/(nx-1.) 
!
   if (ncly==0) dy=yly/ny 
   if (ncly==1 .or. ncly==2) dy=yly/(ny-1.) 
!
   if (nz>1) then 
      if (nclz==0) dz=zlz/nz 
      if (nclz==1 .or. nclz==2) dz=zlz/(nz-1.) 
   else 
      dz=zlz 
   endif
!
     x=0.
     do i=1,nx
        xx(i)=x
        x=x+dx 
     enddo
!
     y=0.
     do i=1,ny
        yy(i)=y
        y=y+dy 
     enddo
!
     z=0.
     do i=1,nz
        zz(i)=z
        z=z+dz 
     enddo
!
  pi=acos(-1.)
!     yinf cest la borne inf du domaine que l'on veut obtenir
  yinf=-yly/2.
  den=2.*beta*yinf
  xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
  alpha=abs(xnum/den)

  xcx=1./beta/alpha
!
     if (alpha.ne.0.) then
        if (istret.eq.1) yp(1)=0.
        if (istret.eq.2) yp(1)=0.
        if (istret.eq.1) yeta(1)=0.
        if (istret.eq.2) yeta(1)=-0.5	
        if (istret.eq.3) yp(1)=0.
        if (istret.eq.3) yeta(1)=-0.5
         do j=2,ny
            if (istret==1) then
               if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)
               if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))
            endif
            if (istret==2) then
               if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)-0.5
               if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))-0.5
            endif
           if (istret==3) then
               if (ncly.eq.0) yeta(j)=((j-1.)*(1./2./ny)-0.5)
               if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=((j-1.)*(1./2./(ny-1.))-0.5)
            endif
            den1=sqrt(alpha*beta+1.)
            xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
            den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
            den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
            den4=2.*alpha*beta-cos(2.*pi*yeta(j))+1.
            xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
            cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
            if (istret==1) then
               if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf
               if (yeta(j).eq.0.5) yp(j)=0.-yinf
               if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
            endif
            if (istret==2) then
               if (yeta(j).lt.0.5) yp(j)=xnum1-cst+yly
               if (yeta(j).eq.0.5) yp(j)=0.+yly
               if (yeta(j).gt.0.5) yp(j)=xnum1+cst+yly
            endif
            if (istret==3) then
               if (yeta(j).lt.0.5) yp(j)=(xnum1-cst+yly)*2.
               if (yeta(j).eq.0.5) yp(j)=(0.+yly)*2.
               if (yeta(j).gt.0.5) yp(j)=(xnum1+cst+yly)*2.
            endif
         enddo
   endif
   if (alpha.eq.0.) then
      yp(1)=-1.e10
      do j=2,ny
         yeta(j)=(j-1.)*(1./ny)
         yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
      enddo
   endif
!
!*****************************************************************
     open (10,file='avs.x',form='unformatted')
     write (10) (xx(i),i=1,nx)
     close (10)
!*****************************************************************
     if (istret.eq.0) then
        open (10,file='avs.y',form='unformatted')
        write (10) (yy(j),j=1,ny)
        close (10)   
     else
        open (10,file='avs.y',form='unformatted')
        write (10) (yp(j),j=1,ny)
        close (10)   
     endif
!*****************************************************************
     open (10,file='avs.z',form='unformatted')
     write (10) (zz(i),i=1,nz)
     close (10)
!
   end program grille
