
!***********************************************************
!
program post
!
!***********************************************************

implicit none

integer,parameter :: nzb=1
integer,parameter :: nx1=64,ny1=129,nz1=48
character(len=3) suffix,suffix1
character(len=20) nfichier,nfichier1,nchamp
integer :: num,i,j,k,nb_fich,nb_proc_total,nb_proc,longueur,long,itime
integer :: nxyz,count,itime1,itime2,itime3

real(8),dimension(nx1,ny1,nz1) :: umean1,vmean1,wmean1,uumean1,vvmean1,wwmean1
real(8),dimension(nx1,ny1,nz1) :: uvmean1,uwmean1,vwmean1
real(4),dimension(ny1) :: yp,ypi
!real(8),dimension(nx1,ny1,nz1) :: umean2,vmean2,wmean2,uumean2,vvmean2,wwmean2
!real(8),dimension(nx1,ny1,nz1) :: uvmean2,uwmean2,vwmean2
real(4) :: u_to,u_to1,u_to2,re,xl2,xl3,xnu
real(8),dimension(ny1,6) :: q_stat


open (15,file='yp.dat',form='formatted',status='unknown')
do j=1,ny1
   read(15,*) yp(j),ypi(j)
!print *,yp(j)
enddo
close(15)
print *,'read yp'

itime1=7000


OPEN(11,FILE='umean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) umean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
   !  print *,k,'UX',umean1(nx1/2,ny1/2,k)/itime
  ENDDO
  CLOSE(11)
OPEN(11,FILE='vmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
    ! print *,k,'UY',vmean1(nx1/2,ny1/2,k)/itime
  ENDDO
  CLOSE(11)
OPEN(11,FILE='wmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) wmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
    ! print *,k,'UZ',wmean1(nx1/2,ny1/2,k)/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uumean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uumean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
    ! print *,k,'UXUX',uumean1(nx1/2,ny1/2,k)/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vvmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
    ! print *,k,'UYUY',vvmean1(nx1/2,ny1/2,k)/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='wwmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) wwmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
    ! print *,k,'UZUZ',wwmean1(nx1/2,ny1/2,k)/itime
  ENDDO
  CLOSE(11)

  print *,'READ DATA DONE 1'


umean1=umean1/itime1
vmean1=vmean1/itime1
wmean1=wmean1/itime1
uumean1=uumean1/itime1
vvmean1=vvmean1/itime1
wwmean1=wwmean1/itime1



!DO AN AVERAGE IN X AND Z

q_stat=0.
do j=1,ny1
   do k=1,nz1
   do i=1,nx1
      q_stat(j,1)=q_stat(j,1)+umean1(i,j,k)
      q_stat(j,2)=q_stat(j,2)+vmean1(i,j,k)
      q_stat(j,3)=q_stat(j,3)+wmean1(i,j,k)
      ! Direct Reynolds stresses
      q_stat(j,4)=q_stat(j,4)+uumean1(i,j,k)
      q_stat(j,5)=q_stat(j,5)+vvmean1(i,j,k)
      q_stat(j,6)=q_stat(j,6)+wwmean1(i,j,k)
   enddo
   enddo
!   print *,j,umean1(nx1/2,j,nz1/2)
!   print 100,q_stat(j,1)/nx1/nz1,q_stat(j,4)/nx1/nz1,q_stat(j,4)/nx1/nz1-q_stat(j,1)/nx1/nz1*q_stat(j,1)/nx1/nz1
enddo
q_stat(:,:)=q_stat(:,:)/nx1/nz1




xnu=1./10322.629
re=10322.629
   xl3=xnu*(q_stat(2,1)-q_stat(1,1))/(yp(2)-yp(1))
   xl2=sqrt(xl3)/xnu
   print *,'hi+ = ',xl2
   u_to1=xl2
   xl3=xnu*(q_stat(ny1-1,1)-q_stat(ny1,1))/(yp(ny1)-yp(ny1-1))
   xl2=sqrt(xl3)/xnu
   print *,'hs+ = ',xl2
   u_to2=xl2

   print *, 'ecart',abs(u_to1-u_to2) 

   u_to=(u_to1+u_to2)/2.
   xl2=u_to
   print *,'re_to Moy',u_to
   
   u_to1=((u_to1+u_to2)/2.)
   u_to=u_to/re

do j=1,ny1
!   print 100,q_stat(j,1),q_stat(j,4),q_stat(j,4)-q_stat(j,1)*q_stat(j,1)
   q_stat(j,4)=sqrt((q_stat(j,4)-q_stat(j,1)*q_stat(j,1)))/u_to
   q_stat(j,5)=sqrt((q_stat(j,5)-q_stat(j,2)*q_stat(j,2)))/u_to
   q_stat(j,6)=sqrt((q_stat(j,6)-q_stat(j,3)*q_stat(j,3)))/u_to
enddo

100 format(2F12.6,2F12.6,2F12.6,2F12.6,2F12.6,2F12.6)

  open (144,file='test1.dat',form='formatted',status='unknown')
   do j=1,ny1/2+1
      write(144,100) yp(j)+1.,yp(j)*xl2+xl2,((q_stat(j,1)/u_to)+(q_stat(j,1)/u_to))/2.,&
           (q_stat(j,4)+q_stat(j,4))/2.,(q_stat(j,5)+q_stat(j,5))/2.,&
           (q_stat(j,6)+q_stat(j,6))/2.
   enddo
   close(144)

 open (144,file='test2.dat',form='formatted',status='unknown')
   do j=1,ny1/2+1
      write(144,100) yp(j)+1.,yp(j)*xl2+xl2,((q_stat(ny1-j+1,1)/u_to)+(q_stat(ny1-j+1,1)/u_to))/2.,&
           (q_stat(ny1-j+1,4)+q_stat(ny1-j+1,4))/2.,(q_stat(ny1-j+1,5)+q_stat(ny1-j+1,5))/2.,&
           (q_stat(ny1-j+1,6)+q_stat(ny1-j+1,6))/2.
   enddo
   close(144)


  open (144,file='test.dat',form='formatted',status='unknown')
   do j=1,ny1/2+1
      write(144,100) yp(j)+1.,yp(j)*xl2+xl2,((q_stat(ny1-j+1,1)/u_to)+(q_stat(j,1)/u_to))/2.,&
           (q_stat(ny1-j+1,4)+q_stat(j,4))/2.,(q_stat(j,5)+q_stat(ny1-j+1,5))/2.,&
           (q_stat(j,6)+q_stat(ny1-j+1,6))/2.
   enddo
   close(144)
   
  


end program post
!******************************************************************
!
subroutine numcar (num,car)
!
!******************************************************************!

character(len=3) car

if (num.ge.100) then
   write (car,1) num
1  format (i3)
else
   if (num.ge.10) then
      write (car,2) num
2     format ('0',i2)
   else
      write (car,3) num
3     format ('00',i1)
   endif
endif

return
end subroutine numcar
