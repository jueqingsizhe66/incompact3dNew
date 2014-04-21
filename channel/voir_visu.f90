!BE CAREFUL WITH FORTRAN COMPILER

!ifort -shared-intel -mcmodel=large -assume byterecl voir_visu.f90!

program visu

  implicit none

  integer, parameter :: nx=128 , ny=129, nz=84
  real(8),dimension(nx,ny,nz) :: ux
  real(4),dimension(nx,ny,nz) :: uy
  integer :: i,j,k,count,nfil
  real(4),dimension(nx) :: y1
  real(4),dimension(ny) :: y2,yp
  real(4),dimension(nz) :: y3
  real(4) :: pi,xa

pi=acos(-1.)

  OPEN(10,FILE='ux037',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           READ(10,REC=COUNT) ux(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     print *,k,ux(nx/2,ny/2,k)
  ENDDO
  CLOSE(10)

open(12,file='yp.dat',form='formatted',status='unknown')
do j=1,ny
   read(12,*) yp(j),xa
enddo
close(12)


do i=1,nx
   y1(i)=(i-1)*4.*pi/128.
enddo


do j=1,ny
   y2(j)=yp(j)+1.
enddo


do k=1,nz
   y3(k)=(k-1)*4./3.*pi/84.
enddo

uy=0.;uy=ux
  nfil=41
  open(nfil,file='tampon22.vtr')
  write(nfil,*)'<VTKFile type="RectilinearGrid" version="0.1"',' byte_order="LittleEndian">'
  write(nfil,*)'  <RectilinearGrid WholeExtent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(nfil,*)'    <Piece Extent=','"1 ',nx,' 1 ',ny,' 1 ',nz,'">'
  write(nfil,*)'      <Coordinates>'
  write(nfil,*)'        <DataArray type="Float32"',' Name="X_COORDINATES"',' NumberOfComponents="1">'
  write(nfil,*) (y1(i),i=1,nx)
  write(nfil,*)'        </DataArray>'
  write(nfil,*)'        <DataArray type="Float32"',' Name="Y_COORDINATES"',' NumberOfComponents="1">'
  write(nfil,*) (y2(j),j=1,ny)
  write(nfil,*)'        </DataArray>'
  write(nfil,*)'        <DataArray type="Float32"',' Name="Z_COORDINATES"',' NumberOfComponents="1">'
  write(nfil,*) (y3(k),k=1,nz)
  write(nfil,*)'        </DataArray>'
  write(nfil,*)'      </Coordinates>'
  write(nfil,*)'      <PointData Scalars="scalar">'
  write(nfil,*)'        <DataArray Name="test"',' type="Float32"',' NumberOfComponents="1"',' format="ascii">'
  write(nfil,*) (((uy(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  write(nfil,*)'        </DataArray>'
  write(nfil,*)'      </PointData>'
  write(nfil,*)'    </Piece>'
  write(nfil,*)'  </RectilinearGrid>'
  write(nfil,*)'</VTKFile>'
  close(nfil)



    end program visu
