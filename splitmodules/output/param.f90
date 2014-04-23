module param

use decomp_2d, only : mytype

  integer :: nclx,ncly,nclz
  integer :: ifft, ivirt,istret,iforc_entree,iturb
  integer :: itype, iskew, iin, nscheme, ifirst, ilast, iles
  integer :: isave,ilit,idebmod, imodulo, idemarre, icommence, irecord
  integer :: iscalar
  integer :: nxboite, istat,iread,iadvance_time
  real(mytype) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2
  real(mytype) :: dt,xnu,noise,noise1,pi,twopi,u1,u2,sc
  real(mytype) :: t,xxk1,xxk2
  integer :: itr,itime
  character :: filesauve*80, filenoise*80, &
       nchamp*80,filepath*80, fileturb*80, filevisu*80
  real(mytype), dimension(5) :: adt,bdt,cdt,gdt
end module param
