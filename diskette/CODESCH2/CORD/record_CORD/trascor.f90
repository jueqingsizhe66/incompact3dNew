program trascor
!      include 'param.f'
      use tramod
      implicit none
      character(len=60) :: dummy

  301 format(a4)
      open(15,file='trascor.d')
      read(15,301) dummy
      read(15,*)n1
      n1m=n1-1
      dx=1./float(n1-1)                                                 
      dx=1./dx                              
      dxq=dx*dx
!
!
      call coordi 

contains
      real function dgal(al,zitaf)   ! 求那一陀三角函数的倒数关于al
!
!   analytical derivatives of the function f 
!   to compare the discrete derivatives
!
  !   use tramod,only :al
      implicit none
      real,intent(IN) :: al,zitaf
      real ::alzitf,dg1,dg2,dg3


      alzitf=al*zitaf
      dg1=1./(sinh(alzitf)*cosh(alzitf))
      dg2=-alzitf/sinh(alzitf)**2
      dg3=-alzitf/cosh(alzitf)**2
      dgal=dg1+dg2+dg3
      return
      end function dgal
      

!
      real function gal(al,zitaf)
!
!   function f for which the derivatives are evaluated
!
      implicit none
      real,intent(IN) :: al,zitaf
      real ::alzitf

      alzitf=al*zitaf
      gal=al/(cosh(alzitf)*sinh(alzitf))
      return
      end function gal

      subroutine coordi
!                                                                       
!     In this routine the grid in the physical space are
!     calculated and these are functions of the grid points
!     in the computational space.
!                                                                       
!      include 'param.f'
      use tramod
      implicit none
      real,dimension(m1)  :: xt1,xt2
      
      real :: gam,csi,eta,zita,zitaf
      real :: etac,etaf,beta,pi
      character(len=60) ::  namfile,dummy
      integer ::i,iter
      integer :: nc
      real    :: xc,xf,alin
      real    :: hbig,hsma
      real    :: al1,al2
      real    :: csic,val,res,tang,dal
      real    :: siet,azit,alet

      pi=2.*asin(1.)
  110 format(2f12.5)
!                                                                       
!     Transformation 2.4.a                                              
!                                                                       
      namfile='coord1a.plo'
      open(18,file=namfile)
  301 format(a4)
      read(15,301) dummy
      read(15,*) al
      gam=tanh(al*0.5)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      eta=al*(csi-0.5)
      xt1(i)=0.5*(1.+tanh(eta)/gam)
      y(i)=(-0.5+xt1(i))*2.
      write(18,110)csi,y(i)
      enddo
!                                                                       
!     Transformation 2.4.b      
!
      namfile='coord1b.plo'
      open(18,file=namfile) 
      read(15,301) dummy
      read(15,*)  nc,xc,xf,alin
      csic=(nc-1)/float(n1m)*0.5+0.5
      etac=abs(csic-0.5)
      etaf=0.5
      zitaf=etaf-etac
      al=alin
      beta=xc/etac
      iter=1
      val=beta/(xf-xc)
!     write(6,*)val,zitaf
!  !3 continue
      do 
      res=gal(al,zitaf)-val
      if(abs(res).lt..1e-04) then
          exit
      else
      tang=dgal(al,zitaf)   ! newton-rapson数学迭代方案
      dal=res/tang
      write(6,*)iter,res,al,tang,dal
      al=al-dal
      iter=iter+1
      if(iter.gt.21) then
          exit 
      endif
      endif
      enddo
      write(6,*)'convergence',iter,res,al
      gam=tanh(al*zitaf)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      eta=csi-0.5
      if(abs(eta).le.etac)  then
      y(i)=eta*xc/etac
                            else 
      siet=eta/abs(eta) 
      if(siet.ge.0) then
      zita=(eta-etac)     ! 有疑问地方2 应该修改为-eta/etac
      !zita=(-eta/etac)    ! 修改之后直接NAN  无法继续 
      azit=al*(zita-zitaf)
      y(i)=(xf+tanh(azit)/gam*(xf-xc))
                            else 
      zita=(eta+etac)     ! 有疑问地方3  应该修改为+eta/etac
      !zita=(eta/etac)    
      azit=al*(zita+zitaf)
      y(i)=(-xf+tanh(azit)/gam*(xf-xc))
                            endif                                   
                            endif                                   
      write(18,110)csi,y(i)
      enddo
!                                                                       
!     Transformation 2.4.c                                              
!                                                                       
      namfile='coord1c.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  al      
      gam=tanh(al)                                                      
      alet=2.*al
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      if(csi.le.0.5) then
      eta=alet*csi
      y(i)=-1+tanh(eta)/gam                                   
                     else
      eta=alet*(csi-1.)
      y(i)=1.+tanh(eta)/gam                                   
                     endif
      write(18,110)csi,y(i)
      enddo
!                                                                       
!                                                                       
!     Transformation 2.4.d                                              
!                                                                       
      namfile='coord1d.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  al      
      gam=atanh(al*0.25)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      if(csi.le.0.5) then
      eta=(0.25-csi)
      alet=al*eta
      y(i)=0.5*(-1-atanh(alet)/gam)                                  
                     else
      eta=(csi-0.75)
      alet=al*eta
      y(i)=0.5*(+1+atanh(alet)/gam)                                  
                     endif
      write(18,110)csi,y(i)
      enddo
!                                                                       
!                                                                       
!     Transformation 2.4.e                                              
!                                                                       
      namfile='coord1e.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  al,hbig,hsma
      gam=tanh(al*0.25)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      if(csi.le.0.5) then
      eta=(0.25-csi)
      alet=al*eta
      y(i)=0.5*(+1.-tanh(alet)/gam)*(hbig-hsma)
                     else
      eta=(csi-0.75)
      alet=al*eta
      y(i)=hbig-hsma+hsma*(1.+tanh(alet)/gam)                                  
                     endif
      write(18,110)csi,y(i)
      enddo
!                                                                       
!     Transformation 2.4.f                                              
!                                                                       
      namfile='coord1f.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  al1,al2,xc,xf,nc    !al1  al2为选定值 是经验值
      csic=(nc-1)/float(n1m)  !csir 的含义  在中心区
      gam=tanh(al1*csic)                                                      
      do i=1,n1
      csi=(i-1)/float(n1m)
      xt1(i)=xc/xf*tanh(al1*csi)/gam   ! xr=xc/xf  r在这边的意思就是中心圈内
      enddo
      write(6,*)'2.4.f  ',csic,xt1(n1),xt1(1)
      do i=1,n1
      csi=(i-1)/float(n1m)
      xt2(i)=1./xt1(n1)+(1.-1./xt1(n1))*tanh(al2*(csi-1.))/tanh(al2*(csic-1.))
      !xt2(i)=1./xt1(1)+(1.-1./xt1(1))*tanh(al2*(csi-1.))/tanh(al2*(csic-1.))
      ! 若是按照书本，则会是错误的结果
      !修改从  n1   -> 1
      !xt1(1) 变为了xt1(n1) n1是62 也就是总的点数
      !所以这边也是有问题的！n1/n1==1才对
      y(i)=xt1(i)*xt2(i)*xf
      write(18,110)csi,y(i)
      enddo
      write(6,*)'2.4.f  ',xt2(1),xt2(n1)
!
      return
      end subroutine coordi

end program  trascor
