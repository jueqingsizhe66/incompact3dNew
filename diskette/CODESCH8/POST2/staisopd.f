c************************************************************
c****************  stfunn   **********************************
c************************************************************
c***********************************************************
c
c  evaluates the PDF for the calculation of the structure functions 
c           of velocity   or vorticity
c     this is a very long calculation because all the power are done
c
      subroutine stfunn(q,linqu,sfuma)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3)
      common/stfrms/riip(mlr),rkkp(mlr),rjjp(mlr)
      dimension qiip(22,mlr),qkkp(22,mlr),qjjp(22,mlr)
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/alpds/alasfo(201),alaqo(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
      common/nsfto/ntsfl(9,mlr),ntsfu(9,mlr)
      dimension zit1ma(mlr),zit2ma(mlr),zit3ma(mlr)
      dimension zit1mi(mlr),zit2mi(mlr),zit3mi(mlr)
      character*2 npse
      character*60 namfi3
c
c   linqu=0  for velocity
c   linqu=3  for vorticity
c
      write(6,*)'   struct funct linqu=',linqu
      nfil0=60
      do np=1,21
      write(npse,178) np
  178 format(i2.2)
      if(linqu.eq.0) namfi3='stfvel'//npse//'.out'
      if(linqu.eq.3) namfi3='stfvor'//npse//'.out'
      nfil=nfil0+np
      open(nfil,file=namfi3)      
      enddo
c
c    here the rms are calculated
c
        do k=1,nlr,nlrju
          rkkp(k)=0.
          rjjp(k)=0.
          riip(k)=0.
            do np=1,21
          qkkp(np,k)=0.
          qjjp(np,k)=0.
          qiip(np,k)=0.
        enddo
        enddo
c
c rms of struct funct in 1 direction for q_1
c
      do j=1,n2m
         do k=1,n3m
            do i= 1,n1m
        do lr= 1,nlr,nlrju
               ii=i+lr
               if(ii.gt.n1m) ii=ii-n1m
               qip= (q(1,ii,j,k)-q(1,i,j,k))
            do np=1,20
               qiip(np,lr)= qiip(np,lr)+qip**np
            enddo
               qiip(21,lr)= qiip(21,lr)+abs(qip)**3
        end do
c
c rms of struct funct in 2 direction for q_2
c
        do lr= 1,nlr,nlrju
               jj=j+lr
               if(jj.gt.n2m) jj=jj-n2m
               qjp= (q(2,i,jj,k)-q(2,i,j,k))
            do np=1,20
               qjjp(np,lr)= qjjp(np,lr)+qjp**np
            enddo
               qjjp(21,lr)= qjjp(21,lr)+abs(qjp)**3
        end do
c
c rms of struct funct in 3 direction for q_3
c
        do lr= 1,nlr,nlrju
               kk=k+lr
               if(kk.gt.n3m) kk=kk-n3m
               qkp= (q(3,i,j,kk)-q(3,i,j,k))
            do np=1,20
               qkkp(np,lr)= qkkp(np,lr)+qkp**np
            enddo
               qkkp(21,lr)= qkkp(21,lr)+abs(qkp)**3
        end do
            end do
        end do
      end do
        averg=1./float(n1m*n2m*n3m)
        do k=1,nlr,nlrju
            do np=1,21
          qkkp(np,k)=qkkp(np,k)*averg
          qjjp(np,k)=qjjp(np,k)*averg
          qiip(np,k)=qiip(np,k)*averg
            enddo
          rkkp(k)=sqrt(qkkp(2,k))
          rjjp(k)=sqrt(qjjp(2,k))
          riip(k)=sqrt(qiip(2,k))
            do np=3,20
          qkkp(np,k)=qkkp(np,k)/rkkp(k)**np
          qjjp(np,k)=qjjp(np,k)/rjjp(k)**np
          qiip(np,k)=qiip(np,k)/riip(k)**np
            enddo
          qkkp(21,k)=qkkp(21,k)/rkkp(k)**3
          qjjp(21,k)=qjjp(21,k)/rjjp(k)**3
          qiip(21,k)=qiip(21,k)/riip(k)**3
        enddo
      nfil0=60
            do np=1,21
      nfil=nfil0+np
        do k=1,nlr,nlrju
      write(nfil,*) k,qiip(np,k),qjjp(np,k),qkkp(np,k)
        enddo
      close(nfil)
            enddo
c
c PDF of normalised structure funct in 1 direction for q_1
c
              do lr= 1,nlr,nlrju
      zit1ma(lr)=-100.
      zit1mi(lr)=+100.
      zit2ma(lr)=-100.
      zit2mi(lr)=+100.
      zit3ma(lr)=-100.
      zit3mi(lr)=+100.
              enddo
      do j=1,n2m
         do k=1,n3m
          do i= 1,n1m
      lsf=1+linqu
              do lr= 1,nlr,nlrju
               ii=i+lr
               if(ii.gt.n1m) ii=ii-n1m
               zitalr=(q(1,ii,j,k)-q(1,i,j,k))/riip(lr)
          zit1ma(lr)=max(zitalr,zit1ma(lr))  
          zit1mi(lr)=min(zitalr,zit1mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
c PDF of normalised structure funct in 2 direction for q_2
c
      lsf=2+linqu
              do lr= 1,nlr,nlrju
         jj=j+lr
         if(jj.gt.n2m) jj=jj-n2m
               zitalr=(q(2,i,jj,k)-q(2,i,j,k))/rjjp(lr)
          zit2ma(lr)=max(zitalr,zit2ma(lr))  
          zit2mi(lr)=min(zitalr,zit2mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
c PDF of normalised structure funct in 3 direction for q_3
c
      lsf=3+linqu
              do lr= 1,nlr,nlrju
         kk=k+lr
         if(kk.gt.n3m) kk=kk-n3m
               zitalr=(q(3,i,j,kk)-q(3,i,j,k))/rkkp(lr)
          zit3ma(lr)=max(zitalr,zit3ma(lr))  
          zit3mi(lr)=min(zitalr,zit3mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
           end do
         end do
      end do
          do n=1,3
      lsf=n+linqu
              do lr= 1,nlr,nlrju
      ntsfu(lsf,lr)=ntsfu(lsf,lr)+nousf(lsf,lr)
      do ll=1,lipdf
      ntsfl(lsf,lr)=ntsfl(lsf,lr)+nllsf(lsf,lr,ll)
      enddo
              enddo
          enddo

              do lr= 1,nlr,nlrju
      lp=lr
      write(6,133) lp, zit1ma(lr), zit1mi(lr), zit2ma(lr), 
     1                 zit2mi(lr), zit3ma(lr), zit3mi(lr),
     1          ntsfl(1+linqu,lr),ntsfl(2+linqu,lr),ntsfl(3+linqu,lr)
  133 format(3x,i4,6e12.5,3x,3i8)
              enddo
      return
      end
c************************************************************
c****************  stfun   **********************************
c************************************************************
c***********************************************************
c
c  evaluates the PDF for the calculation of the structure functions 
c           of velocity   or vorticity
c
      subroutine stfun(q,linqu,sfuma)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3)
      dimension riip(mlr),rkkp(mlr),rjjp(mlr)
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/alpds/alasfo(201),alaqo(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
      common/nsfto/ntsfl(9,mlr),ntsfu(9,mlr)
      dimension zit1ma(mlr),zit2ma(mlr),zit3ma(mlr)
      dimension zit1mi(mlr),zit2mi(mlr),zit3mi(mlr)
c
c PDF of structure funct in 1 direction for q_1
c
              do lr= 1,nlr,nlrju
      zit1ma(lr)=-100.
      zit1mi(lr)=+100.
      zit2ma(lr)=-100.
      zit2mi(lr)=+100.
      zit3ma(lr)=-100.
      zit3mi(lr)=+100.
              enddo
      do j=1,n2m
         do k=1,n3m
          do i= 1,n1m
      lsf=1+linqu
              do lr= 1,nlr,nlrju
               ii=i+lr
               if(ii.gt.n1m) ii=ii-n1m
               zitalr=(q(1,ii,j,k)-q(1,i,j,k))
          zit1ma(lr)=max(zitalr,zit1ma(lr))  
          zit1mi(lr)=min(zitalr,zit1mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
c PDF of structure funct in 2 direction for q_2
c
      lsf=2+linqu
              do lr= 1,nlr,nlrju
         jj=j+lr
         if(jj.gt.n2m) jj=jj-n2m
               zitalr=(q(2,i,jj,k)-q(2,i,j,k))
          zit2ma(lr)=max(zitalr,zit2ma(lr))  
          zit2mi(lr)=min(zitalr,zit2mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
c PDF of  structure funct in 3 direction for q_3
c
      lsf=3+linqu
              do lr= 1,nlr,nlrju
         kk=k+lr
         if(kk.gt.n3m) kk=kk-n3m
               zitalr=(q(3,i,j,kk)-q(3,i,j,k))
          zit3ma(lr)=max(zitalr,zit3ma(lr))  
          zit3mi(lr)=min(zitalr,zit3mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
           end do
         end do
      end do
          do n=1,3
      lsf=n+linqu
              do lr= 1,nlr,nlrju
      ntsfu(lsf,lr)=ntsfu(lsf,lr)+nousf(lsf,lr)
      do ll=1,lipdf
      ntsfl(lsf,lr)=ntsfl(lsf,lr)+nllsf(lsf,lr,ll)
      enddo
              enddo
          enddo

              do lr= 1,nlr,nlrju
      lp=lr
      write(6,133) lp, zit1ma(lr), zit1mi(lr), zit2ma(lr), 
     1                 zit2mi(lr), zit3ma(lr), zit3mi(lr),
     1          ntsfl(1+linqu,lr),ntsfl(2+linqu,lr),ntsfl(3+linqu,lr)
  133 format(3x,i4,6e12.5,3x,3i8)
              enddo
      return
      end
c************************************************************
c****************  pdfqua   **********************************
c************************************************************
c***********************************************************
c
c  evaluates the PDF for the calculation of the structure functions 
c           of velocity   or vorticity
c
      subroutine pdfqua(q,qnorm,nllqu,nouqu,qmax)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/paqdf/sfuma,qquma,sfoma,qqoma
      dimension q(m1,m2,m3)
      dimension nllqu(201)
c
c PDF of normalised general quantity 
c
      zitama=-100.
      zitami=+100.
      do j=1,n2m
         do k=1,n3m
          do i= 1,n1m
               zitalr=q(i,j,k)/qnorm
          zitama=max(zitalr,zitama)  
          zitami=min(zitalr,zitami)  
               if(abs(zitalr).le.qmax) then
      all=lipdh*(qmax+zitalr)/qmax+1.5
      ll=all
      nllqu(ll)=nllqu(ll)+1
                                        else
      nouqu=nouqu+1
                                        endif
           end do
         end do
      end do
      write(6,*)' in pdfqua  ',qmax,zitama,zitami
      return
      end
c************************************************************
c****************  pdfqlo   **********************************
c************************************************************
c***********************************************************
c
c  evaluates the PDF for the calculation of the structure functions 
c           of velocity   or vorticity
c
      subroutine pdfqlo(q,qnorm,nllqu,nouqu,qmax)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/paqdf/sfuma,qquma,sfoma,qqoma
      dimension q(m1,m2,m3)
      dimension nllqu(201)
c
c PDF of normalised general quantity 
c
      zitama=-100.
      zitami=+100.
      do j=1,n2m
         do k=1,n3m
          do i= 1,n1m
               zitalr=q(i,j,k)/qnorm
          zitama=max(zitalr,zitama)  
          zitami=min(zitalr,zitami)  
c     alo=lipdn*zitalr/qmax+1.
      alo=lipdn*(1.+(alog(zitalr)/alog(10.)-1.5)/3.)+1.
      ll=alo
      if(ll.ge.1.and.ll.le.lipdn)        then
      nllqu(ll)=nllqu(ll)+1
                                        else
      nouqu=nouqu+1
                                        endif
           end do
         end do
      end do
      write(6,*)' in pdfqlo  ',qmax,zitama,zitami
      return
      end
c************************************************************
c****************  pdfpdf   **********************************
c************************************************************
c***********************************************************
c
c
      subroutine pdfpdf(q,qcap,ln,qnorm,qmax,linqu)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension qcap(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/paqdf/sfuma,qquma,sfoma,qqoma
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/alpds/alasfo(201),alaqo(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
      dimension nllq(201)
      common/nouto/ntqll(7),ntqou(7),ntsll(7),ntsou(7)
      common/pdfqdf/pdfsf(9,mlr,201),pdfqu(7,201),pdfsc(7,201)
c
c
c    PDF of ln component of a vector quantity 
c       linqu=0   for velocity   linqu=3  for vorticity
c
      n=ln+linqu
      ntqll(n)=0
      ntqou(n)=0
      do ll=1,lipdf
      nllq(ll)=0
      nouq=0
      enddo
      vmax=0.
      do i=1,n1m
      do j=1,n2m
      do k=1,n3m
      qcap(i,j,k)=q(ln,i,j,k)
      vmax=max(abs(qcap(i,j,k)),vmax)
      enddo
      enddo
      enddo
      call pdfqua(qcap,qnorm,nllq,nouq,qmax)
      do ll=1,lipdf
      nllqu(n,ll)=nllq(ll)
      nouqu(n)=nouq
      ntqll(n)=ntqll(n)+nllq(ll)
      ntqou(n)=ntqou(n)+nouq
      enddo
      write(18,131)ln,linqu,n,ntqll(n),ntqou(n),vmax
      do ll=1,lipdf
      pdfqu(n,ll)=nllqu(n,ll)/float(ntqll(n))
c     write(18,133)ll,nllqu(n,ll),pdfqu(n,ll)
  133 format(3x,i4,3x,i10,3x,e12.5)
      enddo
      write(6,131)ln,linqu,n,ntqll(n),ntqou(n),vmax,lipdf,lipdn,lipdh
  131 format(3x,'in pdfpdf',3x,5i10,2x,e12.5,3x,3i4)
      return
      end
c
c  ****************************** subrout pdfini **********************
c in this routine the abscissa of the PDF are calculated
c depending on the interval between minimum and maximum
c
      subroutine pdfini   
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/paqdf/sfuma,qquma,sfoma,qqoma
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/alpds/alasfo(201),alaqo(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
c************************************************************************
      lipdh=(lipdf-1)/2
      lipdn=lipdf-1
      open(18,file='pdf.out')
      write(6,*)' in pdfini lipdf,lipdh,lipdn',lipdf,lipdh,lipdn
c
c   for structure function of u_i  om_i and pr
c
      do ll=1,lipdf
      alasfu(ll)=(ll-1-lipdh)/float(lipdh)*sfuma
      alasfo(ll)=(ll-1-lipdh)/float(lipdh)*sfoma
      enddo
         do ll=1,lipdf
            do lr=1,nlr,nlrju
               do np=1,9
      nllsf(np,lr,ll)=0
      nousf(np,lr)=0
               enddo
            enddo
         enddo
c
c   for pdf of u_i  om_i    pr
c
      do ll=1,lipdn
      alaqu(ll)=(ll-0.5-lipdh)/float(lipdh)*qquma
      alaqo(ll)=(ll-0.5-lipdh)/float(lipdh)*qqoma
      enddo
               do np=1,7
      nouqu(np)=0
         do ll=1,lipdf
      nllqu(np,ll)=0
         enddo
               enddo
c
c   for pdf of helicity and dissipation
c
      do ll=1,lipdn
      alahe(ll)=(ll-0.5-lipdh)/float(lipdh)
      enddo
      do ll=1,lipdn
c     aladi(ll)=(ll-0.5)/float(lipdn)*sfuma
      aladi(ll)=(1.5-3.*(1.-(ll-0.5)/float(lipdn)))
      enddo
               do np=1,2
      nousc(np)=0
         do ll=1,lipdf
      nllsc(np,ll)=0
         enddo
               enddo
      return
      end
c
c  ****************************** subrout pdffiq **********************
c  the pdf of the statistics are printed
c
c
      subroutine pdffiq(time)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/paqdf/sfuma,qquma,sfoma,qqoma
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/alpds/alasfo(201),alaqo(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
      common/nouto/ntqll(7),ntqou(7),ntsll(7),ntsou(7)
      common/nsfto/ntsfl(9,mlr),ntsfu(9,mlr)
      common/pdfqdf/pdfsf(9,mlr,201),pdfqu(7,201),pdfsc(7,201)
      dimension pdfsft(9,mlr),pdfqut(7),pdfsct(7)
      common/pdfcog/pdfco(6,201),pdfcot(6)
c      
      character*4 pntim
      character*60 namfi3
      character*3 njpse
      itim=nint(time)
      write(pntim,67) itim
   67 format(i4.4)
  177 format(i3.3)
  133 format(15x,10i10)
  198 format(2x,e12.4,3x,6e12.5)
      do n=1,7
      pdfqut(n)=0.
      enddo
      write(18,133) (ntqll(n),n=1,7)
      write(18,133) (ntqou(n),n=1,7)
      namfi3='pdfvor'//'.'//pntim
      open(25,file=namfi3)
      namfi3='pdfvel'//'.'//pntim
      open(27,file=namfi3)
      namfi3='pdfpre'//'.'//pntim
      open(28,file=namfi3)
                 do ll=1,lipdn
      write(27,198)alaqu(ll),(pdfqu(n,ll),n=1,3)
      write(25,198)alaqo(ll),(pdfqu(n,ll),n=4,6)
      write(28,198)alaqo(ll),pdfqu(7,ll)
      do n=1,7
      pdfqut(n)=pdfqut(n)+pdfqu(n,ll)
      enddo
                  enddo
      write(18,197) (pdfqut(n),n=1,7)
      close(28)
      close(27)
      close(25)
      do n=1,2
      pdfsct(n)=0.
      enddo
      write(18,133) (ntsll(n),n=1,2)
      write(18,133) (ntsou(n),n=1,2)
      namfi3='pdfhel'//'.'//pntim
      open(25,file=namfi3)
      namfi3='pdfdis'//'.'//pntim
      open(27,file=namfi3)
                 do ll=1,lipdn
      alevc=alahe(ll)
      write(25,198)alahe(ll),pdfsc(1,ll)
      write(27,198)aladi(ll),pdfsc(2,ll)
      do n=1,6
      pdfsct(n)=pdfsct(n)+pdfsc(n,ll)
      enddo
                  enddo
      write(18,197) (pdfsct(n),n=1,2)
  197 format(2x,3x,8e12.5)
      close(27)
      close(25)
      namfi3='pdfheco'//'.'//pntim
      open(25,file=namfi3)
      namfi3='pdflaco'//'.'//pntim
      open(26,file=namfi3)
                 do ll=1,lipdn
      write(25,198)alahe(ll),pdfco(1,ll),pdfco(2,ll),pdfco(3,ll)
      write(26,198)alahe(ll),pdfco(4,ll),pdfco(5,ll),pdfco(6,ll)
      do n=1,6
      pdfcot(n)=pdfcot(n)+pdfco(n,ll)
      enddo
                  enddo
      write(18,197) (pdfcot(n),n=1,6)
      close(26)
      close(25)
      return
      end

c
c  ****************************** subrout pdffis **********************
c   the PDF of the velocity differences are printed
c
c
      subroutine pdffis(time)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/paqdf/sfuma,qquma,sfoma,qqoma
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/alpds/alasfo(201),alaqo(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
      common/nouto/ntqll(7),ntqou(7),ntsll(7),ntsou(7)
      common/nsfto/ntsfl(9,mlr),ntsfu(9,mlr)
      common/pdfqdf/pdfsf(9,mlr,201),pdfqu(7,201),pdfsc(7,201)
      dimension pdfsft(9,mlr),pdfqut(7),pdfsct(7)
c      
      character*4 pntim
      character*60 namfi3
      character*3 njpse
      itim=nint(time)
      write(pntim,67) itim
   67 format(i4.4)
      write(18,*) ' pdf structure funct '
           do jc=1,nlr,nlrju
      do n=1,9
      pdfsft(n,jc)=0.
      enddo
      write(njpse,177) jc
  177 format(i3.3)
  132 format(i5,5x,10i10)
  133 format(15x,10i10)
      write(18,132) jc,(ntsfl(n,jc),n=1,9)
      write(18,132) jc,(ntsfu(n,jc),n=1,9)
      namfi3='stfuvor'//njpse//'.'//pntim
      open(25,file=namfi3)
      namfi3='stfuvel'//njpse//'.'//pntim
      open(27,file=namfi3)
      namfi3='stfupre'//njpse//'.'//pntim
      open(28,file=namfi3)
                 do ll=1,lipdn
      write(27,*)alasfu(ll),(pdfsf(n,jc,ll),n=1,3)
      write(25,*)alasfo(ll),(pdfsf(n,jc,ll),n=4,6)
      write(28,*)alasfo(ll),(pdfsf(n,jc,ll),n=7,9)
      do n=1,9
      pdfsft(n,jc)=pdfsft(n,jc)+pdfsf(n,jc,ll)
      enddo
                 enddo
      close(25)
      close(27)
      rc=jc/dx1
      write(18,198) rc,(pdfsft(n,jc),n=1,9)
  198 format(2x,e12.4,3x,9e12.5)
           enddo
      return
      end
c************************************************************
c****************  sfprnn   **********************************
c************************************************************
c***********************************************************
c
c  evaluates the PDF for the calculation of the structure functions 
c           of pressure   
c
      subroutine sfprnn(pr,sfuma)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      common/stfrms/riip(mlr),rkkp(mlr),rjjp(mlr)
      dimension qiip(20,mlr),qkkp(20,mlr),qjjp(20,mlr)
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/alpds/alasfo(201),alaqo(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
      common/nsfto/ntsfl(9,mlr),ntsfu(9,mlr)
      dimension zit1ma(mlr),zit2ma(mlr),zit3ma(mlr)
      dimension zit1mi(mlr),zit2mi(mlr),zit3mi(mlr)
      character*2 npse
      character*60 namfi3
c
      write(6,*)'   struct funct pressure'
      nfil0=60
      do np=1,20
      write(npse,178) np
  178 format(i2.2)
      namfi3='stfpre'//npse//'.out'
      nfil=nfil0+np
      open(nfil,file=namfi3)      
      enddo
        do k=1,nlr,nlrju
          rkkp(k)=0.
          rjjp(k)=0.
          riip(k)=0.
            do np=1,20
          qkkp(np,k)=0.
          qjjp(np,k)=0.
          qiip(np,k)=0.
        enddo
        enddo
c
c rms of struct funct in 1 direction for pr
c
      do j=1,n2m
         do k=1,n3m
            do i= 1,n1m
        do lr= 1,nlr,nlrju
               ii=i+lr
               if(ii.gt.n1m) ii=ii-n1m
               qip= (pr(ii,j,k)-pr(i,j,k))
            do np=1,20
               qiip(np,lr)= qiip(np,lr)+qip**np
            enddo
        end do
c
c rms of struct funct in 2 direction for pr
c
        do lr= 1,nlr,nlrju
               jj=j+lr
               if(jj.gt.n2m) jj=jj-n2m
               qjp= (pr(i,jj,k)-pr(i,j,k))
            do np=1,20
               qjjp(np,lr)= qjjp(np,lr)+qjp**np
            enddo
        end do
c
c rms of struct funct in 3 direction for q_3
c
        do lr= 1,nlr,nlrju
               kk=k+lr
               if(kk.gt.n3m) kk=kk-n3m
               qkp= (pr(i,j,kk)-pr(i,j,k))
            do np=1,20
               qkkp(np,lr)= qkkp(np,lr)+qkp**np
            enddo
        end do
            end do
        end do
      end do
        averg=1./float(n1m*n2m*n3m)
        do k=1,nlr,nlrju
            do np=1,20
          qkkp(np,k)=qkkp(np,k)*averg
          qjjp(np,k)=qjjp(np,k)*averg
          qiip(np,k)=qiip(np,k)*averg
            enddo
          rkkp(k)=sqrt(qkkp(2,k))
          rjjp(k)=sqrt(qjjp(2,k))
          riip(k)=sqrt(qiip(2,k))
            do np=3,20
          qkkp(np,k)=qkkp(np,k)/rkkp(k)**np
          qjjp(np,k)=qjjp(np,k)/rjjp(k)**np
          qiip(np,k)=qiip(np,k)/riip(k)**np
            enddo
        enddo
      nfil0=60
            do np=1,20
      nfil=nfil0+np
        do k=1,nlr,nlrju
      write(nfil,133) k,qiip(np,k),qjjp(np,k),qkkp(np,k)
        enddo
      close(nfil)
            enddo
c
c PDF of normalised structure funct in 1 direction for pr
c
      linqu=6
              do lr= 1,nlr,nlrju
      zit1ma(lr)=-100.
      zit1mi(lr)=+100.
      zit2ma(lr)=-100.
      zit2mi(lr)=+100.
      zit3ma(lr)=-100.
      zit3mi(lr)=+100.
              enddo
      do j=1,n2m
         do k=1,n3m
          do i= 1,n1m
      lsf=linqu+1
              do lr= 1,nlr,nlrju
               ii=i+lr
               if(ii.gt.n1m) ii=ii-n1m
               zitalr=(pr(ii,j,k)-pr(i,j,k))/riip(lr)
          zit1ma(lr)=max(zitalr,zit1ma(lr))  
          zit1mi(lr)=min(zitalr,zit1mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
c PDF of normalised structure funct in 2 direction for pr
c
      lsf=linqu+2
              do lr= 1,nlr,nlrju
         jj=j+lr
         if(jj.gt.n2m) jj=jj-n2m
               zitalr=(pr(i,jj,k)-pr(i,j,k))/rjjp(lr)
          zit2ma(lr)=max(zitalr,zit2ma(lr))  
          zit2mi(lr)=min(zitalr,zit2mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
c PDF of normalised structure funct in 3 direction for pr
c
      lsf=linqu+3
              do lr= 1,nlr,nlrju
         kk=k+lr
         if(kk.gt.n3m) kk=kk-n3m
               zitalr=(pr(i,j,kk)-pr(i,j,k))/rkkp(lr)
          zit3ma(lr)=max(zitalr,zit3ma(lr))  
          zit3mi(lr)=min(zitalr,zit3mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
           end do
         end do
      end do
          do n=1,3
      lsf=n+linqu
              do lr= 1,nlr,nlrju
      ntsfu(lsf,lr)=ntsfu(lsf,lr)+nousf(lsf,lr)
      do ll=1,lipdf
      ntsfl(lsf,lr)=ntsfl(lsf,lr)+nllsf(lsf,lr,ll)
      enddo
              enddo
          enddo

              do lr= 1,nlr,nlrju
      lp=lr
      write(6,133) lp, zit1ma(lr), zit1mi(lr), zit2ma(lr), 
     1                 zit2mi(lr), zit3ma(lr), zit3mi(lr),
     1          ntsfl(1+linqu,lr),ntsfl(2+linqu,lr),ntsfl(3+linqu,lr)
  133 format(3x,i4,6e12.5,3x,3i8)
              enddo
      return
      end
c************************************************************
c****************  sfprn   **********************************
c************************************************************
c***********************************************************
c
c  evaluates the PDF for the calculation of the structure functions 
c           of pressure 
c
      subroutine sfprn(pr,sfuma)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension riip(mlr),rkkp(mlr),rjjp(mlr)
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/alpds/alasfo(201),alaqo(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
      common/nsfto/ntsfl(9,mlr),ntsfu(9,mlr)
      dimension zit1ma(mlr),zit2ma(mlr),zit3ma(mlr)
      dimension zit1mi(mlr),zit2mi(mlr),zit3mi(mlr)
c
c PDF of structure funct in 1 direction for pr
c
      linqu=6
              do lr= 1,nlr,nlrju
      zit1ma(lr)=-100.
      zit1mi(lr)=+100.
      zit2ma(lr)=-100.
      zit2mi(lr)=+100.
      zit3ma(lr)=-100.
      zit3mi(lr)=+100.
              enddo
      do j=1,n2m
         do k=1,n3m
          do i= 1,n1m
      lsf=1+linqu
              do lr= 1,nlr,nlrju
               ii=i+lr
               if(ii.gt.n1m) ii=ii-n1m
               zitalr=(pr(ii,j,k)-pr(i,j,k))
          zit1ma(lr)=max(zitalr,zit1ma(lr))  
          zit1mi(lr)=min(zitalr,zit1mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
c PDF of structure funct in 2 direction for pr
c
      lsf=2+linqu
              do lr= 1,nlr,nlrju
         jj=j+lr
         if(jj.gt.n2m) jj=jj-n2m
               zitalr=(pr(i,jj,k)-pr(i,j,k))
          zit2ma(lr)=max(zitalr,zit2ma(lr))  
          zit2mi(lr)=min(zitalr,zit2mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
c PDF of structure funct in 3 direction for pr
c
      lsf=3+linqu
              do lr= 1,nlr,nlrju
         kk=k+lr
         if(kk.gt.n3m) kk=kk-n3m
               zitalr=(pr(i,j,kk)-pr(i,j,k))
          zit3ma(lr)=max(zitalr,zit3ma(lr))  
          zit3mi(lr)=min(zitalr,zit3mi(lr))  
               if(abs(zitalr).le.sfuma) then
      all=lipdh*(sfuma+zitalr)/sfuma+1.5
      ll=all
      nllsf(lsf,lr,ll)=nllsf(lsf,lr,ll)+1
                                        else
      nousf(lsf,lr)=nousf(lsf,lr)+1
                                        endif
              end do
c
           end do
         end do
      end do
          do n=1,3
      lsf=n+linqu
              do lr= 1,nlr,nlrju
      ntsfu(lsf,lr)=ntsfu(lsf,lr)+nousf(lsf,lr)
      do ll=1,lipdf
      ntsfl(lsf,lr)=ntsfl(lsf,lr)+nllsf(lsf,lr,ll)
      enddo
              enddo
          enddo

              do lr= 1,nlr,nlrju
      lp=lr
      write(6,133) lp, zit1ma(lr), zit1mi(lr), zit2ma(lr), 
     1                 zit2mi(lr), zit3ma(lr), zit3mi(lr),
     1          ntsfl(1+linqu,lr),ntsfl(2+linqu,lr),ntsfl(3+linqu,lr)
  133 format(3x,i4,6e12.5,3x,3i8)
              enddo
      return
      end
