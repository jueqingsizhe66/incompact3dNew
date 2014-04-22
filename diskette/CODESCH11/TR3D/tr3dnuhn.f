c***********************************************************************
c
c     h term for the q1 momentum equation at i,j+1/2,k+1/2
c
c***********************************************************************
      subroutine hdnl1
      include 'param.f'
      do 10 kc=1,n3m
      udx3 = dx3/g3rm(kc)*0.25
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 10 jc=2,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      jp=jc+1
      jm=jc-1
      udx1 = dx1/rm(jc)*0.25
      udx1vis=2.*dx1/rm(jc)**2/ren
      udx2=dx2/(g2rm(jc)*rm(jc))*0.25
      do 10 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c    q1 q1 term
c
c
c             1   d  q_t q_t 
c            --- -----------
c             r   d   t      
c
      h11=( (q1(ip,jc,kc)+q1(ic,jc,kc))*(q1(ip,jc,kc)+q1(ic,jc,kc))
     %     -(q1(im,jc,kc)+q1(ic,jc,kc))*(q1(im,jc,kc)+q1(ic,jc,kc))
     %    )*udx1
c
c    q1 q2 term
c
c             1   d  q_t q_r     q_t q_r
c            --- -----------  +  --------    jc > 2
c             r   d   r            r^2
c
c
      q1nn=(q1(ic,jpp,kc)+q1(ic,jc,kc))
      q1ss=(q1(ic,jc,kc)+q1(ic,jmm,kc))
      h12d=( (q2(ic,jp,kc)+q2(im,jp,kc))*q1nn
     %      -(q2(ic,jc,kc)+q2(im,jc,kc))*q1ss
     %     )*udx2
      h12n=  q1(ic,jc,kc)/rm(jc) * (
     %                  (q2(ic,jp,kc)+q2(im,jp,kc))/rc(jp)
     %                 +(q2(ic,jc,kc)+q2(im,jc,kc))/rc(jc)  )*0.25
       h12=h12d+h12n
c
c   q1 q3 term
c
c
c                 d  q_t q_x 
c                -----------
c                 d   x      
c
      h13=((q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kpp)+q1(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm))
     %    )*udx3
      hq1=h11+h12+h13
c
c   first derivative of q2 with respect to x1
c
c
c              2      d  q_r / r
c            ------  -----------
c            Re r^2   d   t      
c
      q2e=(q2(ic,jp,kc)/rc(jp)+q2(ic,jc,kc)/rc(jc))*0.5
      q2w=(q2(im,jp,kc)/rc(jp)+q2(im,jc,kc)/rc(jc))*0.5
      d11q2e=(q2e-q2w)*udx1vis
c
c   coriolis term
c
c
c              1      q_r 
c           - -----  -----
c              Ro      r      
c
      coriol= - ((q2(ic,jc,kc)+q2(im,jc,kc))/rc(jc) +
     %           (q2(ic,jp,kc)+q2(im,jp,kc))/rc(jp)  )*0.25*ros
c
      dq(ic,jc,kc)=d11q2e-hq1+coriol
   10 continue
c
c     JC = 1   SPECIAL TREATMENT IS NEEDED FOR RADIAL AND VISCOUS 
c     DERIVATIVES AT THE AXIS
c
      jc = 1
      jp = jc + 1
      jmm=jmv(jc)
      jpp=jpv(jc)
      udx1 = dx1/rm(jc)*0.25
      udx1vis=2.*dx1/rm(jc)**2/ren
      udx2=dx2/(g2rm(jc)*rm(jc))*0.25
      do 11 kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 11 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      insy = isym(ic)
      insym = isym(im)
c
c    q1 q1 term
c
c
c             1   d  q_t q_t 
c            --- -----------
c             r   d   t      
c
      h11=( (q1(ip,jc,kc)+q1(ic,jc,kc))*(q1(ip,jc,kc)+q1(ic,jc,kc))
     %     -(q1(im,jc,kc)+q1(ic,jc,kc))*(q1(im,jc,kc)+q1(ic,jc,kc))
     %    )*udx1
c
c    q1 q2 term
c
c             1   d  q_t q_r     q_t q_r
c            --- -----------  +  --------    jc > 2
c             r   d   r            r^2
c
c
      q1nn=(q1(ic,jpp,kc)+q1(ic,jc,kc))
      q1ss=(q1(ic,jc,kc)+q1(ic,jmm,kc))
      h12d=( (q2(ic,jp,kc)+q2(im,jp,kc))*q1nn
     %      -(q2(ic,jc,kc)+q2(im,jc,kc))*q1ss
     %     )*udx2
       q2s1 = (q2(ic,jp,kc) - q2(insy,jp,kc))*0.5/rc(jp)
       q2s2 = (q2(im,jp,kc) - q2(insym,jp,kc))*0.5/rc(jp)
      h12n=  q1(ic,jc,kc) / rm(jc) * (
     %        ((q2(ic,jp,kc)+q2(im,jp,kc))/rc(jp)+q2s1 + q2s2)
     %      )*0.25
       h12=h12d+h12n
c
c   q1 q3 term
c
c
c                 d  q_t q_x 
c                -----------
c                 d   x      
c
      h13=((q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kpp)+q1(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm))
     %    )*udx3
      hq1=h11+h12+h13
c
c   first derivative of q2 with respect to x1
c
c
c              2      d  q_r / r
c            ------  -----------
c            Re r^2   d   t      
c
      q2e=(q2(ic,jp,kc)/rc(jp)+q2s1)*0.5
      q2w=(q2(im,jp,kc)/rc(jp)+q2s2)*0.5
      d11q2e=(q2e-q2w)*udx1vis
c
c
c   coriolis term
c
c
c              1      q_r 
c           - -----  -----
c              Ro      r      
c
      coriol= - ( q2s1        +q2s2 +
     %           (q2(ic,jp,kc)+q2(im,jp,kc))/rc(jp)  )*0.25*ros
c
      dq(ic,jc,kc)=d11q2e-hq1+coriol
   11 continue
      return
      end
c
c***********************************************************************
c
c     h term for the q2 momentum equation at i+1/2,j,k+1/2
c
c***********************************************************************
      subroutine hdnl2
      include 'param.f'
      do 20 kc=1,n3m
      udx3 = dx3/g3rm(kc)*0.25
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 20 jc=3,n2m
      jm=jc-1
      jp=jc+1
      udx1 = dx1*0.25/rc(jc)
      udx2 = 1./g2rc(jc)*dx2*0.25
      udx1vis = dx1/(rc(jc)*ren)
      do 20 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
c
c     q2 q1 term
c
c
c             1   d  q_t q_r 
c            --- -----------
c             r   d   t      
c
      h21=((q1(ipp,jc,kc)+q1(ipp,jm,kc))
     %     *(q2(ipp,jc,kc)+q2(ic,jc,kc))
     %     -(q1(ic,jc,kc)+q1(ic,jm,kc))
     %     *(q2(ic,jc,kc)+q2(imm,jc,kc)))*udx1
c
c     q2 q2 term
c
c
c                 d  q_r q_r / r
c                ---------------
c                 d   r      
c
      h22=( (q2(ic,jp,kc)/rc(jp)+q2(ic,jc,kc)/rc(jc))
     %     *(q2(ic,jp,kc)+q2(ic,jc,kc))
     %     -(q2(ic,jc,kc)/rc(jc)+q2(ic,jm,kc)/rc(jm))
     %     *(q2(ic,jc,kc)+q2(ic,jm,kc))
     %    )*udx2
c
c   q1 q1 term  not differentiated
c
c             ( q_t )^2
c
      h11n=((q1(ipp,jc,kc)+q1(ipp,jm,kc)
     %       +q1(ic,jc,kc)+q1(ic,jm,kc))*0.25)**2
c
c     q2 q3 term
c
c
c                 d  q_x q_r 
c                -----------
c                 d   x      
c
      h23=((q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kpp)+q2(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm))
     %    )*udx3
      hq2=h21+h22+h23-h11n
c
c   second derivative of q1 with respect to x1
c
c
c              2      d  q_t    
c          - ------  --------
c            Re r     d   t      
c
      q1e=(q1(ipp,jc,kc)+q1(ipp,jm,kc))
      q1w=(q1(ic,jc,kc)+q1(ic,jm,kc))
      d11q2e=-(q1e-q1w)*udx1vis
c
c
c   coriolis term
c
c
c              1            
c             -----   r q_t
c              Ro             
c
      coriol=(q1(ic,jc,kc)+q1(ipp,jc,kc)+q1(ic,jm,kc)+
     %        q1(ipp,jm,kc))*0.25*rc(jc)*ros
c
      dph(ic,jc,kc)=d11q2e-hq2+coriol
c
   20 continue
c
c     JC = 2   SPECIAL TREATMENT IS NEEDED FOR RADIAL 
c     DERIVATIVE AT THE AXIS
c
      jc=2
      jm=jc-1
      jp=jc+1
      udx1 = dx1*0.25/rc(jc)
      udx2 = 1./g2rc(jc)*dx2*0.25
      udx1vis = dx1/(rc(jc)*ren)
      do 21 kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 21 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
      insy = isym(ic)
c
c     q2 q1 term
c
c
c             1   d  q_t q_r 
c            --- -----------
c             r   d   t      
c
      h21=((q1(ipp,jc,kc)+q1(ipp,jm,kc))
     %     *(q2(ipp,jc,kc)+q2(ic,jc,kc))
     %     -(q1(ic,jc,kc)+q1(ic,jm,kc))
     %     *(q2(ic,jc,kc)+q2(imm,jc,kc)))*udx1
c
c     q2 q2 term
c
c
c                 d  q_r q_r / r
c                ---------------
c                 d   r      
c
       q2s1 = (q2(ic,jc,kc) - q2(insy,jc,kc))*0.5/rc(jc)
      h22=( (q2(ic,jp,kc)/rc(jp)+q2(ic,jc,kc)/rc(jc))
     %     *(q2(ic,jp,kc)+q2(ic,jc,kc))
     %     -(q2(ic,jc,kc)/rc(jc)+q2s1 )
     %     *(q2(ic,jc,kc)             )
     %    )*udx2
c
c   q1 q1 term  not differentiated
c
c             ( q_t )^2
c
      h11n=((q1(ipp,jc,kc)+q1(ipp,jm,kc)
     %       +q1(ic,jc,kc)+q1(ic,jm,kc))*0.25)**2
c
c     q2 q3 term
c
c
c                 d  q_x q_r 
c                -----------
c                 d   x      
c
      h23=((q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kpp)+q2(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm))
     %    )*udx3
      hq2=h21+h22+h23-h11n
c
c   second derivative of q1 with respect to x1
c
c
c              2      d  q_t    
c          - ------  --------
c            Re r     d   t      
c
      q1e=(q1(ipp,jc,kc)+q1(ipp,jm,kc))
      q1w=(q1(ic,jc,kc)+q1(ic,jm,kc))
      d11q2e=-(q1e-q1w)*udx1vis
c
c
c   coriolis term
c
c
c              1            
c             -----   r q_t
c              Ro             
c
      coriol=(q1(ic,jc,kc)+q1(ipp,jc,kc)+q1(ic,jm,kc)+
     %        q1(ipp,jm,kc))*0.25*rc(jc)*ros
c
      dph(ic,jc,kc)=d11q2e-hq2+coriol
c
   21 continue
      return
      end
c
c***********************************************************************
c
c     h term for the q3 momentum equation at i+1/2,j+1/2,k
c
c***********************************************************************
      subroutine hdnl3
      include 'param.f'
      do 30 kc=2,n3m
      udx3 = dx3/g3rc(kc)*0.25
      km=kmv(kc)
      kp=kc+1
      do 30 jc=1,n2m
      jp=jc+1
      jmm=jmv(jc)
      jpp=jpv(jc)
      udx2=dx2/(g2rm(jc)*rm(jc))*0.25
      udx1=dx1/rm(jc)*0.25
      do 30 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
c
c    q3 q1 term
c
c
c            1    d  q_x q_t 
c           ---  -----------
c            r    d   t      
c
      h31=((q1(ipp,jc,kc)+q1(ipp,jc,km))
     %        *(q3(ipp,jc,kc)+q3(ic,jc,kc))
     %        -(q1(ic,jc,kc)+q1(ic,jc,km))
     %        *(q3(ic,jc,kc)+q3(imm,jc,kc)))*udx1
c
c    q3 q2 term
c
c
c            1    d  q_x q_r 
c           ---  -----------
c            r    d   r      
c
      h32=(((q2(ic,jp,kc)+q2(ic,jp,km))
     %     *(q3(ic,jpp,kc)+q3(ic,jc,kc)))
     %    -((q2(ic,jc,kc)+q2(ic,jc,km))
     %     *(q3(ic,jc,kc)+q3(ic,jmm,kc))))*udx2
c
c    q3 q3 term
c
c
c                 d  q_x q_x 
c                -----------
c                 d   x      
c
      h33=((q3(ic,jc,kp)+q3(ic,jc,kc))*(q3(ic,jc,kp)+q3(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jc,km))
     %    )*udx3
c
      qcap(ic,jc,kc)=-(h31+h32+h33)
   30 continue
      return
      end
c
c***********************************************************************
c
c     h term for the passive scalar at i+1/2,j+1/2,k+1/2
c
c***********************************************************************
      subroutine hdnlps
      include 'param.f'
      do 30 kc=1,n3m
      udx3 = dx3/g3rm(kc)*0.5
      km=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 30 jc=1,n2m
      jp=jc+1
      jmm=jmv(jc)
      jpp=jpv(jc)
      udx2=dx2/(g2rm(jc)*rm(jc))*0.5
      udx1=dx1/rm(jc)*0.5
      do 30 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
c
c    rho q1 term
c
c
c            1    d  rho q_t 
c           ---  -----------
c            r    d   t      
c
      h31=(q1(ipp,jc,kc)*(psc(ipp,jc,kc)+psc(ic,jc,kc))-
     %     q1(ic,jc,kc)*(psc(ic,jc,kc)+psc(imm,jc,kc))
     %    )*udx1
c
c    rho q2 term
c
c
c            1    d  rho q_r 
c           ---  -----------
c            r    d   r      
c
      h32=(q2(ic,jp,kc)*(psc(ic,jpp,kc)+psc(ic,jc,kc))-
     %     q2(ic,jc,kc)*(psc(ic,jc,kc)+psc(ic,jmm,kc))
     %    )*udx2
c
c    rho q3 term
c
c
c                 d  rho q_x 
c                -----------
c                 d   x      
c
      h33=(q3(ic,jc,kp)*(psc(ic,jc,kpp)+psc(ic,jc,kc))-
     %     q3(ic,jc,kc)*(psc(ic,jc,kc)+psc(ic,jc,km))
     %    )*udx3
      qcap(ic,jc,kc)=-(h31+h32+h33)
   30 continue
      return
      end
