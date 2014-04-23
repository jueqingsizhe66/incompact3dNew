diff --git a/channel/poisson.f90 b/channel/poisson.f90
index 1fad3f7..6309b94 100755
--- a/channel/poisson.f90
+++ b/channel/poisson.f90
@@ -43,6 +43,7 @@ module decomp_2d_poisson
   private        ! Make everything private unless declared public
 
 !  real(mytype), private, parameter :: PI = 3.14159265358979323846_mytype
+!                                   This 0.0_mytype!!!!!!!!!!!!!1111
 
 #ifdef DOUBLE_PREC
   real(mytype), parameter :: epsilon = 1.e-16
@@ -91,6 +92,9 @@ contains
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Initialise Poisson solver for given boundary conditions
+  !                               given
+  !                               given
+  !   just for init
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine decomp_2d_poisson_init(bcx1, bcy1, bcz1)
 
@@ -129,6 +133,7 @@ contains
        allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
        allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
        allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
+       ! (000)  ^c  ^kxyz ^a
     else if (bcx==1 .and. bcy==0 .and. bcz==0) then
        allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                  sp%xst(3):sp%xen(3)))
@@ -145,6 +150,7 @@ contains
        allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
        allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
        allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
+       !(100)      consist :: cw1 cw1b   rw1   rw1b   rw2    ^a
     else if (bcx==0 .and. bcy==1 .and. bcz==0) then
        allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
             ph%yst(3):ph%yen(3)))
@@ -165,6 +171,7 @@ contains
        allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
        allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
        allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
+       !(010)  rw2 rw2b cw2 cw22 cw2b cw2c kxyz  ^a
     else if (bcx==1 .and. bcy==1) then
        allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                  sp%xst(3):sp%xen(3)))
@@ -186,6 +193,7 @@ contains
             ph%yst(3):ph%yen(3)))
        allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
             ph%yst(3):ph%yen(3)))
+        !(11*) cw1 cw1b cw2 cw2b  cw2c  rw1 rw1b rw2b
        if (bcz==1) then  
           allocate(rw3(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
        end if
@@ -194,6 +202,7 @@ contains
        allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
        allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
        allocate(a3(sp%yst(1):sp%yen(1),nym,sp%yst(3):sp%yen(3),5))      
+       !(111) rw3  kxyz  ^a
     end if 
 
     call waves()
@@ -221,14 +230,19 @@ contains
 
     if (bcx==0 .and. bcy==0 .and. bcz==0) then
        deallocate(cw1)
+       ! (000) cw1
     else if (bcx==1 .and. bcy==0 .and. bcz==0) then
        deallocate(cw1,cw1b,rw1,rw1b,rw2)
+       ! (100)   cw1 cw1b rw1 rw1b rw2
     else if (bcx==0 .and. bcy==1 .and. bcz==0) then
        deallocate(cw1,cw2,cw2b,rw2,rw2b)
+       ! (010)   cw1 cw2  cw2b  rw2 rw2b
     else if (bcx==1 .and. bcy==1) then
        deallocate(cw1,cw1b,cw2,cw2b,rw1,rw1b,rw2,rw2b)
+       ! (11*)   cw1,cw1b,cw2,cw2b,rw1,rw1b,rw2,rw2b
        if (bcz==1) then
           deallocate(rw3)
+        !(111)    rw3
        end if
     end if
 
@@ -265,6 +279,8 @@ contains
 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Solving 3D Poisson equation with periodic B.C in all 3 dimensions
+  !                                                      3
+  !                                                      3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine poisson_000(rhs)
 
@@ -297,10 +313,13 @@ contains
        fft_initialised = .true.
     end if
 
-    ! compute r2c transform 
+    ! compute r2c transform   r2c:real to complex by zhaoliang
+    ! cw1 is 3dim data structrue
     call decomp_2d_fft_3d(rhs,cw1)
+    
 
     ! normalisation
+    !  why this below is called normalisation ??????  just divide 3dim grid number
     cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
          / real(nz, kind=mytype)
 
@@ -311,17 +330,20 @@ contains
              ! post-processing in spectral space
 
              ! POST PROCESSING IN Z
-             tmp1 = real(cw1(i,j,k), kind=mytype)
-             tmp2 = aimag(cw1(i,j,k))
+             tmp1 = real(cw1(i,j,k), kind=mytype)   !  get the real number of cw1
+             tmp2 = aimag(cw1(i,j,k))               !  get the imagenumber of cw1
              cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
-                  tmp2*bz(k)-tmp1*az(k), kind=mytype)
+                  tmp2*bz(k)-tmp1*az(k), kind=mytype)   ! modify the cw1 in the spectral space Z
+                                                    !  bz  az
+                                                    !  by  ay
+                                                    !  bx  ax
 
              ! POST PROCESSING IN Y
              tmp1 = real(cw1(i,j,k), kind=mytype)
              tmp2 = aimag(cw1(i,j,k))
              cw1(i,j,k) = cmplx(tmp1*by(j)+tmp2*ay(j), &
                   tmp2*by(j)-tmp1*ay(j), kind=mytype)
-             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)
+             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)   ! why should be axisymmetry!
 
              ! POST PROCESSING IN X
              tmp1 = real(cw1(i,j,k), kind=mytype)
@@ -334,7 +356,8 @@ contains
              tmp1=real(kxyz(i,j,k), kind=mytype)
              tmp2=aimag(kxyz(i,j,k))
              ! CANNOT DO A DIVISION BY ZERO
-             if ((tmp1.lt.epsilon).or.(tmp2.lt.epsilon)) then
+             !  Yes ! division  by zero is impossible!!!-----------------------------------<
+             if ((tmp1.lt.epsilon).or.(tmp2.lt.epsilon)) then   !epsilon?  what does it mean?
                 cw1(i,j,k)=0._mytype
 !                print *,'DIV 0',i,j,k,epsilon
              else
@@ -354,6 +377,9 @@ contains
              tmp2 = aimag(cw1(i,j,k))
              cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                   -tmp2*bz(k)-tmp1*az(k), kind=mytype)
+              !                                          bz    az
+              !                                          by    ay
+              !                                          bx    ax
 
              ! POST PROCESSING IN Y
              tmp1 = real(cw1(i,j,k), kind=mytype)
@@ -374,7 +400,7 @@ contains
     end do
              
     ! compute c2r transform
-    call decomp_2d_fft_3d(cw1,rhs)
+    call decomp_2d_fft_3d(cw1,rhs)    !  from complex  to real!
     
  !   call decomp_2d_fft_finalize
 
@@ -401,6 +427,7 @@ contains
     nz = nz_global
 
     ! rhs is in Z-pencil but requires global operations in X
+    !   two steps to change from z to x
     call transpose_z_to_y(rhs,rw2,ph)
     call transpose_y_to_x(rw2,rw1,ph)
     do k=ph%xst(3),ph%xen(3)
