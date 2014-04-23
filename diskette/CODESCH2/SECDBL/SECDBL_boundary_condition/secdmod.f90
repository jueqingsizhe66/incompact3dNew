module secd
      
      implicit none
      integer,parameter :: m1 = 2001
      real,dimension(m1) :: y,ys  !cord  cordinate
      real,dimension(m1) :: sd1,sd2,sd3  !nusol numerical solution
      real,dimension(m1) :: sda  ! ansol :analysit solution
      real,dimension(m1) :: u
      real,dimension(m1) :: dcy,dmy,ddcy  !metri  度量量
      real,dimension(m1) :: dcay,dmay     !metra 额外度量量
      real   :: dx,dxq    ! mesh网格
      integer  :: n1,n1m    !dim尺寸
      integer  :: imor     !icas 迭代继续吗
      real     :: al,alfu  !parco
      real     :: yprer    !ypr 物理表示
      real     :: uup,uwall  !ubou u的边界量

      
end module secd
