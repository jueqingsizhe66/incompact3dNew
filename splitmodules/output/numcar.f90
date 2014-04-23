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
