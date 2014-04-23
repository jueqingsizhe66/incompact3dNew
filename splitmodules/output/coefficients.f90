subroutine coefficients
!
!*********************************************************************

USE parfiX
USE parfiY
USE parfiZ

implicit none

fialx=0.45
fibex=(3. - 2.*fialx)/10.
fiaix=(2. + 3.*fialx)/4.
fibix=(6. + 7.*fialx)/8.
ficix=(6. + fialx)/20.
fidix=(2. - 3.*fialx)/40.
fibix=fibix*0.5
ficix=ficix*0.5
fidix=fidix*0.5
fialy=0.45
fibey=(3. - 2.*fialy)/10.
fiaiy=(2. + 3.*fialy)/4.
fibiy=(6. + 7.*fialy)/8.
ficiy=(6. + fialy)/20.
fidiy=(2. - 3.*fialy)/40.
fibiy=fibiy*0.5
ficiy=ficiy*0.5
fidiy=fidiy*0.5
fialz=0.45
fibez=(3. - 2.*fialz)/10.
fiaiz=(2. + 3.*fialz)/4.
fibiz=(6. + 7.*fialz)/8.
ficiz=(6. + fialz)/20.
fidiz=(2. - 3.*fialz)/40.
fibiz=fibiz*0.5
ficiz=ficiz*0.5
fidiz=fidiz*0.5
fia1x=15./16.
fib1x=1./4.
fic1x=-3./8.
fid1x=1./4.
fie1x=-1./16.

fia1y=15./16.
fib1y=1./4.
fic1y=-3./8.
fid1y=1./4.
fie1y=-1./16.

fia1z=15./16.
fib1z=1./4.
fic1z=-3./8.
fid1z=1./4.
fie1z=-1./16.

fia2x=3./4.
fib2x=1./16.
fic2x=3./8.
fid2x=-1./4.
fie2x=1./16.

fia2y=3./4.
fib2y=1./16.
fic2y=3./8.
fid2y=-1./4.
fie2y=1./16.

fia2z=3./4.
fib2z=1./16.
fic2z=3./8.
fid2z=-1./4.
fie2z=1./16.

fia3x=5./8.
fib3x=-1./16.
fic3x=1./4.
fid3x=4./16.
fie3x=-1./16.

fia3y=5./8.
fib3y=-1./16.
fic3y=1./4.
fid3y=4./16.
fie3y=-1./16.

fia3z=5./8.
fib3z=-1./16.
fic3z=1./4.
fid3z=4./16.
fie3z=-1./16.

fianx=fia1x
fibnx=fib1x
ficnx=fic1x
fidnx=fid1x
fienx=fie1x

fiany=fia1y
fibny=fib1y
ficny=fic1y
fidny=fid1y
fieny=fie1y

fianz=fia1z
fibnz=fib1z
ficnz=fic1z
fidnz=fid1z
fienz=fie1z

fiamx=fia2x
fibmx=fib2x
ficmx=fic2x
fidmx=fid2x
fiemx=fie2x

fiamy=fia2y
fibmy=fib2y
ficmy=fic2y
fidmy=fid2y
fiemy=fie2y

fiamz=fia2z
fibmz=fib2z
ficmz=fic2z
fidmz=fid2z
fiemz=fie2z

fiapx=fia3x
fibpx=fib3x
ficpx=fic3x
fidpx=fid3x
fiepx=fie3x

fiapy=fia3y
fibpy=fib3y
ficpy=fic3y
fidpy=fid3y
fiepy=fie3y

fiapz=fia3z
fibpz=fib3z
ficpz=fic3z
fidpz=fid3z
fiepz=fie3z

return
end subroutine coefficients
