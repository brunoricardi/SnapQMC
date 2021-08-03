program cluster
Implicit none
Integer, parameter:: N_snaps=5000,N_sites=58,Lx=7,Ly=7
Integer:: i,j,k
Real(8), Dimension(N_snaps,N_sites):: x,y
Integer, Dimension(N_snaps,N_sites):: n
Real(8), Dimension(N_sites):: chist
Integer :: csize

open (80,file='snapcfg_hex.in')

do i=1,N_snaps
    do j=1,N_sites
        read(80,*) x(i,j), y(i,j), n(i,j)
    end do
end do

chist=0._8
do i = 1, N_snaps
do j = 1, N_sites
if (n(i,j).eq.0) then
        csize=0
        n(i,j)=-1
        csize=csize+1
        call count_hole(j,n(i,:),x(i,:),y(i,:),csize)
        chist(csize)=chist(csize)+1._8
        if (i.eq.1) write(*,*) j, csize
end if
end do
end do

chist=chist/float(N_snaps)


open (13,file='chist.out')

do i=1,N_sites
write(13,*) i, chist(i)
end do

end program

recursive subroutine count_hole(il,nl,xl,yl,c)
Implicit none

Integer, parameter :: N_sites = 58
Integer, Dimension(N_sites) :: nl
Real(8), Dimension(N_sites) :: xl, yl
Real(8) :: r2
Integer:: il,c,is,js


do js =1, N_sites
    r2=(xl(il)-xl(js))**2+(yl(il)-yl(js))**2
    if (r2.gt.0.99_8.AND.r2.lt.1.01_8) then
        if (nl(js).eq.0._8) then
            c=c+1
            nl(js)=-1
            is=js
    end if
    end if
end do


end

