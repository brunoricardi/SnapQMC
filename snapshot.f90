

      Module functions
      Implicit none
      Integer, parameter :: kd=8
      Integer, Dimension(4), save :: mm=(/ 502,1521,4071,2107/), ll=(/   0,   0,   0,   1/), &
              &                      irn=(/   0,   0,   0,   1/)  
      Contains

        Real(kd) function rannyu()
        Real(kd), Parameter :: ooto12=1.0_kd/4096.0_kd
        Integer, Parameter :: itwo12=4096
        Integer :: i1, i2, i3, i4
        i1=ll(1)*mm(4)+ll(2)*mm(3)+ll(3)*mm(2)+ll(4)*mm(1)
        i2=ll(2)*mm(4)+ll(3)*mm(3)+ll(4)*mm(2)
        i3=ll(3)*mm(4)+ll(4)*mm(3)
        i4=ll(4)*mm(4)
        ll(4)=mod(i4,itwo12)
        i3=i3+i4/itwo12
        ll(3)=mod(i3,itwo12)
        i2=i2+i3/itwo12
        ll(2)=mod(i2,itwo12)
        ll(1)=mod(i1+i2/itwo12,itwo12)
        rannyu=ooto12*(float(ll(1)) + ooto12*(float(ll(2)) + ooto12*(float(ll(3)) + ooto12*(float(ll(4))))))
        end function rannyu

      subroutine setrn(iseed)
        Integer, Dimension(4) :: iseed
        Integer :: i, j
        Integer :: isn, ipe, ipd, id
        do j = 1, 4
          isn = 0
          do i = 1, 4
            ipe = 4 - i
            ipd = 10 ** ipe
            id = iseed(j) / ipd
            isn = isn + id * 8 ** ipe
            iseed(j) = iseed(j) - id * ipd
          end do
          iseed(j) = isn
        end do
        ll=iseed
        ll(4)=2*(ll(4)/2)+1
      end subroutine setrn

      subroutine savern(iseed)
        Integer, Dimension(4) :: iseed
        Integer :: i, j
        Integer :: ipe, ipo, id, isn
        iseed=ll
        do j = 1,4
          isn = 0
          do i = 1,4
            ipe = 4 - i
            ipo = 8 ** ipe
            id = iseed(j) / ipo
            isn = isn + id * 10 ** ipe
            iseed(j) = iseed(j) - ipo * id
          end do
          iseed(j) = isn
       end do
      end subroutine savern

    FUNCTION FindDet(matrix, n)
    Complex(8) :: FindDet        
    Complex(8), DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    Complex(8) :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == cmplx(0._8,0._8)) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= cmplx(0._8,0._8)) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = cmplx(0._8,0._8)
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
   
END FUNCTION FindDet

end module


program snapshot
use functions
Implicit none
Real(8), parameter:: zero=0._8, one=1._8,half=0.5_8
Complex(8), Dimension(:,:), Allocatable:: psi_up, psi_dwn, psi_uph, psi_dwnh
Integer, Dimension(:), Allocatable:: sup_old,sup_new,sdwn_old,sdwn_new,nconf,sconf
Real(8), Dimension(:), Allocatable:: n_ave, n_err, n_sum, n_ac, n_ac2
Complex(8), Dimension(:,:), Allocatable:: Iup, Idwn, oup, odwn
Complex(8) :: detoup, detodwn
Real(8) :: pup_old, pup_new,pdwn_old,pdwn_new,pratio
Real(8), Dimension(:,:), Allocatable:: psi, sites
Real(8) :: dy, xb, yb, p0, acc, c_stp
Integer:: N_up, N_dwn, N_sites, N_snap, N_metrop, N_eq, N_spb
Integer:: Lx, Ly
Integer:: i,j, iblock, imetrop, ipart, istep, is


call setrn(irn)

Lx=4
Ly=4
N_up=8
N_dwn=8
N_sites=Lx*Ly
N_eq=10
N_snap=50
N_spb=400
N_metrop=5

! set lattice positions

Allocate(sites(1:N_sites,1:2))

dy=half*sqrt(3._8)

do i=1,N_sites-Lx+1,Lx

xb = -float((i-1)/Lx)*half
yb = -float((i-1)/Lx)*dy

do j=0,Lx-1

is = i+j

sites(is,1)=xb+float(j)
sites(is,2)=yb

end do

end do

open (5,file='psi.in')

Allocate(sconf(1:N_sites))
Allocate(nconf(1:N_sites))
Allocate(n_ave(1:N_sites))
Allocate(n_err(1:N_sites))
Allocate(n_sum(1:N_sites))
Allocate(n_ac(1:N_sites))
Allocate(n_ac2(1:N_sites))
Allocate(psi_up(1:N_sites,1:N_up))
Allocate(psi_dwn(1:N_sites,1:N_dwn))
Allocate(psi(1:N_sites,1:N_sites))

! reading eigenstates of the free-particle Hamiltonian
! The eigenstates were pre-computed in octave

do i=1,N_up
do j=1,N_sites

read(5,*) psi(j,i)
end do
end do

! define the spin-up (down) sectors of the wave function

do i=1,N_sites
do j=1,N_up
psi_up(i,j)=cmplx(psi(i,j),zero)
end do
do j=1,N_dwn
psi_dwn(i,j)=cmplx(psi(i,j),zero)
end do
end do

! define the conjugate-transpose of psi

Allocate(psi_uph(1:N_up,1:N_sites))
Allocate(psi_dwnh(1:N_dwn,1:N_sites))

do i=1, N_sites
do j=1, N_up
psi_uph(j,i)=conjg(psi_up(i,j))
end do
do j=1, N_dwn
psi_dwnh(j,i)=conjg(psi_dwn(i,j))
end do
end do

! seting initial configurations

Allocate(sup_old(1:N_up))
Allocate(sup_new(1:N_up))
Allocate(sdwn_old(1:N_dwn))
Allocate(sdwn_new(1:N_dwn))

do i=1,N_up
sup_old(i)=2*i
sdwn_old(i)=2*i-1
end do


! creating the Slater Determinants Iup and Idwn

Allocate(Iup(1:N_sites,1:N_up))
Allocate(Idwn(1:N_sites,1:N_dwn))

Iup = cmplx(zero,zero)
Idwn= cmplx(zero,zero)

do i=1, N_up
Iup(sup_old(i),i)=cmplx(one,zero)
end do

do i=1, N_dwn
Idwn(sdwn_old(i),i)=cmplx(one,zero)
end do

! calculating the overlap with psi

Allocate(oup(1:N_up,1:N_up))
Allocate(odwn(1:N_dwn,1:N_dwn))

oup=matmul(psi_uph,Iup)
odwn=matmul(psi_dwnh,Idwn)

! evaluating the probabilities

!detoup = FindDet(oup,N_up)
!detodwn = FindDet(odwn,N_dwn)
call matinv(oup,N_up,detoup)
call matinv(odwn,N_dwn,detodwn)
pup_old=real(detoup)
pdwn_old=real(detodwn)

open(6,file='snapshot.out')
open(8,file='snapcfg.in')
p0=pup_old*pdwn_old

!..........................
!Equilibration phase

do iblock=1,N_eq

do istep=1, N_spb

do imetrop=1, N_metrop
! sample the overlap probability with Metropolis

do ipart=1,N_up

sup_new=sup_old

sup_new(ipart) = int(N_sites*rannyu())+1

Iup = cmplx(zero,zero)

do i=1, N_up
Iup(sup_new(i),i)=cmplx(one,zero)
end do

! calculating the overlap with psi

oup=matmul(psi_uph,Iup)

! evaluating the probabilities

!detoup = FindDet(oup,N_up)
call matinv(oup,N_up,detoup)
pup_new=real(detoup)

pratio=(pup_new/pup_old)**2
if (pratio .gt. rannyu()) then
sup_old(ipart)=sup_new(ipart)
pup_old=pup_new
end if

end do

do ipart=1,N_dwn

sdwn_new=sdwn_old

sdwn_new(ipart) = int(N_sites*rannyu())+1

Idwn = cmplx(zero,zero)

do i=1, N_dwn
Idwn(sdwn_new(i),i)=cmplx(one,zero)
end do

! calculating the overlap with psi

odwn=matmul(psi_dwnh,Idwn)

! evaluating the probabilities

!detodwn = FindDet(odwn,N_dwn)
call matinv(odwn,N_dwn,detodwn)
pdwn_new=real(detodwn)

pratio=(pdwn_new/pdwn_old)**2
if (pratio .gt. rannyu()) then
sdwn_old(ipart)=sdwn_new(ipart)
pdwn_old=pdwn_new
end if

end do

end do

end do

end do

!...........................

! start to create snapshots

n_ac=zero
n_ac2=zero
c_stp=zero
acc=zero

do iblock=1,N_snap

n_sum=zero

do istep=1,N_spb

do imetrop=1, N_metrop
! sample the overlap probability with Metropolis

do ipart=1,N_up

sup_new=sup_old

sup_new(ipart) = int(N_sites*rannyu())+1

Iup = cmplx(zero,zero)

do i=1, N_up
Iup(sup_new(i),i)=cmplx(one,zero)
end do

! calculating the overlap with psi

oup=matmul(psi_uph,Iup)

! evaluating the probabilities

!detoup = FindDet(oup,N_up)
call matinv(oup,N_up,detoup)

pup_new=real(detoup)

pratio=(pup_new/pup_old)**2
c_stp=c_stp+one
if (pratio .gt. rannyu()) then
acc=acc+one
sup_old(ipart)=sup_new(ipart)
pup_old=pup_new
end if

end do

do ipart=1,N_dwn

sdwn_new=sdwn_old

sdwn_new(ipart) = int(N_sites*rannyu())+1

Idwn = cmplx(zero,zero)

do i=1, N_dwn
Idwn(sdwn_new(i),i)=cmplx(one,zero)
end do

! calculating the overlap with psi

odwn=matmul(psi_dwnh,Idwn)

! evaluating the probabilities

!detodwn = FindDet(odwn,N_dwn)
call matinv(odwn,N_dwn,detodwn)

pdwn_new=real(detodwn)

pratio=(pdwn_new/pdwn_old)**2
c_stp=c_stp+one
if (pratio .gt. rannyu()) then
acc=acc+one
sdwn_old(ipart)=sdwn_new(ipart)
pdwn_old=pdwn_new
end if

end do

end do

nconf=0
sconf=0

do i=1,N_up
nconf(sup_old(i))=nconf(sup_old(i))+1
sconf(sup_old(i))=sconf(sup_old(i))+1
end do
do i=1,N_dwn
nconf(sdwn_old(i))=nconf(sdwn_old(i))+1
sconf(sdwn_old(i))=sconf(sdwn_old(i))-1
end do

n_sum=n_sum+nconf

end do
write(6,*) 'prob', (pup_old*pdwn_old/p0)**2 
do i=1,N_sites
write(6,*) sites(i,1), sites(i,2), nconf(i), sconf(i)
write(8,*) sites(i,1), sites(i,2), nconf(i), sconf(i)
end do

n_sum=n_sum/float(N_spb)
n_ac = n_ac + n_sum
n_ac2 = n_ac2 + n_sum**2

end do

acc=acc/c_stp
write(6,*) 'final results, acc =', acc
do i=1, N_sites
n_ave(i) = n_ac(i)/float(N_snap)
n_err(i) = sqrt(abs(n_ave(i)**2-(n_ac2(i)/float(N_snap)))/float(N_snap-1))
write(6,*) sites(i,1), sites(i,2), n_ave(i), n_err(i)
end do


close(6)


end program

     subroutine matinv(aa,nsub,det)
      Implicit none
      Real(8), parameter:: one=1._8, zero=0._8
      Integer, parameter :: nmax=65536
      Integer, parameter :: nbmx=256
      Integer, parameter :: nbmxs=nbmx*nbmx

!
! routine to calculate inverse and determinant of matrix a
! assumed to be dimensioned a(nsub,nsub).
! the matrix a is replaced by its inverse.
!
      Complex(8) :: aa(nbmxs), atemp(nmax), adiag, adiagi
      Complex(8) :: t, determ, det
      Integer :: ipivot(nmax), idiag
      Integer :: n, nsub, iclm, itemp
      Integer :: i, j, k, jn

      n=nsub
      do 5  i=1,n
  5   ipivot(i)=i
!
! initialize determinant
!
      determ=one
!
! loop through columns
!
      iclm=-n
      do 10 i=1,n
      iclm=iclm+n
!
! loop through rows and select row with largest element
!
      adiag=aa(ipivot(i)+iclm)
      idiag=i
      do 15 k=i+1,n
      if (abs(aa(ipivot(k)+iclm)).gt.abs(adiag)) then
                                          adiag=aa(ipivot(k)+iclm)
                                          idiag=k
                                          endif
   15 continue
!
! interchange pointers if row different from
! original is selected and change sign of determinant because
! of interchange
!
      if (idiag.ne.i) then
                      determ=-determ
                      itemp=ipivot(i)
                      ipivot(i)=ipivot(idiag)
                      ipivot(idiag)=itemp
                      endif
!
! update determinant
!
      determ=adiag*determ
      aa(ipivot(i)+iclm)=cmplx(one,zero)
      adiagi=cmplx(one,zero)/adiag
!
! scale row by inverse diagonal
!
      call sscal(n,adiagi,aa(ipivot(i)),n)
!
! loop through other rows
! if not current row, then row reduce
!
      do 20 j=1,n
      if (j.ne.ipivot(i)) then
                  t=-aa(j+iclm)
                  aa(j+iclm)=cmplx(zero,zero)
                  call thissaxpy(n,t,aa(ipivot(i)),n,aa(j),n)
                  endif
   20 continue
   10 continue
!
! interchange elements to unpivot inverse matrix
! the following is equivalent to:
!      anew(i,ipivot(j))=aold(ipivot(i),j)
!
      jn=-n
      do 30 j=1,n
      jn=jn+n
      do 40 i=1,n
   40 atemp(i)=aa(i+jn)
      do 50 i=1,n
   50 aa(i+jn)=atemp(ipivot(i))
   30 continue
      do 55 j=1,n
   55 ipivot(j)=(ipivot(j)-1)*n
      do 60 i=1,n
      jn=-n
      do 70 j=1,n
      jn=jn+n
   70 atemp(j)=aa(i+jn)
      do 80 j=1,n
   80 aa(i+ipivot(j))=atemp(j)
   60 continue
      det=determ
      return
      end
     subroutine sscal(n,scalor,x,nskip)
! routine to scale x by scalor
! n=number of elements in x to be scaled
! nskip=stride

      Implicit none
      !Integer, parameter :: nbmxs=nbmx*nbmx
      Complex(8) :: x(1), scalor
      Integer :: n, nskip
      Integer :: i, ix

      ix=1-nskip
      do 10 i=1,n
      ix=ix+nskip
   10 x(ix)=scalor*x(ix)
      return
      end

      subroutine thissaxpy(n,aa,x,nxskip,y,nyskip)
!
! routine to calculate y=a*x+y where
! x and y are arrays, and a is a scalar
! n=number of elements in x and y to calculate
! nxskip=x stride
! nyskip=y stride
!
      Implicit none
      !Integer, parameter :: nbmxs=nbmx*nbmx
      Complex(8) :: x(1),y(1),aa
      Integer :: n, nxskip, nyskip
      Integer :: i, ix, iy

      ix=1-nxskip
      iy=1-nyskip
      do 10 i=1,n
      ix=ix+nxskip
      iy=iy+nyskip
   10 y(iy)=aa*x(ix)+y(iy)
      return
      end

