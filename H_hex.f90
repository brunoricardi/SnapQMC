!!!!
!! File: H_hex.f90
!! Description: Create and diagonalize Hamiltonian for Fermi-Hubbard model on a triangular lattice
!! Author: Bruno R. de Abreu  |  babreu at illinois dot edu
!! National Center for Supercomputing Applications (NCSA)
!!  
!! Creation Date: Friday, 6th August 2021, 11:25:48 am
!! Last Modified: Friday, 6th August 2021, 11:27:01 am
!!  
!! Copyright (c) 2021, Bruno R. de Abreu, National Center for Supercomputing Applications.
!! All rights reserved.
!! License: This program and the accompanying materials are made available to any individual
!!          under the citation condition that follows: On the event that the software is
!!          used to generate data that is used implicitly or explicitly for research
!!          purposes, proper acknowledgment must be provided in the citations section of
!!          publications. This includes both the author's name and the National Center
!!          for Supercomputing Applications. If you are uncertain about how to do
!!          so, please check this page: https://github.com/babreu-ncsa/cite-me.
!!          This software cannot be used for commercial purposes in any way whatsoever.
!!          Omitting this license when redistributing the code is strongly disencouraged.
!!          The software is provided without warranty of any kind. In no event shall the
!!          author or copyright holders be liable for any kind of claim in connection to
!!          the software and its usage.
!!!!
program H_hex
    implicit none
    integer :: Lx, Ly  !! sides of the lattice
    integer :: N_sites, N_up, N_down   !! number of sites, spins up, spins down
    integer :: i, j, k, l, lmax, my   !! integer helpers
    real*8 :: dx, dy    !! lattice parameters
    real*8 :: x, y, r2  !! double helpers
    real*8 :: global_dens, local_dens   !! densities
    complex*16 :: diag       !! diagonal elements
    real*8, allocatable :: sites(:,:)   !! matrix with lattice sites positions
    real*8, allocatable :: eigvals(:) !! eigenvalues of the Hamiltonian
    complex*16, allocatable :: H_0(:,:), eigvecs(:,:) !! hamiltonian, its eigenvectors 
    
    Lx=10
    Ly=9
    N_sites=106
    N_up=53
    N_down=53
    my = (Ly+1) / 2
    dx=0.5
    dy = sqrt(3.0)/2.0

    allocate(sites(0:N_sites-1,2))
    sites = 0.0
    allocate(H_0(N_sites,N_sites))
    H_0 = (0.0, 0.0)
    allocate(eigvecs(N_sites,N_sites))
!    eigvecs = (0.0, 0.0)
    allocate(eigvals(N_sites))
!    eigvals = 0.0


!! start giving lattice sites their positions (this needs to be re-thought, lattice only needs to be created once)
    k=1
    do i = 0, my-1
        x = -i*dx
        y = -i*dy
        do j = 1, Lx+i
            sites(k,1) = x + j - 1;
            sites(k,2) = y
            k = k + 1
        enddo
    enddo
    
    lmax = Lx + my + 1
    l=1
    do i = my, Ly-1
        x = x + dx
        y = -i*dy
        do j = 1, lmax - l
            sites(k,1) = x + j - 1
            sites(k,2) = y
            k = k + 1
        enddo
    enddo


    !! construct Hamiltonian
    do i = 1, N_sites
        do j = 1, N_sites
            r2 = (sites(i,1) - sites(j,1))**2.0 + (sites(i,2) - sites(j,2)**2.0)
            if (r2 <= 1.5) then
                if(r2 >= 0.5) then
                    H_0(i,j) = (-1.0,0.0)
                endif
            endif
        enddo
    enddo

    ! diagonalize hamiltonian
    call hermitean_diagonalization(H_0,N_sites,eigvals,eigvecs)

    ! write eigen vector to file
    call write_Umatrix_to_file(N_sites,eigvecs)
    
    ! calculate diagonal densities ?? IMPLEMENT

    deallocate(sites)
    deallocate(H_0)
    deallocate(eigvecs)
    deallocate(eigvals)
end program

subroutine hermitean_diagonalization(matrix,order,eigvals,eigvecs)
    !! This uses LAPACK's ZHEEV to diagonalize HERMITEAN matrices
    IMPLICIT NONE
    COMPLEX*16, intent(in) :: matrix(order,order)
    INTEGER, intent(in) :: order
    REAL*8, intent(out) :: eigvals(order)
    COMPLEX*16, intent(out) :: eigvecs(order,order)
    INTEGER LWMAX
    PARAMETER (LWMAX=1000)
    INTEGER INFO, LWORK
    REAL*8 RWORK(3*order-2)
    COMPLEX*16, allocatable :: WORK(:)

    allocate(WORK(LWMAX))
    eigvecs = matrix
    LWORK=-1
    call ZHEEV('V', 'U', order, eigvecs, order, eigvals, WORK, LWORK, RWORK, INFO)
    LWORK = min(LWMAX, int(WORK(1)))
    deallocate(WORK)
    allocate(WORK(LWORK))
    call ZHEEV('V', 'U', order, eigvecs, order, eigvals, WORK, LWORK, RWORK, INFO)
    if(INFO == 0) then
        write(*,*) 'Diagonalization performed.'
    else
        write(*,*) 'Diagonalization failed.'
        return
    endif

end subroutine hermitean_diagonalization

subroutine write_Umatrix_to_file(order,matrix)
    !! This writes the UPPER TRIANGULAR part of a matrix to a file
    implicit none
    integer, intent(in) :: order
    complex*16, intent(in) :: matrix(order,order)
    integer :: i, j     ! helpers
    character*20 :: filename

    filename = 'psif90_hex.in'
    open(1, file=filename, status='new')

    do j=1,order
        do i=j,order
            write(1,*) matrix(i,j)
        enddo
    enddo

    close(1)
end subroutine write_Umatrix_to_file

