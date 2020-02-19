module gausJordaMod
contains

  ! ========================================================== !
  ! === Gauss Jordan Matrix Inversion                      === !
  ! ========================================================== !
  subroutine gaussJordan( Amat, Bmat, nSize )
    implicit none
    integer         , intent(in)  :: nSize
    double precision, intent(in)  :: Amat(nSize,nSize)
    double precision, intent(out) :: Bmat(nSize,nSize)
    integer                       :: i, j, k
    double precision              :: Dinv
    double precision, parameter   :: eps = 1.d-10
    double precision, allocatable :: Umat(:,:)

    ! ----------------------------------------- !
    ! --- [1] Preparation                   --- !
    ! ----------------------------------------- !
    !  -- [1-1] allocate Umat               --  !
    allocate( Umat(nSize,nSize) )
    !  -- [1-2] Copy AMat & bvec            --  !
    do j=1, nSize
       do i=1, nSize
          Umat(i,j) = Amat(i,j)
       enddo
    enddo
    
    ! ----------------------------------------- !
    ! --- [2] Upper Triangular Matrization  --- !
    ! ----------------------------------------- !
    do k=1, nSize
       
       !  -- [2-1] Diagonal Component       --  !
       if ( abs( Umat(k,k) ).lt.eps ) then
          write(6,*) '[gaussElimin] Amat != Regular Matrix :: No Solution End :: @ k= ', k
          stop
       endif
       Dinv      = 1.d0 / Umat(k,k)
       Umat(k,k) = 1.d0

       !  -- [2-2] Non-Diagonal Component   --  !
       if ( k.eq.nSize ) then
          ! -- [    Last Row :: k == nSize ] -- !
          Bmat(k,:) = Dinv * Bmat(k,:)
       else
          ! -- [Not Last Row :: k != nSize ] -- !
          !  - Division    -  !
          Umat(k,k+1:nSize) = Dinv * Umat(k,k+1:nSize)
          Bmat(k,:)         = Dinv * Bmat(k,:)
          !  - subtraction -  !
          do j=k+1,nSize
             Umat(j,k+1:nSize) = Umat(j,k+1:nSize) - Umat(j,k) * Umat(k,k+1:nSize)
             Bmat(j,:)         = Bmat(j,:)         - Umat(j,k) * Bmat(k,:)
             Umat(j,k)         = 0.d0
          enddo
       endif
       
    end do

    ! ----------------------------------------- !
    ! --- [3] Backward Substituition        --- !
    ! ----------------------------------------- !
    do k=nSize-1, 1, -1
       do i=nSize, k+1, -1
          Bmat(k,:) = Bmat(k,:) - Umat(k,i)*Bmat(i,:)
       enddo
    enddo

    return
  end subroutine gaussJordan

end module gausJordaMod
