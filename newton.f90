module newton

    ! module parameters:
    implicit none
    integer, parameter :: maxiter = 1000
    real(kind=8), parameter :: tol = 1.d-10

contains

subroutine solve(f, fp, Ta0, Ta, iters, debug)

    ! Estimate the zero of f(Ta) using Newton's method. 
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   Ta0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate Ta satisfying f(Ta)=0 (assumes Newton converged!) 
    !   the number of iterations iters
     
    implicit none
    real(kind=8), intent(in) :: Ta0
    real(kind=8), external :: f, fp
    logical, intent(in) :: debug
    real(kind=8), intent(out) :: Ta
    integer, intent(out) :: iters

    ! Declare any local variables:
    real(kind=8) :: deltaTa, fTa, fTaprime
    integer :: k


    ! initial guess
    Ta = Ta0

    if (debug) then
        print 1, Ta
 1     format('Initial guess: Ta = ', e22.15)
        endif

    ! Newton iteration to find a zero of f(Ta) 

    do k=1,maxiter

        ! evaluate function and its derivative:
        fTa = f(Ta)
        fTaprime = fp(Ta)

        if (abs(fTa) < tol) then
            exit  ! jump out of do loop
            endif

        ! compute Newton increment Ta:
        deltaTa = fTa/fTaprime

        ! update Ta:
        Ta = Ta - deltaTa

        if (debug) then
            print 2, k,Ta
2          format('After', i3, ' iterations, Ta = ', e22.15)
            endif

        enddo

    if (k > maxiter) then
        ! might not have converged

        fTa = f(Ta)
        if (abs(fTa) > tol) then
            print *, '*** Warning: has not yet converged'
            endif
        endif 

    ! number of iterations taken:
    iters = k-1


end subroutine solve

end module newton