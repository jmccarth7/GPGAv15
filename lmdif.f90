subroutine lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
  diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf ) !, iunit ) ! 20131209

!*****************************************************************************80
!
!! LMDIF minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!
!      subroutine fcn ( m, n, x, fvec, iflag )
!      integer ( kind = 4 ) m
!      integer ( kind = 4 ) n
!
!      real fvec(m)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.  N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Input, real ( kind = 8 ) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  Therefore, XTOL
!    measures the relative error desired in the approximate solution.  XTOL
!    should be nonnegative.
!
!    Input, real ( kind = 8 ) GTOL. termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable step length for
!    the forward-difference approximation.  This approximation assumes that
!    the relative errors in the functions are of the order of EPSFCN.
!    If EPSFCN is less than the machine precision, it is assumed that the
!    relative errors in the functions are of the order of the machine
!    precision.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  This bound is
!    set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of iterates if it
!    is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See the description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of squares
!       are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN has reached or exceeded MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of the array FJAC.
!    LDFJAC must be at least M.
!
!    Output, integer ( kind = 4 ) IPVT(N), defines a permutation matrix P such that
!    JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column IPVT(J)
!    of the identity matrix.
!
!    Output, real ( kind = 8 ) QTF(N), the first N elements of Q'*FVEC.
!
  implicit none

  !integer(kind=4),intent(in) ::  iunit   ! jjm

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) actred
  real ( kind = 8 ) delta
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) dirder
  real ( kind = 8 ) enorm
  real ( kind = 8 ) epsfcn
  real ( kind = 8 ) epsmch
  real ( kind = 8 ) factor
  external fcn
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fnorm
  real ( kind = 8 ) fnorm1
  real ( kind = 8 ) ftol
  real ( kind = 8 ) fvec(m)
  real ( kind = 8 ) gnorm
  real ( kind = 8 ) gtol
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) nprint
  real ( kind = 8 ) par
  real ( kind = 8 ) pnorm
  real ( kind = 8 ) prered
  real ( kind = 8 ) qtf(n)
  real ( kind = 8 ) ratio
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) wa1(n)
  real ( kind = 8 ) wa2(n)
  real ( kind = 8 ) wa3(n)
  real ( kind = 8 ) wa4(m)

  real ( kind = 8 ) x(n)

  real ( kind = 8 ) xnorm
  real ( kind = 8 ) xtol

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0

  if ( n <= 0 ) then
    go to 300
  else if ( m < n ) then
    go to 300
  else if ( ldfjac < m ) then
    go to 300
  else if ( ftol < 0.0D+00 ) then
    go to 300
  else if ( xtol < 0.0D+00 ) then
    go to 300
  else if ( gtol < 0.0D+00 ) then
    go to 300
  else if ( maxfev <= 0 ) then
    go to 300
  else if ( factor <= 0.0D+00 ) then
    go to 300
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 300
      end if
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  !if( iunit > 0 ) write(6,'(A)') 'lmdif: call fcn '

  call fcn ( m, n, x, fvec, iflag )

  !if( iunit > 0 ) write(6,'(A,1x,I6)') 'lmdif: aft call fcn iflag = ', iflag

  nfev = 1

  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = enorm ( m, fvec )
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0D+00
  iter = 1
!
!  Beginning of the outer loop.
!
30 continue
!
!  Calculate the jacobian matrix.
!
  iflag = 2

  !if( iunit > 0 ) write(6,'(A)') 'lmdif: call fdjac2 '

  call fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn )

  !if( iunit > 0 ) write(6,'(A,1x,I6)') 'lmdif: aft call fdjac2 iflag = ', iflag

  nfev = nfev + n

  !if( iunit > 0 ) write(6,'(A,1x,I10)') 'lmdif: aft call fdjac2  nfev = ', nfev

  if ( iflag < 0 ) then
    go to 300
  end if
!
!  If requested, call FCN to enable printing of iterates.
!
     if ( 0 < nprint ) then
       iflag = 0
       if ( mod ( iter-1, nprint ) == 0 ) then
         call fcn ( m, n, x, fvec, iflag )
       end if
       if ( iflag < 0 ) then
         go to 300
       end if
     end if
!
!  Compute the QR factorization of the jacobian.
!
     call qrfac ( m, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )
!
!  On the first iteration and if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     if ( iter == 1 ) then

       if ( mode /= 2 ) then
         diag(1:n) = wa2(1:n)
         do j = 1, n
           if ( wa2(j) == 0.0D+00 ) then
             diag(j) = 1.0D+00
           end if
         end do
       end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm ( n, wa3 )
       delta = factor * xnorm
       if ( delta == 0.0D+00 ) then
         delta = factor
       end if
     end if
!
!  Form Q' * FVEC and store the first N components in QTF.
!
     wa4(1:m) = fvec(1:m)

     do j = 1, n

       if ( fjac(j,j) /= 0.0D+00 ) then
         sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
         temp = - sum2 / fjac(j,j)
         wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
       end if

       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)

     end do
!
!  Compute the norm of the scaled gradient.
!
     gnorm = 0.0D+00

     if ( fnorm /= 0.0D+00 ) then

       do j = 1, n

         l = ipvt(j)

         if ( wa2(l) /= 0.0D+00 ) then
           sum2 = 0.0D+00
           do i = 1, j
             sum2 = sum2 + fjac(i,j) * ( qtf(i) / fnorm )
           end do
           gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
         end if

       end do

     end if
!
!  Test for convergence of the gradient norm.
!
     if ( gnorm <= gtol ) then
       info = 4
       go to 300
     end if
!
!  Rescale if necessary.
!
     if ( mode /= 2 ) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  Beginning of the inner loop.
!
200  continue
!
!  Determine the Levenberg-Marquardt parameter.
!
        call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min ( delta, pnorm )
        end if
!
!  Evaluate the function at X + P and calculate its norm.
!
        iflag = 1
        call fcn ( m, n, wa2, wa4, iflag )

        !if( iunit > 0 ) write(6,'(A,1x,I6)') 'lmdif:2 aft call fcn    iflag = ', iflag

        nfev = nfev + 1

        !if( iunit > 0 ) write(6,'(A,1x,I10)') 'lmdif:2 aft call fcn     nfev = ', nfev

        if ( iflag < 0 ) then
          go to 300
        end if
        fnorm1 = enorm ( m, wa4 )
!
!  Compute the scaled actual reduction.
!
        if ( 0.1D+00 * fnorm1 < fnorm ) then
          actred = 1.0D+00 - ( fnorm1 / fnorm )**2
        else
          actred = -1.0D+00
        end if
!
!  Compute the scaled predicted reduction and the scaled directional derivative.
!
        do j = 1, n
          wa3(j) = 0.0D+00
          l = ipvt(j)
          temp = wa1(l)
          wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        end do

        temp1 = enorm ( n, wa3 ) / fnorm
        temp2 = ( sqrt ( par ) * pnorm ) / fnorm
        prered = temp1**2 + temp2**2 / 0.5D+00
        dirder = - ( temp1**2 + temp2**2 )
!
!  Compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0D+00
        if ( prered /= 0.0D+00 ) then
          ratio = actred / prered
        end if
!
!  Update the step bound.
!
        if ( ratio <= 0.25D+00 ) then

           if ( actred >= 0.0D+00 ) then
             temp = 0.5D+00
           endif

           if ( actred < 0.0D+00 ) then
             temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
           end if

           if ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) then
             temp = 0.1D+00
           end if

           delta = temp * min ( delta, pnorm / 0.1D+00  )
           par = par / temp

        else

           if ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) then
             delta = 2.0D+00 * pnorm
             par = 0.5D+00 * par
           end if

        end if
!
!  Test for successful iteration.
!

!
!  Successful iteration. update X, FVEC, and their norms.
!
        if ( 0.0001D+00 <= ratio ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:m) = wa4(1:m)
          xnorm = enorm ( n, wa2 )
          fnorm = fnorm1
          iter = iter + 1
        end if
!
!  Tests for convergence.
!
        if ( abs ( actred) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) then
          info = 1
        end if

        if ( delta <= xtol * xnorm ) then
          info = 2
        end if

        if ( abs ( actred) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) info = 3
        if ( info /= 0 ) go to 300
!
!  Tests for termination and stringent tolerances.
!
        !if( iunit > 0 ) write(6,'(A,2(1x,I10))') 'lmdif:3 nfev, maxfev ', nfev , maxfev

        if ( nfev >= maxfev ) then
          info = 5
        end if

        if ( abs ( actred) <= epsmch .and. prered <= epsmch &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) info = 6
        if ( delta <= epsmch * xnorm ) info = 7
        if ( gnorm <= epsmch ) info = 8

        if ( info /= 0 ) then
          go to 300
        end if
!
!  End of the inner loop.  Repeat if iteration unsuccessful.
!
        if ( ratio < 0.0001D+00 ) then
          go to 200
        end if
!
!  End of the outer loop.
!
     go to 30

300 continue

!if( iunit > 0 ) write(6,'(A,2(1x,I10))') 'lmdif:4 aft 300 iflag', iflag
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( nprint > 0 ) then
    call fcn ( m, n, x, fvec, iflag )
  end if

  return
end
