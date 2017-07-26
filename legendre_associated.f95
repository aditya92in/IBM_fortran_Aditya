subroutine legendre_associated ( n, m, x, cx )

!*****************************************************************************80
!
!! LEGENDRE_ASSOCIATED evaluates the associated Legendre functions.
!
!  Differential equation:
!
!    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N,X) )
!
!    P00 =    1
!    P10 =    1 X
!    P20 = (  3 X^2 -     1)/2
!    P30 = (  5 X^3 -   3 X)/2
!    P40 = ( 35 X^4 -  30 X^2 +     3)/8
!    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
!    P60 = (231 X^6 - 315 X^4 + 105 X^2 -    5)/16
!    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
!
!    M = 1
!
!    P01 =   0
!    P11 =   1 * SQRT(1-X*X)
!    P21 =   3 * SQRT(1-X*X) * X
!    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
!    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
!
!    M = 2
!
!    P02 =   0
!    P12 =   0
!    P22 =   3 * (1-X*X)
!    P32 =  15 * (1-X*X) * X
!    P42 = 7.5 * (1-X*X) * (7*X*X-1)
!
!    M = 3
!
!    P03 =   0
!    P13 =   0
!    P23 =   0
!    P33 =  15 * (1-X*X)^1.5
!    P43 = 105 * (1-X*X)^1.5 * X
!
!    M = 4
!
!    P04 =   0
!    P14 =   0
!    P24 =   0
!    P34 =   0
!    P44 = 105 * (1-X*X)^2
!
!  Recursion:
!
!    if N < M:
!      P(N,M) = 0
!    if N = M:
!      P(N,M) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      P(N,M) = X*(2*M+1)*P(M,M)
!    if M+1 < N:
!      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
!
!  Special values:
!
!    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
!    function of the first kind equals the Legendre polynomial of the
!    first kind.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to be
!    evaluated.  X must satisfy -1 <= X <= 1.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 functions.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(0:n)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  real ( kind = 8 ) somx2
  real ( kind = 8 ) x

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop 1
  end if
 
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    stop 1
  end if
 
  if ( x < -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no less than -1.'
    stop 1
  end if
 
  if ( 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no more than 1.'
    stop 1
  end if
  
  cx(0:m-1) = 0.0D+00

  cx(m) = 1.0D+00
  somx2 = sqrt ( 1.0D+00 - x * x )
 
  fact = 1.0D+00
  do i = 1, m
    cx(m) = -cx(m) * fact * somx2
    fact = fact + 2.0D+00
  end do
 
  if ( m + 1 <= n ) then
    cx(m+1) = x * real ( 2 * m + 1, kind = 8 ) * cx(m)
  end if

  do i = m + 2, n
    cx(i) = ( real ( 2 * i     - 1, kind = 8 ) * x * cx(i-1) &
            + real (   - i - m + 1, kind = 8 ) *     cx(i-2) ) &  
            / real (     i - m,     kind = 8 )
  end do

  return
end
