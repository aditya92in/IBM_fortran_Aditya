!This file is part of Notus 0.2.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 24-07-2015, antoine.lemoine@bordeaux-inp.fr

!This software is a computer program whose purpose is to simulate fluid flows.

!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software.  You can  use,
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info".

!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software's author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability.

!In this respect, the user's attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software's suitability as regards their
!requirements in conditions enabling the security of their systems and/or
!data to be ensured and,  more generally, to use and operate it in the
!same conditions as regards security.

!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.



module mod_cg2_line_polygon
use mod_cg2_polygon
use mod_cg2_line_segment
implicit none

contains

!-------------------------
!   Detection functions
!-------------------------

logical pure function cg2_is_line_intersect_polygon(polygon, l0, l1) result(intersect)
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(in) :: l0, l1

integer :: minloc_dist
double precision :: minval_dist

intersect = .false.

call cg2_minloc_line_point_signed_distance(polygon, l0, l1, polygon%point(:,1), minval_dist, minloc_dist)

if (minval_dist <= 0.0d0) intersect = .true.
end function cg2_is_line_intersect_polygon

!---------------------------
!   Computation functions
!---------------------------

pure subroutine cg2_line_polygon_clipping(polygon, l0, l1, clipped_polygon)
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(in) :: l0, l1
type(polygon_2d), intent(inout) :: clipped_polygon

integer :: i
double precision, dimension(2) :: s0, s1, intersection_point
logical :: is_intersection

clipped_polygon%nb_points = 0

s0 = polygon%point(:,polygon%nb_points)
do i = 1, polygon%nb_points
s1 = polygon%point(:,i)

if (cg2_is_point_left_of_or_on_line(l0, l1, s0)) then
! s0 belongs to the clipping region
if (cg2_is_point_left_of_or_on_line(l0, l1, s1)) then
! s1 belongs to the clipping region
clipped_polygon%nb_points = clipped_polygon%nb_points + 1
clipped_polygon%point(:,clipped_polygon%nb_points) = s1
else
! s1 is outside the clipping region
call cg2_line_segment_inclusive_intersection(l0, l1, s0, s1, is_intersection, intersection_point)

! if no intersection if found, suppose that s1 belongs to the line
if (.not. is_intersection) intersection_point = s1

clipped_polygon%nb_points = clipped_polygon%nb_points + 1
clipped_polygon%point(:,clipped_polygon%nb_points) = intersection_point
end if
else
! s0 is outside the clipping region
if (cg2_is_point_left_of_or_on_line(l0, l1, s1)) then
! s1 belongs to the clipping region
call cg2_line_segment_inclusive_intersection(l0, l1, s0, s1, is_intersection, intersection_point)

! if no intersection if found, suppose that s0 belongs to the line
if (.not. is_intersection) intersection_point = s0

clipped_polygon%nb_points = clipped_polygon%nb_points + 1
clipped_polygon%point(:,clipped_polygon%nb_points) = intersection_point

if (norm2(intersection_point - s1) > 0.0d0) then ! Avoid zero-length edge
clipped_polygon%nb_points = clipped_polygon%nb_points + 1
clipped_polygon%point(:,clipped_polygon%nb_points) = s1
end if
end if
end if

s0 = s1
end do
end subroutine cg2_line_polygon_clipping

pure subroutine cg2_line_polygon_boundary_intersection(p0, p1, polygon, s0, s1, i, intersection_found)
double precision, dimension(2), intent(in) :: p0, p1
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(out) :: s0, s1
integer, intent(out) :: i
logical, intent(out) :: intersection_found

intersection_found = .false.
s0 = 0.0d0; s1 = 0.0d0;

do i = 1, polygon%nb_points - 1
call cg2_segments_overlap(p0, p1, polygon%point(:,i+1), polygon%point(:,i), s0, s1, intersection_found)

if (intersection_found) return
end do

! Last edge of the polygon
i = polygon%nb_points

call cg2_segments_overlap(p0, p1, polygon%point(:,1), polygon%point(:,i), s0, s1, intersection_found)
end subroutine cg2_line_polygon_boundary_intersection

pure subroutine cg2_split_polygon(polygon, p0, p1, i0, i1, polygon_left, polygon_right)
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(in) :: p0, p1
integer, intent(in) :: i0, i1
type(polygon_2d), intent(out) :: polygon_left, polygon_right

integer :: nb_left, nb_right, i0p1, i1p1, i

i0p1 = modulo(i0, polygon%nb_points) + 1
i1p1 = modulo(i1, polygon%nb_points) + 1

! Construct the polygon on the left side of the line
if (i0p1 <= i1) then
nb_left = i1 - i0p1 + 3
call cg2_initialize_polygon(polygon_left, nb_left)

polygon_left%point(:,1) = p0

do i = 0, nb_left - 3
polygon_left%point(:,i+2) = polygon%point(:,i+i0p1)
end do

polygon_left%point(:, nb_left) = p1
else
nb_left = polygon%nb_points - i0p1 + i1 + 3
call cg2_initialize_polygon(polygon_left, nb_left)

polygon_left%point(:,1) = p0

do i = i0p1, polygon%nb_points
polygon_left%point(:,i-i0p1+2) = polygon%point(:,i)
end do

do i = 1, i1
polygon_left%point(:,polygon%nb_points-i0p1+2+i) = polygon%point(:,i)
end do

polygon_left%point(:, nb_left) = p1
end if

! Construct the polygon on the right side of the line
if (i1p1 <= i0) then
nb_right = i0 - i1p1 + 3
call cg2_initialize_polygon(polygon_right, nb_right)

polygon_right%point(:,1) = p1

do i = 0, nb_right - 3
polygon_right%point(:,i+2) = polygon%point(:,i+i1p1)
end do

polygon_right%point(:, nb_right) = p0
else
nb_right = polygon%nb_points - i1p1 + i0 + 3
call cg2_initialize_polygon(polygon_right, nb_right)

polygon_right%point(:,1) = p1

do i = i1p1, polygon%nb_points
polygon_right%point(:,i-i1p1+2) = polygon%point(:,i)
end do

do i = 1, i0
polygon_right%point(:,polygon%nb_points-i1p1+2+i) = polygon%point(:,i)
end do

polygon_right%point(:, nb_right) = p0
end if
end subroutine cg2_split_polygon

pure subroutine cg2_split_polygon_with_line(polygon, l0, l1, is_intersection, polygon_left, polygon_right)
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(in) :: l0, l1
logical, intent(out) :: is_intersection
type(polygon_2d), intent(out) :: polygon_left, polygon_right

! Numerical experiments show that optimized algorithm is faster than brute force algorithm if the number of points exceeds 100.
if (polygon%nb_points <= 100) then
call cg2_brute_force_split_polygon_with_line(polygon, l0, l1, is_intersection, polygon_left, polygon_right)
else
call cg2_optimized_split_polygon_with_line(polygon, l0, l1, is_intersection, polygon_left, polygon_right)
end if
end subroutine cg2_split_polygon_with_line

pure subroutine cg2_brute_force_split_polygon_with_line(polygon, l0, l1, is_intersection, polygon_left, polygon_right)
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(in) :: l0, l1
logical, intent(out) :: is_intersection
type(polygon_2d), intent(out) :: polygon_left, polygon_right

integer :: i, i0, i1, i1_start
double precision, dimension(2) :: left_point, right_point

i0 = huge(1)
i1 = huge(1)

! Look for first intersection
call cg2_line_segment_intersection(l0, l1, polygon%point(:,polygon%nb_points), polygon &
%point(:,1), is_intersection, left_point)

if (is_intersection) then
i0 = polygon%nb_points
i1_start = 2
else
do i = 2, polygon%nb_points
call cg2_line_segment_intersection(l0, l1, polygon%point(:,i-1), polygon &
%point(:,i), is_intersection, left_point)

if (is_intersection) then
i0 = i-1
exit
end if
end do

i1_start = i0+2
end if

if (.not. is_intersection) return

! Look for second intersection
do i = i1_start, polygon%nb_points
call cg2_line_segment_intersection(l0, l1, polygon%point(:,i-1), polygon &
%point(:,i), is_intersection, right_point)

if (is_intersection) then
i1 = i-1
exit
end if
end do

if (.not. is_intersection) return

! Split the polygon
call cg2_split_polygon(polygon, left_point, right_point, i0, i1, polygon_left, polygon_right)
end subroutine cg2_brute_force_split_polygon_with_line

pure subroutine cg2_optimized_split_polygon_with_line(polygon, l0, l1, is_intersection, polygon_left, polygon_right)
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(in) :: l0, l1
logical, intent(out) :: is_intersection
type(polygon_2d), intent(out) :: polygon_left, polygon_right

integer :: minloc_dist, maxloc_dist, tmp
double precision :: minval_dist, maxval_dist

integer :: inf, sup, half, left, right, left_next, right_next
double precision :: dist, dist_inf_left, dist_inf_right
double precision, dimension(2) :: left_point, right_point

call cg2_minloc_line_point_signed_distance(polygon, l0, l1, polygon%point(:,1), minval_dist,  minloc_dist)
call cg2_maxloc_line_point_signed_distance(polygon, l0, l1, polygon%point(:,1), maxval_dist,  maxloc_dist)

! Intersection? minval_dist == 0 means intersection with a point
if (minval_dist < -epsilon(1.0d0)) then
is_intersection = .true.
else
is_intersection = .false.
return
end if

if (maxval_dist > epsilon(1.0d0)) then
is_intersection = .true.
else
is_intersection = .false.
return
end if

! Assume minloc_dist < maxloc_dist
if (minloc_dist > maxloc_dist) then
tmp = minloc_dist
minloc_dist = maxloc_dist
maxloc_dist = tmp
end if

! Binary search for the first point
inf = minloc_dist
sup = maxloc_dist
half = (sup + inf)/2
dist = cg2_line_point_signed_distance(l0, l1, polygon%point(:,half), polygon%point(:,1))
dist_inf_left = cg2_line_point_signed_distance(l0, l1, polygon%point(:,inf), polygon%point(:,1))

do while (sup - inf > 1)
if (dist_inf_left == 0.0d0) exit

if (dist_inf_left*dist < 0.0d0) then
sup = half
else
inf = half
end if

half = (sup + inf)/2
dist = cg2_line_point_signed_distance(l0, l1, polygon%point(:,half), polygon%point(:,1)) 
dist_inf_left = cg2_line_point_signed_distance(l0, l1, polygon%point(:,inf), &
polygon%point(:,1))
end do

left = inf

left_next = modulo(left, polygon%nb_points) + 1 
if (dist_inf_left*cg2_line_point_signed_distance(l0, l1, polygon%point(:,left_next), &
polygon%point(:,1)) > 0.0d0) then 
left_next = modulo(left - 2, polygon%nb_points) + 1
end if

! Binary search for the second point
inf = maxloc_dist
sup = polygon%nb_points + minloc_dist

half = (sup + inf)/2
dist = cg2_line_point_signed_distance(l0, l1, polygon%point(:, &
modulo(half-1,polygon%nb_points)+1), polygon%point(:,1))
dist_inf_right = cg2_line_point_signed_distance(                                 &
              l0, l1, polygon%point(:, modulo(inf-1,polygon%nb_points)+1), &
              polygon%point(:,1)                                           &
           )

do while (sup - inf > 1)
if (dist_inf_right == 0.0d0) exit

if (dist_inf_right*dist < 0.0d0) then
sup = half
else
inf = half
end if

half = (sup + inf)/2
dist = cg2_line_point_signed_distance(l0, l1, polygon%point(:, &
modulo(half-1,polygon%nb_points)+1), polygon%point(:,1))
dist_inf_right = cg2_line_point_signed_distance(                                 &
              l0, l1, polygon%point(:, modulo(inf-1,polygon%nb_points)+1), &
              polygon%point(:,1)                                           &
           )
end do

right = modulo(inf-1, polygon%nb_points)+1

right_next = modulo(right, polygon%nb_points) + 1
if (dist_inf_right*cg2_line_point_signed_distance(l0, l1, polygon%point(:,right_next),&
polygon%point(:,1)) > 0.0d0) then
right_next = modulo(right - 2, polygon%nb_points) + 1
end if

! Compute intersection points
call cg2_line_segment_intersection(l0, l1, polygon%point(:,left), &
polygon%point(:,left_next), is_intersection, left_point)

if (.not. is_intersection) return

call cg2_line_segment_intersection(l0, l1, polygon%point(:,right), &
polygon%point(:,right_next), is_intersection, right_point)

if (.not. is_intersection) return

call cg2_split_polygon(polygon, left_point, right_point, left, right,  polygon_left,polygon_right)
end subroutine cg2_optimized_split_polygon_with_line

pure subroutine cg2_minloc_line_point_signed_distance(polygon, p0, p1, pref,&
minimum_distance, minloc_distance)
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(in) :: p0, p1, pref
double precision, intent(out) :: minimum_distance
integer, intent(out) :: minloc_distance

integer, dimension(47), parameter :: FIB = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, &
233, 377, 610, 987, &
1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, &
317811, 514229, 832040,  &
1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, &
 102334155,     &
165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, huge(minloc_distance)]

integer :: n, l, r, lk, rk, a, b, i
double precision :: gl, gr, gl0, gr0, tmp

n = 1
tmp = 1

! A0
do while (fib(n) < polygon%nb_points)
n = n + 1
end do

! A1
if (n == 1) then
minloc_distance = 1
minimum_distance = cg2_line_point_signed_distance(p0, p1, polygon%point(:,1), pref)
return
end if

! A2
l = 1
r = l + fib(n-1)

! A3
gl = cg2_line_point_signed_distance(p0, p1, polygon%point(:,l), pref)
gr = cg2_line_point_signed_distance(p0, p1, polygon%point(:,r), pref)
gl0 = gl
gr0 = gr

! A4
if (gl0 >= gr0) then
a = l
b = l + fib(n)
l = a + fib(n-2)
gl = cg2_line_point_signed_distance(p0, p1, polygon%point(:,&
modulo(l - 1, polygon%nb_points) + 1), pref)
else
b = r
a = r - fib(n)
r = b - fib(n-2)
gr = cg2_line_point_signed_distance(p0, p1, polygon%point(:,&
modulo(r - 1, polygon%nb_points) + 1), pref)
end if

! A5
i = 2
do while (i > 1)
if (gl >= gr) then
tmp = gr
lk = l
l = r
r = 2*lk - a
a = lk
gr = cg2_line_point_signed_distance(p0, p1, &
polygon%point(:,modulo(r - 1, polygon%nb_points) + 1), pref)
i = b - lk
gl = tmp
else
tmp = gl
rk = r
r = l
l = 2*rk - b
b = rk
gl = cg2_line_point_signed_distance(p0, p1,&
polygon%point(:,modulo(l - 1, polygon%nb_points) + 1), pref)
i = rk - a
gr = tmp
end if
end do

if (gl <= gr) then
minloc_distance = modulo(l - 1, polygon%nb_points) + 1
else
minloc_distance = modulo(r - 1, polygon%nb_points) + 1
end if

minimum_distance = tmp
end subroutine cg2_minloc_line_point_signed_distance

pure subroutine cg2_maxloc_line_point_signed_distance(polygon, p0, p1, pref,&
maximum_distance, maxloc_distance)
type(polygon_2d), intent(in) :: polygon
double precision, dimension(2), intent(in) :: p0, p1, pref
double precision, intent(out) :: maximum_distance
integer, intent(out) :: maxloc_distance

integer, dimension(47), parameter :: FIB = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, & 
144, 233, 377, 610, 987,& 
1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, &
196418, 317811, 514229, 832040,& 
1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986,& 
102334155, &
165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, huge(maxloc_distance)]

integer :: n, l, r, lk, rk, a, b, i
double precision :: gl, gr, gl0, gr0, tmp

n = 1
tmp = 1

! A0
do while (fib(n) < polygon%nb_points)
n = n + 1
end do

! A1
if (n == 1) then
maxloc_distance = 1
maximum_distance = cg2_line_point_signed_distance(p0, p1, polygon%point(:,1), pref)
return
end if

! A2
l = 1
r = l + fib(n-1)

! A3
gl = cg2_line_point_signed_distance(p0, p1, polygon%point(:,l), pref)
gr = cg2_line_point_signed_distance(p0, p1, polygon%point(:,r), pref)
gl0 = gl
gr0 = gr

! A4
if (gl0 <= gr0) then
a = l
b = l + fib(n)
l = a + fib(n-2)
gl = cg2_line_point_signed_distance(p0, p1,&
polygon%point(:,modulo(l - 1, polygon%nb_points) + 1), pref)
else
b = r
a = r - fib(n)
r = b - fib(n-2)
gr = cg2_line_point_signed_distance(p0, p1,&
polygon%point(:,modulo(r - 1, polygon%nb_points) + 1), pref)
end if

! A5
i = 2
do while (i > 1)
if (gl <= gr) then
tmp = gr
lk = l
l = r
r = 2*lk - a
a = lk
gr = cg2_line_point_signed_distance(p0, p1,&
polygon%point(:,modulo(r - 1, polygon%nb_points) + 1), pref)
i = b - lk
gl = tmp
else
tmp = gl
rk = r
r = l
l = 2*rk - b
b = rk
gl = cg2_line_point_signed_distance(p0, p1,&
polygon%point(:,modulo(l - 1, polygon%nb_points) + 1), pref)
i = rk - a
gr = tmp
end if
end do

if (gl >= gr) then
maxloc_distance = modulo(l - 1, polygon%nb_points) + 1
else
maxloc_distance = modulo(r - 1, polygon%nb_points) + 1
end if

maximum_distance = tmp
end subroutine cg2_maxloc_line_point_signed_distance

end module mod_cg2_line_polygon
