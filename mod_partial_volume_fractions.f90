module mod_partial_volume_fractions
   !use variables_mpi
   !use variables_mpif
   use mod_cg_transformation
   use mod_cg2_polygon
   use mod_cg2_polygon_polygon
   use mod_polyarea
   implicit none
contains 
subroutine set_phi(x,y,ilow,ihigh,jlow,jhigh,dx,dy,phi,x_general,y_general,n)
implicit none 
type(polygon_2d) :: polygon1, polygon2, intersection
logical :: is_intersection;
real,allocatable,dimension(:):: a,b;
real start,finish;
real,allocatable,dimension(:)::x,y,x_general,y_general;
real, allocatable, dimension(:,:)::phi;
integer i,j
integer ilow,ihigh,jlow,jhigh,n;
real dx,dy;
!real x_c,y_c,area_circle;
!REAL, PARAMETER :: pi = 3.14159271
call cpu_time(start);
call initialize(polygon1,n)
polygon1%point(1,:)=x_general;!x_c+r*cos(i*h);
polygon1%point(2,:)=y_general;!y_c+r*sin(i*h);
call initialize(polygon2, 4)
do j=jlow,jhigh
 do i=ilow,ihigh
      polygon2%point(1,:)=[x(j)-dx/2, x(j+1)-dx/2, x(j+1)-dx/2 ,x(j)-dx/2];
      polygon2%point(2,:)= [y(i)-dy/2, y(i)-dy/2, y(i+1)-dy/2, y(i+1)-dy/2];
call cg2_convex_polygon_clipping(polygon1, polygon2, intersection,is_intersection)
!print*,'intersecting',is_intersection;
if(is_intersection.eqv.(.true.)) then
a=intersection%point(1,:);
b=intersection%point(2,:);
phi(i,j)=polyarea(a,b,size(a))/(dx*dy);
else
phi(i,j)=0.0;
end if
end do 
end do 
call cpu_time(finish);
print*,'total time=',finish-start;
end subroutine set_phi
end module mod_partial_volume_fractions



















