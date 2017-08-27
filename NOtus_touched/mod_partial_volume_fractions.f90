module mod_partial_volume_fractions
   !use variables_mpi
   !use variables_mpif
   use mod_cg_transformation
   use mod_cg2_polygon
   use mod_cg2_polygon_polygon
   use mod_polyarea
   implicit none
contains 
subroutine set_phi(z,r,ilow,ihigh,jlow,jhigh,dz,dr,phi,z_general,r_general,n)
implicit none 
type(polygon_2d) :: polygon1, polygon2, intersection
logical :: is_intersection;
real,allocatable,dimension(:):: a,b;
real start,finish;
real,allocatable,dimension(:)::z,r,z_general,r_general;
real, allocatable, dimension(:,:)::phi;
double precision, dimension(2)::centroid_1,centroid_2
integer i,j
integer ilow,ihigh,jlow,jhigh,n;
real dz,dr
double precision vol_temp;
real c_z,c_r!area_circle;
!REAL, PARAMETER :: pi = 3.14159271
call cpu_time(start);
open(unit=15,file='z_polygon_1.txt');
open(unit=16,file='r_polygon_1.txt');
open(unit=77,file='r_polygon_2.txt');
open(unit=78,file='z_polygon_2.txt');
!print*,'enter subroutine set phi';
call initialize(polygon1,n)
polygon1%point(2,:)=z_general;!x_c+r*cos(i*h);
polygon1%point(1,:)=r_general;!y_c+r*sin(i*h);
write(15,*) polygon1%point(2,:)
write(16,*) polygon1%point(1,:)
call initialize(polygon2, 4)
!if(ilow==0) then
!ilow=ilow+1
!endif
do j= jlow,jhigh
 do i= ilow,ihigh
      polygon2%point(2,:)=[z(i), z(i),z(i+1), z(i+1)];
      polygon2%point(1,:)= [r(j), r(j+1), r(j+1), r(j) ];
write(77,*) polygon2%point(2,:)
write(78,*) polygon2%point(1,:)
vol_temp=cg2_polygon_volume(polygon2);
!print*,'vol_temp=',vol_temp;
call cg2_polygon_centroid(polygon2,centroid_2,vol_temp);
!print*,'centroid',centroid_2
call cg2_convex_polygon_clipping(polygon1, polygon2, intersection,is_intersection)
!print*,'intersecting',is_intersection;
if(is_intersection.eqv.(.true.)) then
!print*,'it is interesecting'
allocate(a(size(intersection%point(1,:))),b(size(intersection%point(2,:))))
a=intersection%point(1,:);
b=intersection%point(2,:);
vol_temp=cg2_polygon_volume(intersection);
!print*,'vol_temp',vol_temp;
call cg2_polygon_centroid(intersection,centroid_1,vol_temp)

!print*,'phi_a',(polyarea(a,b,size(a))/(dr*dz));
!print*,'centroid_2',centroid_2;
!print*,'centroid_1',centroid_1;
phi(i,j)=(abs(polyarea(a,b,size(a)))/(dr*dz)) 
!*(2*4*atan(1.0)*(abs(r(1)-centroid_1(1))))/ &
!(2*4*atan(1.0)*(abs(r(1)-centroid_2(1))));
!print*,'phi(i,j)',phi(i,j);
deallocate(a,b);
else
phi(i,j)=0.0;
end if
end do 
end do 
call cpu_time(finish);
print*,'total time=',finish-start;
end subroutine set_phi
end module mod_partial_volume_fractions



















