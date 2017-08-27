module mod_partial_volume_fractions
   !use variables_mpi
   !use variables_mpif
   use mod_cg_transformation
   use mod_cg2_polygon
   use mod_cg2_polygon_polygon
   use mod_polyarea
   implicit none
contains 
subroutine set_phi(z,r,ilow,ihigh,jlow,jhigh,dz,dr,phi,z_general,r_general,n,c_z,c_r)
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
REAL, PARAMETER :: pi = 3.14159271
call cpu_time(start);
print*,'enter subroutine set phi';
call initialize(polygon1,n)
polygon1%point(1,:)=z_general;!x_c+r*cos(i*h);
polygon1%point(2,:)=r_general;!y_c+r*sin(i*h);
call initialize(polygon2, 4)
!if(ilow==0) then
!ilow=ilow+1
!endif
do j= jlow,jhigh
 do i= ilow,ihigh
      polygon2%point(1,:)=[z(j)-dz/2, z(j+1)-dz/2, z(j+1)-dz/2 ,z(j)-dz/2];
      polygon2%point(2,:)= [r(i)-dr/2, r(i)-dr/2, r(i+1)-dr/2, r(i+1)-dr/2];
vol_temp=cg2_polygon_volume(polygon2);
!print*,'vol_temp=',vol_temp;
call cg2_polygon_centroid(polygon2,centroid_2,vol_temp);
!print*,'centroid',centroid_2
call cg2_convex_polygon_clipping(polygon1, polygon2, intersection,is_intersection)
!print*,'intersecting',is_intersection;
if(is_intersection.eqv.(.true.)) then
!print*,'it is interesecting'
a=intersection%point(1,:);
b=intersection%point(2,:);
vol_temp=cg2_polygon_volume(intersection);
!print*,'vol_temp',vol_temp;
call cg2_polygon_centroid(intersection,centroid_1,vol_temp)

!print*,'phi_a',(polyarea(a,b,size(a))/(dr*dz));
!print*,'centroid_2',centroid_2;
!print*,'centroid_1',centroid_1;
phi(i,j)=(polyarea(a,b,size(a))/(dr*dz)) &
*(2*pi*(abs(r(1)-centroid_1(2))))/ &
(2*pi*(abs(r(1)-centroid_2(2))));
!print*,'phi(i,j)',phi(i,j);
else
phi(i,j)=0.0;
end if
end do 
end do 
call cpu_time(finish);
print*,'total time=',finish-start;
end subroutine set_phi
end module mod_partial_volume_fractions



















