module mod_partial_volume_fractions
   !use variables_mpi
   !use variables_mpif
   use mod_cg_transformation
   use mod_cg2_polygon
   use mod_cg2_polygon_polygon
   use mod_polyarea
   implicit none
contains 
subroutine set_phi(x,y,ilow,ihigh,jlow,jhigh,dx,dy,phi_a,phi_v,x_general,y_general,n)
implicit none 
type(polygon_2d) :: polygon1, polygon2, intersection
logical :: is_intersection;
real,allocatable,dimension(:):: a,b;
real start,finish;
real,allocatable,dimension(:)::x,y,x_general,y_general;
real, allocatable, dimension(:,:)::phi_a,phi_v;
integer i,j
integer ilow,ihigh,jlow,jhigh,n;
real dx,dy;
real(kind=8) vol_temp;
double precision, dimension(2)::centroid_1,centroid_2
vol_temp=0.0;
call cpu_time(start);
call initialize(polygon1,n)
polygon1%point(1,:)=x_general;
polygon1%point(2,:)=y_general;
call initialize(polygon2, 4)
do j=jlow,jhigh
!if(j==jlow) then
  do i=ilow,ihigh
     vol_temp=0.0;
       polygon2%point(1,:)=[x(j), x(j+1), x(j+1) ,x(j)];
       polygon2%point(2,:)= [y(i), y(i), y(i+1), y(i+1)];
       vol_temp=cg2_polygon_volume(polygon2);
       call cg2_polygon_centroid(polygon2,centroid_2,vol_temp);
       call cg2_convex_polygon_clipping(polygon1, polygon2, intersection,is_intersection)
        if(is_intersection.eqv.(.true.)) then
         !print*,'is intersecting in jlow=1';
         a=intersection%point(1,:);
         b=intersection%point(2,:);
         phi_a(i,j)=polyarea(a,b,size(a))/(dx*dy);
         vol_temp=cg2_polygon_volume(intersection);
         call cg2_polygon_centroid(intersection,centroid_1,vol_temp)
          phi_v(i,j)=phi_a(i,j)*(2*4*atan(1.0)*(abs(x(1)-centroid_1(1))))/ &
               (2*4*atan(1.0)*(abs(x(1)-centroid_2(1))));
        else
        phi_a(i,j)=0.0;
        phi_v(i,j)=0.0;
        end if 
end do
!else
 !do i=ilow,ihigh
  !    vol_temp=0.0;
   !   polygon2%point(1,:)=[x(j), x(j+1), x(j+1) ,x(j)];
    !  polygon2%point(2,:)= [y(i), y(i), y(i+1), y(i+1)];
     ! vol_temp=cg2_polygon_volume(polygon2);
      !call cg2_polygon_centroid(polygon2,centroid_2,vol_temp);
      !call cg2_convex_polygon_clipping(polygon1, polygon2, intersection,is_intersection)
      !if(is_intersection.eqv.(.true.)) then
      !print*,'is intersecting in jlow!=1';
      !a=intersection%point(1,:);
      !b=intersection%point(2,:);
      !phi_a(i,j)=polyarea(a,b,size(a))/(dx*dy);
      !vol_temp=cg2_polygon_volume(intersection); 
      !call cg2_polygon_centroid(intersection,centroid_1,vol_temp)
      !phi_v(i,j)=phi_a(i,j)*(2*4*atan(1.0)*(abs(x(1)-centroid_1(1))))/ &
      !         (2*4*atan(1.0)*(abs(x(1)-centroid_2(1))));
       !else
       !phi_a(i,j)=0.0;
       !phi_v(i,j)=0.0;
       !end if
  !end do
  !end if 
end do
call cpu_time(finish);
end subroutine set_phi
end module mod_partial_volume_fractions
