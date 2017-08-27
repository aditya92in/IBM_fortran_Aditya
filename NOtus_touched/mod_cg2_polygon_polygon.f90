!This file is part of Notus 0.2.0
 
 !Copyright Bordeaux-INP, Université de Bordeaux, CNRS
 
 !Contributors:
 !Antoine Lemoine, 24-07-2015, antoine.lemoine@bordeax-inp.fr
 
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
 
 
 
 module mod_cg2_polygon_polygon
    use mod_cg2_points
    use mod_cg2_line_segment
    use mod_cg2_line_polygon
    use mod_cg2_polygon
    implicit none
 
 contains
 
    pure subroutine cg2_convex_polygon_intersection(polygon1, polygon2, polygon_intersection, intersection_found)
       type(polygon_2d), intent(in) :: polygon1, polygon2
       type(polygon_2d), intent(out) :: polygon_intersection
       logical, intent(out) :: intersection_found
 
       integer, parameter :: INSIDE_POLYGON1 = 1
       integer, parameter :: INSIDE_POLYGON2 = 2
       integer, parameter :: UNDEFINED = -1
 
       integer :: m, n, vertex_number
       integer :: a, b, a1, b1, aa, ba
       integer :: in_flag
       double precision, dimension(2, polygon1%nb_points+polygon2%nb_points+1) :: vertex_list
       double precision, dimension(2) :: vec_a, vec_b, intersection_point
       double precision :: signed_area
       logical :: first_point, is_intersection
       logical :: a_in_half_vec_b, b_in_half_vec_a
       logical :: polygon1_in_polygon2, polygon2_in_polygon1
 
       in_flag = undefined ! 0: inside, 1: outside, other: undefined
       first_point = .true.
       intersection_found = .false.
 
       polygon1_in_polygon2 = .true.
       polygon2_in_polygon1 = .true.
 
       if (polygon1%nb_points < 3 .or. polygon2%nb_points < 3) return
 
       n = polygon1%nb_points
       m = polygon2%nb_points
       vertex_number = 0
       a = 1
       b = 1
       aa = 1
       ba = 1
 
       do while((aa <= n .or. ba <= m) .and. aa <= 2*n .and. ba <= 2*m)
          ! Compute previous points indices
          a1 = modulo(a - 2, n) + 1
          b1 = modulo(b - 2, m) + 1
 
          ! If vec_a and vec_b intersect
          vec_a = polygon1%point(:,a) - polygon1%point(:,a1)
          vec_b = polygon2%point(:,b) - polygon2%point(:,b1)
 
          signed_area = cg2_triangle_double_area([0.0d0, 0.0d0], vec_a, vec_b)
          b_in_half_vec_a = cg2_is_point_left_of_line(polygon1%point(:,a1), polygon1%point(:,a), polygon2%point(:,b))
          a_in_half_vec_b = cg2_is_point_left_of_line(polygon2%point(:,b1), polygon2%point(:,b), polygon1%point(:,a))
 
          polygon2_in_polygon1 = polygon2_in_polygon1 .and. b_in_half_vec_a
          polygon1_in_polygon2 = polygon1_in_polygon2 .and. a_in_half_vec_b
 
          call cg2_segments_intersection(polygon1%point(:,a1), polygon1%point(:,a), polygon2%point(:,b1), polygon2%point(:,b), &
             &                         is_intersection, intersection_point)
 
          if (is_intersection) then
             if (in_flag < 0 .and. first_point) then
                aa = 1
                ba = 1
                first_point = .false.
             end if
 
             vertex_number = vertex_number + 1
             if (vertex_number > size(vertex_list,2)) return
             vertex_list(:,vertex_number) = intersection_point
 
             if (a_in_half_vec_b) then
                in_flag = inside_polygon1
             else if (b_in_half_vec_a) then
                in_flag = inside_polygon2
             end if
          end if
 
          ! Advance rules
          if (signed_area == 0.0d0 .and. .not. a_in_half_vec_b .and. .not. b_in_half_vec_a) then
             if (in_flag == inside_polygon1) then
                ba = ba + 1
                b = modulo(b, m) + 1
             else
                aa = aa + 1
                a = modulo(a, n) + 1
             end if
          else if (signed_area >= 0.0d0) then
             if (b_in_half_vec_a) then
                if (in_flag == inside_polygon1) then
                   vertex_number = vertex_number + 1
                   if (vertex_number > size(vertex_list,2)) return
                   vertex_list(:,vertex_number) = polygon1%point(:,a)
                end if
                aa = aa + 1
                a = modulo(a, n) + 1
             else
                if (in_flag == inside_polygon2) then
                   vertex_number = vertex_number + 1
                   if (vertex_number > size(vertex_list,2)) return
                   vertex_list(:,vertex_number) = polygon2%point(:,b)
                end if
                ba = ba + 1
                b = modulo(b, m) + 1
             end if
          else
             if (a_in_half_vec_b) then
                if (in_flag == inside_polygon2) then
                   vertex_number = vertex_number + 1
                   if (vertex_number > size(vertex_list,2)) return
                   vertex_list(:,vertex_number) = polygon2%point(:,b)
                end if
                ba = ba + 1
                b = modulo(b, m) + 1
             else
                if (in_flag == inside_polygon1) then
                   vertex_number = vertex_number + 1
                   if (vertex_number > size(vertex_list,2)) return
                   vertex_list(:,vertex_number) = polygon1%point(:,a)
                end if
                aa = aa + 1
                a = modulo(a, n) + 1
             end if
          end if
       end do
 
       ! Create the intersection polygon if an intersection was found
       if (in_flag > 0) then
          call cg2_initialize_polygon(polygon_intersection, vertex_number)
          polygon_intersection%point = vertex_list(:,:vertex_number)
          intersection_found = .true.
       else if (polygon1_in_polygon2) then
          call cg2_initialize_polygon(polygon_intersection, polygon1%nb_points)
          polygon_intersection%point = polygon1%point
 
          intersection_found = .true.
       else if (polygon2_in_polygon1) then
          call cg2_initialize_polygon(polygon_intersection, polygon2%nb_points)
          polygon_intersection%point = polygon2%point
 
          intersection_found = .true.
       end if
    end subroutine cg2_convex_polygon_intersection
 
    pure subroutine cg2_convex_polygon_clipping(reference_polygon, clip_polygon, clipped_polygon, intersection_found)
       type(polygon_2d), intent(in) :: reference_polygon, clip_polygon
       type(polygon_2d), intent(out) :: clipped_polygon
       logical, intent(out) :: intersection_found
 
       type(polygon_2d) :: tmp_polygon
       double precision, dimension(2) :: l0, l1
       integer :: i
 
       intersection_found = .false.
 
       call initialize(tmp_polygon, reference_polygon%nb_points + clip_polygon%nb_points)
 
       tmp_polygon%nb_points = clip_polygon%nb_points
       tmp_polygon%point(:,:tmp_polygon%nb_points) = clip_polygon%point(:,:tmp_polygon%nb_points)
 
       call initialize(clipped_polygon, reference_polygon%nb_points + clip_polygon%nb_points)
 
       clipped_polygon%nb_points = clip_polygon%nb_points
       clipped_polygon%point(:,:clipped_polygon%nb_points) = clip_polygon%point(:,:clipped_polygon%nb_points)
 
       l0 = reference_polygon%point(:,reference_polygon%nb_points)
       do i = 1, reference_polygon%nb_points
          l1 = reference_polygon%point(:,i)
 
          call cg2_line_polygon_clipping(clipped_polygon, l0, l1, tmp_polygon)
 
          if (tmp_polygon%nb_points < 3) then
             call unalloc(tmp_polygon)
             return
          end if
 
          if (i == reference_polygon%nb_points) exit
 
          clipped_polygon%nb_points = tmp_polygon%nb_points
          clipped_polygon%point(:,:tmp_polygon%nb_points) = tmp_polygon%point(:,:tmp_polygon%nb_points)
 
          l0 = l1
       end do
 
       call initialize(clipped_polygon, tmp_polygon%nb_points)
       clipped_polygon%point = tmp_polygon%point(:,:tmp_polygon%nb_points)
 
       call unalloc(tmp_polygon)
       intersection_found = .true.
    end subroutine cg2_convex_polygon_clipping
 
    pure subroutine cg2_polygon_boundary_intersection(polygon1, polygon2, s0, s1, i, j, intersection_found)
       type(polygon_2d), intent(in) :: polygon1, polygon2
       double precision, dimension(2), intent(out) :: s0, s1
       integer, intent(out) :: i, j
       logical, intent(out) :: intersection_found
 
       intersection_found = .false.
       s0 = 0.0d0; s1 = 0.0d0;
 
       do i = 1, polygon1%nb_points-1
          do j = 1, polygon2%nb_points-1
             call cg2_segments_overlap(polygon1%point(:,i), polygon1%point(:,i+1), &
                &                    polygon2%point(:,j), polygon2%point(:,j+1), &
                &                    s0, s1, intersection_found)
 
             if (intersection_found) return
          end do
 
          ! Last edge of the second polygon
          j = polygon2%nb_points
 
          call cg2_segments_overlap(polygon1%point(:,i), polygon1%point(:,i+1), &
             &                    polygon2%point(:,j), polygon2%point(:,1), &
             &                    s0, s1, intersection_found)
 
          if (intersection_found) return
       end do
 
       ! Last edge of the first polygon
       i = polygon1%nb_points
 
       do j = 1, polygon2%nb_points-1
          call cg2_segments_overlap(polygon1%point(:,i), polygon1%point(:,1), &
             &                    polygon2%point(:,j), polygon2%point(:,j+1), &
             &                    s0, s1, intersection_found)
 
          if (intersection_found) return
       end do
 
       ! Last edge of the second polygon
       j = polygon2%nb_points
 
       call cg2_segments_overlap(polygon1%point(:,i), polygon1%point(:,1), &
          &                    polygon2%point(:,j), polygon2%point(:,1), &
          &                    s0, s1, intersection_found)
    end subroutine cg2_polygon_boundary_intersection
 
 end module mod_cg2_polygon_polygon
