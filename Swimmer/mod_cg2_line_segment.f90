!This file is part of Notus 0.2.0
 
 !Copyright Bordeaux-INP, Université de Bordeaux, CNRS
 !
 !Contributors:
 !Antoine Lemoine, 23-07-2015, antoine.lemoine@bordeaux-inp.fr
 
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
 
 
 
 module mod_cg2_line_segment
    use mod_cg2_points
    implicit none
 
 contains
 
    !-------------------------
    !   Detection functions
    !-------------------------
 
    logical pure function cg2_does_segments_intersect(p0, p1, q0, q1) result(r)
       double precision, dimension(2), intent(in) :: p0, p1, q0, q1
 
       ! Check if an extreme point of one line segment is on the other line segment
       if (cg2_is_point_on_segment(p0, p1, q0)) then
          r = .true.
          return
       else if (cg2_is_point_on_segment(p0, p1, q1)) then
          r = .true.
          return
       else if (cg2_is_point_on_segment(q0, q1, p0)) then
          r = .true.
          return
       else if (cg2_is_point_on_segment(q0, q1, p1)) then
          r = .true.
          return
       end if
 
       ! Check whether an extreme point of one line segment is equal to an exteme point of the other line segment
       if (all(p0 == q0) .or. all(p0 == q1) .or. all(p1 == q0) .or. all(p1 == q1)) then
          r = .true.
          return
       end if
 
       ! Check if two line segments intersect [@p p0, @p p1] ∩ [@p q0, @p q1]
       r    =       (cg2_are_counterclockwise(p0, q0, q1) .neqv. cg2_are_counterclockwise(p1, q0, q1)) &
          & .and. (cg2_are_counterclockwise(p0, p1, q0) .neqv. cg2_are_counterclockwise(p0, p1, q1))
    end function cg2_does_segments_intersect
 
    logical pure function cg2_does_line_intersect_segment(l0, l1, s0, s1) result(r)
       double precision, dimension(2), intent(in) :: l0, l1, s0, s1
 
       if (all(l0 == s0) .or. all(l0 == s1) .or. all(l1 == s0) .or. all(l1 == s1)) then
          r = .true.
          return
       end if
 
       r = cg2_are_counterclockwise(l0, l1, s1) .neqv. cg2_are_counterclockwise(l0, l1, s0)
    end function cg2_does_line_intersect_segment
 
    !---------------------------
    !   Computation functions
    !---------------------------
 
    double precision pure function cg2_line_point_distance(l0, l1, p) result(distance)
       double precision, dimension(2), intent(in) :: l0
       double precision, dimension(2), intent(in) :: l1
       double precision, dimension(2), intent(in) :: p
 
       double precision, dimension(2) :: w
 
       w = l1 - l0
 
       distance = abs(((p(1)-l0(1))*w(2) - (p(2)-l0(2))*w(1))/norm2(w))
    end function cg2_line_point_distance
 
    double precision pure function cg2_line_point_signed_distance(l0, l1, p, pref) result(distance)
       double precision, dimension(2), intent(in) :: l0
       double precision, dimension(2), intent(in) :: l1
       double precision, dimension(2), intent(in) :: p
       double precision, dimension(2), intent(in) :: pref
 
       double precision, dimension(2) :: w
 
       w = l1 - l0
 
       distance = ((p(1)-l0(1))*w(2) - (p(2)-l0(2))*w(1)) &
          &   * sign(1.0d0, (pref(1)-l0(1))*w(2) - (pref(2)-l0(2))*w(1)) &
          &   / norm2(w)
    end function cg2_line_point_signed_distance
 
    pure subroutine cg2_point_line_segment_closest_point(s0, s1, p, point)
       double precision, dimension(2), intent(in) :: s0, s1, p
       double precision, dimension(2), intent(out) :: point
 
       double precision, dimension(2) :: u, v
       double precision :: length
 
       u = s1 - s0
       v = p - s0
 
       length = norm2(u)
 
       if (length < tiny(1.0d0)) then
          point = s0
          return
       end if
 
       u = u/length
 
       ! Orthogonal projection
       point = s0 + dot_product(v, u)*u
 
       ! Test if the point belongs to the segment
       if (abs(s1(1) - s0(1)) > abs(s1(2) - s0(2))) then
          if (s0(1) < s1(1)) then
             if (point(1) <= s0(1)) then
                point = s0
             else if (point(1) >= s1(1)) then
                point = s1
             end if
          else
             if (point(1) <= s1(1)) then
                point = s1
             else if (point(1) >= s0(1)) then
                point = s0
             end if
          end if
       else
          if (s0(2) < s1(2)) then
             if (point(2) <= s0(2)) then
                point = s0
             else if (point(2) >= s1(2)) then
                point = s1
             end if
          else
             if (point(2) <= s1(2)) then
                point = s1
             else if (point(2) >= s0(2)) then
                point = s0
             end if
          end if
       end if
    end subroutine cg2_point_line_segment_closest_point
 
    pure subroutine cg2_line_segment_intersection(l0, l1, s0, s1, is_intersection, intersection_point)
       double precision, dimension(2), intent(in) :: l0, l1, s0, s1
       logical, intent(out) :: is_intersection
       double precision, dimension(2), intent(out) :: intersection_point
 
       double precision, dimension(2) :: u, v, w
       double precision :: alpha_i, perp, lu, lv, epsperp
 
       is_intersection = cg2_does_line_intersect_segment(l0, l1, s0, s1)
       intersection_point = 0.0d0
 
       if (.not. is_intersection) return
 
       u = s1 - s0
       v = l1 - l0
       w = s0 - l0
       lu = norm2(u)
       lv = norm2(v)
 
       perp = v(1)*u(2) - v(2)*u(1)
 
       epsperp = max(lu, lv)*epsilon(1.0d0)
       if (abs(perp) < epsperp) then
          is_intersection = .false.
          return
       end if
 
       alpha_i = (v(2)*w(1) - v(1)*w(2))/perp
       intersection_point = s0 + alpha_i*u
 
       if (alpha_i < 0.0d0) then
          intersection_point = s0
       else if (alpha_i > 1.0d0) then
          intersection_point = s1
       end if
    end subroutine cg2_line_segment_intersection
 
    pure subroutine cg2_line_segment_intersection2(origin, v, s0, s1, is_intersection, intersection_point)
       double precision, dimension(2), intent(in) :: origin, v, s0, s1
       logical, intent(out) :: is_intersection
       double precision, dimension(2), intent(out) :: intersection_point
 
       double precision, dimension(2) :: u, w
       double precision :: alpha_i, perp, lu, lv, epsperp
 
       is_intersection = cg2_does_line_intersect_segment(origin, origin+v, s0, s1)
       intersection_point = 0.0d0
 
       if (.not. is_intersection) return
 
       u = s1 - s0
       w = s0 - origin
       lu = norm2(u)
       lv = norm2(v)
 
       perp = v(1)*u(2) - v(2)*u(1)
 
       epsperp = max(lu, lv)*epsilon(1.0d0)
       if (abs(perp) < epsperp) then
          is_intersection = .false.
          return
       end if
 
       alpha_i = (v(2)*w(1) - v(1)*w(2))/perp
       intersection_point = s0 + alpha_i*u
 
       if (alpha_i < 0.0d0) then
          intersection_point = s0
       else if (alpha_i > 1.0d0) then
          intersection_point = s1
       end if
    end subroutine cg2_line_segment_intersection2
 
    pure subroutine cg2_line_segment_inclusive_intersection(l0, l1, s0, s1, is_intersection, intersection_point)
       double precision, dimension(2), intent(in) :: l0, l1, s0, s1
       logical, intent(out) :: is_intersection
       double precision, dimension(2), intent(out) :: intersection_point
 
       double precision, dimension(2) :: u, v, w
       double precision :: lu, lv, perp_uv, perp_vw, alpha_i
 
       is_intersection = .false.
 
       u = s1 - s0
       v = l1 - l0
       w = s0 - l0
       lu = norm2(u)
       lv = norm2(v)
 
       ! Check if u and v are collinear
       perp_uv = v(1)*u(2) - v(2)*u(1)
 
       ! Check if s0 belongs to the line (l0, l1)
       perp_vw = v(1)*w(2) - v(2)*w(1)
 
       if (perp_uv == 0.0d0) then
          ! If u and v are collinear and s0 does not belong to (l0, l1), no intersection
          if (abs(perp_vw) > 0.0d0) return
 
          ! Return the second point of the segment
          intersection_point = s1
          is_intersection = .true.
          return
       end if
 
       is_intersection = .true.
 
       ! Compute the intersection
       alpha_i = -perp_vw/perp_uv
       intersection_point = s0 + alpha_i*u
 
       if (alpha_i < 0.0d0) then
          intersection_point = s0
       else if (alpha_i > 1.0d0) then
          intersection_point = s1
       end if
    end subroutine cg2_line_segment_inclusive_intersection
 
    pure subroutine cg2_segments_intersection(p0, p1, q0, q1, is_intersection, intersection_point)
       double precision, dimension(2), intent(in) :: p0, p1, q0, q1
       logical, intent(out) :: is_intersection
       double precision, dimension(2), intent(out) :: intersection_point
 
       double precision, dimension(2) :: u, v, w
       double precision :: alpha_i, perp, lu, lv, epsperp
 
       is_intersection = cg2_does_segments_intersect(p0, p1, q0, q1)
       intersection_point = 0.0d0
 
       if (.not. is_intersection) return
 
       u = q1 - q0
       v = p1 - p0
       w = q0 - p0
       lu = norm2(u)
       lv = norm2(v)
 
       perp = v(1)*u(2) - v(2)*u(1)
 
       epsperp = tiny(0.0d0)
       if (abs(perp) < epsperp) then
          is_intersection = .false.
          return
       end if
 
       alpha_i = (v(2)*w(1) - v(1)*w(2))/perp
       intersection_point = q0 + alpha_i*u
 
       if (alpha_i < 0.0d0) then
          intersection_point = q0
       else if (alpha_i > 1.0d0) then
          intersection_point = q1
       end if
    end subroutine cg2_segments_intersection
 
    pure subroutine cg2_segments_overlap(p0, p1, q0, q1, s0, s1, are_overlapping)
       double precision, dimension(2), intent(in) :: p0, p1, q0, q1
       double precision, dimension(2), intent(out) :: s0, s1
       logical, intent(out) :: are_overlapping
 
       double precision :: perp_product, distance, length_p, length_q, alpha
       double precision, dimension(2) :: vector_p, vector_q
 
       s0 = 0.0d0; s1 = 0.0d0
       are_overlapping = .false.
 
       ! Compute vectors from the edge points
       vector_p = p1 - p0
       vector_q = q1 - q0
 
       ! Check the length of each vectors
       length_p = norm2(vector_p)
       if (length_p <= tiny(0.0d0)) return
 
       length_q = norm2(vector_q)
       if (length_q <= tiny(0.0d0)) return
 
       ! Normalize vector_p and vector_q
       vector_p = vector_p/length_p
       vector_q = vector_q/length_q
 
       ! Compute the perpendicular product of the two current line segments
       perp_product = vector_p(1)*vector_q(2) - vector_p(2)*vector_q(1)
 
       ! Check if the two line segments are **almost** parallel
       if (abs(perp_product) > 1.0d2*epsilon(1.0d0)) return
 
       ! Compute the distance between the two line segments. Since they are considered parallel, just compute the distance
       ! from one point of the second line segment to the first line segment.
       distance = cg2_line_point_distance(p0, p1, q0)
 
       ! Check if the two segments are close enough to each other
       if (distance > 1.0d2*epsilon(1.0d0)) return
 
       ! Compute the first point of the final segment
       alpha = dot_product(vector_p, q0 - p0)
       if (alpha <= epsilon(1.0d0)) then
          s0 = p0
       else if (alpha/length_p - 1.0d0 >= epsilon(1.0d0)) then
          s0 = p1
       else
          s0 = p0 + alpha*vector_p
       end if
 
       ! Compute the second point of the final segment
       alpha = dot_product(vector_p, q1 - p0)
       if (alpha <= epsilon(1.0d0)) then
          s1 = p0
       else if (alpha/length_p - 1.0d0 >= epsilon(1.0d0)) then
          s1 = p1
       else
          s1 = p0 + alpha*vector_p
       end if
 
       ! Check the distance between the two points of the segment
       if (norm2(s1 - s0) <= epsilon(1.0d0)) return
 
       are_overlapping = .true.
    end subroutine cg2_segments_overlap
 
    subroutine cg2_line_write_vtk_file(s0, s1, filename)
       double precision, dimension(2) :: s0, s1
       character(len=*), intent(in) :: filename
 
       integer :: gridunit
 
       open(newunit=gridunit, file=trim(filename), status="unknown")
 
       write(gridunit,'("# vtk DataFile Version 3.0")')
       write(gridunit,'(a)') trim(filename)
       write(gridunit,'("ASCII")')
       write(gridunit,'("DATASET POLYDATA")')
       write(gridunit,*)
       write(gridunit,'("POINTS ",i8," float")') 2
 
       write(gridunit,*) sngl(s0), sngl(0.0d0)
       write(gridunit,*) sngl(s1), sngl(0.0d0)
 
       write(gridunit,*)
       write(gridunit,'("LINES 1 ",i8)') 3
 
       write(gridunit,'("2 0 1")')
 
       close(gridunit)
    end subroutine cg2_line_write_vtk_file
 
 end module mod_cg2_line_segment
