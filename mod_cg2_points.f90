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
 
 
 
 module mod_cg2_points
    implicit none
 
 contains
 
    !-------------------------
    !   Detection functions
    !-------------------------
 
    logical pure function cg2_are_counterclockwise(p0, p1, p2) result(r)
       double precision, dimension(2), intent(in) :: p0, p1, p2
 
       r = (p2(2)-p0(2))*(p1(1)-p0(1)) > (p1(2)-p0(2))*(p2(1)-p0(1))
    end function cg2_are_counterclockwise
 
    logical pure function cg2_is_point_left_of_line(l0, l1, p) result(r)
       double precision, dimension(2), intent(in) :: l0, l1, p
 
       r = (l1(1)-l0(1))*(p(2)-l0(2)) > (p(1)-l0(1))*(l1(2)-l0(2))
    end function cg2_is_point_left_of_line
 
    logical pure function cg2_is_point_left_of_or_on_line(l0, l1, p) result(r)
       double precision, dimension(2), intent(in) :: l0, l1, p
 
       r = (l1(1)-l0(1))*(p(2)-l0(2)) >= (p(1)-l0(1))*(l1(2)-l0(2))
    end function cg2_is_point_left_of_or_on_line
 
    logical pure function cg2_are_points_collinear(p0, p1, p2) result(r)
       double precision, dimension(2), intent(in) :: p0, p1, p2
 
       r = (p1(1)-p0(1))*(p2(2)-p0(2)) == (p2(1)-p0(1))*(p1(2)-p0(2))
    end function cg2_are_points_collinear
 
    logical pure function cg2_is_point_on_segment(s0, s1, p) result(r)
       double precision, dimension(2), intent(in) :: s0, s1, p
 
       r = cg2_are_points_collinear(s0, s1, p)
 
       if (.not. r) return
 
       if (s0(1) /= s1(1)) then
          r =         ((s0(1) <= p(1)) .and. (p(1) <= s1(1))) &
             & .or. ((s0(1) >= p(1)) .and. (p(1) >= s1(1)))
       else
          r =         ((s0(2) <= p(2)) .and. (p(2) <= s1(2))) &
             & .or. ((s0(2) >= p(2)) .and. (p(2) >= s1(2)))
       end if
    end function cg2_is_point_on_segment
 
    !---------------------------
    !   Computation functions
    !---------------------------
 
    double precision pure function cg2_perp_product(p0, p1) result(r)
       double precision, dimension(2), intent(in) :: p0, p1
 
       r = p0(1)*p1(2) - p0(2)*p1(1)
    end function cg2_perp_product
 
    double precision pure function cg2_triangle_double_area(p0, p1, p2) result(r)
       double precision, dimension(2), intent(in) :: p0, p1, p2
 
       r = (p1(1)-p0(1))*(p2(2)-p0(2)) - (p2(1)-p0(1))*(p1(2)-p0(2))
    end function cg2_triangle_double_area
 
 end module mod_cg2_points
