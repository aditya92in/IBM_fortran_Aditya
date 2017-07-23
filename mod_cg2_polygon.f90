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
 
 
 
 module mod_cg2_polygon
    use mod_cg2_points
    implicit none
 
    type polygon_2d
       double precision, dimension(:,:), allocatable :: point
       integer :: nb_points = 0
    end type polygon_2d
 
    interface initialize
       module procedure cg2_initialize_polygon
    end interface initialize
 
    interface unalloc
       module procedure cg2_finalize_polygon
    end interface unalloc
 
 contains
 
    pure subroutine cg2_initialize_polygon(polygon, n)
       type(polygon_2d), intent(inout) :: polygon
       integer, intent(in) :: n
 
       if (allocated(polygon%point)) then
          if (polygon%nb_points /= n) then
             deallocate(polygon%point)
             allocate(polygon%point(2,n))
          end if
       else
          allocate(polygon%point(2,n))
       end if
 
       polygon%point = 0.0d0
       polygon%nb_points = n
    end subroutine cg2_initialize_polygon
 
    elemental subroutine cg2_finalize_polygon(polygon)
       type(polygon_2d), intent(inout) :: polygon
 
       if (allocated(polygon%point)) deallocate(polygon%point)
       polygon%nb_points = 0
    end subroutine cg2_finalize_polygon
 
    double precision pure function cg2_polygon_volume(polygon) result(volume)
       type(polygon_2d), intent(in) :: polygon
 
       integer :: i
       double precision :: x0, y0
 
       x0 = polygon%point(1,1)
       y0 = polygon%point(2,1)
 
       volume = 0.0d0
 
       do i = 2, polygon%nb_points - 1
          volume = volume + (polygon%point(1,i) - x0)*(polygon%point(2,i+1) - y0) &
             &          - (polygon%point(2,i) - y0)*(polygon%point(1,i+1) - x0)
       end do
 
       volume = 0.5d0*abs(volume)
    end function cg2_polygon_volume
 
    pure subroutine cg2_polygon_centroid(polygon, centroid, volume)
       type(polygon_2d), intent(in) :: polygon
       double precision, dimension(2), intent(out) :: centroid
       double precision, intent(out), optional :: volume
 
       integer :: i
       double precision, dimension(2) :: pref
       double precision :: cross_product, polygon_volume
 
       pref = polygon%point(:,1)
 
       polygon_volume = 0.0d0
       centroid = 0.0d0
 
       do i = 2, polygon%nb_points - 1
          cross_product = (polygon%point(1,i) - pref(1))*(polygon%point(2,i+1) - pref(2)) &
             &        - (polygon%point(2,i) - pref(2))*(polygon%point(1,i+1) - pref(1))
 
          polygon_volume = polygon_volume + cross_product
 
          centroid = centroid + (polygon%point(:,i) - pref + polygon%point(:,i+1) - pref)*cross_product
       end do
 
       polygon_volume = 0.5d0*abs(polygon_volume)
 
       if (polygon_volume > 0.0d0) then
          centroid = pref + centroid/(6.0d0*polygon_volume)
       end if
 
       if (present(volume)) volume = polygon_volume
    end subroutine cg2_polygon_centroid
 
    integer pure function cg2_polygon_winding_number(polygon, point) result(winding_number)
       type(polygon_2d), intent(in) :: polygon
       double precision, dimension(2), intent(in) :: point
 
       integer :: i
 
       winding_number = 0
 
       do i = 1, polygon%nb_points-1
          if (polygon%point(2,i) <= point(2)) then
             if (polygon%point(2,i+1) <= point(2)) cycle
             ! Check if the point is on the left of the segment
             if (.not. cg2_is_point_left_of_line(polygon%point(:,i), polygon%point(:,i+1), point)) cycle
 
             winding_number = winding_number + 1
          else
             if (polygon%point(2,i+1) > point(2)) cycle
             ! Check if the point is on the right of the segment
             if (cg2_is_point_left_of_line(polygon%point(:,i), polygon%point(:,i+1), point)) cycle
 
             winding_number = winding_number - 1
          end if
       end do
 
       ! Last segment [n+1, 1]
       if (polygon%point(2,polygon%nb_points) <= point(2)) then
          if (polygon%point(2,1) <= point(2)) return
          ! Check if the point is on the left of the segment
          if (.not. cg2_is_point_left_of_line(polygon%point(:,polygon%nb_points), polygon%point(:,1), point)) return
 
          winding_number = winding_number + 1
       else
          if (polygon%point(2,1) > point(2)) return
          ! Check if the point is on the right of the segment
          if (cg2_is_point_left_of_line(polygon%point(:,polygon%nb_points), polygon%point(:,1), point)) return
 
          winding_number = winding_number - 1
       end if
    end function cg2_polygon_winding_number
 
    elemental subroutine cg2_move_alloc_polygon(polygon1, polygon2)
       type(polygon_2d), intent(inout) :: polygon1, polygon2
 
       call move_alloc(polygon1%point, polygon2%point)
 
       polygon2%nb_points = polygon1%nb_points
       polygon1%nb_points = 0
    end subroutine cg2_move_alloc_polygon
 
    subroutine cg2_polygon_write_vtk_file(polygon, filename)
       type(polygon_2d), intent(in) :: polygon
       character(len=*), intent(in) :: filename
 
       integer :: gridunit, i
 
       open(newunit=gridunit, file=trim(filename), status="unknown")
 
       write(gridunit,'("# vtk DataFile Version 3.0")')
       write(gridunit,'(a)') trim(filename)
       write(gridunit,'("ASCII")')
       write(gridunit,'("DATASET POLYDATA")')
       write(gridunit,*)
       write(gridunit,'("POINTS ",i8," float")') polygon%nb_points
 
       do i = 1, polygon%nb_points
          write(gridunit,*) sngl(polygon%point(:,i)), sngl(0.0d0)
       end do
 
       write(gridunit,*)
       write(gridunit,'("POLYGONS 1 ",i8)') polygon%nb_points+1
 
       write(gridunit,'(i8)',advance="no") polygon%nb_points
       do i = 1, polygon%nb_points
          write(gridunit,'(i8)',advance="no") i-1
       end do
 
       close(gridunit)
    end subroutine cg2_polygon_write_vtk_file
 
 end module mod_cg2_polygon
