module grids
use global_values_function
contains
  !subroutine for constricting equispaced grid
  subroutine equispaced_grid(left, right, n_grid, grid_space, grid_array)
    implicit none
    real(8), intent(in)                      :: left! left (lower) bound of grid
    real(8), intent(in)                      :: right! right (upper) bound of grid
    integer, intent(in)                      :: n_grid! number of grid points
    real(8), intent(out)                     :: grid_space! gap between two grid points
    real(8),  dimension(n_grid), intent(out) :: grid_array
!
    integer                                  :: i
!
    grid_space = (right-left)/real(n_grid*1.0d0-1.0d0)
    do i = 1,n_grid
      grid_array(i) = left + (i-1)*grid_space
    end do
  end subroutine equispaced_grid
!
  !function for talbot parametrization
  function talbot_parametrization(talbot_var_value, prefactor)result(lap_var_value)
    implicit none
    real(8), intent(in)       :: talbot_var_value
    real(8), intent(in)       :: prefactor
    complex(8)                :: lap_var_value
!
    real(8)                   :: cotangent
!
    cotangent = 1.0d0/tan(talbot_var_value)
    lap_var_value = prefactor*talbot_var_value*(cotangent + iota)
  end function talbot_parametrization
!
  !subroutine for constructing talbot variable grid
  subroutine talbot_grid(t_left, t_right, t_n_grid, t_prefactor, lap_var)
    implicit none
    real(8), intent(in)                          :: t_left! left bound of talbot variable
    real(8), intent(in)                          :: t_right! right bound of talbot variable
    integer, intent(in)                          :: t_n_grid! number of talbot variable grid points
    real(8), intent(in)                          :: t_prefactor! prefactor in the parametrization
    complex(8), dimension(t_n_grid-2), intent(out) :: lap_var! laplace space variable array
  !
    real(8), dimension(t_n_grid)                 :: talbot_var! Talbot variable array
    real(8)                                      :: t_grid_space! gap between talbot grid points
    integer                                      :: i
!
    call equispaced_grid(t_left, t_right, t_n_grid, t_grid_space, talbot_var)
    do i = 1,(t_n_grid-2)
      lap_var(i) = talbot_parametrization(talbot_var(i+1), t_prefactor)
    end do
  end subroutine talbot_grid
end module grids
