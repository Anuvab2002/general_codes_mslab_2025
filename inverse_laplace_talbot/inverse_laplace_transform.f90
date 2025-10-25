module inv_lap_trans_module
use global_values_function
use grids
use quadratures
  implicit none
contains
  !Function for calculating the weight factor term (due to Jacobian conversion)
  ! the term is = (1+i*\sigma(\theta));
  ! \sigma(\theta) = \cot(\theta)(\theta \cot(\theta)-1)
  function calculate_weight_factor(talbot_var_value)result(weight_factor)
    implicit none
    real(8), intent(in)    :: talbot_var_value
    complex(8)             :: weight_factor
!
    real(8)                :: sigma
    real(8)                :: cotangent
!
    cotangent = 1.0d0/tan(talbot_var_value)
    sigma = talbot_var_value + cotangent*(talbot_var_value*cotangent - 1.0d0)
!
    weight_factor = 1.0d0 + iota*sigma
  end function calculate_weight_factor
!
  !the main inverse laplace code
  subroutine talbot_inverse_laplace(lap_fn, lap_var, inv_lap_var, inv_lap_value)
    !lap_fn = function in the laplace space
    !lap_var = variable of the laplace space
    !inv_lap_var = variable of the inverse laplace function
    !inv_laplce_val = value of the inverse laplace function at particular point
    implicit none
    complex(8), intent(in)      :: lap_fn(:)! array of the laplace space function
    complex(8), intent(in)      :: lap_var(:)! array of the laplace space variable
    real(8), intent(in)         :: inv_lap_var! inverse laplace variable
    real(8), intent(out)        :: inv_lap_value! inverse laplace function value
!
    real(8)                     :: t_left
    real(8)                     :: t_right
    integer                     :: t_n_grid! number of talbot variable grid points
    real(8)                     :: t_grid_space
    real(8), allocatable        :: talbot_var(:)
    integer                     :: k
    real(8)                     :: t_prefactor
    complex(8), allocatable     :: weight_factor(:)
    real(8)                     :: exponential_term
    real(8), allocatable        :: integrand(:)
    logical                     :: status
!
    t_left = -pi
    t_right = pi
    t_n_grid = 1000
!
    allocate(talbot_var(t_n_grid))
    call equispaced_grid(t_left, t_right, t_n_grid, t_grid_space, talbot_var)
    t_prefactor = (2.0d0*t_n_grid)/(5.0d0*inv_lap_var)
!
    allocate(weight_factor(t_n_grid-2))
    do k = 1,(t_n_grid-2)
      weight_factor(k) = calculate_weight_factor(talbot_var(k+1))
      exponential_term = exp(inv_lap_var*lap_var(k))
      integrand(k) = real(exponential_term*lap_fn(k)*weight_factor(k))
    end do
    call simpson_onethird(integrand, t_grid_space, t_right, t_left, status, inv_lap_value)
    inv_lap_value = inv_lap_value * t_prefactor/(2.0d0*pi)
  end subroutine talbot_inverse_laplace
!
end module inv_lap_trans_module
