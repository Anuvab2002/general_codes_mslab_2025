program test_inverse_laplace_transform
!use global_values_function
!use inv_lap_trans_module
  implicit none
  call inverse_laplace_test_01()
contains
  subroutine inverse_laplace_test_01()
    implicit none
    print*, "Test 1"
  end subroutine inverse_laplace_test_01
end program test_inverse_laplace_transform
