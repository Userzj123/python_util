program write_array_to_file
  implicit none
  integer, parameter :: nx=1, ny=2, nz=3
  integer :: i, j, k, t
  real(kind=8) :: data(nx,ny,nz)

  t = 0
  ! Fill the array with values from 0 to 5
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        data(i,j,k) = t
        t = t + 1
      end do
    end do
  end do

  ! Open the binary file for writing
  open(unit=10, file="data.bin", form="unformatted", access="stream")

  ! Write the array to the binary file in column-major order
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        write(unit=10) data(i,j,k)
      end do
    end do
  end do

  ! Close the binary file
  close(unit=10)

end program write_array_to_file