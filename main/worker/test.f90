subroutine test (x)
  implicit none
  real(8), dimension(1:16) :: x

  print *, x(1:16)
  return
end subroutine

program a
  implicit none
  real(8), dimension(0:3,0:3) :: y

  y(:,0) = (/1, 2, 3, 4/)
  y(:,1) = (/5, 6, 7, 8/)
  y(:,2) = (/9,10,11,12/)
  y(:,3) = (/13,14,15,16/)
  call test(y)

end program
