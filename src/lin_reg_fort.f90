module lin_reg_fort

contains

!================================================================
subroutine calc_lin_reg(input_file,x_start, x_end)
 implicit none
!=====
 character(len=100),intent(in) :: input_file
 real*8,intent(in)             :: x_start, x_end
!=====
 character(len=100) :: ch_x_start
 character(len=100) :: ch_x_end
 integer            :: nrow, io, irow, unit_file,niter,iiter
 integer            :: unit_output
 integer            :: num
 logical            :: file_exists
 real*8             :: x_sum, y_sum, xy_sum, x2_sum
 real*8             :: slope, mean_x, mean_y, mean_xy, mean_x2
 real*8,allocatable :: points_arg(:), points_val(:),weight(:)
!=====
 
 num=0
 x_sum=0.d0
 y_sum=0.d0
 xy_sum=0.d0
 x2_sum=0.d0

 inquire(file=input_file,exist=file_exists)
 if(.NOT. file_exists) then
   write(*,'(a)') 'STOP: input file not found'
   stop
 endif

!get number of strings
 nrow=get_number_of_rows(input_file)
 allocate(points_arg(nrow))
 allocate(points_val(nrow))
 allocate(weight(nrow))

 open(newunit=unit_file, file = input_file)
 do irow=1,nrow
   read(unit_file,*)points_arg(irow),points_val(irow) 
 end do
 close(unit_file)

 do irow=1, nrow
   if(points_arg(irow)>=x_start .AND. points_arg(irow)<=x_end) then
     x_sum=x_sum+points_arg(irow)
     y_sum=y_sum+points_val(irow)
     xy_sum=xy_sum+points_arg(irow)*points_val(irow) 
     x2_sum=x2_sum+points_arg(irow)**2.d0
     num=num+1
  end if
 end do 

 mean_x = x_sum / num
 mean_y = y_sum / num
 mean_xy = xy_sum / num
 mean_x2 = x2_sum / num
 slope = (mean_xy - (mean_x*mean_y)) / (mean_x2 - (mean_x*mean_x))

 write(*,*) slope
! open(newunit=unit_output,file="integral.dat")
!   write(unit_int_val,*) aver_int, int_error, av_error,aver_impact
! close(unit_output)

end subroutine calc_lin_reg
!=======================
function get_number_of_rows(input_file) result(num)
 implicit none
 character(len=*),intent(in)  :: input_file
 integer                      :: num
!=====
 integer                      :: unit_file,io
 character(len=100)             :: cur_string
!=====
 num=0
 open(newunit=unit_file, file = input_file)
 do
   read(unit_file,'(A)',iostat=io)cur_string
   if (io/=0) exit
   num = num + 1
 end do
 close (unit_file)
end function get_number_of_rows 

end module lin_reg_fort
!===============================================================
