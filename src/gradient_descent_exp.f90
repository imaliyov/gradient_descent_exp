!=========================================================================
! This is the code gradient_descent_exp which allows to fit the data over 
! the regular grid (corresponding to the bcc latice for physical 
! applications) using the adapted gradient descent algorithm
! with the fitting function of on each center of the form:
!                       f_i(x)=A_i*exp(-B_i*x)+C_i, 
! where A,B and C are the fitting parameters of the center "i". 
! Then, the final fitting function has the form:
!                        F(x)=sum_i{f_i(x)}
!
! For the moment the code works for 5 centers with given coordinates 
! which are (ratom(i,:)) but it can be easily generalized to any number 
! and any geometrical configuration of parameters. The code allows precise
! and simple fitting parameters optimization for different particular
! cases. The code also performs the integration of the fitted function
! over the selected region of 2D space (however, this part is still 
! in the developement phase).
!=========================================================================
program gradient_descent_exp
 use rotate_integrate
 use lin_reg_fort
 implicit none
!=====
 character(len=100) :: input_file
 integer            :: nrow, io, irow, unit_file,niter,iiter
 integer            :: unit_output,unit_param,unit_surf,unit_error
 integer            :: unit_int_val, unit_int_err, unit_av_err
 logical            :: file_exists
 real*8             :: gradient(3),var_x,var_y
 real*8             :: par_a,par_alpha,par_b
 real*8             :: step_a,step_alpha,step_b,learning_rate
 real*8             :: prev_update(3),momentum
 real*8,allocatable :: points_arg(:), points_val(:),weight(:)
 integer            :: time1, time2
 real*8             :: calc_time
 integer            :: count_rate,count_max,natom
 real*8             :: ratom(5,2)
 real*8             :: alat
 real*8             :: aver_int,aver_impact
 real*8             :: var_a,var_alpha,var_b,delta_a,delta_alpha,delta_b,da,dalpha,db
 real*8             :: gauss_val1,gauss_val2 
 real*8             :: rad1,rad2,int_error,cur_error,fit_val,av_error
 real*8             :: delta1,delta2
 real*8,parameter   :: pi=3.14159265358979323d0
 real*8,parameter   :: pi2=pi**2
 real*8             :: aver_int_tri
!=====

 alat=3.49d0
 natom=5
 ratom(1,:)=(/  0.d0,       0.d0 /)
 ratom(2,:)=(/ -alat/2.d0,  alat/2.d0 /)
 ratom(3,:)=(/  alat/2.d0,  alat/2.d0 /)
 ratom(4,:)=(/ -alat/2.d0, -alat/2.d0 /)
 ratom(5,:)=(/  alat/2.d0, -alat/2.d0 /)

 call GET_COMMAND_ARGUMENT(1,input_file)

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

 do irow=1,nrow-1
   weight(irow)=points_arg(irow)**2
 end do
 weight(nrow)=(alat**2-pi*points_arg(nrow)**2)
 weight(:)=weight/SUM(points_arg(:))

!parameters initialisation
 par_b=points_val(nrow)
! par_a=points_val(1)-par_b
 par_a=-(points_val(2)-points_val(1))*points_arg(1)/(points_arg(2)-points_arg(1))/2+points_val(1)-par_b
 par_alpha=-1.d0/points_arg(2)**2*LOG((points_val(2)-par_b)/par_a)
! par_b=points_val(nrow)-fit_func_1_gauss(points_arg(nrow),par_a,par_alpha,par_b)
 par_b=points_val(nrow)

! par_a=3.17989
! par_alpha=4231.73
! par_b=0.127919

 step_a=0.1d0
 step_alpha=10000.0d0 ! step_alpha=10000.d0
 step_b=0.1d0
 
 momentum=1.0d0
 prev_update(:)=0.d0

 niter=100000 ! niter=20000 

 write(*,*) par_a,par_alpha,par_b,error_func_1_gauss(nrow,par_a,par_alpha,par_b,points_arg,weight,points_val)

!==========MAIN LOOP====================
 open(newunit=unit_param,file="param.dat")
 call system_clock(COUNT_RATE=count_rate,COUNT_MAX=count_max)
 call system_clock(COUNT=time1)

 do iiter=1,niter
   call calc_gradient_1_gauss(nrow,par_a,par_alpha,par_b,points_arg,points_val,weight,gradient)  
!   call calc_gradient_5_gauss(nrow,ratom,natom,par_a,par_alpha,par_b,points_arg,points_val,gradient)

   par_a     =par_a     - step_a*gradient(1)     - prev_update(1)*momentum 
   par_alpha =par_alpha - step_alpha*gradient(2) - prev_update(2)*momentum  
   par_b     =par_b     - step_b*gradient(3)     - prev_update(3)*momentum

   prev_update(1)=step_a*gradient(1)
   prev_update(2)=step_alpha*gradient(2)
   prev_update(3)=step_b*gradient(3)

   write(unit_param,*) iiter, par_a,par_alpha,par_b
 end do
 call system_clock(COUNT=time2)
 calc_time=MODULO( time2 - time1, count_max) / REAL(count_rate)
 write(*,*) "Time: ", calc_time
 close(unit_param)
!============END MAIN LOOP=================

 write(*,*) par_a,par_alpha,par_b,  & 
            & error_func_1_gauss(nrow,par_a,par_alpha,par_b,points_arg,weight,points_val) !, &
!            & error_func_5_gauss(nrow,ratom,natom,par_a,par_alpha,par_b,points_arg,points_val)

 open(newunit=unit_output,file="fitted.dat")
   var_x=0.d0
   do while( var_x <=1.7d0  )
     write(unit_output,*) var_x, fit_func_1_gauss(var_x,par_a,par_alpha,par_b) !, &
!          & fit_func_5_gauss(ratom,natom,(/ var_x, 0.d0 /),par_a,par_alpha,par_b)
     var_x=var_x+0.001d0
   end do
 close(unit_output)


! write(*,*) "Surface Integral: ", integrate_5_gauss(ratom,natom,par_a,par_alpha,par_b,alat)
! write(*,*) "Estimated integral: ", 2.d0*par_a/(par_alpha*alat**2.d0)*pi+par_b

 gauss_val1=integrate_gaussian(par_alpha,alat,alat/2.d0)
 gauss_val2=integrate_gaussian(par_alpha,alat,alat) 

 aver_int=4.d0*par_a/(alat**2.d0)*(gauss_val1**2.d0+gauss_val2**2.d0)+par_b
 aver_impact=SQRT(1.d0/par_alpha*LOG(par_a/(aver_int-par_b)))
! av_error=ABS(aver_int-fit_func_1_gauss(0.048017279d0,par_a,par_alpha,par_b))/aver_int*100.d0
 av_error=ABS(aver_int-points_val(nrow))/aver_int*100.d0


!Error estimation

! open(newunit=unit_error,file="int_error.dat")
 int_error=0.d0

 do irow=1,nrow-1

   fit_val=fit_func_1_gauss(points_arg(irow),par_a,par_alpha,par_b)
   cur_error=pi*(points_arg(irow+1)**2-points_arg(irow)**2)*(points_val(irow)-fit_val)
   int_error=int_error+cur_error
!   write(unit_error,*) irow,cur_error,int_error
 
 end do

 fit_val=fit_func_1_gauss(points_arg(nrow),par_a,par_alpha,par_b)
 cur_error=(alat**2-pi*points_arg(nrow)**2)*(points_val(nrow)-fit_val)
 int_error=int_error+cur_error
! int_error=ABS(int_error)/aver_int*100.d0

! write(unit_error,*) nrow,cur_error,int_error
! close(unit_error)

 int_error=int_error/alat**2

 open(newunit=unit_int_val,file="integral.dat")
   write(unit_int_val,*) aver_int, int_error, av_error,aver_impact,par_b
 close(unit_int_val)

 !=====INTEGRATION OVER THE TRIANGLE REGION=======
 call calculate_integral(input_file,alat,aver_int_tri)
 write(*,*) "Integration over the triangle: ", aver_int_tri

 !======LINEAR REGRESSION OF THE INPUT FILE======
 call calc_lin_reg(input_file,0.d0,3.d0)

contains

!========================
function fit_func_1_gauss(var_x,par_a,par_alpha,par_b) result(res)
 implicit none
 real*8   :: res,var_x, par_a,par_alpha,par_b
 res=par_a*EXP(-par_alpha*var_x**2.d0)+par_b
end function fit_func_1_gauss

!========================
function fit_func_5_gauss(ratom,natom,rpoint,par_a,par_alpha,par_b) result(res)
 implicit none
 integer  :: natom
 real*8   :: res,par_a,par_alpha,par_b
 real*8   :: ratom(natom,2)
 real*8   :: rpoint(2)
!=====
 integer  :: iatom
!===== 

 res=0.d0 

 do iatom=1,natom
   res=res+EXP(-par_alpha*NORM2(rpoint(:)-ratom(iatom,:))**2.d0)
 end do
 res=res*par_a
 res=res+par_b

end function fit_func_5_gauss

!=======================
function integrate_gaussian(par_alpha,alat,lim) result(res)
  implicit none
  real*8  :: par_alpha,res,alat,lim
!=====
  real*8  :: var_x,dx
!=====
  dx=0.001
  res=0.d0
  var_x=0.d0
  do while(var_x<=lim)
    res=res+EXP(-par_alpha*var_x**2.d0)*dx
    var_x=var_x+dx
  end do
end function

!========================
function integrate_5_gauss(ratom,natom,par_a,par_alpha,par_b,alat) result(res)
 implicit none
 integer  :: natom
 real*8   :: res,par_a,par_alpha,par_b,alat
 real*8   :: ratom(natom,2)
!=====
 real*8   :: var_x,var_y,dx,dy,xmin,xmax,ymin,ymax
 integer  :: iatom
!===== 

 dx=alat/2000.d0
 dy=alat/2000.d0

 xmin=0.d0
 ymin=0.d0
 xmax=alat/2.d0
 ymax=alat/2.d0

 res=0.d0
 var_x=xmin
 do while (var_x<xmax)
   var_y=ymin
   do while (var_y<ymax)
     res=res+fit_func_5_gauss(ratom,natom,(/ var_x, var_y /),par_a,par_alpha,par_b)*dx*dy
     var_y=var_y+dy
   end do
   var_x=var_x+dx
 end do

 res=4.d0*res/(alat**2.d0)

end function integrate_5_gauss

!========================
function error_func_1_gauss(nrow,par_a,par_alpha,par_b,points_arg,weight,points_val) result(res)
 implicit none
 integer,intent(in)      :: nrow
 real*8,intent(in)       :: par_a,par_alpha,par_b
 real*8, intent(in)      :: points_arg(nrow),points_val(nrow),weight(nrow)
!=====
 integer  :: irow
 real*8   :: res
!=====
 res=0.d0
 do irow=1,nrow
   res=res+weight(irow)*(points_val(irow)-fit_func_1_gauss(points_arg(irow),par_a,par_alpha,par_b))**2.d0
 end do
end function error_func_1_gauss

!========================
function error_func_5_gauss(nrow,ratom,natom,par_a,par_alpha,par_b,points_y,points_val) result(res)
 implicit none
 integer,intent(in)      :: nrow,natom
 real*8,intent(in)       :: par_a,par_alpha,par_b
 real*8, intent(in)      :: points_y(nrow),points_val(nrow)
 real*8,intent(in)       :: ratom(natom,2)
!=====
 integer  :: irow
 real*8   :: res
 real*8   :: rpoints(nrow,2)
!=====

 rpoints(:,2)=points_y(:)
 rpoints(:,1)=0.d0

 res=0.d0
 do irow=1,nrow
   res=res+(points_val(irow)-fit_func_5_gauss(ratom,natom,rpoints(irow,:),par_a,par_alpha,par_b))**2.d0
 end do
end function error_func_5_gauss

!======================
subroutine calc_gradient_1_gauss(nrow,par_a,par_alpha,par_b,points_arg,points_val,weight,gradient) 
 implicit none
 integer,intent(in)      :: nrow
 real*8,intent(in)       :: par_a,par_alpha,par_b
 real*8, intent(in)      :: points_arg(nrow),points_val(nrow),weight(nrow)
 real*8,intent(inout)    :: gradient(3)
 
!=====
 integer  :: irow
 real*8   :: common_part
 
 gradient(:)=0.d0
 common_part=0.d0
 do irow=1,nrow
   common_part=2.d0*weight(irow)*(points_val(irow)-par_a*EXP(-par_alpha*points_arg(irow)**2.d0)-par_b)/nrow
   gradient(1)=gradient(1)+common_part*(-EXP(-par_alpha*points_arg(irow)**2.d0))
   gradient(2)=gradient(2)+common_part*par_a*points_arg(irow)**2.d0*EXP(-par_alpha*points_arg(irow)**2.d0)
   gradient(3)=gradient(3)+common_part*(-1)
 end do

end subroutine calc_gradient_1_gauss

!======================
subroutine calc_gradient_5_gauss(nrow,ratom,natom,par_a,par_alpha,par_b,points_y,points_val,gradient) 
 implicit none
 integer,intent(in)      :: nrow,natom
 real*8,intent(in)       :: par_a,par_alpha,par_b
 real*8,intent(in)       :: points_y(nrow),points_val(nrow)
 real*8,intent(in)       :: ratom(natom,2)
 real*8,intent(inout)    :: gradient(3)
!=====
 integer  :: irow,iatom
 real*8   :: common_part,term1,term2,tmp
 real*8   :: rpoints(nrow,2)

 rpoints(:,2)=points_y(:)
 rpoints(:,1)=0.d0


 gradient(:)=0.d0
 common_part=0.d0

 do irow=1,nrow
   term1=0.d0
   term2=0.d0
   common_part=0.d0
   do iatom=1,natom
     tmp=EXP(-par_alpha*NORM2(rpoints(irow,:)-ratom(iatom,:))**2.d0)
     term1=term1-tmp
     term2=term2+par_a*NORM2(rpoints(irow,:)-ratom(iatom,:))**2.d0*tmp
     common_part=common_part+par_a*tmp
   end do
   common_part=points_val(irow)-common_part-par_b
   gradient(1)=gradient(1)+common_part*term1
   gradient(2)=gradient(2)+common_part*term2
   gradient(3)=gradient(3)+common_part*(-1)

 end do
 
 gradient(:)=gradient(:)*2.d0/nrow

end subroutine calc_gradient_5_gauss

end program gradient_descent_exp

