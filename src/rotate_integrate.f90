module rotate_integrate
 use lin_reg_fort
contains

!=========================================================================
subroutine calculate_integral(input_file,alat,aver)
 implicit none
!=====
 character(len=100),intent(in) :: input_file
 real*8,intent(in)             :: alat
 real*8,intent(inout)          :: aver
!=====
 integer            :: iray,nray,nrow,irow,unit_file
 logical            :: file_exists
 real*8             :: area,dphi,phi,stop_inter,s0
 real*8,allocatable :: coor(:),stopping(:)
 real*8,parameter   :: pi    = 3.14159265358979323d0

 nray=10000
 ! alat=3.49d0
 area=alat**2/8.d0
 dphi=pi/2.d0/nray

 inquire(file=input_file,exist=file_exists)
 if(.NOT. file_exists) then
   write(*,'(a)') 'STOP: input file not found'
   stop
 endif

 nrow=get_number_of_rows(input_file)+1
 allocate(stopping(nrow))
 allocate(coor(nrow)) 

 open(newunit=unit_file, file = input_file)
 do irow=2,nrow
   read(unit_file,*) coor(irow),stopping(irow)
 end do
 close(unit_file)

 coor(1)=0.d0
 stopping(1)=stopping(nrow)

 s0=stopping(nrow)
 
 do irow=1,nrow
   stopping(irow)=stopping(irow)-stopping(nrow)
 end do

 do iray=0,nray-1
   phi=iray*dphi
   do irow=2,nrow
     if(coor(irow)<=rmax(phi,alat)) then
       aver=aver+(stopping(irow)+stopping(irow-1))/2.d0*(coor(irow)-coor(irow-1))*(coor(irow)+coor(irow-1))/2.d0
     else
       stop_inter=(rmax(phi,alat)-coor(irow-1))*(stopping(irow)-stopping(irow-1))/(coor(irow)-coor(irow-1))+stopping(irow-1)
       aver=aver+(stop_inter+stopping(irow-1))/2.d0*(rmax(phi,alat)-coor(irow-1))*(coor(irow)+coor(irow-1))/2.d0
       exit
     end if
   end do
 end do

 aver=aver*dphi/area+s0

 write(*,*) aver

end subroutine calculate_integral
 
!=======================
function rmax(phi,alat)
 implicit none
 real*8   :: phi,alat,rmax
!===
 rmax=alat/2.d0/(COS(phi)+SIN(phi)) 
end function rmax

end module rotate_integrate
!=========================================================================
