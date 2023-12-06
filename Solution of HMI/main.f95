program main 
   implicit none !declaring the variable for the program 
   
   double precision , allocatable :: x(:) , si(:) , si1(:) ! here we use for both radial and angular part 
   integer :: num_of_intervals , j                         ! by just reallocate the array
   double precision :: x_upper , x_lower , h , energy, normalizer , c2 , d , c_s 
   logical :: f1 , f2 , angular = .false.
   
   d = 2.0 ! it is the inter nuclear distance 
   
   print*, "Enter the value of the energy : " 
   read*, energy
   
   num_of_intervals = 10000
   
   c2 = d*d*energy/4.0  ! c_s can be taken by both value it depend upon the what type of system it is.
   c_s = 0.6043499*c2 - 0.006188*c2*c2 - 0.0000589*c2*c2*c2 +2
   !c_s = 0.3127477*c2 - 0.0231669*c2*c2 - 0.0005110*(c2**3) - 0.0000045*(c2**4) 
   
   ! for radial part solution 
   
   inquire(file="radial.dat",exist=f1) ! opening the file for the store the data for plot 
   if(f1) then
   	open(1,file="radial.dat",status="replace")
   else
   	open(1,file="radial.dat",status="new",action="write")
   endif
   

   x_upper = 10         ! This is the upper limit we can extent it also for higher energies
   x_lower = 1 + 1.0e-7 ! we can't use the 1 because the function is not defined at 1
   h = (x_upper - x_lower)/real(num_of_intervals)
   
   allocate(x(num_of_intervals+1),si(num_of_intervals+1),si1(num_of_intervals+1))
   
   ! here x = e , si = R(e) , si1 = R'(e)

   x(1) = x_lower
   si(1) = 1
   si1(1) = -(2*d + c2 - c_s)/2.0 ! this is the boundary condition we derive with m=0
   normalizer = 0  ! normalizer to use to normalize the function it use Trapezoidal method integration

   do j = 2,num_of_intervals
       call rkmethod(j)       ! here we call the rk method to solve for next value 
       normalizer = normalizer + 0.5*h*(si(j)*si(j) + si(j-1)*si(j-1))
   enddo
   
   normalizer = sqrt(2*normalizer) ! we multiply it with 2 because we find only one side of the radial fxn
                                   ! And for we getting exact value when we do for both side 
  
   do j = 1, num_of_intervals
       si(j) = si(j)/normalizer
       write(1,*)x(j),si(j),si1(j),si(j)**2,si(j)*x(j)
   enddo
   
   deallocate(x,si,si1) ! here we deallocate the array for use in angular part 
   
   close(1)
   
   ! for Angular wave function solution 
   
   angular = .true.  ! this change the func into angular equation
   
   inquire(file='angular.dat',exist=f2)
   if(f2) then 
   	open(2,file="angular.dat",status="replace")
   else
   	open(2,file="angular.dat",status="new",action="write")
   endif
   ! x = n , si = S(n) , si1 = S'(n)
   allocate(x(num_of_intervals+1),si(num_of_intervals+1),si1(num_of_intervals+1))
   
   x_upper = 0    ! here we done for half of the and generate the other value using these values
   x_lower = -1 + 1.0e-4  
   h = (x_upper - x_lower)/real(num_of_intervals)
   
   x(1) = x_lower
   si(1) = 1
   si1(1) = (c2-c_s)/2  ! here m = 0 and we get this equation in Boundary condition section 
   normalizer = 0       ! normalizer is same as the radial part but we not using it.
 
   do j = 2,num_of_intervals
       call rkmethod(j) 
       normalizer = normalizer + 0.5*h*(si(j)*si(j) + si(j-1)*si(j-1))
   enddo
   
   normalizer = sqrt(2*normalizer)
   do j = 1,num_of_intervals       ! this gives value between -1 to 0 
      !si(j) = si(j)/normalizer
       write(2,*)x(j),si(j),si1(j)/normalizer,si(j)/x(j)
   enddo
  do j = num_of_intervals,1,-1    ! here we generate the value b/w 0 to 1 
      !si(j) = si(j)/normalizer   ! using the other values 
      write(2,*)-x(j),si(j),si1(j)/normalizer,si(j)/x(j)
   enddo
   
   deallocate(x,si,si1) ! here we finally deallocate the memory
   close(2)
   	  
   print*,'simulation is end'
   stop
   
   contains  ! func is the wave equations

   	    double precision function func(val_x,val_y,val_y1)  
            double precision , intent(in) :: val_x,val_y,val_y1
            double precision :: result , mid_1 , mid_2 
            if(angular) then    
            	mid_1 = (c2*val_x*val_x - c_s)*val_y ! angular function 
                mid_2 = (2*val_x*val_y1)
                result = (mid_1 + mid_2)/(1 - val_x*val_x)
                else
                mid_1 = ((2*d*val_x + c2*val_x*val_x - c_s)*val_y)
                mid_2 = 2*val_x*val_y1                         ! radial function 
                result = (-mid_1 - mid_2)/(val_x*val_x - 1)
                endif
                    func = result
            return 
        end function func
                                ! here we calculate the value by the rk method and this is for 2nd order equation
        subroutine rkmethod(i)    
            integer , intent(in) :: i
            double precision :: k1,k2,k3,k4,k11,k12,k13,k14,result  
            k11 = h*si1(i-1)                    
            k1 = h*func(x(i-1),si(i-1),si1(i-1))  ! here we use 4th order rk method for 2nd order ODE
            k12 = (si1(i-1) + k1/2.0)*h
            k2 = h*func((x(i-1)+h/2.0),(si(i-1)+k11/2.0),(si1(i-1)+k1/2.0))
            k13 = (si1(i-1) + k2/2.0)*h
            k3 = h*func((x(i-1)+h/2.0),(si(i-1)+k13/2.0),(si1(i-1)+k2/2.0))
            k14 = (si1(i-1) + k3)*h
            k4 = h*func((x(i-1)+h),(si(i-1)+k14),(si1(i-1)+k3))
            x(i) = x(i-1) + h
            si(i) = si(i-1) + (k11+2*k12+2*k13+k14)/6.0
            si1(i) = si1(i-1) + (k1+2*k2+2*k3+k4)/6.0
            return 
        end subroutine rkmethod
end program 
   
