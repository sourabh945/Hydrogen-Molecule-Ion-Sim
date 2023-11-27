program main 
   implicit none
   
   double precision , allocatable :: x(:) , si(:) , si1(:) , potential_energy(:)
   integer :: num_of_intervals , j , last_point
   double precision :: x_upper , x_lower , h , energy , normalizer , c2 , d , c_s
  
   logical :: f1 , peak = .false. 
   
   print*, "Enter the value of the energy : " 
   read*, energy
   print*, "Enter the value of D: "
   read*, d
   print*, "Enter the value of the c_s : "
   read*, c_s
   
   c2 = energy*d*d/4.0
   
   inquire(file="angularpoints.dat",exist=f1) ! opening the file for the store the data for plot 
   if(f1) then
   open(1,file="angularpoints.dat",status="replace")
   else
   open(1,file="angularpoints.dat",status="new",action="write")
   endif
   
   num_of_intervals = 5000
   x_upper = 0
   x_lower = -1 + 1.0e-7
   h = (x_upper - x_lower)/real(num_of_intervals)
   
   allocate(x(num_of_intervals+1),si(num_of_intervals+1),si1(num_of_intervals+1))
   
   x(1) = x_lower + 1
   si(1) = -1
   si1(1) = -1
   normalizer = 0 
   last_point = num_of_intervals
   do j = 2,num_of_intervals
       call rkmethod(j) 
       normalizer = normalizer + 0.5*h*(si(j)*si(j) + si(j-1)*si(j-1))
   enddo
   
   normalizer = sqrt(normalizer)
   do j = 1, num_of_intervals
       !si(j) = si(j)/normalizer
       write(1,*)x(j)-1,si(j),si1(j)/normalizer,x(j)
   enddo
    do j = num_of_intervals,1,-1
       !si(j) = si(j)/normalizer
       write(1,*)-1*x(j)+1,-si(j),si1(j)/normalizer,-x(j)
   enddo
   close(1)
   	
   print*,'simulation is end'
   stop
   
   contains
   	
   	double precision function func(val_x,val_y,val_y1) ! here the ODE we wanted to solve 
            double precision , intent(in) :: val_x,val_y,val_y1
            double precision :: result , mid_1 , mid_2
            mid_1 = ((c2*val_x*((val_x -1)**2)*(val_x -2) - c_s*val_x*(val_x-2))*val_y )
            mid_2 = ((2*val_x*(val_x-1)*(val_x-2))*val_y1)
            result = (-mid_1 - mid_2)/((val_x**2)*((val_x - 2)**2))
            func = result 
            return 
        end function func

        subroutine rkmethod(i)   ! here we calculate the value by the rk method and this is for 2nd order equation 
            integer , intent(in) :: i
            double precision :: k1,k2,k3,k4,k11,k12,k13,k14,result  ! here we use 4th order rk method for 2nd order ODE
            k11 = h*si1(i-1)
            k1 = h*func(x(i-1),si(i-1),si1(i-1))
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
   
