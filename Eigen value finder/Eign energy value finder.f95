program main 
   implicit none
   
   double precision , allocatable :: x(:) , si(:) , si1(:) , realsi1(:)
   integer :: num_of_intervals , j , last_point , k 
   double precision :: x_upper , x_lower , h , energy, step_size , energy_vibrator , normalizer , pot_energy , max_si , c2 , d , c_s
   double precision :: energy_upper , energy_step
   logical :: f1 , peak = .false. , zero = .false. , angular=.false.
   print*, "Enter the distance " 
   read*, d
   
   print*, "Enter the value of the energy : " 
   read*, energy_upper 
   
   energy_step = energy_upper/200
   
   inquire(file="points.dat",exist=f1) ! opening the file for the store the data for plot 
   if(f1) then
   open(1,file="points.dat",status="replace")
   else
   open(1,file="points.dat",status="new",action="write")
   endif
   
   do k = 1,200
   
   energy = -energy_step*k
   
   
   c2 = d*d*energy/4.0
   c_s = 0.6043499*c2 - 0.006188*c2*c2 - 0.0000589*c2*c2*c2 +2 
   !c_s = 0.3127477*c2 - 0.0231669*c2*c2 - 0.0005110*(c2**3) - 0.0000045*(c2**4) 
   print*, c_s
   
   
   num_of_intervals = 10000
   x_upper = 10
   x_lower = 1 + 1.0e-7
   h = (x_upper - x_lower)/real(num_of_intervals)
   
   allocate(x(num_of_intervals+1),si(num_of_intervals+1),si1(num_of_intervals+1))
   
   x(1) = x_lower
   si(1) = 1
   si1(1) = -(2*d + c2 - c_s)/2.0
   normalizer = 0 
   last_point = num_of_intervals
   do j = 2,num_of_intervals
       call rkmethod(j) 
       normalizer = normalizer + 0.5*h*(si(j)*si(j) + si(j-1)*si(j-1))
       
   enddo
   
   normalizer = sqrt(normalizer)
   print*, normalizer
   do j = 1, num_of_intervals
       si(j) = si(j)/normalizer
       if((abs(si(j)) < 1.0e-4 ) .and. (abs(si(j-10)) < 1.0e-4)) zero=.true. 
       write(1,*)x(j),si(j),si1(j),si(j)**2 , energy
       
   enddo
   write(1,*)""
   if (zero) then 
   	print*, "The solution of the equation is the founded at " , energy , " and at ", k , "point"
   	deallocate(x,si,si1)
   	stop
   endif
   deallocate(x,si,si1)
   
   enddo
   
   close(1)
   	
   print*,'simulation is end'
   stop
   
   contains
   	subroutine potential_energy_cal(po_energy,val_x)
   	    double precision , intent(out) :: po_energy 
   	    double precision, intent(in) :: val_x
   	    po_energy = (0.765625e-2)*(val_x**2)
   	    return 
   	end subroutine potential_energy_cal
   	
   	double precision function func(val_x,val_y,val_y1) ! here the ODE we wanted to solve 
            double precision , intent(in) :: val_x,val_y,val_y1
            double precision :: result , mid_1 , mid_2
            if(angular) then
            	mid_1 = (c2*val_x*val_x - c_s)*val_y
		mid_2 = (2*val_x*val_y1)
		result = (mid_1 + mid_2)/(1 - val_x*val_x)
	    else
		 mid_1 = ((2*d*val_x + c2*val_x*val_x - c_s)*val_y)
		 mid_2 = 2*val_x*val_y1
		 result = (-mid_1 - mid_2)/(val_x*val_x - 1)
	    endif
            !mid_1 = ((2*d*val_x + c2*val_x*val_x - c_s)*val_y)
            !mid_2 = 2*val_x*val_y1
            !result = (-mid_1 - mid_2)/(val_x*val_x - 1)
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
   
