program main 
   implicit none
   
   double precision , allocatable :: x(:) , si(:) , si1(:) 
   integer :: num_of_intervals , j ,count , points,i,value
   double precision :: x_upper , x_lower , h , energy, normalizer , c2 , d , c_s , sum_all,sum_si=0 ,x1,y1,t,r
   double precision , dimension(10) :: prob
   real , parameter :: pi = 3.141592653589793
   logical :: f1 
   
   d = 2.0
   
   print*, "Enter the value of the energy : " 
   read*, energy
   
   
   c2 = d*d*energy/4.0
   c_s = 0.6043499*c2 - 0.006188*c2*c2 - 0.0000589*c2*c2*c2 +2 
   !c_s = 0.3127477*c2 - 0.0231669*c2*c2 - 0.0005110*(c2**3) - 0.0000045*(c2**4) 
  
   inquire(file="points.dat",exist=f1) ! opening the file for the store the data for plot 
   if(f1) then
   open(1,file="points.dat",status="replace")
   else
   open(1,file="points.dat",status="new",action="write")
   endif
   
   num_of_intervals = 9000
   x_upper = 10
   x_lower = 1 + 1.0e-4
   h = (x_upper - x_lower)/real(num_of_intervals)
   
   allocate(x(num_of_intervals+1),si(num_of_intervals+1),si1(num_of_intervals+1))
   
   x(1) = x_lower
   si(1) = 1
   si1(1) = -(2*d + c2 - c_s)/2.0
   normalizer = 0 

   do j = 2,num_of_intervals
       call rkmethod(j) 
       normalizer = normalizer + 0.5*h*(si(j)*si(j) + si(j-1)*si(j-1))
   enddo
   
   normalizer = sqrt(2*normalizer)
   count = 0
   do j = 1, num_of_intervals
       if(j == 1000*count+1) then
       		prob(count+1) = si(j)
       		count = count+1
       	endif
       	sum_all = sum_all + si(j)
   enddo
   do j = 1,5
   	points = prob(j)*100
   	do i = 1, 1000
   		if (points < i) then 
	   		call random_number(t)
	   		call random_number(r)
	   		r = r*real(j)
	   		x1 = r*cos(t*2*pi) + 1
	   		y1 = r*sin(2*pi*t) 
	   		write(1,*)x1,y1
	   		x1 = r*cos(t*2*pi) - 1
	   		y1 = r*sin(2*pi*t) 
	   		write(1,*)x1,y1
   		endif
   	enddo
   	write(1,*)""
   enddo
   
   deallocate(x,si,si1)
   
   close(1)
   
   	  
   print*,'simulation is end'
   stop
   
   contains
   	    	
   	double precision function func(val_x,val_y,val_y1) ! here the ODE we wanted to solve 
            double precision , intent(in) :: val_x,val_y,val_y1
            double precision :: result , mid_1 , mid_2 , val
            val = (val_x)
            mid_1 = ((2*d*val+ c2*val*val - c_s)*val_y)
            mid_2 = 2*val*val_y1
            result = (-mid_1 - mid_2)/(val*val - 1)
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
   
