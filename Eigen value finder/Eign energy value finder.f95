program main 
   implicit none

   ! here we are only find the energy eigenvalue by using the shooting method and for this we 
   ! we go from a range of energy and solve the radial equation for all the value of the energy and 
   ! when the solution of radial equation behave like we wanted then we print this energy value and 
   ! use the value to getting a very good radial and angular solution and also the electron cloud 
   
    double precision , allocatable :: x(:) , si(:) , si1(:) 
    integer :: num_of_intervals , j , k 
    double precision :: x_upper , x_lower , h , energy, step_size  , normalizer , c2 , d , c_s
    double precision :: energy_upper , energy_step
    logical :: f1 , zero = .false. 

    print*, "Enter the distance "  ! here we getting the inter nuclear distance 
    read*, d
    
    print*, "Enter the value of the energy : "  ! upper limit of the energy 
    read*, energy_upper 
                                        ! here we divide energy into 200 steps ,u can increase it 
    energy_step = energy_upper/200    ! if u have a very good computer
    
    inquire(file="points.dat",exist=f1) ! opening the file for the store the data for plot 
    if(f1) then
    open(1,file="points.dat",status="replace")
    else
    open(1,file="points.dat",status="new",action="write")
    endif
   
    do k = 1,200  ! this is our main loop for the energy eigenvalue
   
        energy = -energy_step*k ! calculating the energy c2 and c_s
        
        c2 = d*d*energy/4.0
        c_s = 0.6043499*c2 - 0.006188*c2*c2 - 0.0000589*c2*c2*c2 +2 
        !c_s = 0.3127477*c2 - 0.0231669*c2*c2 - 0.0005110*(c2**3) - 0.0000045*(c2**4) 

        num_of_intervals = 10000 ! it for the radial equation 
        x_upper = 10
        x_lower = 1 + 1.0e-7
        h = (x_upper - x_lower)/real(num_of_intervals)
        
        allocate(x(num_of_intervals+1),si(num_of_intervals+1),si1(num_of_intervals+1))
        
        x(1) = x_lower  ! all the initial values and boundary condition
        si(1) = 1
        si1(1) = -(2*d + c2 - c_s)/2.0
        normalizer = 0   ! here we also normalize the function 
        do j = 2,num_of_intervals
            call rkmethod(j)       ! here we calculate the radial function for given energy
            normalizer = normalizer + 0.5*h*(si(j)*si(j) + si(j-1)*si(j-1))
        enddo    ! we use Trapezoidal method for integration 
        
        normalizer = sqrt(normalizer)
        do j = 1, num_of_intervals
            si(j) = si(j)/normalizer   ! in this condition we calculate that the function is 
            if((abs(si(j)) < 1.0e-4 ) .and. (abs(si(j-10)) < 1.0e-4)) zero=.true. 
            write(1,*)x(j),si(j),si1(j),si(j)**2 , energy  ! dying or not when go away 
        enddo                          ! if the function is dying then we get the energy value 
        write(1,*)""                   ! we print it's energy and if it's not then repeat for another value
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
   	
   	    double precision function func(val_x,val_y,val_y1) ! here the radial equation 
            double precision , intent(in) :: val_x,val_y,val_y1
            double precision :: result , mid_1 , mid_2
            mid_1 = ((2*d*val_x + c2*val_x*val_x - c_s)*val_y)
            mid_2 = 2*val_x*val_y1
            result = (-mid_1 - mid_2)/(val_x*val_x - 1)
            func = result 
            return 
        end function func
 
        ! here we calculate the value by the rk method and this is for 2nd order equation

        subroutine rkmethod(i)    ! here we use 4th order rk method for 2nd order ODE
            integer , intent(in) :: i
            double precision :: k1,k2,k3,k4,k11,k12,k13,k14,result  
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