!-----------------------------------------------------------------------
! Blowout Program
! By Olaf Palmer Val Pinheiro 
! Part 01 - Single-Chamber Decompression
! 
!-----------------------------------------------------------------------

program blowout_calculation

! Loading modules
use io_manage
use single_chamber

implicit none

! Variable declaration
!-----------------------------------------------------------------------

character*60            :: file_input
type(input_type)        :: inp

integer               ::  time_steps, i
real*8                ::  rho_ext, press, temp      ! initial conditions
real*8, allocatable   ::  P(:,:), rho(:,:), T(:,:)  ! state parameters
!-----------------------------------------------------------------------


write(*,*) '==========================================='
write(*,*) '          SINGLE CHAMBER DECOMPRESSION'
write(*,*) '==========================================='
write(*,*)

!-----------------------------------------------------------------------
! 1) READ INPUT FILE AND PRE PROCESSING
!-----------------------------------------------------------------------
write(*,*) '1) Reading inputs and initializing state variables...'
file_input = "blowout_input.txt"
inp = get_inputs(file_input)

time_steps = floor(inp%time_end/inp%dt) + 1

allocate(P(time_steps,3))
allocate(T(time_steps,3))
allocate(rho(time_steps,3))

! Initializing
P(1) = inp%p_int 
T(1) = inp%t_int 
rho(1) = P(1) / (R * T(1))
rho_ext = inp%p_ext / (R * inp%t_ext)

!-----------------------------------------------------------------------
! 2) DECOMPRESSION CALCULATION
!-----------------------------------------------------------------------
write(*,*) '2) Decompression Calculation '
do i = 2,time_steps
    
    rho(i) = rk4_step(inp%dt,inp%area,inp%comp_vol,inp%Cd, T(i-1), P(i-1),inp%p_ext,rho(i-1),rho_ext, rho(1), P(1), T(1))   

    ! Updating the parameters     
    call update_param(P(i),T(i),rho(i),rho(1), P(1), T(1))  
    
    write(*,'(A,I6,A,I6,F10.1)') 'Iteration: ', i-1, ' of ', time_steps-1, P(i-1) 
end do    

!-----------------------------------------------------------------------
! 3) OUTPUT PLOTTING
!-----------------------------------------------------------------------
write(*,*) '3) Writing Outputs '

call write_output(inp%dt,time_steps, P, T, rho)
read(*,*)

!-----------------------------------------------------------------------


deallocate(P)
deallocate(T)
deallocate(rho)


end program blowout_calculation
