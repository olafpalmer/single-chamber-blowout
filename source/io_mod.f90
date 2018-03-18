module io_manage

    type input_type
        real*8   :: dt, p_int, p_ext, t_int, t_ext
        real*8   :: area, comp_vol, cd, time_end
        integer  :: ncomp
    end type input_type
    
    real*8, parameter  :: T0 = 273.15
    
    contains
        
        function get_inputs(input_file) result(inp)
            
            character*60, intent(in)  :: input_file
            character*60              :: line, aux
            integer              :: stat
            type(input_type)     :: inp        
            
            open(unit=10,file=input_file,status='unknown')
    
            do 
                
                read(10,'(A)',iostat=stat) line            
                if (is_iostat_end(stat)) exit
                
                select case (line(1:8)) 
                 
                    case('TIMESTEP')   ! Time Step
                        read(line,*) line(1:8), inp%dt, aux
                    case('EXTPRESS')   ! External Pressure
                        read(line,*) line(1:8), inp%p_ext, aux
                    case('INTPRESS')   ! Internal Pressure
                        read(line,*) line(1:8), inp%p_int, aux
                    case('EXT_TEMP')   ! External Temperature
                        read(line,*) line(1:8), inp%t_ext, aux
                        inp%t_ext = inp%t_ext + T0
                    case('INT_TEMP')   ! Internal Temperature
                        read(line,*) line(1:8), inp%t_int, aux
                        inp%t_int = inp%t_int + T0
                    case('COMP_VOL')   ! Compartment volume
                        read(line,*) line(1:8), inp%comp_vol, aux        
                    case('OPENAREA')   
                        read(line,*) line(1:8), inp%area, aux
                    case('DISCOEFF')   ! Discharge coefficient
                        read(line,*) line(1:8), inp%cd, aux     
						case('TIME_END')   ! Time evaluation
                        read(line,*) line(1:8), inp%time_end, aux   
                              
                endselect 
                
            enddo
            close(10)
            
        endfunction get_inputs
    
        !subroutine init_parameters(rho,temp,press_in,press_out)
        
        !endsubroutine
        
        subroutine write_output(dt,n_steps,P,T,rho)
        
            integer, intent(in)   :: n_steps
            integer               :: i
            real*8, intent (in)   :: dt             
            real*8, intent (in), allocatable   ::  P(:), T(:), rho(:)
            character*60          :: file_output
        
            file_output = "single_chamber_results.txt"
            open(unit=11,file=file_output,status='unknown')
            
            write(11,*)  "DECOMPRESSION HISTORY"
            write(11,*)  "TIME (S)    PRESSURE (Pa)   TEMPERATURE (Celsius)   DENSITY (kg/m^3)"
            
            do i = 1,n_steps
            
                write(11,'(F6.3,F16.1,F19.1,F22.3)') (i-1)*dt, P(i), T(i) - T0, rho(i) 
            
            enddo 
            close(11)
            
        endsubroutine write_output
        
end module io_manage
