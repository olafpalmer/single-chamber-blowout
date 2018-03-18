module single_chamber

real*8, parameter   ::  gamma = 1.4, R = 287.053

contains
    
    !------------------------------------------------------------------
    ! Update Pressure and Temperature 
    !------------------------------------------------------------------
    subroutine update_param(press, temp, rho, rho_0, press_0, temp_0)
    
        real*8, intent(out)  :: press, temp
        real*8, intent(in)  :: rho, rho_0, press_0, temp_0
        
        temp = temp_0 * (rho / rho_0)**(gamma - 1) 
        press = press_0 * (rho / rho_0)**(gamma)
    
    endsubroutine update_param
    
    !------------------------------------------------------------------
    ! Density derivative evaluation 
    !------------------------------------------------------------------
    function drho_dt(Area,Vi,Cd,Ti,Pi,Pj,rhoi,rhoj) 
    
        real*8, intent(in)  :: Area, Vi, Ti, Pi, Pj, rhoi, rhoj, Cd
        real*8              :: drho_dt
        real*8  :: P_chocked, Aeff, mij_dot
        
        P_chocked = Pj*((gamma + 1)/2)**(gamma/(gamma - 1))
        Aeff = Area*Cd
        
        if (Pi >= P_chocked) then
        
            mij_dot = rhoi*((2/(gamma+1))**(1/(gamma-1)))*Aeff*sqrt(2*gamma*R*Ti/(gamma + 1))
            
        elseif (Pj > Pi) then
        
            mij_dot = -Aeff*sqrt(2*Pj*rhoj*(gamma/(gamma - 1))*((Pi/Pj)**(2/gamma))*(1 - (Pi/Pj)**((gamma - 1)/gamma)))
            
        else
        
            mij_dot = Aeff*sqrt(2*Pi*rhoi*(gamma/(gamma - 1))*((Pj/Pi)**(2/gamma))*(1 - (Pj/Pi)**((gamma - 1)/gamma)))
        
        endif 
    
        drho_dt = (1/Vi)*(-mij_dot)
        
    endfunction drho_dt
    
    !------------------------------------------------------------------
    ! Time integration Scheme - 4th order Runge Kutta
    !------------------------------------------------------------------
    function rk4_step(h,Area,Vi,Cd,Ti,Pi,Pj,rhoi,rhoj,rho0,P0,T0)
    
        real*8, intent(in)  ::  Area, Vi, rhoi, rhoj, h, Cd, rho0,P0,T0, Pi, Pj, Ti
        real*8              ::  rk4_step, press, temp
        real*8              ::  k(4)
        
        press = Pi
        temp = Ti
        
        k(1) = drho_dt(Area,Vi,Cd,temp,press,Pj,rhoi,rhoj)
                
        call update_param(press,temp,rhoi + 0.5*h*k(1),rho0, P0, T0)
        k(2) = drho_dt(Area,Vi,Cd,temp,press,Pj,rhoi + 0.5*h*k(1),rhoj)        
                
        call update_param(press,temp,rhoi + 0.5*h*k(2),rho0, P0, T0)
        k(3) = drho_dt(Area,Vi,Cd,temp,press,Pj,rhoi + 0.5*h*k(2),rhoj) 
               
        call update_param(press,temp,rhoi + h*k(3),rho0, P0, T0)
        k(4) = drho_dt(Area,Vi,Cd,temp,press,Pj,rhoi + h*k(3),rhoj)
    
        rk4_step = rhoi + h*(k(1) + 2*k(2) + 2*k(3) + k(4))/6        

         
            
    end function rk4_step

end module single_chamber
