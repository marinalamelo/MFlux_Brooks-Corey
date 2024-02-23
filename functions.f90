module functions

contains

real(kind=8) function f(Ta)
    implicit none
    real(kind=8), intent(in) :: Ta
    real pi
    real a, r0, rx, Kroot, hw, hb, Ks, RLD, z, Ll, raiz3
    real rm, lambda, Ms
    real hs, hl, h0, hx, rho, phi, p, LAMBDA31, AA, BB, PP
    real Mx, M0, Ta0, Taf, S
    character*50 AuxFile
    
        AuxFile = "auxi.qvl"
        open (unit = 2, file = AuxFile , status ='old')
         read (2,*) RLD, a, r0, rx, Kroot, Ll, hw, hl
         read (2,*) hb, p, Ks, hs, z
        close (2) 
        pi = 3.14159265358979                                                                     !adm     constant
        rm = 0.01*sqrt(1.0/pi/RLD)                                                                !m       calculated (plant)
        lambda = (p-1)/3                                                                          !adm     calculated (soil)
        raiz3 = 2.0**(1.0/3.0)                                                                    !adm     calculated (soil)
        LAMBDA31 = -Ks*hb/p                                                                       !m²d-1   calculated (soil)
        Ms = LAMBDA31*((hb/hs)**p-(hb/hw)**p)                                                     !m²d-1   calculated (soil)
        rho = 4/(r0**2.0-(a*rm)**2+2*(rm**2.0+r0**2)*LOG(a*rm/r0))                                !m-2     calculated (plant)
        phi = rho*rm**2*LOG(r0/rx)/2.0/Kroot                                                      !d m-1   calculated (plant)
        AA = phi*LAMBDA31*(-hb)**p                                                                !m^(1+p) calculated (plant)
        Ta0 = 4.37665046487299E-03                                                                !m d-1   constant   (plant)
        hx = hl+Ta0/Ll                                                                            !m       calculated (plant)
        Mx = Kroot*(hx-hw)                                                                        !m²d-1   calculated (plant)
        BB = hx+phi*Ms+AA/((-hw)**p)                                                              !m       calculated (plant)
        PP =(3.0*sqrt(3.0)*sqrt(27.0*AA**2.0-4.0*AA*BB**3.0)+27.0*AA-2.0*BB**3.0)**(1.0/3.0)      !adm     calculated (plant)
        h0 = (BB-PP/raiz3-raiz3*BB**2.0/PP)/3.0                                                   !m       calculated (plant)                                                
        M0 = -Ks*hb/p*((hb/h0)**p-(hb/hw)**p)                                                     !m²d-1   calculated (plant)
        S = rho*(Ms-M0)                                                                           !d-1     calculated (plant)
        Taf = S*z                                                                                 !m d-1   calculated (plant)
    
        f = abs(Taf-Ta)

end function f

real(kind=8) function fprime(Ta)
    implicit none
    real(kind=8), intent(in) :: Ta
    
    fprime = 1.d0 * Ta**(1-1)

end function fprime

end module functions