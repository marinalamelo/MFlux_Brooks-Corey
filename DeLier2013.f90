program DeLier2013

    use newton, only: solve
    use functions, only: f, fprime

    implicit none
    real(kind=8) :: Ta, Ta0, fTa
    integer :: iters
	logical :: debug         ! set to .true. or .false.
    
   ! Variables and parameters de Lier 2013
    real pi
    real a, r0, rx, Kroot, hw, hb, Ks, RLD, z, Ll, raiz3
    real rm, lambda, Ms
    real hs, hl, h0, hx, rho, phi, p, LAMBDA31, AA, BB, PP
    real Mx, M0, Taf, S
    integer i
    character*50 InputFile, AuxFile, OutputFile 
    character*1 dum
      
    !Body of the model
        InputFile = "input.qvl"
        open (unit = 1, file = InputFile , status ='old')
          do i=1,4  
            read (1,*) dum
          enddo
          read (1,*) RLD, a, r0, rx, Kroot, Ll, hw, hl
          do i=1,3  
            read (1,*) dum
          enddo
          read (1,*) hb, p, Ks, hs, z
        AuxFile = "auxi.qvl"
        open (unit = 2, file = AuxFile , status ='replace')
         write (2,1) RLD, a, r0, rx, Kroot, Ll, hw, hl
         write (2,2) hb, p, Ks, hs, z
1        format (8f22.15)
2        format (5f22.15)
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
        Ta = 4.37665046487299E-03                                                                 !m d-1   constant   (plant)
        hx = hl+Ta/Ll                                                                             !m       calculated (plant)
        Mx = Kroot*(hx-hw)                                                                        !m²d-1   calculated (plant)
        BB = hx+phi*Ms+AA/((-hw)**p)                                                              !m       calculated (plant)
        PP =(3.0*sqrt(3.0)*sqrt(27.0*AA**2.0-4.0*AA*BB**3.0)+27.0*AA-2.0*BB**3.0)**(1.0/3.0)      !adm     calculated (plant)
        h0 = (BB-PP/raiz3-raiz3*BB**2.0/PP)/3.0                                                   !m       calculated (plant)                                                
        M0 = -Ks*hb/p*((hb/h0)**p-(hb/hw)**p)                                                     !m²d-1   calculated (plant)
        S = rho*(Ms-M0)                                                                           !d-1     calculated (plant)
        Taf = S*z                                                                                 !m d-1   calculated (plant)
    
    print *, "Test routine for computing zero of f (f = Taf - Ta)"
    debug = .true.

        Ta0 = 1.0
		print *, ' '  ! blank line
        call solve(f, fprime, Ta0, Ta, iters, debug)

        print 3, Ta, iters
3       format('solver returns Ta = ', e22.15, ' after', i3, ' iterations')

        fTa = f(Ta)
        print 4, fTa
4       format('the value of f(Ta) is ', e22.15)

        if (abs(Ta-Taf) > 1d-10) then
            print 5, Ta
5          format('*** Unexpected result: Ta = ', e22.15)
        endif
    
    !Recalculate values
        hx = hl+Ta/Ll                                                                             !m       calculated (plant)
        Mx = Kroot*(hx-hw)                                                                        !m²d-1   calculated (plant)
        BB = hx+phi*Ms+AA/((-hw)**p)                                                              !m       calculated (plant)
        PP =(3.0*sqrt(3.0)*sqrt(27.0*AA**2.0-4.0*AA*BB**3.0)+27.0*AA-2.0*BB**3.0)**(1.0/3.0)      !adm     calculated (plant)
        h0 = (BB-PP/raiz3-raiz3*BB**2.0/PP)/3.0                                                   !m       calculated (plant)                                                
        M0 = -Ks*hb/p*((hb/h0)**p-(hb/hw)**p)                                                     !m²d-1   calculated (plant)
        S = rho*(Ms-M0)                                                                           !d-1     calculated (plant)
        Taf = S*z                                                                                 !m d-1   calculated (plant)
        close (1)
        close (2)
        OutputFile = "Info.out"
        open (unit = 3, file = OutputFile, status ='replace')
         write (3,9)'___________________________________________ Outputs of simple MFlux &
         model ___________________________________________ ' 
         write (3,7)'hx', 'h0', 'Mx', 'M0', 'Ms','S', 'Ta'
         write (3,8) hx, h0, Mx, M0, Ms, S, Taf
         write (3,7)'m', 'm', 'm2 d-1', 'm2 d-1', 'm2 d-1', 'd-1', 'm d-1'
7        format (7a17)
8        format (7f17.10)
9        format (a120)
        close (3)
end program DeLier2013