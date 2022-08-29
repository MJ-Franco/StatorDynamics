*María-José Franco Oñate
************************************************************************
*Subroutine for calculate the traces of recuitment/detachemt of bacterial
*flagellar motor using the model of lattice gas
************************************************************************
*OCTOBER 2020
************************************************************************
       
      subroutine tracesFM(eJ,eMu0,ntot,n0,n_v,n_sim,n_st,phi_v,r_v,
     & av_h,var_h,sd_h,avr_h,sdr_h)
      
      !eJ =  Value of J (cooperativity)
      !eMu =  total chemical potential of the motor (mu+ep)
      !eMu0 = constant of reservoir chemical potential 
      !ntot = total number of stators in the system
      !n_v = Vector size
      !n_sim = Number of simulations
      !n_st = Number of steps in each simulation
      !n0 =Initial number of stators in the rotor
      
      implicit real*8 (a-h,o-z)
      external ran2
      external Average
      real*8 Average
      external var
      real*8 var
      
      dimension s(n_v),t_v(n_sim,n_st),r_v(n_sim,n_st),av_h(n_sim),
     & sd_h(n_sim),var_h(n_sim),avr_h(n_sim),varr_h(n_sim),sdr_h(n_sim),
     & phi_v(n_sim,n_st),t2h_v(n_sim,n_st),phi2_v(n_sim,n_st),
     & av2_h(n_sim)
      
      !s(n_v) variable that denotes the state of occupation of a place in the rotor (0,1) 
      !t_v(n_sim,n_st) total number of stators in simulation m at step n
      !r_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
      !av_h (n_sim) final average of each history
      !sd_h (n_sim) final sd of each history
      
!      tot=real(ntot)
      
!      eMu0 =eMu-dlog(tot)
      
!      write(*,*) 'mu0=',eMu0
******************************
*Parameters for random numbers
****************************** 
      iseed = -123456789 
      
      an= 1.0d0
      bn= real(n_v)+1.0d0
      
*****************************
*File opening
*****************************      
      if(n0.eq.0) then  
        open(10,file='Traces_Res_All_Variances.dat')
        open(20,file='Traces_Res_All_Averages.dat')
        open(30,file='Traces_Res.dat')
        open(40,file='Reservoir_Res.dat')
      else
        open(10,file='Traces_Stall_All_Variances.dat')
        open(20,file='Traces_Stall_All_Averages.dat')
        open(30,file='Traces_Stall.dat')
        open(40,file='Reservoir_Stall.dat')
      end if
      
      open(50,file='emu_emur.dat')
      

******************************
*Simulation
******************************

************************************************************************      
      do j=1,n_sim  !Do for several histories
************************************************************************        
******************************
*Vector initialisation
******************************

        !Fill vectors with 0
        do i=1,n_v
            s(i)=0.0d0
        end do
        
        
        !If n0=0 then it keepts going
        if(n0.eq.0) then
            continue
        !If n0=\0 it fills the vectors with n0 1s randomly    
        else    
            do k=1,n0
        
5               t=ran2(iseed)
                tt=an+(bn-an)*t
                i=int(tt)
                
                if(s(i).eq.0.0d0) then  
                    s(i)=1.0d0
                else
                    goto 5
                end if
            end do
         end if
         
         
*******************************
*mur initialisation (reservoir potential)
*******************************
      if(ntot.lt.n0) then
        write(*,*)'Stators in the system less than stators in the motor'
        stop
      end if
      
      eMur = eMu0 + dlog(real(ntot)-real(n0)+1.0d0)
      
      phi0=real(n0)/real(n_v)
      
      B=dexp(-eJ)
      A=2.0d0*phi0-1.0d0
            
      eMu = 2.0d0*dlog(A*dsqrt(B)+dsqrt(1.0d0+A**2.0d0*(B-1)))-
     & dlog(1.0d0-A**2.0d0) - eJ  
      
************************************************************************
        do l=1,n_st   !The simulation for one story begins
************************************************************************         
           
*****************************
*Random site choosing
*****************************           
            r=ran2(iseed)
            rr=an+(bn-an)*r
            i=int(rr)
            
************************************************************************            
            do k=1,n_v   !Montecarlo step
************************************************************************            
            
******************************
*Boundary conditions
******************************
                if((i+1).gt.n_v) then 
                    suma = s(i-1)+s(1)    
                elseif ((i-1).lt.1) then
                    suma = s(i+1) + s(n_v)
                else
                    suma = s(i+1)+s(i-1)
                end if
                
******************************
*Energy state
****************************** 
                if(suma.eq.2.0d0) then
                    H0 = -eMur 
                    H1 = -2.0d0*eJ - 2.0d0*eMur   
                else if(suma.eq.1.0d0) then
                    H0 = -eMur/2.0d0
                    H1 = -eJ-1.5d0*eMur
                else if(suma.eq.0.0d0) then
                    H0 = 0.0d0
                    H1 = -eMur
                end if


*******************************
*Metropolis algorithm
*******************************
                !Possible change from 0 to 1
                if(s(i).eq.0.0d0) then
        
                    if(H1.gt.H0) then
                        p0=dexp(-(H1-H0))
                        r0=ran2(iseed)
                        if(p0.gt.r0) then
                            s(i)=1.0d0
                        else
                            continue
                        end if
                    else
                        s(i) = 1.0d0
                    end if
                !Possible change from 1 t0 0
                else if(s(i).eq.1.0d0) then
                
                    if(H0.gt.H1) then
                        p1=dexp(-(H0-H1))
                        r1=ran2(iseed)
                        if(p1.gt.r1) then
                            s(i)=0.0d0
                        else
                            continue
                        end if
                    else
                        s(i)=0.0d0
                    end if
                end if
        
************************************************************************
            end do !End of Montecarlo step
************************************************************************

*********************************
*Total number of stators in the rotor
*********************************
          
            t_v(j,l)=0.0d0
            do i=1,n_v
                t_v(j,l) = t_v(j,l) + s(i)
            end do
            
            phi_v(j,l)=t_v(j,l)/real(n_v)
            
            !We write in a file this number
            write(30,*) j, l, phi_v(j,l)
            
            t2h_v(j,l)=0.0d0
            do i=1,n_v
                do k=1,n_v
                    t2h_v(j,l)= t2h_v(j,l) + s(i)*s(k)
                end do
            end do

            phi2_v(j,l) = t2h_v(j,l)/real(n_v)**2.0d0

        
**********************************
*Recalculation of mur (reservoir chemical potential)
**********************************
            
            r_v(j,l)=real(ntot) - t_v(j,l)
            
            if(r_v(j,l).le.-1.0d0) then
                eMur=-100.0d0
            else
                eMur = eMu0 + dlog(real(ntot)-t_v(j,l)+1.0d0)
            end if
            
            
            B=dexp(-eJ)
            A=2.0d0*phi_v(j,l)-1.0d0
            
            eMu = 2.0d0*dlog(A*dsqrt(B)+dsqrt(1.0d0+A**2.0d0*(B-1)))-
     &      dlog(1.0d0-A**2.0d0) - eJ   
         
         
            tol=0.001d0
            
            
            
            if(abs(eMu-eMur).le.tol) then
                write(50,*) j,l, 'equilibrium'
            end if
            
            
            write(40,*) j, l , r_v(j,l)
            

            
************************************************************************      
         end do !End of simulation step
************************************************************************
        

**********************************
*Average and sd of each history
**********************************
      !Stators in the rotor
*10    
      sav=0.0d0
      sav2=0.0d0
      do l=1,n_st
        sav=sav+phi_v(j,l)
        sav2=sav2 + phi2_v(j,l)
      end do
      
      
      av_h(j) = sav/real(n_st)
      av2_h(j) = sav2/real(n_st)
      var_h(j) = av2_h(j)-av_h(j)**2.0d0
      sd_h(j)=dsqrt(av2_h(j)-av_h(j)**2.0d0)
      
      write(10,*) j, var_h(j)
      write(20,*) j, av_h(j), sd_h(j)
      
      write(10,*) j, var_h(j)
      write(20,*) j, av_h(j), sd_h(j)
      
      !Stators in the reservoir
*      savr=0.0d0
*      svarr=0.0d0
*      do l=1,n_stop
*        savr=savr+r_v(j,l)
*        svarr=svarr+r_v(j,l)**2.0d0
*      end do
      
*      avr_h(j) = savr/real(n_stop)
*      varr_h(j) = svarr/real(n_stop)
*      sdr_h(j) = dsqrt(abs(varr_h(j)-avr_h(j)**2.0d0))
      
      
************************************************************************
      end do !End of history
************************************************************************        
     
     
      
******************************
*File closing
******************************      
      close(10)
      close(20)
      close(30)
      close(40)
      close(50)

      
      return
      end
      


************************************************************************      
      FUNCTION ran2(idum)
************************************************************************      
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.0d0/IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.0d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END FUNCTION

************************************************************************
      real*8 function Average(it_min,it_max,dat)
************************************************************************
      implicit real*8(a-h,o-z)
      dimension dat(int(it_max-it_min))
        
        sm=0.0d0
        do j=it_min,it_max
            sm=sm+dat(j)
        end do
        
        Average=sm/real(it_max-it_min)
        
      return  
      end 

************************************************************************      
      real*8 function Av_2(it_min,it_max,dat)
************************************************************************
      implicit real*8(a-h,o-z)
      dimension dat(int(it_max-it_min))
        
        sm=0.0d0
        do j=it_min,it_max
            sm=sm+dat(j)**2.0d0
        end do
        
        Av_2=sm/real(it_max-it_min)
        
      return  
      end 
