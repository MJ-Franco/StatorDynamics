
!Maria Jose Franco Onate

!Subroutine for calculate the traces of recuitment/detachemt of bacterial
!flagellar motor using the model of lattice gas

!OCTOBER 2020

module TRACES

contains
       
    subroutine tracesFM(eJ,eMu,n0,n_v,n_sim,n_st,phi_v,av_h,var_h,sd_h)
    
      
      !eJ =  Value of J (cooperativity)
      !eMu =  total chemical potential of the motor (mu+ep)
      !eMu0 = constant of reservoir chemical potential 
      !ntot = total number of stators in the system
      !n_v = Vector size
      !n_sim = Number of simulations
      !n_st = Number of steps in each simulation
      !n0 =Initial number of stators in the rotor
      
      implicit real*8 (a-h,o-z)

      
      real*8 :: s(n_v),t_v(n_sim,n_st),av_h(n_sim),var_h(n_sim), &
      av2_h(n_sim), phi2_v(n_sim,n_st), t2h_v(n_sim,n_st), &
      phi_v(n_sim,n_st), sd_h(n_sim)
      
      !s(n_v) variable that denotes the state of occupation of a place in the rotor (0,1) 
      !t_v(n_sim,n_st) total number of stators in simulation m at step n
      !r_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
      !av_h (n_sim) final average of each history
      !sd_h (n_sim) final sd of each history
      

!Parameters for random numbers

      iseed = -123456789 
      
      an= 1.0d0
      bn= real(n_v)+1.0d0
      

!File opening
     
      if(n0.eq.0) then  
        open(10,file='Traces_Res_All_Variances.dat')
        open(20,file='Traces_Res_All_Averages.dat')
        open(30,file='Traces_Res.dat')
      else
        open(10,file='Traces_Stall_All_Variances.dat')
        open(20,file='Traces_Stall_All_Averages.dat')
        open(30,file='Traces_Stall.dat')
      end if
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      do j=1,n_sim  !Do for several histories
              

!Vector initialisation


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
         


        do l=1,n_st   !The simulation for one story begins
         
           

!Random site choosing
           
            r=ran2(iseed)
            rr=an+(bn-an)*r
            i=int(rr)
            
            
            do k=1,n_v   !Montecarlo step
            
            
!Boundary conditions

                if((i+1).gt.n_v) then 
                    suma = s(i-1)+s(1)    
                elseif ((i-1).lt.1) then
                    suma = s(i+1) + s(n_v)
                else
                    suma = s(i+1)+s(i-1)
                end if
                

!Energy state
 
                if(suma.eq.2.0d0) then
                    H0 = -eMu 
                    H1 = -2.0d0*eJ - 2.0d0*eMu   
                else if(suma.eq.1.0d0) then
                    H0 = -eMu/2.0d0
                    H1 = -eJ-1.5d0*eMu
                else if(suma.eq.0.0d0) then
                    H0 = 0.0d0
                    H1 = -eMu
                end if



!Metropolis algorithm

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
            

            end do !End of Montecarlo step



!Total number of stators in the rotor

          
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
            
            
            
      
         end do !End of simulation step



!Average and sd of each history

      !Stators in the rotor
      sav=0.0d0
      sav2=0.0d0
      do l=1,n_st
        sav=sav+phi_v(j,l)
        sav2=sav2 + phi2_v(j,l)
      end do
      
      
      av_h(j) = sav/real(n_st)
      av2_h(j) = sav2/real(n_st)
      var_h(j) = av2_h(j)-av_h(j)**2.0d0
      sd_h(j) = dsqrt(var_h(j))
      
      write(10,*) j, var_h(j)
      write(20,*) j, av_h(j)

        
      
      

      end do !End of history
        
        
      

!File closing
   
      close(10)
      close(20)
      close(30)

      
      return
      end
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      FUNCTION ran2(idum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      IMPLICIT REAL*8 (A-H,O-Z)
!      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
!      REAL*8 ran2,AM,EPS,RNMX
      INTEGER,PARAMETER :: IM1=2147483563,IM2=2147483399,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB
      REAL*8, PARAMETER :: AM=1.0d0/IM1,EPS=1.2d-7,RNMX=1.0d0-EPS
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 function Average(it_min,it_max,dat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit real*8(a-h,o-z)
      dimension dat(int(it_max-it_min))
        
        sm=0.0d0
        do j=it_min,it_max
            sm=sm+dat(j)
        end do
        
        Average=sm/real(it_max-it_min)
        
      return  
      end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      real*8 function Av_2(it_min,it_max,dat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit real*8(a-h,o-z)
      dimension dat(int(it_max-it_min))
        
        sm=0.0d0
        do j=it_min,it_max
            sm=sm+dat(j)**2.0d0
        end do
        
        Av_2=sm/real(it_max-it_min)
        
      return  
      end
      


end module