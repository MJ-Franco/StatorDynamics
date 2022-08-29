*María-José Franco Oñate
************************************************************************
*Program for obtain the dependence of variance on Nmax
************************************************************************
*DECEMBER 2020
************************************************************************ 

      include 'SubTraces_w_Depletion.f'
      implicit real*8(a-h,o-z)
      
      external Average
      real*8 Average
      
!      parameter(n_v=12)    !Vector size
      parameter(n_sim=500) !Number of simulations
      parameter(n_st=5000) !Number of steps in each simulation
      parameter(n0=0)      !Initial number of stators in the rotor
      parameter(n=10)
      character*3 m
      character*3 mn
      
      dimension t_v(n_sim,n_st),av_histories(n), av_h(n_sim),
     & var_histories(n), sd_var(n),var_h(n_sim),r_v(n_sim,n_st),
     & sd_h(n_sim), avr_h(n_sim),sd_histories(n)
     
     
       !t_v(n_sim,n_st) total number of stators in simulation m in step n
       !av_h (n_sim) final average of each history
       !sd_h (n_sim) final sd of each history
       !av_step(n_st) average number of stators in each step
       !sd_step(n_st) standard deviation in each step
       !av_histories(n) average of all histories
       !var_histories(n)_ average var of all histories
       
       
****************************************     
*Constants
****************************************
      eJ=2.0d0     !Cooperativity term 
      eMu0 =-2.0d0  !Total chemical potential of the motor (mu-ep)
      ntot=11

*****************************************
*Opening file
*****************************************
      Mu=eMu0
      
      J=eJ
      
      write(m,'(i3)') Mu
      write(mn,'(i3)') J
      
      open(70,file='Var_Dependence_J='//mn//'mu0='//m//'.dat')
      
****************************************
*Variance dependence on Nmax
****************************************      
      
      do i=1,n
      
        n_v=10*i
        
        call tracesFM(eJ,eMu0,ntot,n0,n_v,n_sim,n_st,t_v,r_v,
     & av_h,var_h,sd_h,avr_h,sdr_h)
        
        av_histories(i)=Average(1,n_sim,av_h)

        var_histories(i)=Average(1,n_sim,var_h)
        
        sd_histories(i)=Average(1,n_sim,sd_h)
        
        sdv=0.0d0
        do j=1,n_sim
            sdv = sdv+(var_h(j)-var_histories(i))**2.0d0
        end do
        
        sd_var(i)=dsqrt(sdv/real(n_sim))
        
        
        write(70,*) i, av_histories(i), sd_histories(i),  
     & var_histories(i), sd_var(i)
        
        write(*,*) i, av_histories(i), var_histories(i)
        
      
      end do
      
      
      write(*,*)
      
      
      
      
      write(*,*) 'Program finished'
      end program
