!María-José Franco Oñate

!Program for obtain the dependence of Nss on ntot

!OCTOBER 2020

include 'Glauber_Traces_w_Depletion.f90'

use GLAUBER
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13,n_sim=10000,n_st=3000, n0=10, n=20


real*8 :: t_v(n_sim,n_st),av_h(n_sim),var_h(n_sim),av_step(n_st), &
std_step(n_st),sd_h(n_sim),r_v(n_sim,n_st),avr_h(n_sim),sdr_h(n_sim), &
av_res(n_st),sd_phi(n),real_av_total(n)

character*6 m
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

eJ = 0.00d0    !Value of J (cooperativity)
eMu0 = 0.47d0  !Total chemical potential of the motor (mu-ep)
      
!0 0.7    !0.25 0.37  !0.5 0.05 !0.75 -0.25     !0 -0.69 !0.25 -0.85 !0.5 -1.03 !0.75 -1.23
!1 -0.57  !1.25 -0.87 !1.50 -1.15 !1.75 -1.45   !1 -1.41 !1.25 -1.62 !1.50 -1.83 !1.75 -2.04
!2 -1.74  !2.25 -2.01 !2.50 -2.29 !2.75 -2.57   !2 -2.25 !2.25 -2.47 !2.50 -2.7 !2.75 -2.92
!3 -2.84  !3.25 -3.11 !3.50 -3.37 !3.75 -3.64   !3 -3.15 !3.25 -3.39 !3.50 -3.62 !3.75 -3.86
!4 -3.87                                        !4 -4.1
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!File opening

J=eJ

write(m,'(i3)') J
      
open(70,file='ntot_dependence_J='//m//'.dat')
      
      
do i=1,n
      
	ntot=10*i   !total number of stators in each simulation
      
    call GTraces(eJ,eMu0,ntot,n0,n_v,n_sim,n_st,t_v,r_v,av_h, &
		var_h,sd_h,avr_h,sdr_h)
     
     do l=1,n_st
        ssp=0.0d0
        do j=1,n_sim
          ssp=ssp+real(t_v(j,l))
        end do
          av_step(l)=ssp/real(n_sim)

     end do
        
        !Nss
     real_av_total(i)=Average(1,n_st,av_step)
        
     sdv=0.0d0
     do l=1,n_st
        sdv = sdv+(av_step(l)-real_av_total(i))**2.0d0
     end do
      
     sd_phi(i)=dsqrt(sdv/real(n_st))
        
        
     write(70,*) ntot, real_av_total(i), sd_phi(i)
        
          
end do
      
close(70)
      
write(*,*) 'Program finished'
end
