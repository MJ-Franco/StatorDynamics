!Program for calculate the curves of recruitment and stall from the
!bacterial flagellar motor with depletion and Glauber dynamics

!04 march 2021

!María-José Franco Oñate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include 'Glauber_Traces_w_Depletion.f90'

use GLAUBER
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=10000, n_st=7000, ntot=13


real*8 :: t_v(n_sim,n_st),av_h(n_sim),var_h(n_sim),av_step(n_st), &
std_step(n_st),sd_h(n_sim),r_v(n_sim,n_st),avr_h(n_sim),sdr_h(n_sim), &
av_res(n_st),n0(n_sim)

	
character(len=6) :: FMT
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


eJ = 5.00d0    !Value of J (cooperativity)
eMu0 = -4.15d0  !Total chemical potential of the motor (mu-ep)

!0 1.55   0.25 1.25   0.50 0.90   0.75 0.60       !0 0.55   0.25 0.30   0.50 0.05   0.75 -0.21    !0 -0.40  0.25 -0.57  0.50 -0.75  0.75 -0.95
!1 0.30   1.25 0.00   1.50 -0.29  1.75 -0.57      !1 -0.45  1.25 -0.7   1.50 -0.95  1.75 -1.20    !1 -1.14  1.25 -1.33  1.50 -1.52  1.75 -1.73
!2 -0.87  2.25 -1.15  2.50 -1.42  2.75 -1.70      !2 -1.45  2.25 -1.70  2.50 -1.95  2.75 -2.20    !2 -1.90  2.25 -2.13  2.50 -2.32  2.75 -2.56 
!3 -1.97  3.25 -2.24  3.50 -2.50  3.75 -2.77      !3 -2.45  3.25 -2.65  3.50 -2.87  3.75 -3.10    !3 -2.78  3.25 -2.98  3.50 -3.22  3.75 -3.44
!4 -3.02                                          !4 -3.35                                        !4 -3.66
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(FMT,'(F4.2)') eJ

open(70,file='Initial_conditions_J='//FMT//'.dat')

do j=1,n_sim
  !read(70,*) n0(j)
  
  n0(j)=0.0d0/real(n_v)
end do

write(*,*) n0(10)

call GTraces(eJ,eMu0,ntot,n0,n_v,n_sim,n_st,t_v,r_v,av_h, &
var_h,sd_h,avr_h,sdr_h)


!File opening

      write(FMT,'(F4.2)') eJ
      
      if(n0(10).eq.0.0d0) then
        open(50, file='Traces_Res_J='//FMT//'.dat')
      else
        open(50, file='Traces_Stall_J='//FMT//'.dat')
      end if
      
   
!Average and sd of all histories 

      !Stators in the rotor 
      av_histories=Average(1,n_sim,av_h)
      
      Av_N= av_histories*real(n_v)
      
      write(*,*) 'Average  phi= ', av_histories

      var_histories=Average(1,n_sim,var_h)
      
      var_N=var_histories*real(n_v)
      
      sd_histories=sqrt(var_histories)
      
      write(*,*) 'Average sd phi = ', sd_histories
      
      write(*,*)
      
      write(*,*) 'Nss = ', Av_N
      
      write(*,*) 'var= ', var_N
      
      write(*,*)
      

!Average and sd of each step


      do l=1,n_st
        ssp=0.0d0
        do j=1,n_sim
            ssp=ssp+t_v(j,l)
        end do
        av_step(l)=ssp/real(n_sim)
      end do
      
      !Real average of the curve
      real_av_total=Average(100,n_st,av_step)
      
      !Standard deviation of average curve
      do l=1,n_st
        ssd=0.0d0
        do j=1,n_sim
            ssd=ssd+(t_v(j,l)-real_av_total)**2.0d0
        end do
        std_step(l)=sqrt(ssd/real(n_sim))
        
        write(50,*) l, av_step(l), std_step(l)
      end do  
      
      real_std_total=Average(100,n_st,std_step)

      write(*,*)
      
      write(*,*) 'Nss=', real_av_total
      write(*,*) 'sd=', real_sd_total
      
      write(*,*) 
      

!File closing
     
      close(50)


      
      write(*,*) 'Program finished'
      end 
