!Program for calculate the curves of recruitment and stall from the
!bacterial flagellar motor with Kinetic Monte Carlo

!29 march 2021

!María-José Franco Oñate

include 'Berg_CoopKMC.f90'


use Berg_CoopKMC
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=10000, n_st=20000

!integer :: ntt(n_sim)

real*8 :: av_h(n_sim),rn0(n_sim)


real*8, allocatable :: t(:,:), av_step(:),std_step(:),t_step(:),t_v(:,:)
	
character(len=4) :: FMT1
character(len=6) :: FMT
character(len=5) :: fm
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

eJ=5.00d0
eMu=-5.092d0

rk0=0.00575d0
alpha=30.44d0

write(FMT1,'(F4.2)') eJ !F4.2
write(FMT,'(F5.2)') eMu !F4.2
write(fm,'(F5.2)') alpha


open(70,file='Initial_conditions_J='//FMT1//'_a='//fm// &
   '_N=6.dat')

do j=1,n_sim
  !read(70,*) rn0(j)
  
  rn0(j)=0.0d0
end do


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*) rn0(1)

call BergCoopKMC(eJ,eMu,rk0,alpha,rn0,n_v,n_sim,n_st,t_v,t,ntt)


allocate(av_step(ntt))
allocate(std_step(ntt))
allocate(t_step(ntt))

write(*,*) ntt

!File opening


if(rn0(1)==0.0d0) then
   open(50,file='Traces_Res_KMC_Coop_Berg_J='//FMT1//'_a='//fm// &
   '_N=4.dat')
else
   open(50,file='Traces_Stall_KMC_Coop_Berg_J='//FMT1//'_a='//fm// &
   '_N=4.dat')
end if

open(60,file='Av_KMC_Coop_Berg_J='//FMT1//'_mu='//FMT//'_&
a='//fm//'.dat')


do i=1,ntt

  !Average value and time of each step

  sum_step=0.0d0
  sum_t=0.0d0
  
  
  
  do j=1,n_sim
    
    
    sum_step = sum_step + t_v(j,i)
    sum_t = sum_t + t(j,i)
    
  end do
  
  av_step(i) = sum_step/real(n_sim)
  t_step(i) = sum_t/real(n_sim)
  
 
  
  !Standard deviation of each step
  
  sum_sd=0.0d0
  
  do j=1,n_sim
  
    sum_sd = sum_sd + (t_v(j,i)-av_step(i))**2.0d0
  
  end do

  std_step(i) = dsqrt(sum_sd/real(n_sim))
  
  !File filling
  
  write(50,*) t_step(i), av_step(i), std_step(i)
  
end do

total_av = Average(ntt-3000,ntt,av_step)
total_sd = Average(ntt-3000,ntt,std_step)
!av_time = Average(1,ntt,t_step)

write(60,*) eJ, total_av, total_sd

write(*,*) '<N> = ', total_av, ' +- ', total_sd
!write(*,*) '<t> = ', av_time
  




!File closing
     
close(50)
close(60)

write(*,*)      
write(*,*) 'Program finished'
end 
