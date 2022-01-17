!Program for calculate the curves of recruitment and stall from the
!bacterial flagellar motor with Kinetic Monte Carlo

!29 march 2021

!María-José Franco Oñate

include 'KMC_V6.f90'


use KMC_T
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=10000, n_st=5000

!integer :: ntt(n_sim)

real*8 :: t_v(n_sim,n_st),av_h(n_sim),rn0(n_sim)

integer, allocatable :: nt_v(:,:)
real*8, allocatable :: t(:,:), av_step(:),std_step(:),t_step(:)
	
character(len=6) :: FMT1
character(len=6) :: FMT
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rkon = 2.18d-3

rkoff = 3.96d-3

!eJ = 0.00d0    !Value of J (cooperativity)
!eMu =-0.846d0 !Total chemical potential of the motor (mu-ep)

!ek0=0.00575d0
!alpha=3.06d0

write(FMT1,'(F4.2)') eJ !F4.2
write(FMT,'(F5.2)') eMu !F4.2

open(70,file='Initial_conditions_J='//FMT1//'.dat')

do j=1,n_sim
  read(70,*) rn0(j)
  
  !rn0(j)=0.0d0
end do


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!call KMC_Traces(rkon,rkoff,rn0,n_v,n_sim,n_st,t_v,t,av_h)

call KMC_Traces(rkon,rkoff,rn0,n_v,n_sim,n_st,nt_v,t,ntt)

allocate(av_step(ntt))
allocate(std_step(ntt))
allocate(t_step(ntt))


!File opening


if(rn0(10)==0.0d0) then
   open(50,file='Traces_Res_KMC.dat')
else
   open(50,file='Traces_Stall_KMC.dat')
end if


do i=1,ntt

  !Average value and time of each step

  sum_step=0.0d0
  sum_t=0.0d0
  
  do j=1,n_sim
  
    sum_step = sum_step + nt_v(j,i)
    sum_t  = sum_t + t(j,i)
  
  end do
  
  av_step(i) = sum_step/real(n_sim)
  t_step(i) = sum_t/real(n_sim)
  
  
  !Standard deviation of each step
  
  sum_sd=0.0d0
  
  do j=1,n_sim
  
    sum_sd = sum_sd + (nt_v(j,i)-av_step(i))**2.0d0
  
  end do

  std_step(i) = dsqrt(sum_sd/real(n_sim))
  
  !File filling
  
  write(50,*) t_step(i), av_step(i), std_step(i)
  
end do

total_av = Average(100,ntt,av_step)
total_sd = Average(100,ntt,std_step)
av_time = Average(1,ntt,t_step)

write(*,*) '<N> = ', total_av, ' +- ', total_sd
write(*,*) '<t> = ', av_time
  




!File closing
     
close(50)


write(*,*)      
write(*,*) 'Program finished'
end 
