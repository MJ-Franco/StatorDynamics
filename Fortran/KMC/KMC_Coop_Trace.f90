!Program for calculate the curves of recruitment and stall from the
!bacterial flagellar motor with Kinetic Monte Carlo

!29 march 2021

!María-José Franco Oñate

include 'CoopKMC_v2.f90'


use CoopKMC
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=10000, n_st=5000

!integer :: ntt(n_sim)

real*8 :: av_h(n_sim),rn0(n_sim)


real*8, allocatable :: t(:,:), av_step(:),std_step(:),t_step(:),t_v(:,:)
	
character(len=4) :: FMT
character(len=4) :: FMT1
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step
  
integer, parameter :: nPROCESS=0
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

eJ = 4.00d0
eMu= -4.00d0

!initial conditions
eJi = 2.00d0 
eMui = -1.48d0

sumn0=0.0d0
do j=1,n_sim

   rn0(j)=0.00d0
   
   sumn0 = sumn0 + rn0(j)
  
end do

!if (nPROCESS==1) then
    
!    call KMCoop(eJi,eMui,rn0,n_v,n_sim,n_st,t0_v,t0,ntt0)
    
!    sumn0=0.0d0
!    do j=1,n_sim
!      rn0(j)=t0_v(j,ntt0)
      
!      sumn0 = sumn0 + rn0(j)
!    end do

!end if
  

rmrn0=sumn0/real(n_sim)   
   
exN_inf =  6.50d0

write(*,*)
write(*,*)

write(*,*) 'J=', eJ
write(*,*) 'N_inf=', exN_inf
write(*,*) '<phi_0>=', rmrn0

write(*,*)

write(FMT1,'(F4.2)') eJ !F4.2
write(FMT,'(F4.2)') exN_inf !F4.2

!File opening
if(nPROCESS==0) then
  open(50, file='Traces_Res_J='//FMT1//'_N='//FMT//'.dat')
  
        
  write(*,*) 'Resurrection'
  write(*,*)
  
else if(nPROCESS==1) then
  open(50, file='Traces_Stall_J='//FMT1//'_N='//FMT//'.dat')
  !open(70, file='All_Traces_Stall_J='//FMT1//'_N='//FMT//'.dat')
        
  write(*,*) 'Release'
  write(*,*)
  
else
  open(50, file='Traces_Stall_J='//FMT1//'_N='//FMT//'_fixed.dat')
  !open(70, file='All_Traces_Stall_J='//FMT1//'_N='//FMT//'_fixed.dat')
        
  write(*,*) 'Fixed release'
  write(*,*)
end if

open(70, file='All_Traces_Res_J='//FMT1//'_N='//FMT//'.dat')
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



call KMCoop(eJ,eMu,rn0,n_v,n_sim,n_st,t_v,t,ntt)

do i=1,n_sim
    do j = 1,ntt
      write(70,*) i, t(i,j), t_v(i,j)
    end do
  end do
   


allocate(av_step(ntt))
allocate(std_step(ntt))
allocate(t_step(ntt))

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
