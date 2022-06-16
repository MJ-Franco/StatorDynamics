!Program for calculate the curves of recruitment and stall from the
!bacterial flagellar motor with Kinetic Monte Carlo

!29 march 2021

!María-José Franco Oñate

include 'SubNBerg_KMC.f90'


use KMC_New_Berg
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=5000, n_st=1000

!integer :: ntt(n_sim)

real*8 :: rn0(n_sim)

integer, allocatable :: nt_v(:,:)
real*8, allocatable :: t(:,:), av_step(:),std_step(:),t_step(:)
	
!character(len=6) :: FMT1
!character(len=6) :: FMT
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step
    
integer, parameter :: nPROCESS=1
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rkon = 0.0037d0
rkt = 0.13d0
rkl = 0.0017d0
rkoff = 0.057d0



! character(len=4) :: FMT
! character(len=4) :: FMT1

! open(70,file='Initial_conditions_J='//FMT1//'.dat')

if(nPROCESS==0) then

  sumn0 = 0.0d0
  do j=1,n_sim

      rn0(j)=0.0d0
    
      sumn0 = sumn0 + rn0(j)
    
  end do

end if

if (nPROCESS==1) then

  do j=1,n_sim

    rn0(j)=real(n_v)
    
    sumn0 = sumn0 + rn0(j)
    
  end do
    
!    call KMC_NBerg(rkon,rkt,rkl,rkoff,rn0,n_v,n_sim,n_st,nt_v,t,ntt)
    
!    sumn0=0.0d0
!    do j=1,n_sim
!      rn0(j)=nt_v(j,ntt0)
      
!      sumn0 = sumn0 + rn0(j)
!    end do

  end if

rmrn0=sumn0/real(n_sim)
    
write(*,*) '<phi_0>=', rmrn0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call KMC_NBerg(rkon,rkt,rkl,rkoff,rn0,n_v,n_sim,n_st,nt_v,t,ntt)

allocate(av_step(ntt))
allocate(std_step(ntt))
allocate(t_step(ntt))


!File opening


if(nPROCESS==0) then
   open(50,file='Traces_Res_KMC_NBerg.dat')
   
   open(70, file='All_Traces_Res_KMC_NBerg.dat')

else
   open(50,file='Traces_Stall_KMC_NBerg.dat')
   
   open(70, file='All_Traces_Stall_KMC_NBerg.dat')
end if

do i=1,nsim
  do j = 1,ntt
    write(70,*) i, t(i,j),nt_v(i,j) 
  end do
end do


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
  
deallocate(nt_v)
deallocate(t)



!File closing
     
close(50)


write(*,*)      
write(*,*) 'Program finished'
end 
