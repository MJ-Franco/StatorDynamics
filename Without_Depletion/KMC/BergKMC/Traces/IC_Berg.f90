!Program for calculate the curves of recruitment and stall from the
!bacterial flagellar motor.

!04 march 2021

!María-José Franco Oñate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include 'BergKMC.f90'


use KMC_T
!use Glauber
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=10000, n_st=5000

!integer :: ntt(n_sim)

real*8 :: av_h(n_sim),rn0(n_sim)

integer, allocatable :: nt_v(:,:)
real*8, allocatable :: t(:,:), av_step(:),std_step(:),t_step(:)
	
character(len=6) :: FMT1
character(len=6) :: FMT

	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each stepJ
     
     
write(*,*) 'Initial conditions simulation'
write(*,*)     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


rk0=0.00575
alpha=30.44

eMu=1.19d0

!300nm alpha=3.06, emu=-0.200
!500nm alpha =6.15, eMu=0.800
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j=1,n_sim
  rn0(j)=0.0d0
end do

write(FMT,'(F4.2)') eJ

open(80,file='Initial_conditions_J='//FMT//'.dat')


call KMC_Traces(rk0,alpha,eMu,rn0,n_v,n_sim,n_st,nt_v,t,ntt)

allocate(av_step(ntt))
allocate(std_step(ntt))
allocate(t_step(ntt))


do j=1,n_sim
    write(80,*) nt_v(j,n_st)
end do
      
av_histories=Average(1,n_sim,av_h)
      
Av_N = av_histories
      
write(*,*) 'Average  phi= ', Av_N


!Average of the distribution
sav=0.0d0
do j=1,n_sim
    sav = sav + nt_v(j,n_st)
end do

av_c = sav/(real(n_sim)-1.0d0)

!Standard deviation of the distribution
ssd=0.0d0
do j=1,n_sim
    ssd = ssd + (nt_v(j,n_st)-av_c)**2.0d0
end do

sd_c = dsqrt(ssd/real(n_sim))

!Skewness and kurtosis of the distrubution
ssk=0.0d0
sku=0.0d0
do j=1,n_sim
    ssk = ssk + (nt_v(j,n_st)-av_c)**3.0d0
    sku = sku + (nt_v(j,n_st)-av_c)**4.0d0
end do

sk_c = ssk/(real(n_sim)*sd_c**3.0d0)
rku_c = sku/(real(n_sim)*sd_c**4.0d0)


      
      open(50,file='Statistics_Initial_contidions_J='//FMT//'.dat')
      
      write(50,*) av_c, sd_c, sk_c, rku_c
      
     
close(70)
close(50)	  
	  	  
end program
