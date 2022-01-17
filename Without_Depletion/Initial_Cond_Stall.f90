!Program for calculate the curves of recruitment and stall from the
!bacterial flagellar motor.

!04 march 2021

!María-José Franco Oñate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!include 'KMC.f90'
include 'GlauberTraces.f90'


!use KMC_T
use Glauber
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=10000, n_st=5000


real*8 :: t_v(n_sim,n_st),av_h(n_sim),av_step(n_st), &
std_step(n_st),sd_h(n_sim),t(n_sim,n_st),t_step(n_st), rn0(n_sim),&
var_h(n_sim)

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

eJ = 0.00d0    !Value of J (cooperativity)
eMu = 0.9d0 !Total chemical potential of the motor (mu-ep)

!ek0=0.00575
!alpha=3.06d0  !Total chemical potential of the motor (mu-ep)

!For Nmax = 12 
      
!0 0.7    !0.25 0.37  !0.5 0.05 !0.75 -0.25     !0 -0.69 !0.25 -0.85 !0.5 -1.03 !0.75 -1.23    !0 1.1    !0.25 0.73   !0.50 0.38   !0.75 0.03
!1 -0.57  !1.25 -0.87 !1.50 -1.15 !1.75 -1.45   !1 -1.41 !1.25 -1.62 !1.50 -1.83 !1.75 -2.04   !1 -0.31  !1.25 -0.64  !1.50 -0.96  !1.75 -1.27
!2 -1.74  !2.25 -2.01 !2.50 -2.29 !2.75 -2.57   !2 -2.25 !2.25 -2.47 !2.50 -2.7 !2.75 -2.92    !2 -1.57  !2.25 -1.87  !2.50 -2.17  !2.75 -2.45
!3 -2.84  !3.25 -3.11 !3.50 -3.37 !3.75 -3.64   !3 -3.15 !3.25 -3.39 !3.50 -3.62 !3.75 -3.86   !3 -2.74  !3.25 -3.02  !3.50 -3.29  !3.75 -3.56
!4 -3.89                                        !4 -4.1                                        !4 -3.83
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j=1,n_sim
  rn0=0.0d0
end do

write(FMT,'(F4.2)') eJ

open(80,file='Initial_conditions_J='//FMT//'.dat')


call GTraces(eJ,eMu,rn0,n_v,n_sim,n_st,t_v,av_h,var_h,sd_h)
!call KMC_Traces(eJ,eMu,ek0,alpha,rn0,n_v,n_sim,n_st,t_v,av_h,sd_h,t)


do j=1,n_sim
    write(80,*) t_v(j,n_st)
end do
      
av_histories=Average(1,n_sim,av_h)
      
Av_N = av_histories*real(n_v)
      
write(*,*) 'Average  phi= ', Av_N


!Average of the distribution
sav=0.0d0
do j=1,n_sim
    sav = sav + t_v(j,n_st)
end do

av_c = sav/(real(n_sim)-1.0d0)

!Standard deviation of the distribution
ssd=0.0d0
do j=1,n_sim
    ssd = ssd + (t_v(j,n_st)-av_c)**2.0d0
end do

sd_c = dsqrt(ssd/real(n_sim))

!Skewness and kurtosis of the distrubution
ssk=0.0d0
sku=0.0d0
do j=1,n_sim
    ssk = ssk + (t_v(j,n_st)-av_c)**3.0d0
    sku = sku + (t_v(j,n_st)-av_c)**4.0d0
end do

sk_c = ssk/(real(n_sim)*sd_c**3.0d0)
rku_c = sku/(real(n_sim)*sd_c**4.0d0)


      
      open(50,file='Statistics_Initial_contidions_J='//FMT//'.dat')
      
      write(50,*) av_c, sd_c, sk_c, rku_c
      
     
close(70)
close(50)	  
	  	  
end program
