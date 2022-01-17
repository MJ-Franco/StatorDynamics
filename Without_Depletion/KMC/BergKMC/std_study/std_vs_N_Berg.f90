!Maria José Franco Oñate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program to obtain the dependence of standard deviation on steady state
!using Berg dynamics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!30-11-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include 'BergKMC.f90'

use KMC_t
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=1000, n_st=10000, n=30

real*8 :: av_sim(n_sim), sd_sim(n_sim),rn0(n_sim),av_histories(n), &
sd_histories(n),sd_var(n)

integer, allocatable :: nt_v(:,:)
real*8, allocatable :: t(:,:), av_step(:),std_step(:),t_step(:)

write(*,*) 'Standard deviation vs Nss simulation'

rk0=0.00575

do j=1,n_sim
  rn0(j)=0.0d0
end do

!!!!!!!!!!!!!!!!
!File opening
!!!!!!!!!!!!!!!!

      


if(rn0(1)==0.0d0) then      
    open(70,file='sdt_vs_Nss_alpha=30.44.dat')
    open(80,file='std_vs_Nss_Statistics_Res_dat')
else
    open(70,file='std_vs_Nss_Stall.dat')
    open(80,file='std_vs_Nss_Statistics_Stall.dat')
end if


eMu = 1.19
alpha=0.0d0


do i=1,n
    
        
    alpha = alpha + 1.0d0
            
            
    call KMC_Traces(rk0,alpha,eMu,rn0,n_v,n_sim,n_st,nt_v,t,ntt,av_sim,&
    sd_sim)

    av_histories(i)=Average(1,n_sim,av_sim)
        
    sd_histories(i)=Average(1,n_sim,sd_sim)       
            
    !standard deviation of the standard deviation    
    sdv=0.0d0
    do j=1,n_sim
        sdv = sdv+(sd_sim(j)-sd_histories(i))**2.0d0
    end do
        
    sd_var(i)=dsqrt(sdv/real(n_sim))
        
    write(70,*) alpha, av_histories(i), sd_histories(i), sd_var(i)
        
    
end do

!Average of standard deviation
sav=0.0d0
do i=1,n
    sav = sav + sd_histories(i)
end do

sd_av = sav/(real(n)-1.0d0)

!Standard deviation of the average standard deviation
ssd = 0.0d0
do i=1,n
    ssd = ssd + (sd_histories(i)-sd_av)**2.0d0
end do

sd_sd = dsqrt(ssd/(real(n)-1.0d0))

!Skewness and kurtosis of the standard deviation
ssk = 0.0d0
sku = 0.0d0
do i=1,n
    ssk = ssk + (sd_histories(i)-sd_av)**3.0d0
    sku = sku + (sd_histories(i)-sd_av)**4.0d0
end do

sd_sk = ssk/((real(n)-1.0d0)*sd_sd**3.0d0)
sd_ku = sku/((real(n)-1.0d0)*sd_sd**4.0d0)

write(80,*) sd_av,sd_sd,sd_sk,sd_ku

write(*,*) 'Average std = ', sd_av
write(*,*) 'Standard deviation av_std = ', sd_sd
write(*,*) 'Skewness of av_std = ', sd_sk
write(*,*) 'Kurtosis of av_std = ', sd_ku

write(*,*)
write(*,*) 'Program finished'

close(70)
close(80)


end program	
