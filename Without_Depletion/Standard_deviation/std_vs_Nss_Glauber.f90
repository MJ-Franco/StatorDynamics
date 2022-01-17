!Maria José Franco Oñate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program to obtain the dependence of standard deviation on steady state
!using Glauber dynamics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!12-04-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include 'GlauberTraces.f90'

use GLAUBER
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=13, n_sim=10000, n_st=3000, n=14

character(len=1) m

real*8 :: phi_v(n_sim,n_st),av_h(n_sim),var_h(n_sim),av_step(n_st), &
std_step(n_st),sd_h(n_sim),n0(n_sim),av_histories(n),sd_histories(n), &
sd_var(n)

write(*,*) 'Standard deviation vs Nss simulation'

eJ=0.0d0

do j=1,n_sim
  n0(j)=0.0d0/real(n_v)
end do

!!!!!!!!!!!!!!!!
!File opening
!!!!!!!!!!!!!!!!
J=eJ
      
write(m,'(i1)') J


if(n0(1)==0.0d0) then      
    open(70,file='sdt_vs_Nss_J='//m//'_Res.dat')
    open(80,file='std_vs_Nss_J='//m//'_Statistics_Res_dat')
else
    open(70,file='std_vs_Nss_J='//m//'_Stall.dat')
    open(80,file='std_vs_Nss_J='//m//'_Statistics_Stall.dat')
end if


eMu = -6.0d0

do i=1,n
        
    eMu = eMu + 1.0d0
            
    call GTraces(eJ,eMu,n0,n_v,n_sim,n_st,phi_v,av_h,var_h,sd_h)
            
            
    av_histories(i)=Average(1,n_sim,av_h)
        
    sd_histories(i)=Average(1,n_sim,sd_h)
    
    !standard deviation of the standard deviation    
    sdv=0.0d0
    do j=1,n_sim
        sdv = sdv+(sd_h(j)-sd_histories(i))**2.0d0
    end do
        
    sd_var(i)=dsqrt(sdv/real(n_sim))
        
    write(70,*) eMu, av_histories(i), sd_histories(i), sd_var(i)
        
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
