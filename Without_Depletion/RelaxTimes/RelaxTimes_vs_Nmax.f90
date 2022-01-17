!Maria-Jose Franco Oñate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program for obtaining the relaxation time of stall and resurrection
!depending on the size of the system in order to observe if there is 
!important size effects.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!12-04-2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


include 'GlauberTraces.f90'


use GLAUBER

implicit real*8 (a-h,o-z)

integer, parameter :: n_sim=5000, n_st=30000, n_v=100

real*8 :: t_v(n_sim,n_st),av_h(n_sim),var_h(n_sim),av_step(n_st), &
std_step(n_st),sd_h(n_sim),n0(n_sim)


!character(len=6) :: FMT
character(len=3) :: nsize
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step
    
write(*,*) 'Relaxation times vs Nmax simulation'
write(*,*)
 
eJ = 4.00d0    !Value of J (cooperativity)
eMu = -3.89d0  !Total chemical potential of the motor (mu-ep)

do j=1,n_sim
  n0(j)=0.0d0/real(n_v)
end do


write(nsize,'(I3)') n_v

if(n0(10).eq.0.0d0) then
    open(50, file='RelaxTime_vs_Nmax='//nsize//'_Res.dat')
else
    open(50, file='RelaxTime_vs_Nmax='//nsize//'_Stall.dat')
end if


call GTraces(eJ,eMu,n0,n_v,n_sim,n_st,t_v,av_h,var_h,sd_h)
    
    
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

write(*,*) 'Nmax = ', n_v

write(*,*)
      
write(*,*) 'Nss =', real_av_total
write(*,*) 'sd =', real_sd_total

write(*,*) 


close(50)

write(*,*)
write(*,*) 'Program finished'

end program
