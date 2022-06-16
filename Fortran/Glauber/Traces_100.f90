!Program for calculate the curves of recruitment and stall from the
!bacterial flagellar motor.

!04 march 2021

!María-José Franco Oñate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

include 'GlauberTraces.f90'
!include 'FredTraces.f90'


use GLAUBER
!use FRED
implicit real*8 (a-h,o-z)

integer, parameter :: nPROCESS=2

integer, parameter :: n_v=13, n_sim=10000, n_st=1000, nrep=100

real*8 :: t_v(n_sim,n_st),av_h(n_sim),var_h(n_sim),av_step(nrep,n_st), &
std_step(nrep,n_st),sd_h(n_sim),rn0(n_sim),rin0(n_sim)
	
	
character(len=4) :: FMT1 
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

eJ = 0.00d0   !Value of J (cooperativity)
eMu = 0.47d0  !Total chemical potential of the motor (mu-ep)

!initial conditions
eJi = 0.00d0 
eMui = 0.82d0

sumn0=0.0d0
do j=1,n_sim
  rn0(j)=9.0d0/real(n_v)
  
  sumn0 = sumn0 + rn0(j) 
end do
  
!0 0.7    !0.25 0.37  !0.5 0.05 !0.75 -0.25     !0 -0.69 !0.25 -0.85 !0.5 -1.03 !0.75 -1.23
!1 -0.57  !1.25 -0.87 !1.50 -1.15 !1.75 -1.45   !1 -1.41 !1.25 -1.62 !1.50 -1.83 !1.75 -2.04
!2 -1.74  !2.25 -2.01 !2.50 -2.29 !2.75 -2.57   !2 -2.25 !2.25 -2.47 !2.50 -2.7 !2.75 -2.92
!3 -2.84  !3.25 -3.11 !3.50 -3.37 !3.75 -3.64   !3 -3.15 !3.25 -3.39 !3.50 -3.62 !3.75 -3.86
!4 -3.89                                        !4 -4.1


write(FMT1,'(F4.2)') eJ !F4.2
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,*)
write(*,*) 'J=', eJ

write(*,*)

write(*,*) 'Repetitions=', nrep
write(*,*)

if(nPROCESS==0) then
  open(50, file='Traces_100_Res_J='//FMT1//'_N=8.dat')
  
  write(*,*) 'Resurrection'
  write(*,*)
else
  open(50, file='Traces_100_Stall_J='//FMT1//'_N=8_fixed.dat')
  
  write(*,*) 'Release'
  write(*,*)
end if


do n=1,nrep
  
  if (nPROCESS==1) then
    
    call GTraces(eJi,eMui,rn0,n_v,n_sim,n_st,t_v,av_h,var_h,sd_h)
    
    sumn0=0.0d0
    do j=1,n_sim
      rn0(j)=t_v(j,n_st)
      
      sumn0 = sumn0 + rn0(j)      
    end do
    
  end if
   
  rmrn0=sumn0/real(n_sim)
  
  write(*,*) '<phi_0>=', rmrn0

  call GTraces(eJ,eMu,rn0,n_v,n_sim,n_st,t_v,av_h,var_h,sd_h)
 
  sumst=0.0d0
  do l=1,n_st
        ssp=0.0d0
        do j=1,n_sim
            ssp=ssp+t_v(j,l)
        end do
        av_step(n,l)=ssp/real(n_sim)
        
        sumst=sumst + av_step(n,l)
      end do
      
      total_av = sumst/real(n_st)
      
      !Real average of the curve
      !real_av_total=Average(100,n_st,av_step)
      
      !Standard deviation of average curve
      do l=1,n_st
        ssd=0.0d0
        do j=1,n_sim
            ssd=ssd+(t_v(j,l)-real_av_total)**2.0d0
        end do
        std_step(n,l)=sqrt(ssd/real(n_sim))
        
        !18*l
        write(50,*) n, l, av_step(n,l), std_step(n,l)
        
      end do
      
      write(*,*) n, '<phi>= ', total_av  

end do



!File closing
     
      close(50)


      write(*,*)
      write(*,*) 'Program finished'
      end 
