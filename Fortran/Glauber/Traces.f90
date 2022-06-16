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

integer, parameter :: nPROCESS=1

integer, parameter :: n_v=13, n_sim=10000, n_st=20000

real*8 :: t_v(n_sim,n_st),av_h(n_sim),var_h(n_sim),av_step(n_st), &
std_step(n_st),sd_h(n_sim),rn0(n_sim)
	
character(len=4) :: FMT
character(len=4) :: FMT1
	!t_v(n_sim,n_st) total number of stators in simulation m in step n
    !nr_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
    !av_h (n_sim) final average of each history
    !sd_h (n_sim) final sd of each history
    !av_step(n_st) average number of stators in each step
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!Constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

eJ = 5.00d0   !Value of J (cooperativity)
eMu = -5.00d0  !Total chemical potential of the motor (mu-ep)

!initial conditions
eJi = 5.00d0 
eMui = -4.91d0
    

exN_inf =  0.50d0

write(*,*)
write(*,*)

write(*,*) 'J=', eJ
write(*,*) 'N_inf=', exN_inf

write(*,*)


!0 0.7    !0.25 0.37  !0.5 0.05 !0.75 -0.25     !0 -0.69 !0.25 -0.85 !0.5 -1.03 !0.75 -1.23
!1 -0.57  !1.25 -0.87 !1.50 -1.15 !1.75 -1.45   !1 -1.41 !1.25 -1.62 !1.50 -1.83 !1.75 -2.04
!2 -1.74  !2.25 -2.01 !2.50 -2.29 !2.75 -2.57   !2 -2.25 !2.25 -2.47 !2.50 -2.7 !2.75 -2.92
!3 -2.84  !3.25 -3.11 !3.50 -3.37 !3.75 -3.64   !3 -3.15 !3.25 -3.39 !3.50 -3.62 !3.75 -3.86
!4 -3.89                                        !4 -4.1

write(FMT1,'(F4.2)') eJ !F4.2
write(FMT,'(F4.2)') exN_inf !F4.2

!File opening
if(nPROCESS==0) then
  open(50, file='Traces_Res_J='//FMT1//'_N='//FMT//'.dat')
  !open(70, file='All_Traces_Res_J='//FMT1//'_N='//FMT//'.dat')
        
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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Trajectories calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!open(70,file='Initial_conditions_J='//FMT1//'_N=13.dat')

do j=1,n_sim
   rn0(j)=0.00d0/13.0d0
end do

if (nPROCESS==1) then
    
    call GTraces(eJi,eMui,rn0,n_v,n_sim,n_st,t_v,av_h,var_h,sd_h)

end if

sumn0=0.0d0
do j=1,n_sim
  rn0(j)=t_v(j,n_st)
      
   sumn0 = sumn0 + rn0(j)      
end do
    
   
rmrn0=sumn0/real(n_sim)
  
write(*,*) '<phi_0>=', rmrn0
  

call GTraces(eJ,eMu,rn0,n_v,n_sim,n_st,t_v,av_h,var_h,sd_h)


do j=1,n_sim
  do l=1,n_st
    write(70,*) j, l, t_v(j,l)
  end do
end do
   
     
!Average and sd of each step

sumst=0.0d0
do l=1,n_st
  ssp=0.0d0
  
  do j=1,n_sim
    ssp=ssp+t_v(j,l)
  end do
  
  av_step(l)=ssp/real(n_sim)
        
  sumst=sumst + av_step(l)
end do
      
total_av = sumst/real(n_st)


!Standard deviation of average curve

sumtd=0.0d0
do l=1,n_st
  ssd=0.0d0
  
  do j=1,n_sim
    ssd=ssd+(t_v(j,l)-real_av_total)**2.0d0
  end do
  
  std_step(l)=sqrt(ssd/real(n_sim))
  
  sumtd=sumtd+std_step(l)
        
  write(50,*) l, av_step(l), std_step(l)
end do  

total_std = sumtd/real(n_st)
      
      
write(*,*) '<phi>= ', total_av, '+-', total_std  
      
write(*,*) 
      

!File closing
close(50)
!close(70)


      
write(*,*) 'Program finished'
 
end 
