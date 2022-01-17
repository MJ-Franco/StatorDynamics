
!Maria-Jose Franco Onate

!Program to obtain the relaxation times for each J depending on mu. 

!FEBRUARY 2021

include 'SubTraces_wo_Depletion.f90'

use TRACES
implicit real*8 (a-h,o-z)

integer, parameter :: n_v=100, n_sim=1000, n_st=2000, n0=0, n=100
!integer :: omp_rank
!Vector size; simulation number; step number; initial state

real*8 :: phi_v(n_sim,n_st),av_step(n,n_st),std_step(n,n_st), &
    av_h(n_sim),var_h(n_sim),sd_h(n_sim)

!t_v total number of stators in simulation m at step n
!av_step average of each step
!std_step standard devition of each step
!av_h final average of each history
!var_h variance of each similation
!sd_h final sd of each history

character*1 m

!Parameters of the system
eJ = 0.1d0
eMu = eJ - (2.0d0*eJ+3.0d0) !initial mu


!Opening files

J=eJ

write(m,'(i1)') J

if(n0==0) then
    open(40, file='Traces_Res_Avg_J='//m//'.dat')
else
    open(40, file='Traces_Stall_Avg_J='//m//'.dat')        
end if


!Loop for change mu
 

do i=1,n

    
	!mu changes one by one
    eMu = eMu + 30.0d0/real(n) 
   
    
    !call of the subroutine that produces the traces
    call tracesFM(eJ,eMu,n0,n_v,n_sim,n_st,phi_v,av_h,var_h,sd_h)
    
    !Average and std of each step
    do l=1,n_st
        
        sum_st=0.0d0
        sum_st2=0.0d0
        
        do j=1,n_sim
        
            sum_st = sum_st + phi_v(j,l)
            sum_st2 = sum_st2 + phi_v(j,l)**2.0d0
            
        end do
        
        av_step(i,l) = sum_st/real(n_sim)
        
        av2_step = sum_st2/real(n_sim)
        std_step(i,l) = sqrt(av2_step - av_step(i,l)**2.0d0)
        
        
    end do
    
end do

    
!    write(40,*) (i, i=1,n)

!write in a file the traces for each differnt mu (columns)
do l = 1,n_st
    
    write(40,*) l, (av_step(i,l),std_step(i,l),i=1,n)
		
end do

!Closing files
close(40)

write(*,*) 'Program finished'
write(*,*)


end program


