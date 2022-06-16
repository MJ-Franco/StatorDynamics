!Subroutine for obtaining the traces of resurrection and stall of BFM 
!with cooperativity using KMC.

!24 january 2022

!María-José Franco Oñate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module Berg_CoopKMC
contains
	
subroutine BergCoopKMC(eJ,eMu,rk0,alpha,rn0,n_v,n_sim,n_st,tv_smpl, &
t_smpl,ntt)
	
	!eJ =  Value of J (cooperativity)
    !eMu =  total chemical potential of the motor (mu+ep)
    !eMu0 = constant of reservoir chemical potential 
    !ntot = total number of stators in the system
    !n_v = Vector size
    !n_sim = Number of simulations
    !n_st = Number of steps in each simulation
    !n0 =Initial number of stators in the rotor
      
implicit real*8 (a-h,o-z)

      
real*8 :: t_v(n_sim,n_st),phi2_v(n_sim,n_st),t2h_v(n_sim,n_st), &
      phi_v(n_sim,n_st),rn0(n_sim), rvec(6), w(6), t(n_sim,n_st)
      
integer :: nv(6,n_v), n(6)
         
real*8 :: s(n_v),av_h(n_sim),var_h(n_sim), av2_h(n_sim), &
sd_h(n_sim)
 
real*8, allocatable :: t_smpl(:,:), tv_smpl(:,:)

!allocatable :: cluster(:)

!character(len=7) :: FMT1

!character(len=6) :: FMT 
      
      !s(n_v) variable that denotes the state of occupation of a place in the rotor (0,1) 
      !t_v(n_sim,n_st) total number of stators in simulation m at step n
      !r_v(n_sim,n_st) total number of stators in the bulk in simulation at step n
      !av_h (n_sim) final average of each history
      !sd_h (n_sim) final sd of each history
      
!Parameters for random numbers

iseed = -123456789 
      
an= 1.0d0
bn= real(n_v)+1.0d0

	
      !File opening
      
!      write(FMT1,'(F5.2)') eMu !F4.2
!      write(FMT,'(F4.2)') eJ
      

      

!Simulation
     

      

do j=1,n_sim  !Do for several histories
     
  !Vector initialisation
		
  rn0_t = rn0(j)

  !Fill vectors with 0
   do i=1,n_v
     s(i)=0.0d0
   end do
        
        
  !If n0=0 then it keepts going
  if(int(rn0_t)==0) then
    continue
    !If n0=\0 it fills the vectors with n0 1s randomly    
    else    
      do k=1,int(rn0_t)    
        
5       tr=ran2(iseed)
        tt=an+(bn-an)*tr
        i=int(tt) 
        if(s(i)==0.0d0) then  
          s(i)=1.0d0
        else
          goto 5
        end if
      end do
    end if
    



!First step        
      t_v(j,1) = rn0(j)
      phi_v(j,1)=t_v(j,1)/real(n_v)
            
     !We write in a file this number
!     write(30,*) j, l, phi_v(j,1)
            
      t2h_v(j,1)=0.0d0
      do i=1,n_v
        do k=1,n_v
          t2h_v(j,1)= t2h_v(j,1) + s(i)*s(k)
        end do
      end do

      phi2_v(j,1) = t2h_v(j,1)/real(n_v)**2.0d0
        
    
    t(j,1) = 0.0d0
    
    
   w0 = 0.5d0*(1+dexp(-eMu))*rk0*(1.0d0 - dexp(-alpha/rn0(j)))
      
   !Rate for each transition        
    w(1) = w0*(1.0d0 - dtanh(-eMu/2.0d0))
    w(2) = w0*(1.0d0 - dtanh(-(eJ+eMu)/2.0d0))
    w(3) = w0*(1.0d0 - dtanh(-eJ-eMu/2.0d0)) 
    w(4) = w0*(1.0d0 - dtanh(eJ+eMu/2.0d0))
    w(5) = w0*(1.0d0 - dtanh((eJ+eMu)/2.0d0)) 
    w(6) = w0*(1.0d0 - dtanh(eMu/2.0d0))
    
    
    do l=2,n_st   !The simulation for one story begins
       
      do i=1,6
          n(i)=0
          do k=1,n_v
            nv(i,k)=0
          end do
      end do      
           
      do k=1,n_v
      

!Neighbors of each site

        if((k+1).gt.n_v) then 
          suma = s(k-1)+s(1)    
        elseif ((k-1).lt.1) then
          suma = s(k+1) + s(n_v)
        else
          suma = s(k+1)+s(k-1)
        end if
        
        
!Possible transitions     

        if(suma==0) then !There are no neighbors
        
          if(s(k)==0.0d0) then !Transition to 1
          
            n(1) = n(1) + 1
            
            nv(1,n(1)) = k 
        
          else  !Transition to 0
        
            n(6) = n(6) + 1
            
            nv(6,n(6)) = k
            
          end if
        
        else if(suma==1) then !There is 1 neighbor
       
          if(s(k)==0.0d0) then   !Transition to 1
        
            n(2) = n(2) + 1
            
            nv(2,n(2)) = k
        
          else   !Transition to 0
          
            n(5) = n(5) + 1
            
            nv(5,n(5)) = k
        
          end if
          
        else if(suma==2) then !There are 2 neighbors
        
          if(s(k)==0.0d0) then !Transition to 1
        
            n(3) = n(3) + 1
            
            nv(3,n(3)) = k
        
          else    !Transition to 0
        
            n(4) = n(4) + 1
            
            nv(4,n(4)) = k 
        
          end if               
      
        end if
        
        
      end do
      

        
      Rtot=0.0d0
      do i=1,6
        Rtot = Rtot + real(n(i))*w(i)
      end do  
      
      
      rvec(1) = real(n(1))*w(1)
      do i=2,6
        rvec(i) = rvec(i-1) + real(n(i))*w(i)
      end do
      
      
      
      r1 = ran2(iseed)*Rtot !random number to choose the process
      
      
      do i = 1,6
      

        if(rvec(i)>=r1) then

          m=i    !the process is chosen
		  
		  exit
		  
		  
		else
		  continue  
        end if
      
      end do
      
      
      if(s(i)==0) then
        
      end if
      
      r2 = ran2(iseed) !random number to choose site
      
      nprocess = idnint(n(m)*r2)
     
      
      if(nprocess==0) then
        nprocess=1
      end if
        
      
      nsite = nv(m,nprocess)
      
      
      
      
      if(m<=3) then
      
        s(nsite) = 1.0d0
        
      else
      
        s(nsite) = 0.0d0
        
      end if
      

     r3 = ran2(iseed) !random number for the time


     t(j,l) = t(j,l-1) - dlog(r3)/Rtot 
     
!!!Total occupation of the system!!!
     t_v(j,l) = 0.0d0
     do i=1,n_v     
       t_v(j,l) = t_v(j,l) + s(i)
     end do
     
     !write(*,*) j, l, t_v(j,l)
     
     phi_v(j,l)=t_v(j,l)/real(n_v)
     
!!!For the standard deviation!!!     
     t2h_v(j,l)=0.0d0
      do i=1,n_v
        do k=1,n_v
          t2h_v(j,l)= t2h_v(j,l) + s(i)*s(k)
        end do
      end do     
     
     
     w0 = 0.5d0*(1+dexp(-eMu))*rk0*(1.0d0 - dexp(-alpha/t_v(j,l)))
      
!Rate for each transition        
     w(1) = w0*(1.0d0 - dtanh(-eMu/2.0d0))
     w(2) = w0*(1.0d0 - dtanh(-(eJ+eMu)/2.0d0))
     w(3) = w0*(1.0d0 - dtanh(-eJ-eMu/2.0d0)) 
     w(4) = w0*(1.0d0 - dtanh(eJ+eMu/2.0d0))
     w(5) = w0*(1.0d0 - dtanh((eJ+eMu)/2.0d0)) 
     w(6) = w0*(1.0d0 - dtanh(eMu/2.0d0))
     
     
     
  end do
  
    
end do  


sumt=0.0d0
do j=1,n_sim
  sumt = sumt + t(j,n_st)
end do

tt=sumt/n_sim

ntt=dint(tt)
   

    
allocate(t_smpl(n_sim,ntt))
allocate(tv_smpl(n_sim,ntt))
   

do j=1,n_sim 
    l=1
    t_smpl(j,1)=0.0d0 
    tv_smpl(j,1)=t_v(j,1)
    do i=2,ntt
 
      t_smpl(j,i) = t_smpl(j,i-1) + 0.1d0
      
      if(t_smpl(j,i)>=t(j,l+1)) then
        l=l+1
      else
        continue
      end if
      
      tv_smpl(j,i) = t_v(j,l)
      
     ! write(*,*) j, i, tv_smpl(j,i)
      
    end do
    
 end do

!deallocate(t_smpl)
!deallocate(tv_smpl)
    
end subroutine
  
	
	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
       REAL*8 FUNCTION ran2(idum)
      
      IMPLICIT REAL*8 (A-H,O-Z)
!      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
!      REAL*8 AM,EPS,RNMX
      integer,parameter :: IM1=2147483563,IM2=2147483399,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, &
      IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB
      real*8,parameter :: AM=1.0d0/IM1,EPS=1.2d-7,RNMX=1.0d0-EPS
      
      INTEGER :: idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 function Average(it_min,it_max,dat)

      implicit real*8(a-h,o-z)
      dimension dat(int(it_max-it_min))
        
        sm=0.0d0
        do j=it_min,it_max
            sm=sm+dat(j)
        end do
        
        Average=sm/real(it_max-it_min)
        
      return  
      end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      real*8 function Av_2(it_min,it_max,dat)

      implicit real*8(a-h,o-z)
      dimension dat(int(it_max-it_min))
        
        sm=0.0d0
        do j=it_min,it_max
            sm=sm+dat(j)**2.0d0
        end do
        
        Av_2=sm/real(it_max-it_min)
        
      return  
      end 


end module
