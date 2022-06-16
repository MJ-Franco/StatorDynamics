module KMC_T
contains


subroutine KMC_Traces(rk0,alpha,eMu,rn0,n_v,n_sim,n_st,N_smpl, &
t_smpl,ntt)




implicit real*8 (a-h,o-z)

integer :: N(n_sim,n_st)

integer, allocatable :: N_smpl(:,:)

real*8 :: t(n_sim,n_st),av_h(n_sim)

real*8 :: rn0(n_sim), rvec(2)

real*8, allocatable :: t_smpl(:,:)


write(*,*) 'Kinetic Monte Carlo'

write(*,*)

!Parameters for the random numbers
iseed = -123456789    
      

do j=1,n_sim
  do l=1,n_st
    N(j,l) = 0
  end do
end do


do j=1,n_sim ! Begin of one story

  !Vector initialisation
     
     
  N(j,1) = dint(rn0(j)) !Initial occupancy of stators
      
   
     
  do l=1,n_st
    t(j,l) = 0.0d0   !Time initialisation
  end do
     
     
  do l = 1,n_st-1
  
    if (N(j,l)==0) then
      rup=rk0*(n_v-N(j,l))
      rdw=rk0*dexp(-eMu)*N(j,l)
    else
      rup = rk0*(1-dexp(-alpha/N(j,l)))*(n_v-N(j,l))
      rdw = rk0*(1-dexp(-alpha/N(j,l)))*dexp(-eMu)*N(j,l)
    end if
      
    rtot = rup + rdw 
    
      !Selection of process
    !rvec=0.0d0
    rvec(1)=rdw/rtot
    rvec(2)=(rdw+rup)/rtot
   
      
    ra1 = ran2(iseed)
      
    
    if(rvec(1)>ra1) then
    
      N_n=-1 
    
    else if(rvec(2)>ra1) then
    
      N_n=1
      
    end if
    
    
     
      N(j,l+1) = N(j,l) + N_n
    
      ra2 = ran2(iseed)
    
      t(j,l+1) = t(j,l) - dlog(ra2)/rtot
      
      
    end do
    
    av=0.0d0
    do l =100,n_st
      av = av + N(j,l)
    end do
    
    av_h(j) = av/n_st
    
    !write(*,*) j, 't_total =', t(j,n_st)
    
   
end do  !End of one story


    av_tot = Average(1,n_sim,av_h)

    write(*,*) 'N=', av_tot
    
    write(*,*)
  
    

!SAMPLED TIME AND SAMPLED N FOR THE MEAN TRACE

sumt=0.0d0
do j=1,n_sim
  sumt = sumt + t(j,n_st)
end do

tt=sumt/n_sim

ntt=dint(tt)
   
    
allocate(t_smpl(n_sim,ntt))
allocate(N_smpl(n_sim,ntt))
   
do j=1,n_sim 
    l=1
    t_smpl(j,1)=0.0d0 
    N_smpl(j,1)=N(j,1)
    do i=2,ntt
 
      t_smpl(j,i) = t_smpl(j,i-1) + 0.5d0
      
      if(t_smpl(j,i)>=t(j,l+1)) then
        l=l+1
      else
        continue
      end if
      
      N_smpl(j,i) = N(j,l)
      
    end do
 end do

    
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
