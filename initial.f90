subroutine writeparameter
  use constants
  integer i,j,jprime
  open(0,file='result.dat',access='append')
  write(0,*) '#T=',T,'K',',NN=' ,NN,'(nm**(-2))'
  write(0,*) '#p=',p
  write(0,*) '#pdthita=',pdthita,',pdphi=',pdphi
  write(0,*) '#B=',B,'(T)',',Bphi=',Bphi
  write(0,*) '#lambdaBR=',lambdaBR,'(1/ps)'
  write(0,*) '#M=',M,',N=',N,',nlo=',nlo
  write(0,*) '#aa=',aa
  write(0,*) '#Ez=',Ez,'(V/nm)'
  write(0,*) '#driftE=',driftE,'(V/nm)'
  write(0,*) '#deltat=',deltat,'(ps)'
  write(0,*) '#length=',len,'(nm)'
  write(0,*) '#xl=',xl
  write(0,*) '#The fermi energy is Ef(ev):',Ef
  write(0,*) '#The chemical potential is(ev):'
  write(0,*) '#u1=',u1,',u2=',u2
  write(0,*) '#u0=',u0
  write(0,*) '#nrlo=',nrlo 
  write(0,*) '#nllo=',nllo
  write(0,*) '#cutoff=',cutoff(1),cutoff(2)
  write(0,*) '#kT/detaR=',cutoff(3)
  write(0,*) '#hartree-fork=',h_f
  write(0,*) '#e_lo_remote1=',e_lo_remote1
  write(0,*) '#e_lo_remote2=',e_lo_remote2
  write(0,*) '#e_lo_local1=',e_lo_local1
  write(0,*) '#e_lo_local2=',e_lo_local2
  write(0,*) '#e_e_inter=',e_e_inter
  write(0,*) '#e_e_intra=',e_e_intra
  write(0,*) '#e_ac=',e_ac
  write(0,*) '#Ni=',Ni
  write(0,*) '#shift=',shift
  write(0,*) "#Starting time (ps):", starttime
  close(0)
!  open(30,file='distribution0.dat',access='append')
!  do i=1,M
!      do jprime=1,N+1
!          j=jprime
!          if(jprime.eq.N+1) then
!              j=1
!          end if
!          write(30,'(4e16.8)')  (i-0.5d0)*detak*dcos(deltathita*j),(i-0.5d0)*detak*dsin(deltathita*j),rou(i,j,1,0),rou(i,j,3,0)
          
!      end do
!      write(30,*)
!  end do
  
!  close(30)
end subroutine writeparameter
      
subroutine init
  use constants
  real*8 Efup,Efdown
  real*8 finefermienergy
  integer i,j
  if(mod(xl+1,CPU).ne.0) then
      if(myid.eq.0) then
          print *,'please choose xl or CPU again'
      end if
      stop
  else
     meanx=(xl+1)/CPU
  end if
  pi=dacos(-1.d0)
  hbarv=hbar*vf
  Ef=hbarv*dsqrt(pi*NN)
  Efup=Ef*dsqrt(1.d0+p)
  Efdown=Ef*dsqrt(1.d0-p)
  nrlo(1)=nlo
  detaR=wrlo(1)/nrlo(1)
  detak=detaR/hbarv
  deltathita=2.d0*pi/N
  deltal=len/xl
  pi_2hbar=2.d0*pi*hbar
  vq0coe=pi_2hbar*vf*rs
  coulombcoe=-detak**2/(8.d0*pi**2*N*hbar**2*vf)*(e_e_inter+e_e_intra)
  screencoe=-2.d0*kT/(pi*hbarv**2)
  s0=0.5d0/pi**2*detak*deltathita
  hartreecoe=s0/hbar
  nrlo(2)=nint(wrlo(2)/detaR)
  nllo(1)=nint(wllo(1)/detaR)
  nllo(2)=nint(wllo(2)/detaR)
  u1=Efup-pi**2*kT**2/(6.d0*Efup)
  u2=Efdown-pi**2*kT**2/(6.d0*Efdown)
  u0=Ef-pi**2*kT**2/(6.d0*Ef)
  u1=finefermienergy(u1,0.5d0*NN*(1.d0+p))
  u2=finefermienergy(u2,0.5d0*NN*(1.d0-p))
  u0=finefermienergy(u0,0.5d0*NN)
  betau0=u0/kT
  betau1=u1/kT
  betau2=u2/kT
  pdx=dsin(pdthita*pi)*dcos(pdphi*pi)
  pdy=dsin(pdthita*pi)*dsin(pdphi*pi)
  pdz=dcos(pdthita*pi)
 
  do j=1,N
      rnu(j)=1.d0-2.d0*deltat*vf/deltal*dabs(dcos(deltathita*j))
      rdo(j)=2.d0-rnu(j)
  end do
  
  shift=1
10 if(1.d0-1.d0/(exp(((shift-0.5d0)*detaR-u0)/kT)+1.d0).lt.1.d-5) then
      shift=shift+1
      goto 10
  end if
  shift=17 
  
  do i=1,M
     R(i)=(shift+i-1.5d0)*detaR
     kr(i)=(shift+i-1.5d0)*detak
  end do
end subroutine init

    
real*8 function finefermienergy(u,N0)
  use constants
  real*8 u,N0,utop,ubot,func,detau,utmp
  detau=50.d0*Ef
  ubot=u-detau
  utop=u+detau

1 if(abs(utop/ubot-1.d0).gt.1.d-14) then
      utmp=0.5d0*(ubot+utop)
      if(func(utmp,N0).lt.0d0) then
          ubot=utmp
      else
          utop=utmp
      end if
      goto 1
  end if
  finefermienergy=ubot
  return
end function finefermienergy


real*8 function func(u,N0)
   use constants
   real*8 u,N0
   integer i
   func=0.d0
   
   do i=1,scale*M
      func=func+1.d0/(exp(((i-0.5d0)*detaR-u)/kT)+1.d0)*(i-0.5d0)*detak
   end do
   func=func*detak/pi-N0
   return
 end function func
            
      
 subroutine rou_init
   use mpi
   use constants
   integer i,j,l,w,IO,head
   real*8 coea,coeb
   
   allocate(rou(M,N,4,-1:meanx))
   allocate(roupre(M,N,4,-1:meanx))
   allocate(scdrou(M,N,4,-1:meanx))
      
   do  i=1,M
       coea=0.5d0*(1.d0/(exp((R(i)-u1)/kT)+1.d0)+1.d0/(exp((R(i)-u2)/kT)+1.d0))
       coeb=0.5d0*(1.d0/(exp((R(i)-u1)/kT)+1.d0)-1.d0/(exp((R(i)-u2)/kT)+1.d0))
       do j=1,N
           do l=0,meanx-1
               rou(i,j,1,l)=1.d0/(exp((R(i)-u0)/kT)+1.d0)
               rou(i,j,2,l)=0.d0
               rou(i,j,3,l)=rou(i,j,1,l)
               rou(i,j,4,l)=0.d0
           end do
           if(myid.eq.0) then
               rou(i,j,1,0)=coea+pdz*coeb
               rou(i,j,2,0)=coeb*pdx
               rou(i,j,3,0)=coea-pdz*coeb
               rou(i,j,4,0)=-coeb*pdy
           end if
       end do
   end do
   
    if(myid.eq.0) then
       cutoff(1)=rou(M,N,1,0)/rou(1,N,1,0)
       cutoff(2)=rou(M,N,3,0)/rou(1,N,3,0)
       cutoff(3)=kT/detaR
   end if
   
   if(myid.eq.0) then
       call writeparameter
   end if

!   if(driftE.ne.0.d0) then
   if(myid.eq.0) then
   call evolution
   end if
       call mpi_barrier(mpi_comm_world,ierr)
       call mpi_bcast(rou(1,1,1,meanx-1),M*N*4,mpi_real8,CPU-1,mpi_comm_world,ierr)
       
       head=0    
       if(myid.eq.0) then
           head=1
       end if
       do l=head,meanx-2
           do i=1,M
               do j=1,N
                   do w=1,4
                       rou(i,j,w,l)=rou(i,j,w,meanx-1)
                   end do
               end do
           end do
       end do
 !  end if
   
   if(myid.eq.0) then
       cutoff(1)=rou(M,N,1,0)/rou(1,N,1,0)
       cutoff(2)=rou(M,N,3,0)/rou(1,N,3,0)
       cutoff(3)=kT/detaR
   end if
   
   if(myid.eq.0) then
       call writeparameter
   end if
   
 end subroutine rou_init
 
 subroutine evolution
  use constants
  integer tt,i,j,l,kk,ll,jprime
  real*8 vx,vxprime,kx,ky,time
  real*8 rou_tmp(M,N,4),localscdrou(M,N,4),roumid(4,M,N,4)
  real*8 rou_up,rou_down
  tt=0
  vx=0.d0
  if(myid.eq.0) then
      ll=0
  else if(myid.eq.CPU-1) then
      ll=meanx-1
  end if
   
100 tt=tt+1
  call generatedens
  call writeplotdat(tt)
  time=tt*deltat
     do i=1,M
        do j=1,N
           do l=1,4
              rou_tmp(i,j,l)=rou(i,j,l,ll)
           end do
        end do
     end do
  
     do kk=1,4
         do i=1,M
             do j=1,N
                 do l=1,4
                     localscdrou(i,j,l)=0.d0
                 end do
             end do
         end do
         call coherent(rou_tmp,localscdrou)
         call scat(rou_tmp,localscdrou)
         if(driftE.ne.0.d0) then
             call drift(rou_tmp,driftE,localscdrou)
         end if
         do i=1,M
           do j=1,N
               do l=1,4
                   roumid(kk,i,j,l)=localscdrou(i,j,l)
                   rou_tmp(i,j,l)=int((kk+1)/2)*0.5d0*deltat*roumid(kk,i,j,l)+rou(i,j,l,ll)
              end do
           end do
       end do
      
     end do
     
         
     do i=1,M
        do j=1,N
           do l=1,4
              rou(i,j,l,ll)=rou(i,j,l,ll)+1.d0/6.d0*deltat*(roumid(1,i,j,l)&
                   +2.d0*roumid(2,i,j,l)+2.d0*roumid(3,i,j,l)+roumid(4,i,j,l))
           end do
        end do
     end do
     
     goto 100

     if(mod(tt-1,20).eq.0) then

     if(myid.eq.0) then
         call generatedens
         open(100,file='relaxation.dat',access='append')
         !      write(100,'(4e20.12)') time,fpm(0),NN,fpola(0)
         write(100,'(4e20.12)') time,NN,fpm(0),fpola(0)
         close(100)
     end if

     vxprime=0.d0
        do i=1,M
            do j=1,N
                          
                rou_up=0.5d0*(rou(i,j,1,ll)+rou(i,j,3,ll)+(pdz*(rou(i,j,1,ll)-rou(i,j,3,ll))&
                    +2.d0*(pdx*rou(i,j,2,ll)-pdy*rou(i,j,4,ll))))
                rou_down=0.5d0*(rou(i,j,1,ll)+rou(i,j,3,ll)-(pdz*(rou(i,j,1,ll)-rou(i,j,3,ll))&
                    +2.d0*(pdx*rou(i,j,2,ll)-pdy*rou(i,j,4,ll))))
                
             if(rou_up.lt.0.d0.or.rou_up.gt.1.d0.or.rou_down.lt.0.d0.or.rou_down.gt.1.d0) then
                 print *,'error'
                pause 0
           end if
           vxprime=vxprime+kr(i)*dcos(deltathita*j)*(rou(i,j,1,ll)+rou(i,j,3,ll))
        end do
     end do
     
!     if(dabs(vx/vxprime-1.d0).lt.1.d-7) then
          if(myid.eq.0) then
              open(unit=30+myid,file='distribution0.dat')
              open(unit=40+myid,file='mu0.dat',access='append')
          else 
              open(unit=30+myid,file='distributionxl.dat')
              open(unit=40+myid,file='muxl.dat',access='append')
          end if
          write(40+myid,'(3e16.8)')  time,vxprime*vf*s0/(driftE*NN)*1.d-2,vxprime*vf*s0/(driftE*NN)*1.d-2*sqrt(pi*NN)*hbarv*0.5d0*1.d-4
!write(40+myid,'(3e16.8)')  time,vxprime*vf*s0,vxprime*vf*s0*1.d-2*sqrt(pi*NN)*hbarv*0.5d0*1.d-4
!'cm^2/Vs'
          do i=1,M
              do jprime=1,N+1
                  j=jprime
                  if(jprime.eq.N+1) then
                     j=1
                 end if
                write(30+myid,'(4e16.8)')  (i-0.5d0)*detak*dcos(deltathita*j),(i-0.5d0)*detak*dsin(deltathita*j),rou(i,j,1,0),rou(i,j,3,0)
             end do
         end do
         write(30+myid,*) 
         close(30+myid)
         close(40+myid)
!         stop
!     else
!         if(myid.eq.0) then
!             print *,myid,tt*deltat
!         else if(myid.eq.CPU-1) then
!             print *,'                              ',myid,tt*deltat
!         end if
         !         vx=vxprime
     end if
         goto 100
!     end if
!     end if
end subroutine evolution 
 
subroutine generatedens
  use constants
  integer l,i,j
  real*8 rou_up,rou_down
  do l=0,0!meanx-1
      fpm(myid*meanx+l)=0.d0
      fpola(myid*meanx+l)=0.d0
      do i=1,M
          do j=1,N
              rou_up=0.5d0*(rou(i,j,1,l)+rou(i,j,3,l)+(pdz*(rou(i,j,1,l)-rou(i,j,3,l))&
                  +2.d0*(pdx*rou(i,j,2,l)-pdy*rou(i,j,4,l))))
              rou_down=0.5d0*(rou(i,j,1,l)+rou(i,j,3,l)-(pdz*(rou(i,j,1,l)-rou(i,j,3,l))&
                  +2.d0*(pdx*rou(i,j,2,l)-pdy*rou(i,j,4,l))))
              if(rou_up.lt.0.d0.or.rou_up.gt.1.d0.or.rou_down.lt.0.d0.or.rou_down.gt.1.d0) then
                  print *,'error'
                  pause 1
              end if
              fpm(myid*meanx+l)=fpm(myid*meanx+l)+(rou_up+rou_down)*kr(i)
              fpola(myid*meanx+l)=fpola(myid*meanx+l)+(rou_up-rou_down)*kr(i)
          end do
      end do
      fpm(myid*meanx+l)=fpm(myid*meanx+l)*s0
      fpola(myid*meanx+l)=fpola(myid*meanx+l)*s0
  end do
end subroutine generatedens


subroutine writeplotdat(tt)
   use constants
   integer tt,l
   real*8 time
   time=tt*deltat
   open(unit=0,file='result.dat',access='append')
   do l=0,0!xl
       !  write(0,'(5e16.8)') time,l*len/xl,fpm(l),fpola(l),El(l)
       write(0,'(3e20.12)') time,fpm(l),fpola(l)!,El(l)
   end do
!   write(0,*)
   close(0)
   print *,time
 end subroutine writeplotdat
  
