subroutine scat(rou_tmp,localscdrou)
  use constants
  implicit none
  integer i,j,l
  real*8 scat_eim(4),scat_erlo(4),scat_ello(4),scat_eac(4),scat_ee(4),m2(M,N,M,N)
  real*8 rou_tmp(M,N,4),localscdrou(M,N,4)
  
  if(e_e_intra+e_e_inter.ne.0) then
     call generatem2(rou_tmp,m2)
  end if
! !$omp parallel do schedule(dynamic) private(scat_eim,scat_erlo,scat_ello,scat_eac,scat_ee)  
  do i=1,M
     do j=1,N
        do l=1,4
           scat_eim(l)=0.d0
           scat_erlo(l)=0.d0
           scat_ello(l)=0.d0
           scat_eac(l)=0.d0
           scat_ee(l)=0.d0
        end do
        if(Ni+Ci.ne.0.d0) then
           call scat_e_im(rou_tmp,i,j,scat_eim)
        end if
        if(e_lo_remote1+e_lo_remote2.ne.0) then
           call scat_e_rlo(rou_tmp,i,j,scat_erlo)
        end if
        if(e_lo_local1+e_lo_local2.ne.0) then
           call scat_e_llo(rou_tmp,i,j,scat_ello)
        end if
        if(e_ac.eq.1) then
           call scat_e_ac(rou_tmp,i,j,scat_eac)
        end if
        if(e_e_inter+e_e_intra.ne.0) then
           call scat_e_e(rou_tmp,m2,i,j,scat_ee)
        end if
        do l=1,4
            localscdrou(i,j,l)=localscdrou(i,j,l)+scat_eim(l)+scat_erlo(l)+scat_ello(l)+scat_eac(l)+scat_ee(l)
        end do
if(j==1)then            
!print*,scat_ello(1),scat_ello(2),scat_ello(3),scat_ello(4),i,j                 
!pause   
end if
     end do
 end do
! !$omp end parallel do
 end subroutine scat

subroutine scat_e_im(rou_tmp,i,j,scat_eim)
  use constants
  implicit none
  integer i,j,j1,jj1,l
  real*8 su2,scat_eim(4),rou_tmp(M,N,4)
  do j1=1,N
     if(abs(j-j1).le.N/2) then
        jj1=abs(j-j1)
     else
        jj1=N-abs(j-j1)
     end if
     su2=uq2(i,jj1)
     do l=1,4
        scat_eim(l)=scat_eim(l)+su2*(rou_tmp(i,j,l)-rou_tmp(i,j1,l))
     end do
  end do
end subroutine scat_e_im



subroutine scat_e_rlo(rou_tmp,i,j,scat_erlo)
  use constants
  implicit none
  real*8 gdu,Nump,scat_erlo(4),rou_tmp(M,N,4)
  integer tlo,i,j,j1,jj1,l,tlo1,tlo2
  real*8 emi_abs 
  !e_lo_remote1+e_lo_remote2.ne.0
  if(e_lo_remote1.eq.1) then
     tlo1=1
  else
     tlo1=2
  end if
  
  if(e_lo_remote2.eq.1) then
     tlo2=2
  else
     tlo2=1
  end if
  
  do tlo=tlo1,tlo2
     Nump=numrlo(tlo)
     do j1=1,N
        if(abs(j-j1).le.N/2) then
           jj1=abs(j-j1)
        else
           jj1=N-abs(j-j1)
        end if
       
        if(i.gt.nrlo(tlo)) then
           gdu=grlo(tlo,i,jj1,1)
           do l=1,4
              scat_erlo(l)=scat_erlo(l)+gdu*emi_abs(rou_tmp,Nump,i,j,i-nrlo(tlo),j1,l)
           end do
        end if
        
        if((i+nrlo(tlo)).le.M) then
           gdu=grlo(tlo,i,jj1,2)
           do l=1,4
              scat_erlo(l)=scat_erlo(l)-gdu*emi_abs(rou_tmp,Nump,i+nrlo(tlo),j1,i,j,l)
           end do
        end if
        
     end do
  end do
end subroutine scat_e_rlo


subroutine scat_e_llo(rou_tmp,i,j,scat_ello)
  use constants
  implicit none
  real*8 gdu,Nump,scat_ello(4),rou_tmp(M,N,4)
  integer tlo,i,j,j1,jj1,l,tlo1,tlo2
  real*8 emi_abs 
  if(e_lo_local1.eq.1) then
     tlo1=1
  else
     tlo1=2
  end if
  
  if(e_lo_local2.eq.1) then
     tlo2=2
  else
     tlo2=1
  end if
  
  do tlo=tlo1,tlo2
     Nump=numllo(tlo)
     do j1=1,N
        if(abs(j-j1).le.N/2) then
           jj1=abs(j-j1)
        else
           jj1=N-abs(j-j1)
        end if
        
        if(i.gt.nllo(tlo)) then
           if(tlo.eq.1) then
              gdu=gllo1(i,1)
           else
              gdu=gllo2(i,jj1,1)
           end if
           do l=1,4
              scat_ello(l)=scat_ello(l)+gdu*emi_abs(rou_tmp,Nump,i,j,i-nllo(tlo),j1,l)
           end do
        end if
        
        if((i+nllo(tlo)).le.M) then
           if(tlo.eq.1) then
              gdu=gllo1(i,2)
           else
              gdu=gllo2(i,jj1,2)
           end if
           do l=1,4
              scat_ello(l)=scat_ello(l)-gdu*emi_abs(rou_tmp,Nump,i+nllo(tlo),j1,i,j,l)
           end do
        end if
        
     end do
  end do
end subroutine scat_e_llo


subroutine scat_e_ac(rou_tmp,i,j,scat_eac)
  use constants
  implicit none
  real*8 gdu,Nump,scat_eac(4),rou_tmp(M,N,4),q
  integer i,j,j1,jj1,i1,l,ll
  real*8 emi_abs 
  
  do j1=1,N
     if(abs(j-j1).le.N/2) then
        jj1=abs(j-j1)
     else
        jj1=N-abs(j-j1)
     end if
     if(jj1.ne.0) then
         i1=kacdown(i,jj1)
         if(i1.ge.1.and.i1.le.M) then
             gdu=gacd(i,jj1)
             q=dsqrt(kr(i)**2+kr(i1)**2-2.d0*kr(i)*kr(i1)*dcos(deltathita*jj1))
             Nump=1.d0/(exp(hbar*vac*q/kT)-1.d0)
             do l=1,4
                 scat_eac(l)=scat_eac(l)+gdu*emi_abs(rou_tmp,Nump,i,j,i1,j1,l)
             end do
         end if
     
        do ll=1,2
           i1=kacup(i,jj1,ll)
           if(i1.ge.1.and.i1.le.M) then
              gdu=gacu(i,jj1,ll)
              q=dsqrt(kr(i)**2+kr(i1)**2-2.d0*kr(i)*kr(i1)*dcos(deltathita*jj1))
              Nump=1.d0/(exp(hbar*vac*q/kT)-1.d0)
              do l=1,4
                 scat_eac(l)=scat_eac(l)-gdu*emi_abs(rou_tmp,Nump,i1,j1,i,j,l)
              end do
           end if
        end do
     end if
  end do
end subroutine scat_e_ac


real*8 function emi_abs(rou_tmp,Nump,i,j,i1,j1,l)
  use constants
  implicit none
  integer i,j,i1,j1,l
  real*8 Nump,Nump1,rou_tmp(M,N,4)
  Nump1=Nump+1.d0
  if(l.eq.1) then
     emi_abs=Nump1*rou_tmp(i,j,1)-Nump*rou_tmp(i1,j1,1)&
          -rou_tmp(i,j,1)*rou_tmp(i1,j1,1)&
          -(rou_tmp(i,j,2)*rou_tmp(i1,j1,2)&
          +rou_tmp(i,j,4)*rou_tmp(i1,j1,4))
  else if(l.eq.2) then
     emi_abs=-(0.5d0*(rou_tmp(i,j,1)+rou_tmp(i,j,3))&
          *rou_tmp(i1,j1,2)+0.5d0*(rou_tmp(i1,j1,1)&
          +rou_tmp(i1,j1,3)-2.d0)*rou_tmp(i,j,2)&
          -Nump*(rou_tmp(i,j,2)-rou_tmp(i1,j1,2)))
  else if(l.eq.3) then
     emi_abs=Nump1*rou_tmp(i,j,3)-Nump*rou_tmp(i1,j1,3)&
          -rou_tmp(i,j,3)*rou_tmp(i1,j1,3)&
          -(rou_tmp(i,j,2)*rou_tmp(i1,j1,2)&
          +rou_tmp(i,j,4)*rou_tmp(i1,j1,4))
  else if (l.eq.4) then
     emi_abs=-(0.5d0*(rou_tmp(i,j,1)+rou_tmp(i,j,3))&
          *rou_tmp(i1,j1,4)+0.5d0*(rou_tmp(i1,j1,1)&
          +rou_tmp(i1,j1,3)-2.d0)*rou_tmp(i,j,4)&
          -Nump*(rou_tmp(i,j,4)-rou_tmp(i1,j1,4)))
  end if
end function emi_abs
      
subroutine scat_e_e(rou_tmp,m2,i,j,scat_ee)
  use constants
  implicit none
  integer i,j,i2,j2,jj2,y1,y2,y31,y32
  real*8 vq2,scat_ee(4),rou_tmp(M,N,4),m2(M,N,M,N)
  do j2=1,N
     if(abs(j-j2).le.N/2) then
        jj2=abs(j-j2)
     else
        jj2=N-abs(j-j2)
     end if
     do i2=1,M
        vq2=(real(v(i,i2,jj2))**2+dimag(v(i,i2,jj2))**2)*0.5d0*(1.d0+dcos(deltathita*jj2))*kr(i2)
        scat_ee(1)=scat_ee(1)+vq2*(m2(i,j,i2,j2)*rou_tmp(i,j,1)-m2(i2,j2,i,j)*rou_tmp(i2,j2,1)&
             +(m2(i2,j2,i,j)-m2(i,j,i2,j2))*(rou_tmp(i,j,1)*rou_tmp(i2,j2,1)&
             +rou_tmp(i,j,2)*rou_tmp(i2,j2,2)+rou_tmp(i,j,4)*rou_tmp(i2,j2,4)))
        scat_ee(2)=scat_ee(2)+vq2*(m2(i,j,i2,j2)*rou_tmp(i,j,2)-m2(i2,j2,i,j)*rou_tmp(i2,j2,2)&
             +(m2(i2,j2,i,j)-m2(i,j,i2,j2))*0.5d0*(rou_tmp(i,j,1)*rou_tmp(i2,j2,2)&
             +rou_tmp(i2,j2,1)*rou_tmp(i,j,2)+rou_tmp(i,j,3)*rou_tmp(i2,j2,2)&
             +rou_tmp(i2,j2,3)*rou_tmp(i,j,2)))
        scat_ee(3)=scat_ee(3)+vq2*(m2(i,j,i2,j2)*rou_tmp(i,j,3)-m2(i2,j2,i,j)*rou_tmp(i2,j2,3)&
             +(m2(i2,j2,i,j)-m2(i,j,i2,j2))*(rou_tmp(i,j,3)*rou_tmp(i2,j2,3)&
             +rou_tmp(i,j,2)*rou_tmp(i2,j2,2)+rou_tmp(i,j,4)*rou_tmp(i2,j2,4)))
        scat_ee(4)=scat_ee(4)+vq2*(m2(i,j,i2,j2)*rou_tmp(i,j,4)-m2(i2,j2,i,j)*rou_tmp(i2,j2,4)&
             +(m2(i2,j2,i,j)-m2(i,j,i2,j2))*0.5d0*(rou_tmp(i,j,1)*rou_tmp(i2,j2,4)&
             +rou_tmp(i2,j2,1)*rou_tmp(i,j,4)+rou_tmp(i,j,3)*rou_tmp(i2,j2,4)&
             +rou_tmp(i2,j2,3)*rou_tmp(i,j,4)))
     end do
  end do
 
  scat_ee(1)=coulombcoe*scat_ee(1)
  scat_ee(2)=coulombcoe*scat_ee(2)
  scat_ee(3)=coulombcoe*scat_ee(3)
  scat_ee(4)=coulombcoe*scat_ee(4)
end subroutine scat_e_e



subroutine generatem2(rou_tmp,m2)
  use constants
  implicit none
  integer i,j,i2,j2,i1,j11,j12,i3,j31,j32,i1min,i1max
  real*8 m2init,i1coeffi,rou_tmp(M,N,4),m2(M,N,M,N)
  do i=1,M
     do i2=1,M
        do j=1,N
           do j2=1,N
              m2(i,j,i2,j2)=0.d0
           end do
        end do
     end do
  end do
  
  do j=1,N
     do j2=1,N
        if(j.eq.j2) then
           goto 2
        end if
        do i=1,M
           do i2=1,M
              i1min=i1min_max(i,j,i2,j2,1)
              i1max=i1min_max(i,j,i2,j2,2)
              m2init=0.d0
              do i1=i1min,i1max
                 i3=i1-(i-i2)
                 j11=j1j3(i,j,i2,j2,i1,1)
                 j12=j1j3(i,j,i2,j2,i1,2)
                 j31=j1j3(i,j,i2,j2,i1,3)
                 j32=j1j3(i,j,i2,j2,i1,4)
                 i1coeffi=i1coe(i,j,i2,j2,i1)
                 m2init=m2init+i1coeffi*&
                      (rou_tmp(i3,j31,1)+rou_tmp(i3,j31,3)-rou_tmp(i1,j11,1)*rou_tmp(i3,j31,1)&
                      -rou_tmp(i1,j11,3)*rou_tmp(i3,j31,3)-2.d0*(rou_tmp(i1,j11,2)*rou_tmp(i3,j31,2)&
                      +rou_tmp(i1,j11,4)*rou_tmp(i3,j31,4))&
                      +rou_tmp(i3,j32,1)+rou_tmp(i3,j32,3)-rou_tmp(i1,j12,1)*rou_tmp(i3,j32,1)&
                      -rou_tmp(i1,j12,3)*rou_tmp(i3,j32,3)-2.d0*(rou_tmp(i1,j12,2)*rou_tmp(i3,j32,2)&
                      +rou_tmp(i1,j12,4)*rou_tmp(i3,j32,4)))
              end do
              m2(i,j,i2,j2)=m2init
           end do
        end do
2       continue
     end do
  end do
end subroutine generatem2
