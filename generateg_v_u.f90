subroutine generateuq2
  use constants
  implicit none
  integer i,j
  real*8 q,impcoe,impci
  !  real*8 qs,qs_vacuum
  complex*16 screen
!  qs=4.d0*Ef/hbarv*rs
  !  qs_vacuum=4.d0*Ef/hbarv*rs_vacuum
  impcoe=-2.d0*pi**2*Ni*vf*rs**2/N
  impci=-2.d0*pi**2*Ci*vf*rs**2/N
  do i=1,M
     do j=0,N/2
        if(j.eq.0) then
           uq2(i,j)=0.d0
        else
           q=2.d0*kr(i)*dsin(0.5d0*j*deltathita)
           uq2(i,j)=(impci*exp(-2.d0*q*dci)+impcoe*exp(-2.d0*q*d))*kr(i)*(1.d0+dcos(j*deltathita))&
               /(q**2*(real(screen(q,0.d0))**2+dimag(screen(q,0.d0))**2))
        end if
     end do
 end do

end subroutine generateuq2


subroutine generaterlog
  use constants
  implicit none
  real*8 q,qs,coeg,costhita,betag,gtmp(2),lattice
  integer tlo,i,j,ii
  qs=4.d0*Ef/hbarv*rs
  coeg=-0.5d0*e**2/(epsilon0*hbar**2*N*vf)                                                   
  
  gtmp(1)=5.4d-3
  gtmp(2)=3.5d-2
  lattice=0.142d0

  do tlo=1,2
     numrlo(tlo)=1.d0/(exp(wrlo(tlo)/kT)-1.d0)
     !betag=beta(tlo)*wrlo(tlo)*coeg
     betag=-gtmp(tlo)*vf/(lattice*N)
     do i=1,M
        do j=0,N/2
           costhita=dcos(j*deltathita)
           if(i.le.nrlo(tlo)) then
              grlo(tlo,i,j,1)=0.d0
           else
              ii=i-nrlo(tlo)
              q=dsqrt(kr(i)**2+kr(ii)**2-2.d0*kr(i)*kr(ii)*costhita) 
              grlo(tlo,i,j,1)=betag*kr(ii)*exp(-2.d0*q*d)/(q+qs)*0.5d0*(1.d0+costhita)
           end if
           
           if((i+nrlo(tlo)).gt.M) then
              grlo(tlo,i,j,2)=0.d0
           else
              ii=i+nrlo(tlo)
              q=dsqrt(kr(i)**2+kr(ii)**2-2.d0*kr(i)*kr(ii)*costhita) 
              grlo(tlo,i,j,2)=betag*kr(ii)*exp(-2.d0*q*d)/(q+qs)*0.5d0*(1.d0+costhita)
           end if
        end do
     end do
  end do
      
end subroutine generaterlog


subroutine generatellog
  use constants
  implicit none
  integer tlo,i,j,ii
  real*8 coeg(2),costhita
  do tlo=1,2
     numllo(tlo)=1.d0/(exp(wllo(tlo)/kT)-1.d0)
     coeg(tlo)=-1.d0/(N*vf)*D2(tlo)/(tlo*density*wllo(tlo))
  end do
  do i=1,M
     if(i.le.nllo(1)) then
        gllo1(i,1)=0.d0
     else
        gllo1(i,1)=coeg(1)*kr(i-nllo(1))
     end if
           
     if((i+nllo(1)).gt.M) then
        gllo1(i,2)=0.d0
     else
        gllo1(i,2)=coeg(1)*kr(i+nllo(1))
     end if
  end do
  
  do j=0,N/2
     costhita=dcos(deltathita*j)
     do i=1,M
        if(i.le.nllo(2)) then
           gllo2(i,j,1)=0.d0
        else
           gllo2(i,j,1)=coeg(2)*kr(i-nllo(2))*(1.d0-costhita)
        end if
           
        if((i+nllo(2)).gt.M) then
           gllo2(i,j,2)=0.d0
        else
           gllo2(i,j,2)=coeg(2)*kr(i+nllo(2))*(1.d0-costhita)
        end if
     end do
  end do
      
end subroutine generatellog

 
subroutine generateacg
  use constants
  implicit none
  integer i,j,i1,i2,l
  real*8 vratio,k1,costhita,q
  vratio=vf/vac
  do i=1,M
     do j=1,N/2
        kacdown(i,j)=0
        kacup(i,j,1)=0
        kacup(i,j,2)=0
        gacd(i,j)=0.d0
        gacu(i,j,1)=0.d0
        gacu(i,j,2)=0.d0
     end do
  end do
  
  do j=1,N/2
      costhita=dcos(deltathita*j)
      do i=1,M
          k1=kr(i)*(vratio**2-costhita-dsqrt((vratio**2-costhita)**2-(vratio**2-1.d0)**2))&
              /(vratio**2-1.d0)
          
          if (int(k1/detak)+1-(shift-1).ge.1) then
              kacdown(i,j)=int(k1/detak)+1-(shift-1)
              k1=kr(kacdown(i,j))
              q=dsqrt(kr(i)**2+k1**2-2.d0*kr(i)*k1*costhita)
              gacd(i,j)=-0.25d0*k1*q*Dac**2/(N*density*vac*hbar)&
                  *1.d0/abs(vf+vac*(k1-kr(i)*costhita)/q)*(1.d0+costhita)
              if(kacup(kacdown(i,j),j,1).eq.0) then
                  kacup(kacdown(i,j),j,1)=i
              else if (kacup(kacdown(i,j),j,2).eq.0) then
                  kacup(kacdown(i,j),j,2)=i
              else 
                  pause 3
              end if
          end if
       
         
      end do
  end do
  
  do j=1,N/2
     costhita=dcos(deltathita*j)
     do i=1,M
        do l=1,2
           if(kacup(i,j,l).ne.0) then
              k1=kr(kacup(i,j,l))
              q=dsqrt(kr(i)**2+k1**2-2.d0*kr(i)*k1*costhita)
              gacu(i,j,l)=gacd(kacup(i,j,l),j)/kr(i)*k1
              ! print *,gacu(i,j,l),-0.25d0*k1*q*Dac**2/(N*density*vac*hbar)&
              !      *1.d0/abs(vf-vac*(k1-kr(i)*costhita)/q)*(1.d0+costhita)
           end if
        end do
     end do
  end do
end subroutine generateacg

subroutine generatev
  use constants
  implicit none
  integer i,j,i2
  real*8 q,vq0,omega
  complex*16 screen
  do i=1,M
     do i2=1,i
        v(i,i2,0)=cmplx(0.d0,0.d0)
        v(i2,i,0)=cmplx(0.d0,0.d0)
        v_omega0(i,i2,0)=0.d0
        v_omega0(i2,i,0)=0.d0
        do j=1,N/2
           q=dsqrt(kr(i)**2+kr(i2)**2-2.d0*kr(i)*kr(i2)*dcos(deltathita*j))
           omega=vf*(kr(i)-kr(i2))
           vq0=vq0coe/q
           v(i,i2,j)=vq0/screen(q,omega)
           
           v(i2,i,j)=conjg(v(i,i2,j))     
           v_omega0(i,i2,j)=vq0/screen(q,0.d0)
           
           v_omega0(i2,i,j)=v_omega0(i,i2,j)                                                       
        end do
     end do
 end do

  
end subroutine generatev


subroutine generatekm_j1j3_i1coe
  use constants
  implicit none
  integer i,j,i2,j2,convertfai,i1
  real*8 qx,qy,q,generatefai,y1,y2,y3,thita,k3x,k3y,k3
 
  do j=1,N
     do j2=1,N
        if(j.eq.j2) then
           goto 1
        end if
        do i=1,M
           do i2=1,M
               q=dsqrt(kr(i)**2+kr(i2)**2-2.d0*kr(i)*kr(i2)*dcos(deltathita*(j-j2)))
               qx=kr(i)*dcos(deltathita*j)-kr(i2)*dcos(deltathita*j2)
               qy=kr(i)*dsin(deltathita*j)-kr(i2)*dsin(deltathita*j2)
              
               y1=q
               y2=kr(i)-kr(i2)
               y3=generatefai(qx,qy,q)
              i1min_max(i,j,i2,j2,1)=min(max(1,int(0.5d0*(y1+y2)/detak+0.5d0)+1-(shift-1)),M)
              i1min_max(i,j,i2,j2,2)=min(M,M+i-i2)

              do i1=i1min_max(i,j,i2,j2,1),i1min_max(i,j,i2,j2,2)
                  if(abs((y1**2-y2**2+2.d0*kr(i1)*y2)/(2.d0*kr(i1)*y1)-1.d0).lt.1.d-10) then
                      thita=0.d0
                  else if(abs((y1**2-y2**2+2.d0*kr(i1)*y2)/(2.d0*kr(i1)*y1)+1.d0).lt.1.d-10) then
                      thita=pi
                  else
                      thita=dacos((y1**2-y2**2+2.d0*kr(i1)*y2)/(2.d0*kr(i1)*y1))
                  end if
                 j1j3(i,j,i2,j2,i1,1)=convertfai(y3+thita)
                 j1j3(i,j,i2,j2,i1,2)=convertfai(y3-thita)
                 k3x=kr(i1)*dcos(y3+thita)-qx
                 k3y=kr(i1)*dsin(y3+thita)-qy
                 k3=dsqrt(k3x**2+k3y**2)
                 j1j3(i,j,i2,j2,i1,3)=convertfai(generatefai(k3x,k3y,k3))
                 k3x=kr(i1)*dcos(y3-thita)-qx
                 k3y=kr(i1)*dsin(y3-thita)-qy
                 k3=dsqrt(k3x**2+k3y**2)
                 j1j3(i,j,i2,j2,i1,4)=convertfai(generatefai(k3x,k3y,k3))
                 i1coe(i,j,i2,j2,i1)=dsqrt(abs((4.d0*kr(i1)**2-4.d0*kr(i1)*y2+y2**2-y1**2)&
                     /(y1**2-y2**2)))
              end do
           end do
        end do
1       continue
     end do
 end do
 
end subroutine generatekm_j1j3_i1coe



real*8 function generatefai(qx,qy,q)
  use constants
  implicit none
  real*8 qx,qy,q,fai
  if(qy**2.ge.qx**2) then
     if(qy.gt.0.d0) then
        fai=dacos(qx/q)
     else 
        fai=2.d0*pi-dacos(qx/q)
     end if
  else
     if(qx.gt.0.d0) then
        fai=dasin(qy/q)
     else
        fai=pi-dasin(qy/q)
     end if
  end if
  generatefai=dmod(fai+2.d0*pi,2.d0*pi)
end function generatefai
      
       
integer function convertfai(fai)
  use constants
  implicit none
  real*8 fai
  integer j
  if(fai.lt.0.d0) then
     j=int(abs(fai)/pi*0.5d0)+1
     fai=fai+j*2.d0*pi
  end if
  fai=dmod(fai,2.d0*pi)
  convertfai=int(fai/deltathita+0.5d0)
  if(convertfai.eq.0) then
     convertfai=N
  end if
end function convertfai


complex*16 function screen(q,omega)
  use constants
  implicit none
  real*8 q,omega,vq0,screenreal,screenimag,scoe,r_integral,i_integral
  integral_q=q
  integral_omega=omega
  vq0=vq0coe/q
  scoe=q**2/pi_2hbar*1.d0/dsqrt(q**2*vf**2-omega**2)
  screenreal=1.d0-vq0*(screencoe*(0.5d0*(log(1.d0+exp(betau1))+log(1.d0+exp(-betau1))&
       +log(1.d0+exp(betau2))+log(1.d0+exp(-betau2))))+scoe*(r_integral()-0.5d0*pi))
  screenimag=-vq0*scoe*i_integral()
  screen=cmplx(screenreal,screenimag)
end function screen
  
real*8 function r_integral()
  use constants
  implicit none
  external r_kernel
  integer NEQ,IFLAG
  parameter(NEQ=1)
  real*8 Y(NEQ),TIN,TOUT,RELERR,ABSERR,WORK(6*NEQ+3),IWORK(5)
  TIN=-1.d0
  TOUT=1.d0
  Y(1)=0.D0
  RELERR=1.E-7
  ABSERR=1.E-7
  IFLAG=1
  call RKF45(r_kernel,NEQ,Y,TIN,TOUT,RELERR,ABSERR,IFLAG,WORK,IWORK)
  r_integral=Y(1)
end function r_integral

real*8 function i_integral()
  use constants
  implicit none
  external i_kernel
  integer NEQ,IFLAG
  real*8 init,integrate
  parameter(NEQ=1)
  real*8 Y(NEQ),TIN,TOUT,RELERR,ABSERR,WORK(6*NEQ+3),IWORK(5)
  TIN=1.d-10
  TOUT=1.d0
  Y(1)=0.D0
  RELERR=1.E-7
  ABSERR=1.E-7
  IFLAG=1
  call RKF45(i_kernel,NEQ,Y,TIN,TOUT,RELERR,ABSERR,IFLAG,WORK,IWORK)
  i_integral=Y(1)
end function i_integral

     
subroutine r_kernel(TIN,Y,YP)
  use constants
  implicit none
  real*8 TIN,Y(1),YP(1)
  YP(1)=0.5d0*(1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q*TIN-hbar*integral_omega)/kT-betau1))&
       +1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q*TIN-hbar*integral_omega)/kT+betau1))&
       +1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q*TIN-hbar*integral_omega)/kT-betau2))&
       +1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q*TIN-hbar*integral_omega)/kT+betau2)))&
       *dsqrt(1.d0-TIN**2)
end subroutine r_kernel
  

subroutine i_kernel(TIN,Y,YP)
  use constants
  implicit none
  real*8 TIN,Y(1),YP(1)
  YP(1)=0.5d0*(1.d0/(1.d0+exp(0.5d0*(hbarv*integral_q/TIN+hbar*integral_omega)/kT-betau1))&
       +1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q/TIN+hbar*integral_omega)/kT+betau1))&
       -1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q/TIN-hbar*integral_omega)/kT-betau1))&
       -1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q/TIN-hbar*integral_omega)/kT+betau1))&
       +1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q/TIN+hbar*integral_omega)/kT-betau2))&
       +1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q/TIN+hbar*integral_omega)/kT+betau2))&
       -1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q/TIN-hbar*integral_omega)/kT-betau2))&
       -1.d0/(1.d0+exp(0.5d0*abs(hbarv*integral_q/TIN-hbar*integral_omega)/kT+betau2)))&
       *dsqrt(TIN**(-2.d0)-1.d0)*TIN**(-2.d0)
end subroutine i_kernel

           
           
