module constants
  implicit none
  real*8, parameter::aa=300.d0 !nm
  real*8, parameter::len=10.d3  !nm   
  integer,parameter::xl=149
  integer, parameter::nlo=13
  integer, parameter::M=43
  integer, parameter::scale=50
  integer, parameter::N=12
  real*8, parameter::T=50.d0       !K
  real*8, parameter::p=0.05d0
  real*8, parameter::Vgate=200.d0/7.2d0  !V
  real*8, parameter::Ez=7.d-2!Vgate/aa  !V/nm
!  real*8, parameter::Curv=0.2d0 !K
  real*8, parameter::alpha=7.2d-4 !/(1/nm**2*V)
  real*8, parameter::NN=alpha*Vgate   !1/nm^2
!  real*8, parameter::NN=7.d-3       !1/nm^2
  real*8, parameter::Ni=0.d-3       !1/nm^2
  real*8, parameter::d=0.4d0       !nm
  real*8, parameter::Ci=0.d-3 !1/nm^2  !40.d-3
  real*8, parameter::dci=0.2d0       !nm
  real*8, parameter::B=0.d0        !T
  real*8, parameter::Bphi=0.d0
  real*8, parameter::deltat=0.025d0 !ps
  real*8, parameter::driftE=0.d-4  !V/nm
  real*8, parameter::coeBR=5.d-6   !eV/(V/nm)
  real*8, parameter::cacoe=0.d-5   !eV
  real*8, parameter::pdthita=0d0!direction of polarization;dimension:pi
  real*8, parameter::pdphi=0d0
  integer, parameter::A=10000000                                                                                 
  integer, parameter::h_f=0
  integer, parameter::e_lo_remote1=0
  integer, parameter::e_lo_remote2=0
  integer, parameter::e_lo_local1=0!intravalley
  integer, parameter::e_lo_local2=0!intervalley
  integer, parameter::e_e_inter=1
  integer, parameter::e_e_intra=1
  integer, parameter::e_ac=0
  
    
  !------------------------------------------------------------------------------  
  integer, parameter::e=1
  real*8, parameter::hbar=6.582122d-4,vf=1.d3,k=8.617385d-5,uB=5.78838263d-5,gfactor=2.d0,&
       rs=0.8d0,epsilon0=8.8541878d-2/1.602177d0,density=7.6d0/1.602177d0,Dac=19.d0,vac=20.d0,kT=k*T
  real*8, parameter::lambdaBR=(coeBR*Ez+cacoe)/hbar
  !  real*8, parameter::lambdaBR=(0.34d0*k)/hbar
  !  real*8, parameter::lambdaBR=Ez*coeBR/hbar
  real*8, parameter::wrlo(2)=[0.059d0,0.155d0],beta(2)=[0.025d0,0.062d0],wllo(2)=[0.196d0,&
       0.161d0],D2(2)=[4560.d0,9205.d0]
  
  real*8 pi,Ef,hbarv,s0,u1,u2,u0,detaR,detak,deltathita,deltal,coulombcoe,vq0coe,screencoe,&
       pi_2hbar,betau0,betau1,betau2,hartreecoe,pdx,pdy,pdz
  integer nrlo(2),nllo(2)
  
  real*8,allocatable::rou(:,:,:,:),roupre(:,:,:,:),scdrou(:,:,:,:)
  
  real*8 R(M),kr(M),grlo(2,M,0:N/2,2),gllo1(M,2),gllo2(M,0:N/2,2),&
       gacd(M,1:N/2),gacu(M,1:N/2,2),gama(N,2),uq2(M,0:N/2),i1coe(M,N,M,N,M),nx(M,N,4),ns(M,4),&
       El(0:xl),Gk(M,N,4),fpm(0:xl),fpola(0:xl),v_omega0(M,M,0:N/2),numrlo(2),numllo(2),rnu(N),rdo(N),cutoff(3)
  integer i1min_max(M,N,M,N,2),j1j3(M,N,M,N,M,4),kacdown(M,1:N/2),kacup(M,1:N/2,2)
  complex*16 v(M,M,0:N/2)
 
  integer myid,ierr,CPU,meanx
 
  real*8 integral_q,integral_omega
  
  integer shift
  
end module constants
