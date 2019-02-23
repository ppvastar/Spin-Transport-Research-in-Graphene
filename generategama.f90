subroutine generategama
  use constants
  implicit none
  integer j
    do j=1,N
     gama(j,1)=-lambdaBR*dsin(deltathita*j)+0.5d0*gfactor*uB*B*dcos(Bphi*pi)/hbar
     gama(j,2)=lambdaBR*dcos(deltathita*j)+0.5d0*gfactor*uB*B*dsin(Bphi*pi)/hbar
  end do
  end subroutine generategama
  
  
  
  

     
  
