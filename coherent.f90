subroutine coherent(rou_tmp,localscdrou)
  use constants
  implicit none
  real*8 hartree1,hartree2,hartree4,gamax,gamay
  real*8 rou_tmp(M,N,4),localscdrou(M,N,4)
  integer i,j
  do j=1,N
     gamax=gama(j,1)
     gamay=gama(j,2)
     do  i=1,M
        
       hartree1=0.d0
       hartree2=0.d0
       hartree4=0.d0
        if(h_f.eq.1) then
           call hartree_fork(rou_tmp,i,j,hartree1,hartree2,hartree4)
        end if
        localscdrou(i,j,1)=localscdrou(i,j,1)-2.d0*(gamax*rou_tmp(i,j,4)+gamay*rou_tmp(i,j,2))+hartree1
        localscdrou(i,j,2)=localscdrou(i,j,2)+gamay*(rou_tmp(i,j,1)-rou_tmp(i,j,3))+hartree2
        localscdrou(i,j,3)=localscdrou(i,j,3)+2.d0*(gamax*rou_tmp(i,j,4)+gamay*rou_tmp(i,j,2))-hartree1
        localscdrou(i,j,4)=localscdrou(i,j,4)+gamax*(rou_tmp(i,j,1)-rou_tmp(i,j,3))+hartree4
     end do
  end do
end subroutine coherent
      
subroutine hartree_fork(rou_tmp,i,j,hartree1,hartree2,hartree4)
  use constants
  implicit none
  integer i,j,ii,jj,jjj
  real*8 hartree1,hartree2,hartree4,vtmp,rou_tmp(M,N,4)
  
  do ii=1,M
      do jj=1,N
          if (abs(j-jj).le.N/2) then
              jjj=abs(j-jj)
          else
              jjj=N-abs(j-jj)
          end if
          vtmp=v_omega0(i,ii,jjj)*0.5d0*(1.d0+dcos(deltathita*jjj))*kr(ii)
          hartree1=hartree1+vtmp*(rou_tmp(ii,jj,2)*rou_tmp(i,j,4)-rou_tmp(ii,jj,4)*rou_tmp(i,j,2))
          hartree2=hartree2-vtmp*((rou_tmp(ii,jj,1)-rou_tmp(ii,jj,3))*rou_tmp(i,j,4)&
              -(rou_tmp(i,j,1)-rou_tmp(i,j,3))*rou_tmp(ii,jj,4))
          hartree4=hartree4+vtmp*((rou_tmp(ii,jj,1)-rou_tmp(ii,jj,3))*rou_tmp(i,j,2)&
              -(rou_tmp(i,j,1)-rou_tmp(i,j,3))*rou_tmp(ii,jj,2))

     end do
  end do
  hartree1=hartree1*hartreecoe
  hartree2=hartree2*0.5d0*hartreecoe
  hartree4=hartree4*0.5d0*hartreecoe

end subroutine hartree_fork
