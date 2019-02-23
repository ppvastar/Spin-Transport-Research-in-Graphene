            
subroutine drift(rou_tmp,localE,localscdrou)
  use constants
  integer i,j,w,i1,j1,ll
  real*8 nxij_i1j1,Eedivh,karea,localE
  real*8 rou_tmp(M,N,4),localscdrou(M,N,4)
  real*8 around_dis
  Eedivh=e*localE/hbar
  do i=1,M
      karea=s0*2.d0*pi**2*kr(i)
      do j=1,N
          do ll=1,4
              do w=1,4              
              nxij_i1j1=nx(i,j,w)
              if(localE*nxij_i1j1.gt.0.d0) then
                 if(w.eq.1) then
                    i1=i
                    j1=j-1
                    if(j1.eq.0) then
                       j1=N
                    end if
                 else if(w.eq.2) then
                    i1=i
                    j1=j+1
                    if(j1.eq.N+1) then
                       j1=1
                    end if
                 else if(w.eq.3) then
                    i1=i-1
                    j1=j
!                   if(i1.eq.0) then
!                      i1=i
!                   end if
                 else if(w.eq.4) then
                    i1=i+1
                    j1=j
                    if(i1.eq.M+1) then
                       i1=i
                    end if
                 end if
              else
                 i1=i
                 j1=j
             end if
             if(i1.ne.0) then
                 around_dis=rou_tmp(i1,j1,ll)
             else
                 if(shift.gt.1) then
                     if(ll.eq.1.or.ll.eq.3) then
                         around_dis=1.d0-1.d-6
                     else
                         around_dis=0.d0
                     end if
                 else
                     around_dis=rou_tmp(1,j1,ll)
                 end if
             end if
             localscdrou(i,j,ll)=localscdrou(i,j,ll)+Eedivh*nxij_i1j1*ns(i,w)*around_dis/karea
         end do
       end do
        
     end do
 end do

end subroutine drift

      
subroutine generatenx
  use constants
  implicit none
  integer i,j,i1,j1,w
  real*8 k1x,k1y,k2x,k2y
  real*8 kr_i1
  do i=1,M
     do j=1,N
        do w=1,4
           if(w.eq.1) then
              i1=i
              j1=j-1
              if(j1.eq.0) then
                 j1=N
              end if
           else if(w.eq.2) then
              i1=i
              j1=j+1
              if(j1.eq.N+1) then
                 j1=1
              end if
           else if(w.eq.3) then
              i1=i-1
              j1=j
!              if(i1.eq.0) then
!                i1=i
!              end if
           else if(w.eq.4) then
              i1=i+1
              j1=j
              if(i1.eq.M+1) then
                 i1=i
              end if
           end if
           if(i1.eq.0) then
              if(shift.gt.1) then
                   kr_i1=(shift-1.5)*detak
              else 
                   i1=1
                   kr_i1=kr(1)
              end if
           else
                  kr_i1=kr(i1)
           end if
           if((i-i1.eq.0).and.(j-j1.eq.0)) then
              nx(i,j,w)=0.d0
           else
              k1x=kr(i)*dcos(deltathita*j)
              k1y=kr(i)*dsin(deltathita*j)
              k2x=kr_i1*dcos(deltathita*j1)
              k2y=kr_i1*dsin(deltathita*j1)
              nx(i,j,w)=(k2x-k1x)/dsqrt((k2x-k1x)**2+(k2y-k1y)**2)
           end if
        end do
        
     end do
  end do
 
end subroutine generatenx
      

subroutine generatens
  use constants
  implicit none
  integer i
  real*8 k1,k2
  do i=1,M
     ns(i,1)=detak
     ns(i,2)=detak
     ns(i,3)=(i+shift-2)*detak*deltathita
     if(i.ne.M) then
        ns(i,4)=(i+shift-1)*detak*deltathita
     else                     
        ns(i,4)=0.d0
     end if
  end do
end subroutine generatens

      
   
       
      
      
      
