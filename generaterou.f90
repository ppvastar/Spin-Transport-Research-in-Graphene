subroutine generaterou
  use mpi
  use constants
  integer tt,i,j,l,w,n1,n2,jprime,revl,ca,cb,stat(mpi_status_size)
  real*8 rou_tmp(M,N,4),localscdrou(M,N,4)
  n1=int(0.25d0*N)+1
  n2=int(0.75d0*N)
  do tt=1,A
     ca=0
     cb=1
     call mpi_barrier(mpi_comm_world,ierr)
     call generatedens
     call mpi_gather(fpm(myid*meanx),meanx,mpi_real8,fpm,meanx,mpi_real8,0,mpi_comm_world,ierr)
     call mpi_gather(fpola(myid*meanx),meanx,mpi_real8,fpola,meanx,mpi_real8,0,mpi_comm_world,ierr)
!     if(mod(tt-1,1000).eq.0) then
!         call mpi_gather(rou(1,1,1,myid*meanx),M*N*4*meanx,mpi_real8,rou,M*N*4*meanx,mpi_real8,0,mpi_comm_world,ierr)
!         if(myid.eq.0) then
!             open(unit=20,file='start.dat',form='unformatted')
!             write(20) rou
!             write(20) (tt-1)*deltat
!             close(20)
!         end if
!     end if
     

     if(myid.eq.0) then
         call generateEl
         if(mod(tt-1,200).eq.0) then
             call writeplotdat(tt-1)
         end if
     end if
      
        
      
     call mpi_bcast(El,xl+1,mpi_real8,0,mpi_comm_world,ierr)
      
    !$omp parallel do schedule(dynamic) private(rou_tmp,localscdrou)
     
     do l=0,meanx-1    
         
         
         
         do i=1,M
             do j=1,N
                 do w=1,4
                     roupre(i,j,w,l)=rou(i,j,w,l)
                     rou_tmp(i,j,w)=roupre(i,j,w,l)
                     localscdrou(i,j,w)=0.d0
                 end do
             end do
         end do
         call coherent(rou_tmp,localscdrou)
         call scat(rou_tmp,localscdrou)
!         call drift(rou_tmp,El(myid*meanx+l),localscdrou)
!         call drift(rou_tmp,driftE,localscdrou)

         
       
         do i=1,M
             do j=1,N
                 do w=1,4
                     scdrou(i,j,w,l)=localscdrou(i,j,w)
                 end do
             end do
         end do
       
     end do
     !$omp end parallel do
     
   
     if(myid.ne.0) then
         call mpi_recv(rou(1,1,1,-1),M*N*4,mpi_real8,myid-1,1,mpi_comm_world,stat,ierr)
         call mpi_recv(roupre(1,1,1,-1),M*N*4,mpi_real8,myid-1,2,mpi_comm_world,stat,ierr)
         call mpi_recv(scdrou(1,1,1,-1),M*N*4,mpi_real8,myid-1,3,mpi_comm_world,stat,ierr)
     else
         ca=1
     end if
     
     
     do l=ca,meanx-1
         !$omp parallel do schedule(dynamic) private(j)
        do i=1,M
           do jprime=n2+1,n2+N/2
              j=mod(jprime,N)
              if(j.eq.0) then
                 j=N
              end if
              do w=1,4
                 rou(i,j,w,l)=-rnu(j)/rdo(j)*rou(i,j,w,l-1)&
                      +1.d0/rdo(j)*(roupre(i,j,w,l-1)+roupre(i,j,w,l))&
                      +deltat/rdo(j)*(scdrou(i,j,w,l-1)+scdrou(i,j,w,l))
              end do
           end do
        end do
        !$omp end parallel do
     end do
  
     if(myid.ne.CPU-1) then
        call mpi_send(rou(1,1,1,meanx-1),M*N*4,mpi_real8,myid+1,1,mpi_comm_world,ierr)
        call mpi_send(roupre(1,1,1,meanx-1),M*N*4,mpi_real8,myid+1,2,mpi_comm_world,ierr)
        call mpi_send(scdrou(1,1,1,meanx-1),M*N*4,mpi_real8,myid+1,3,mpi_comm_world,ierr)
    end if

    


     if(myid.ne.CPU-1) then
        call mpi_recv(rou(1,1,1,meanx),M*N*4,mpi_real8,myid+1,11,mpi_comm_world,stat,ierr)
        call mpi_recv(roupre(1,1,1,meanx),M*N*4,mpi_real8,myid+1,12,mpi_comm_world,stat,ierr)
        call mpi_recv(scdrou(1,1,1,meanx),M*N*4,mpi_real8,myid+1,13,mpi_comm_world,stat,ierr)
     else
        cb=2
     end if
     

     do l=meanx-cb,0,-1
        !$omp parallel do schedule(dynamic)
        do i=1,M
           do j=n1,n2         
              do w=1,4
                 rou(i,j,w,l)=-rnu(j)/rdo(j)*rou(i,j,w,l+1)&
                      +1.d0/rdo(j)*(roupre(i,j,w,l+1)+roupre(i,j,w,l))&
                      +deltat/rdo(j)*(scdrou(i,j,w,l+1)+scdrou(i,j,w,l))
              end do
           end do
        end do
        !$omp end parallel do
     end do
   
     if(myid.ne.0) then
        call mpi_send(rou(1,1,1,0),M*N*4,mpi_real8,myid-1,11,mpi_comm_world,ierr)
        call mpi_send(roupre(1,1,1,0),M*N*4,mpi_real8,myid-1,12,mpi_comm_world,ierr)
        call mpi_send(scdrou(1,1,1,0),M*N*4,mpi_real8,myid-1,13,mpi_comm_world,ierr)
     end if
     
 
  end do
     
end subroutine generaterou

           
