         
program main
  use mpi 
  use constants
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,CPU,ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call init 
  call generategama
  call generatenx
  call generatens
  if(Ni+Ci.ne.0.d0) then
      call generateuq2
  end if
  if(e_e_inter+e_e_intra+h_f.ne.0) then
     call generatev
     call generatekm_j1j3_i1coe
  end if
  if(e_lo_remote1+e_lo_remote2.ne.0) then
     call generaterlog
  end if
  if(e_lo_local1+e_lo_local2.ne.0) then
     call generatellog
  end if
  if(e_ac.eq.1) then
     call generateacg
 end if
  call rou_init
  call generaterou
  call mpi_finalize(ierr)
end program main

     
