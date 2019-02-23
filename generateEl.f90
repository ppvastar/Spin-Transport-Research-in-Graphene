subroutine generateEl
      use  constants
      integer i
      real*8 coea(0:xl),coeb(0:xl),coec(0:xl),appf(0:xl),midy(-1:xl),coev(-1:xl),fai(0:xl),midu,fcoe
      coea(0)=1.d0
      coea(xl)=1.d0
      coeb(0)=0.d0
      coeb(xl)=0.d0
      coec(0)=0.d0
      coec(xl)=0.d0
      do i=1,xl-1
       coea(i)=-2.d0
       coeb(i)=1.d0
       coec(i)=1.d0
      end do
      
      fcoe=deltal**2/aa*4.d0*pi*hbarv*rs
      appf(0)=0.d0
      appf(xl)=-driftE*len
      do i=1,xl-1
       appf(i)=fcoe*(fpm(i)-NN)
      end do
      coev(-1)=0.d0
      midy(-1)=0.d0
      do i=0,xl
       midu=coea(i)-coec(i)*coev(i-1)
       coev(i)=coeb(i)/midu
       midy(i)=(appf(i)-coec(i)*midy(i-1))/midu
      end do
      do i=xl,0,-1
       if(i.eq.xl) then
        fai(i)=midy(i)
       else
        fai(i)=midy(i)-coev(i)*fai(i+1)
       end if
      end do
  
      do i=0,xl
          if(i.eq.0) then
              El(i)=(fai(i)-fai(i+1))/deltal
          else if(i.eq.xl) then
              El(i)=(fai(i-1)-fai(i))/deltal
          else
              El(i)=0.5d0*(fai(i-1)-fai(i+1))/deltal
          end if
      end do
      
    end subroutine generateEl

