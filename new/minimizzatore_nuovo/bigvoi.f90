subroutine bigvoi

USE PARAMETERS
USE MIN_MOD
implicit none
  
  integer*4 :: i,j,ik,jk
  Real(8) :: dij2,xij,yij,zij
   
  do i=1,natom
     nvois(i)=0
  enddo
      


  do i=1,natom-1
     do j=i+1,natom
     ik=itype(i)
     jk=itype(j)
     xij=xin(j)-xin(i)
     yij=yin(j)-yin(i)
     zij=zin(j)-zin(i)
     dij2=xij*xij+yij*yij+zij*zij
     if (dij2.lt.cutoff_end2(ik,jk)) then
        nvois(i) = nvois(i)+1
        ivois(nvois(i),i) = j
        nvois(j) = nvois(j)+1
        ivois(nvois(j),j) = i
        endif
     enddo
!     if (nvois(i).ge.natom) then
!        write (*,*) ' too many neighbours for atom ', i,' they are ', nvois(i)
!        stop
!     endif
  enddo
  
  return
  end subroutine bigvoi
