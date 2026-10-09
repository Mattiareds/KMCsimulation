subroutine voisin

USE PARAMETERS  !uso il modulo di definizione dei parametri
USE MIN_MOD     !uso il modulo dove defin1:natomco variabili e parametri cluster

Implicit None
Integer :: i,jv,j,ik,jk
Real*8 :: dij2,xij,yij,zij
! stop
 
do i=1,natom
   nvois(i)=0
enddo

do i=1,natom-1                                                   
   do jv=1,nv4(i)                                                   
      j=iv4(jv,i)
      ik=itype(i)
      jk=itype(j)                         
      xij=xin(j)-xin(i)
      yij=yin(j)-yin(i)
      zij=zin(j)-zin(i)       
      dij2=xij*xij+yij*yij+zij*zij
      !neighbor list update for force_lj.f90
      if (dij2.lt.cutoff_end2(ik,jk)) then
         nvois(i)=nvois(i)+1
         nvois(j)=nvois(j)+1
         ivois(nvois(i),i)=j
         ivois(nvois(j),j)=i
      endif      
   enddo      !fine ciclo su jv
enddo       !fine ciclo su i

!do i=1,natom
!    write (*,*) 'in voisin ',i, nvois(i)
!enddo


!stop
EndSubroutine voisin
