Subroutine bigvoi_md

  USE PARAMETERS  !uso il modulo di definizione dei parametri
  USE MIN_MOD    !uso il modulo dove defin1:natomco variabili e parametri cluster

  Implicit None
  Integer :: i,j
  Real*8 :: dij2,xij,yij,zij

  nv4(1:natom)=0
  iv4(:,:)=0
    
  do i=1,natom-1                                              
     do j=i+1,natom                                                  
        xij=xin(j)-xin(i)                                          
        yij=yin(j)-yin(i)                                           
        zij=zin(j)-zin(i)             
        dij2=xij*xij+yij*yij+zij*zij
     
      
      if (dij2.lt.bigvoi_cage) then
         nv4(i)=nv4(i)+1 
         iv4(nv4(i),i)=j
      endif
      
!     
   enddo
enddo

EndSubroutine bigvoi_md
