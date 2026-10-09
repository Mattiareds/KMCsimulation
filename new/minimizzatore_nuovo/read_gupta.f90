SUBROUTINE READ_GUPTA
  
  USE PARAMETERS
  USE MIN_MOD
  IMPLICIT NONE

  Real(8) :: ar,br,cr,ab,bb,cb,dik0
  Integer i,j
  Character(2) :: pippo(nelemmax)
  Character(2) :: sym
 
!************************************************************************
  write(22,*)
  write (22,*) 'START READING FILE WITH POTENTIAL PARAMETERS'
  write (22,*)
  write (22,*) 'The parameters are contained in ', param_file
  open(19,file=param_file,status='old')
     read(19,*) nelem
     read(19,*)
     read(19,*) (elem(i),i=1,nelem)
     do i=1,nelem
        pippo(i)=elem(i) 
     enddo
     write(22,*) "System " 
     write(22,*) (pippo(i),i=1,nelem) 
     write(22,*)
     read (19,*)
     read (19,*) (mass(i),i=1,nelem)
     write (22,*) "Masses " 
     write (22,*) '             ', (pippo(j),j=1,nelem)
     write(22,*) ' ', (mass(i),i=1,nelem)
     write(22,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(p(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           p(j,i)=p(i,j)
        enddo
     enddo
     write (22,*) 'p(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*)  elem(i), (p(i,j),j=1,nelem)           
     enddo
     write(22,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(q(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           q(j,i)=q(i,j)
        enddo
     enddo
     write(22,*) 'q(i,j)'
     write(22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (q(i,j),j=1,nelem)           
     enddo       
     write(22,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(a(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           a(j,i)=a(i,j)
        enddo
     enddo
     write (22,*) 'a(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (a(i,j),j=1,nelem)           
     enddo  
     write(22,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(qsi(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           qsi(j,i)=qsi(i,j)
        enddo
     enddo
     write (22,*) 'qsi(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (qsi(i,j),j=1,nelem)           
     enddo  
     write(22,*)
     read(19,*)
     read(19,*)
     read(19,*)
     do i=1,nelem
        read(19,*) sym,(rat(i,j),j=1,i)           
     enddo   
     do i=2,nelem
        do j=1,i-1
           rat(j,i)=rat(i,j)
        enddo
     enddo
     write (22,*) 'rat(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (rat(i,j),j=1,nelem)           
     enddo  
     write(22,*)

     if (cutoff.eq.'si') then
        read(19,*)
        read(19,*)
        read(19,*)
        do i=1,nelem
           read(19,*) sym,(cutoff_start(i,j),j=1,i)           
        enddo   
        do i=2,nelem
           do j=1,i-1
              cutoff_start(j,i)=cutoff_start(i,j)
           enddo
        enddo
        read(19,*)
        read(19,*)
        read(19,*)
        do i=1,nelem
           read(19,*) sym,(cutoff_end(i,j),j=1,i)           
        enddo 
        do i=2,nelem
           do j=1,i-1
              cutoff_end(j,i)=cutoff_end(i,j)
           enddo
        enddo
     else 
        do i=1,nelem
           do j=1,nelem
              cutoff_start(i,j)=2000.d0
              cutoff_end(i,j)=2002.d0
           enddo
        enddo
     endif   
  
     write (22,*) 'cutoff_start(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (cutoff_start(i,j),j=1,nelem)           
     enddo  

     write (22,*)
     write (22,*) 'cutoff_end(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (cutoff_end(i,j),j=1,nelem)           
     enddo  
     
     write(22,*)        
     write(22,*)     
     write(22,*) 'END READING FILE ', param_file
!     stop
  close(19)
! ************************************************************************

   
! nearest neighbour distances and parameters of the cutoff polynomials
do i=1,nelem
   do j=1,nelem
      nn(i,j)=2.d0*rat(i,j)     
      dik0=nn(i,j)
        ar=-a(i,j)*dexp(-p(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))**3 
        br=-(p(i,j)/dik0)*a(i,j)*dexp(-p(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))**2
        cr=-((p(i,j)/dik0)**2) &
             *a(i,j)*dexp(-p(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))
        ab=-qsi(i,j)*dexp(-q(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))**3
        bb=-(q(i,j)/dik0)*qsi(i,j)*dexp(-q(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))**2 
        cb=-((q(i,j)/dik0)**2) &
             *qsi(i,j)*dexp(-q(i,j)*(cutoff_start(i,j)/dik0-1.d0))/(cutoff_end(i,j)-cutoff_start(i,j))
        x5(i,j)=(12.d0*ab-6.d0*bb+cb)/(2.d0*(cutoff_end(i,j)-cutoff_start(i,j))**2)
        x4(i,j)=(15.d0*ab-7.d0*bb+cb)/(cutoff_end(i,j)-cutoff_start(i,j))
        x3(i,j)=(20.d0*ab-8.d0*bb+cb)/2.d0
        a5(i,j)=(12.d0*ar-6.d0*br+cr)/(2.d0*(cutoff_end(i,j)-cutoff_start(i,j))**2)
        a4(i,j)=(15.d0*ar-7.d0*br+cr)/(cutoff_end(i,j)-cutoff_start(i,j))
        a3(i,j)=(20.d0*ar-8.d0*br+cr)/2.d0

   enddo
enddo 

  write (22,*)
  write (22,*) 'Parameters of the cutoff polynomials  a5,a4,a3,x5,x4,x3'
  write (22,*)
  
     write (22,*) 'a5(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (a5(i,j),j=1,nelem)           
     enddo  
     write (22,*)     
     write (22,*) 'a4(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (a4(i,j),j=1,nelem)           
     enddo  
     write (22,*)     
     write (22,*) 'a3(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (a3(i,j),j=1,nelem)           
     enddo  
     write (22,*)     
     write (22,*) 'x5(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (x5(i,j),j=1,nelem)           
     enddo  
     write (22,*)     
     write (22,*) 'x4(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (x4(i,j),j=1,nelem)           
     enddo  
     write (22,*)     
     write (22,*) 'x3(i,j)'
     write (22,*) '             ', (pippo(j),j=1,nelem)
     do i=1,nelem
        write(22,*) elem(i), (x3(i,j),j=1,nelem)           
     enddo  
     write (22,*)     

 
 
  
  
  return
  
  END SUBROUTINE READ_GUPTA
