  SUBROUTINE XYZ_WRITER
  
  USE PARAMETERS
  USE MIN_MOD
  IMPLICIT NONE
  Integer iwk,istr,ielem
  
!  write (*,*) iwk,istr, enerwalker(iwk,istr)
  
  open (44,file=name_of_the_file,status='unknown')
  write (44,*) natom
  write (44,*) (elem(ielem),ielem=1,nelem), (n_atom(ielem), ielem=1,nelem), ' energy = ', ener, ' eV'
  do iatom=1,natom
     write (44,'(a2,2x,3f21.12,2x,i6)') spec(iatom),xin(iatom),yin(iatom),zin(iatom),iatom
   enddo
  close(44)
  return
END SUBROUTINE  XYZ_WRITER

