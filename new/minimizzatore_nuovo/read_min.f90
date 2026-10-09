  SUBROUTINE READ_MIN
  
  USE PARAMETERS
  USE MIN_MOD
  IMPLICIT NONE

  !Local variables
  Integer :: i,j
      
   
  open (22,file='output_min.out',status='unknown')
  write(22,*)
!
!************************************************************************
  write (22,*) 'START READING FILE input_bh.in'
  write(22,*)
  open(19,file='input_min.in',status='old')
  read(19,*)!******Interaction potentials
  read (19,*) cutoff !cutoff or no cutoff for the atom-atom interaction
  read (19,*) verlet_cage,rshell !to use the Verlet cage or not, and to determine the thickness of its shell 
  read (19,*) param_file
  write (22,*) 'is the cutoff applied? ',cutoff
  if(verlet_cage) then
    write(22,*) 'the Verlet cage is used with shell thickness ',rshell, ' Angstrom'
  endif
  if (verlet_cage.and.(cutoff.ne.'si')) then
     write (22,*) 'Inconsistent use of the Verlet cage because there is no cutoff'
     stop
  endif
   write(22,*)
  close(19)


!  stop

  return
END SUBROUTINE READ_MIN








