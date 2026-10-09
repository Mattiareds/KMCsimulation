SUBROUTINE MINIMIZATION_IN_LOOP
  
  USE PARAMETERS
  USE MIN_MOD
  
  implicit none
  
  integer :: ind
  
  minimization=.false.
  
  if (verlet_cage) then
     call bigvoi_md
     call voisin
     call force_gupta  
     write (22,*) 'initial energy = ',ener   
     call minimization_gupta_v 
     write (22,*) 'energy of the minimized configuration = ', ener, ' eV'
  else
     call bigvoi
     call force_gupta     
     write (22,*) 'initial energy = ',ener   
     call minimization_gupta
     write (22,*) 'energy of the minimized configuration = ', ener, ' eV'
  endif        
  
  
  if ((ener .lt. 0) .and. (.not. isnan(ener)))then
     minimization=.true.
  else	       
     if (isnan(ener)) then 
         n_err=n_err + 1
         write (22,*) 'warning, minimization not achieved at step',ind_mc, ' error number ', n_err
     endif
  endif
  
return
 
 
END SUBROUTINE MINIMIZATION_IN_LOOP
