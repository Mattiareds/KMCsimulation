SUBROUTINE READ_INITIAL
  
  USE PARAMETERS
  USE MIN_MOD
  IMPLICIT NONE
  
  Real(8) :: xp1(natmax,nelem),yp1(natmax,nelem),zp1(natmax,nelem)
  Character*2 :: sp1(natmax,nelem),metal(nelem)
  Integer naddatom,checkspecies,nmet,ispec1,ispec2,atom_number,sum_of_atoms,ielem
  Integer n_atom_first(nelem)

  write (22,*) 'START READING FILE ', ind_mc,'  ', nomefilein(ind_mc)
  ! Initial positions are determined here only for the flying walkers. Their subwalakers get the same initial positions after 
  ! the first local minimization 
  
  
  do ispec1=1,nelem
     n_atom(ispec1)=0
  enddo
  open(19,file=nomefilein(ind_mc),status='old')
     read (19,*) natom
     write (22,*)'the cluster in ', nomefilein(ind_mc), ' has ', natom, ' atoms' 
     read (19,*) 
     do iatom=1,natom
        checkspecies=0              
        read(19,*) spec(iatom),xin(iatom),yin(iatom),zin(iatom)              
!               write(*,*) spec(iatom),xin(iatom),yin(iatom),zin(iatom)              
        
           
        do ispec1=1,nelem
           if (spec(iatom).eq.elem(ispec1)) then
              itype(iatom)=ispec1
              checkspecies=checkspecies+1
              n_atom(ispec1)=n_atom(ispec1)+1
           endif
        enddo   
        if (checkspecies.eq.0) then 
            write(22,*) 'Error: metals in ', nomefilein(ind_mc) ,' must be the same', &
            & ' as in the parametrization file of the potential.'
            write(22,*) 'atom',iatom,' is ',spec(iatom),' which is not in the parameter file of the potential'        
            stop
        endif                           
     enddo                            
     write (22,*) 'the numbers of atoms of the different species are'
     do ispec1 =1,nelem
        write (22,*) elem(ispec1),n_atom(ispec1)
     enddo         
  close(19)
  
  !counts the total number of coordinates
  
   
  !counts the number of elements actually present in the clusters
  number_of_species=0
  do ispec=1,nelem
     if ((n_atom(ispec).gt.0)) number_of_species=number_of_species+1
  enddo
  write (22,*) 'the number of atomic species actually present in the cluster is ', number_of_species

  !calculates the squares of the cutoff_end
  do ispec1=1,nelem
     do ispec2=1,nelem
        cutoff_end2(ispec1,ispec2)=cutoff_end(ispec1,ispec2)*cutoff_end(ispec1,ispec2)
 !       write (*,*) ispec1,ispec2,cutoff_end2(ispec1,ispec2)
     enddo
  enddo
!stop
  
  ! calculates the maximum cutoff_end
  cutoff_end_max=0.d0
  do ispec1=1,nelem
     do ispec2=1,nelem
     if ((n_atom(ispec1).gt.0).and.(n_atom(ispec2).gt.0)) then
        if(cutoff_end(ispec1,ispec2).gt.cutoff_end_max) cutoff_end_max=cutoff_end(ispec1,ispec2)
     endif
     enddo
  enddo

  write(22,*) 'the maximum cutoff_end is ', cutoff_end_max, ' Angstrom'        
  bigvoi_cage=(cutoff_end_max+rshell)**2

  number_of_parameters=3*natom

!  stop

  return
  
  END SUBROUTINE READ_INITIAL
