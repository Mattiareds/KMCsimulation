
MODULE MIN_MOD
  
  USE PARAMETERS
  
  Implicit None
  SAVE
  Integer :: nfiles ! number of files
  Integer :: natom    ! number of atoms in the cluster
  Integer :: n_atom(nelemmax) !number of atoms for each species
  Integer :: ind_mc ! total number of call to moves and running index (note that the grouping move has several steps inside)
  Integer :: isucc ! flag for the success of local minimization
  Integer :: iatom,imin,ispec,iunit
  Integer :: nelem !number of elements in the file of the potential  
  Integer :: number_of_parameters ! it is the number of coordinates, 3*natom
  Integer :: n_err !counter of failed minimizations
  Integer :: number_of_species !counts the species actually present in the cluster
  character*2 :: elem(nelemmax) ! simboli chimici delle specie
  Character*2 :: spec(natmax) ! species of the atoms of the cluster
  Character*20 :: param_file
  Character(ncharmax) :: nomefilein(nmaxfiles)  !names of the files to be minimized
  Character(ncharmax+2) :: nomefileout !name of the minimized file
  Character(ncharmax+2) :: name_of_the_file !for xyz_writer 
   
  logical :: minimization 
   
  !for the Gupta potential
  Real*8 :: p(nelemmax,nelemmax),q(nelemmax,nelemmax),a(nelemmax,nelemmax),qsi(nelemmax,nelemmax)
  Real*8 :: rat(nelemmax,nelemmax),nn(nelemmax,nelemmax)
  Real*8 :: cutoff_start(nelemmax,nelemmax),cutoff_end(nelemmax,nelemmax)
  Real*8 :: mass(nelemmax)
  Real*8 :: x5(nelemmax,nelemmax),x4(nelemmax,nelemmax),x3(nelemmax,nelemmax)
  Real*8 :: a5(nelemmax,nelemmax),a4(nelemmax,nelemmax),a3(nelemmax,nelemmax)
   
  ! running coordinates in angstrom
  Real*8 :: xin(natmax),yin(natmax),zin(natmax)
      
  !for the interaction potential
  Character*2 :: cutoff
  Integer :: itype(natmax)
  Real*8 :: fx(natmax),fy(natmax),fz(natmax)
  Real*8 :: dfx(natmax),dfy(natmax),dfz(natmax)
  Real*8 :: ener_atom(natmax)
  Real*8 :: ener
        
  !for bigvoi_md and voisin
  logical :: verlet_cage !to use the Verlet cage
  Real*8 :: bigvoi_cage
  Real*8 :: rshell,delta,firstd,secd,dmax,cutoff_end_max
  Real*8 :: cutoff_end2(nelemmax,nelemmax)
  Integer :: nvois(natmax),ivois(nvmax0,natmax)
  Integer :: nv4(natmax),iv4(nvmax0,natmax)
  
END MODULE MIN_MOD

