PROGRAM MINIMIZER
  
  USE PARAMETERS
  USE MIN_MOD
  
  implicit none
  
  Integer indfile,ieof,reason


  open (22,file='output_min.out',status='unknown')
  write(22,*)
  
! READING THE GENERAL INPUT AND NAMING FILES AND VARIABLES

  call read_min !to read the main settings of the simulation

  call read_gupta !to read the parameters of the potential and calculate the coefficients of the cutoff polynomials
 
  write (22,*)  

    !read the list of files to be analyzed
  indfile=0
  ieof=0
  open(12,file="listanomi.in",status="old")
     do while(ieof.eq.0)
        indfile=indfile+1
        read(12,*,iostat=reason) nomefilein(indfile)
        if (reason > 0)  then
           write (22,*)
           write (22,*) "bad record at line ", indfile
           write (22,*) 'program stops'
           stop
        else if (reason < 0) then
           ieof=1           
        endif
     enddo
  close(12)
  nfiles=indfile-1
  if (nfiles.gt.nmaxfiles) then
     write (22,*) "filess are too many, nfiles must be less than ", nmaxfiles
  endif
  write (22,*) 'number of files to be minimized = ', nfiles,' They are listed below'
  do indfile=1,nfiles
     write (22,*),indfile, '  ', nomefilein(indfile)
  enddo
  



!  stop

! all output files have been opened and all the data of the initial configuration written        
! now the initial settings are ready, the Monte Carlo cycle can start
  
  write (22,*) 'all initial setup is ready, the minimization loop can start'
  write (22,*) 
  open(77,file='minimized_energies.out',status='unknown')
     write (77,*) '# file number, file name, minimzed energy (eV)'
  close(77)   
     
  open(77,file='minimized_energies.out',status='old',access='append')
  do ind_mc=1,nfiles
     open(33,file=nomefilein(ind_mc),status='old')
        call read_initial !to read the coordinates of nomefilein(ind_mc)
!        stop
        call minimization_in_loop
        name_of_the_file(1:2)="q_"
        name_of_the_file(3:ncharmax+2)= nomefilein(ind_mc)(1:ncharmax)
        write(22,*) 'the minimized file will be ', name_of_the_file
        call xyz_writer !to write the coordinates of the minimized file
        write (22,*)
     close(33)
     write (77,*) ind_mc,name_of_the_file,ener
  enddo      
  close(77)
  close(22)          
  END PROGRAM MINIMIZER
  

