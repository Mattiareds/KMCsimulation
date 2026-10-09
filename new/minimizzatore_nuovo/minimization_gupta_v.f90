
! DRIVER PER L'USO DI 'setulb'; ROUTINE DI MINIMIZZAZIONE DA MCS ANL. GOV.
! pacchetto di nome 'Lbfgsb.2.1.tar.gz' scaricabile da
! http://www.ece.nwu.edu/~nocedal/lbfgsb.html
!=========================================================================

Subroutine minimization_gupta_v

USE PARAMETERS
USE MIN_MOD

Implicit None

!Integer, Intent(In) :: number_of_parameters
! Needed by the minimization routine setulb()
Character(60) :: task, csave
Logical :: lsave(4)
Integer :: n, m, iprint, nbd(number_of_parameters), &
           iwa(3*number_of_parameters)
Integer :: isave(44)
Real(8) :: f, factr, pgtol
Real(8) :: dsave(29)
Real(8), Allocatable :: wa(:)
Real(8) :: l(number_of_parameters)
Real(8) :: u(number_of_parameters)
Real(8) :: g(number_of_parameters)
Real(8) :: xdble(number_of_parameters)
Real*8 :: partialx(natom),partialy(natom),partialz(natom)
Real*8 :: dx(natom),dy(natom),dz(natom)
! local 
Integer, Parameter :: itmax=2000 !The maximum number of calls to setulb
Real(8), Parameter :: epsi_lbfgsb=1.e-50_8 !Close to zero, used to check that 
                                           !setulb does not return the same energy
                                           !as in the previous step
Integer :: ind !counter
Real(8) :: gradiente(number_of_parameters)
Real(8) :: rtol,funzval_old,df,funzval
Integer :: iter !It counts the number of calls to Setulb()

  iprint = -1 ! We suppress the default output.
! We suppress both code-supplied stopping tests because the
! user is providing his own stopping criteria.
  factr=1.0d0 !1.0d+7
  pgtol=0.d0 !1.0d-5


! We now specify nbd which defines the bounds on the variables:
!                l   specifies the lower bounds,
!                u   specifies the upper bounds.
!     nbd(ind)=0 if x(ind) is unbounded,
!            1 if x(ind) has only a lower bound,
!            2 if x(ind) has both lower and upper bounds, and
!            3 if x(ind) has only an upper bound.

 nbd(1:number_of_parameters)=0 !unbounded minimization

! settaggio delle dimensioni del problema e delle 'limited memory corrections'
! N.B. valori raccomandati dagli autori per m: 3 <= m <= 20
  n=number_of_parameters
  m=7 ! m is the number of limited memory corrections

! allocations/initializations
  Allocate(wa(2*m*n+4*n+12*m*m+12*m))
  iter=0
  iwa=0 
  wa=0.d0
  isave=0
  dsave=0.d0
  lsave=.False.

! associazione del puntatore con le coordinate del punto in ingresso
  do ind=1,natom
     xdble(3*ind-2)=xin(ind)
     xdble(3*ind-1)=yin(ind)
     xdble(3*ind)=zin(ind)
  enddo

!serve per la routine di minimizzazione che vuole argomenti dbl-prc
funzval=ener
f=funzval

  do ind=1,natom
     gradiente(3*ind-2)=-fx(ind)
     gradiente(3*ind-1)=-fy(ind)
     gradiente(3*ind)=-fz(ind)
  enddo

  g(1:number_of_parameters)=gradiente(1:number_of_parameters)
  funzval_old=funzval


!We start the iteration by initializing task.
task = 'START'
!------- the beginning of the loop ----------

  
 

  do ind=1,natom
     partialx(ind)=0.d0
     partialy(ind)=0.d0
     partialz(ind)=0.d0
  enddo
  firstd=0.d0
  secd=0.d0
  call bigvoi_md

!write (*,*) 'sono in minimization_gupta', ind_mc
111  continue

iter=iter+1
!write (22,*) ' iter = ',iter


  call setulb(n,m,xdble,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,&
           csave,lsave,isave,dsave)
 
  do ind=1,natom
     dx(ind)=xdble(3*ind-2)-xin(ind)
     dy(ind)=xdble(3*ind-1)-yin(ind)
     dz(ind)=xdble(3*ind)-zin(ind)
  enddo

  do ind=1,natom
     xin(ind)=xdble(3*ind-2)
     yin(ind)=xdble(3*ind-1)
     zin(ind)=xdble(3*ind)     
  enddo
     
  
  firstd=0.d0
  secd=0.d0  
  do ind=1,natom
     partialx(ind)=partialx(ind)+dx(ind)
     partialy(ind)=partialy(ind)+dy(ind)
     partialz(ind)=partialz(ind)+dz(ind)
     delta=partialx(ind)**2+partialy(ind)**2+partialz(ind)**2
     if ((delta.gt.secd).and.(delta.le.firstd)) then
        secd=delta
     endif
     if (delta.gt.firstd) then
        secd=firstd
        firstd=delta
     endif
  enddo     
  dmax=sqrt(firstd)+sqrt(secd)
  if (dmax.gt.rshell) then !update is needed
     do ind=1,natom
        partialx(ind)=0.d0
        partialy(ind)=0.d0
        partialz(ind)=0.d0
     enddo
     call bigvoi_md
!       write (*,*) ' chiamo bigvoi ', ind_mc, iter
  endif 
  call voisin
  
  call force_gupta

  IF (task(1:2) .eq. 'FG') THEN
     funzval=ener
     f=funzval
     do ind=1,natom
        gradiente(3*ind-2)=-fx(ind)
        gradiente(3*ind-1)=-fy(ind)
        gradiente(3*ind)=-fz(ind)
     enddo
     g(1:number_of_parameters)=gradiente(1:number_of_parameters)
         !   se il numero massimo di iterazioni e' raggiunto il migliore individuo
         !   trovato sino ad ora viene restituito in P0 e funzval
     If (iter.Ge.itmax) Then
!       	    Deallocate(wa)
        write(22,*)'WARNING: at step ', ind_mc, ' the minimiziation exits with iter > itmax '
        write(22,*)
        Return
     EndIf 
!        go back to the minimization routine
     Goto 111
  ENDIF ! chiude if (task(1:2) .eq. 'FG') then

  IF (task(1:5) .eq. 'NEW_X') THEN
    ! the minimization routine has returned with a new iterate.
    ! At this point have the opportunity of stopping the iteration 
    ! or observing the values of certain parameters
    rtol=2._8*abs(funzval_old-funzval)/(abs(funzval_old)+abs(funzval))
    df=abs(funzval_old-funzval)
    funzval_old=funzval
    !   se la tolleranza e' soddisfatta o il numero massimo
    !   di iterazioni e' raggiunto il migliore individuo
    !   trovato sino ad ora viene restituito in P0 e funzval
    isucc=10
    If ((iter.Ge.itmax).and.(df.gt.epsi_lbfgsb)) Then
       write(22,*)'WARNING for file: ',ind_mc, '  ', nomefilein(ind_mc),' the minimization exits with iter > itmax '       
       isucc=11
       Return
    EndIf	
    If ((rtol.Lt.ftoll).and.(df.gt.epsi_lbfgsb)) Then
       call force_gupta
       funzval=ener
       f=funzval
       do ind=1,natom
          gradiente(3*ind-2)=-fx(ind)
          gradiente(3*ind-1)=-fy(ind)
          gradiente(3*ind)=-fz(ind)
       enddo
       g(1:number_of_parameters)=gradiente(1:number_of_parameters)
       Return
     EndIf ! chiude If(rtol.Lt.ftoll) Then
  ENDIF ! chiude if (task(1:5) .eq. 'NEW_X') then   
  
  Goto 111

End Subroutine minimization_gupta_v


