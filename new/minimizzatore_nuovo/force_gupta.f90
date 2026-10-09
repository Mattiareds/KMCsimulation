subroutine force_gupta
use iso_c_binding!in realta` e` uno statement per fortran 2003 ma se compila non mi lamento
USE PARAMETERS  
USE MIN_MOD    


Implicit None
Integer :: i,j,k
Integer :: itypi,itypk
Real*8 :: den(natom),frx(natom),fry(natom),frz(natom)
Real*8 :: dikm,dikm2,dikm3,dikm4,dikm5
Real*8 :: dik0,dik,xik,yik,zik
Real*8 :: ebi,eneri,eri,for,forsudik,denik
Real*8 :: fb(natom,natom)
Real*8 :: espo,pexpp,qexpp,qsiexpq
Real*8 :: psudik0,dueapsudik0,dueq,qsi2,qsi2qsudik0

ener=0.d0                                                        
do i=1,natom
   dfx(i)=fx(i)
   dfy(i)=fy(i)
   dfz(i)=fz(i)
   den(i)=0.d0                                                       
   frx(i)=0.d0                                                       
   fry(i)=0.d0                                                       
   frz(i)=0.d0                                                       
enddo

do i=1,natom

   ebi=0.d0                                                       
   eri=0.d0                                                          
   eneri=0.d0
   !se l-atomo i non ha vicini, pongo la forza = 0
   do j=1,nvois(i) 
      k=ivois(j,i)  
      if (k.gt.i) then
         itypi=itype(i)
         itypk=itype(k)
     
      dik0=nn(itypi,itypk)
      psudik0=p(itypi,itypk)/dik0
      dueapsudik0=2.d0*a(itypi,itypk)*psudik0
      dueq=2.d0*q(itypi,itypk)
      qsi2=qsi(itypi,itypk)*qsi(itypi,itypk)
      qsi2qsudik0=qsi2*q(itypi,itypk)/dik0     
      xik=xin(k)-xin(i)
      yik=yin(k)-yin(i)
      zik=zin(k)-zin(i) 
      dik=dsqrt(xik*xik+yik*yik+zik*zik)
      if (dik.lt.cutoff_start(itypi,itypk)) then                         
         espo=1.d0-dik/dik0                                
         pexpp=dexp(p(itypi,itypk)*espo)
         qexpp=dexp(dueq*espo)
         qsiexpq=qsi2*qexpp
         for=pexpp*dueapsudik0
         eri=eri+2.d0*a(itypi,itypk)*pexpp
         fb(k,i) = (q(itypi,itypk)/dik0)*qsiexpq
      else
         dikm=dik-cutoff_end(itypi,itypk)
         dikm2=dikm*dikm
         dikm3=dikm2*dikm
         dikm4=dikm3*dikm
         dikm5=dikm4*dikm  
         qsiexpq=x5(itypi,itypk)*dikm5+x4(itypi,itypk)*dikm4+x3(itypi,itypk)*dikm3
         for=-2.d0*(5.d0*a5(itypi,itypk)*dikm4+4.d0*a4(itypi,itypk)*dikm3+3.d0*a3(itypi,itypk)*dikm2)
         eri=eri+2.d0*(a5(itypi,itypk)*dikm5+a4(itypi,itypk)*dikm4+a3(itypi,itypk)*dikm3)       
         fb(k,i) = -((5.d0*x5(itypi,itypk)*dikm4+4.d0*x4(itypi,itypk)*dikm3+&
              & 3.d0*x3(itypi,itypk)*dikm2))*qsiexpq
         qsiexpq = qsiexpq**2
      endif
!      write (44,'(2i6,2e18.8)') i,k,dik0,for

      den(i)=qsiexpq+den(i)
      den(k)=qsiexpq+den(k)

      forsudik=for/dik
      frx(i)=frx(i)-forsudik*xik                                   
      fry(i)=fry(i)-forsudik*yik                                    
      frz(i)=frz(i)-forsudik*zik  

      frx(k)=frx(k)+forsudik*xik                                   
      fry(k)=fry(k)+forsudik*yik                                    
      frz(k)=frz(k)+forsudik*zik  

   endif !su k>i   
   enddo !su j                                                    
   
   ebi=dsqrt(den(i))                                                
   den(i)=1.d0/ebi                                                  
   eneri=eri-ebi
   ener=ener+eneri
!   ener_atom(i)=eneri                                                   

enddo !su i  
  do 30 i=1,natom
     IF(nvois(i) .gt. 0) THEN !exist neighbours
        do 40 j=1,nvois(i)
           k=ivois(j,i)
           if(k.gt.i)then
              xik=xin(k)-xin(i)
              yik=yin(k)-yin(i)
              zik=zin(k)-zin(i) 
              dik=dsqrt(xik*xik+yik*yik+zik*zik)
              denik=fb(k,i)*(den(i)+den(k))/dik
              frx(i)=frx(i)+denik*xik
              fry(i)=fry(i)+denik*yik
              frz(i)=frz(i)+denik*zik
              frx(k)=frx(k)-denik*xik
              fry(k)=fry(k)-denik*yik
              frz(k)=frz(k)-denik*zik
           endif!k>i
40      enddo
     ENDIF ! no neighbours
     fx(i)=frx(i)
     fy(i)=fry(i)
     fz(i)=frz(i)
30 enddo

!   is_nan = .false.
EndSubroutine force_gupta
