!===Simple test driver of BestLex. Permutations of X.
!        Integer x(3),q(3)
!        Logical first
!        Data x/1,2,3/
!        first=.true.
!--Next loop is one more than Factorial(3) in this case to show that
!  Algorithm 323 terminates by setting first back to .true.
!        do i=1,7
!          write(*,*)'x=',x
!          call BestLex(x,3,q,first)
!          if (first) goto 1
!        end do
!1       stop
!        end
        Subroutine BestLex(x,n,q,first)
!===ACM Algorithm 323, Generation of Permutations in Lexicographic
!   Order (G6) by R. J. Ord-Smith, CACM 11 (Feb. 1968):117
!   Original Algorithm modified via Certification by I.M. Leitch,
!   17 March 1969.
! Algol to Fortran 77 by H.D.Knoble <hdk@psu.edu>, May 1995.
        Integer n,k,m,t
        Integer x(n),q(n)
        Logical first
        !write(*,*)'mi hai chiamato?'
        !write(*,*)'xB=',x
        if (first) then
           !write(*,*)'sono nell-if'
          first=.false.
          do m=1,n-1
            q(m)=n
            !write(*,*)'poiche` n vale',n,'q vale',q
          end do
        endif
        if(q(n-1).eq.n) then
           !write(*,*)'ok, q(n-1)=n'
          q(n-1)=n-1
          t=x(n)
          !write(*,*)'tvale',t
          x(n)=x(n-1)
          x(n-1)=t
          return
        endif
        do k=n-1,1,-1
          if(q(k).eq.k) then
            q(k)=n
          else
            go to 1
          endif
        end do
        first=.true.
        k=1
        goto 2
1       m=q(k)
        t=x(m)
        !write(*,*)'tvale3',t
        x(m)=x(k)
        x(k)=t
        q(k)=m-1
        k=k+1
2       m=n
3       t=x(m)
        !write(*,*)'tvale2',t
        x(m)=x(k)
        x(k)=t
        m=m-1
        k=k+1
        if(k.lt.m) goto 3
        return
        

end subroutine BestLex
