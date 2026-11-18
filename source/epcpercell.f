        SUBROUTINE PERCELL(NUMERO,II,JJ,VALUE,CMNFLD,FLAG)
!NG This is a stand-alone subroutine to find and
!NG average observations for same cells.
!NG It takes in a number (NUMERO) of values (VALUE) 
!NG at [II,JJ] locations (possibly not unique), including
!NG a common field (CMNFLD, defined for all [II,JJ] but
!NG should be the same for common ones).
!NG It then averages values at common cells, ignoring flagged
!NG values (=FLAG), and carrying over CMNFLD.
!NG It finally returns the redefined set as:
!NG NUMERO=number of unique cells
!NG [II,JJ]=I,J locations of the unique cells
!NG VALUE=averaged values at the unique cells.
!NG Note that the algorithm is designed to retain all
!NG unique cells, includines ones with values that are all
!NG flagged (and are ignored for averaging).
!NG Nickitas Georgas, 05/06/2008

        integer NUMERO
        real FLAG
	real VALUE(NUMERO),VALOUT(NUMERO)
        real CMNFLDO(NUMERO),CMNFLD(NUMERO)
        integer ICNT(NUMERO),INDX(NUMERO)
        integer II(NUMERO),JJ(NUMERO)
        integer IOUT(NUMERO),JOUT(NUMERO)

        IUNQ=1
        do n=1,NUMERO
         VALOUT(n)=0.0
         INDX(n)=N
         ICNT(n)=1
         IOUT(n)=0
         JOUT(n)=0
        enddo
        IUNQ=NUMERO
        do n=1,NUMERO-1
         do m=n+1,NUMERO
          if (II(m).eq.II(n).and.JJ(m).eq.JJ(n)) then
           ICNT(M)=0
           if (INDX(n).ne.N) then
            INDX(m)=INDX(n)
           else
            ICNT(N)=ICNT(N)+1
            INDX(M)=N
            IUNQ=IUNQ-1
           endif
          endif
         enddo
        enddo
        do n=1,NUMERO
         if (VALUE(n).eq.FLAG) then
          ICNT(INDX(n))=ICNT(INDX(n))-1 ! Do not average flagged data
         endif
         iout(INDX(n))=ii(n)
         jout(INDX(n))=jj(n)
         CMNFLDO(INDX(n))=CMNFLD(n)
        enddo
        do n=1,NUMERO
         if (VALUE(n).ne.FLAG)
     +    VALOUT(INDX(n))=VALOUT(INDX(n))+VALUE(n)/float(ICNT(INDX(n)))
        enddo
!NG        write (*,'(2(f10.2,2i5))')
!NG     +  (VALUE(n),ICNT(n),INDX(n),VALOUT(n),iout(n),jout(n),
!NG     +n=1,NUMERO)
        m=0
        do n=1,NUMERO
         if (ICNT(n).eq.0) then 
          VALOUT(n)=FLAG
         endif
         m=m+1
         if (IOUT(n).eq.0.and.JOUT(n).eq.0) then
          m=m-1
         else
          VALOUT(m)=VALOUT(n)
          iout(m)=iout(n)
          jout(m)=jout(n)
          CMNFLDO(m)=CMNFLD(n)
         endif
        enddo
!NG        write (*,*) '---'
!NG        write (*,'(2f10.2,2i5)')
!NG     +   (VALOUT(n),CMNFLDO(n),iout(n),jout(n),n=1,m)

        NUMERO=IUNQ
        do n=1,NUMERO
         VALUE(n)=VALOUT(n)
         ii(n)=iout(n)
         jj(n)=jout(n)
         CMNFLD(n)=CMNFLDO(n)
        enddo
        
        return
        end
