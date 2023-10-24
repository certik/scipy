!       INTEGER FUNCTION I1MACH(I)
!       ! use iso_c_binding, only: c_loc, c_f_pointer
!       INTEGER I
!       INTEGER IMACH(16), OUTPUT, SC, SMALL(2)
!       ! INTEGER SC
!       ! SAVE SC
!       SAVE IMACH, SC
!       REAL RMACH
!       EQUIVALENCE (IMACH(4),OUTPUT), (RMACH,SMALL(1))
!       INTEGER I3, J, K, T3E(3)
!       ! INTEGER, target :: IMACH(16), SMALL(2)
!       ! INTEGER, pointer :: OUTPUT
!       ! REAL, pointer :: RMACH
!       ! call c_f_pointer(c_loc(imach(4)), output)
!       ! call c_f_pointer(c_loc(small(1)), rmach)
!       DATA T3E(1) / 9777664 /
!       DATA T3E(2) / 5323660 /
!       DATA T3E(3) / 46980 /
!       DATA SC/0/
!       IF (SC .NE. 987) THEN
!          SMALL(2) = 0
!          RMACH = 1E13
!          IF (SMALL(2) .NE. 0) THEN
!             IF (      (SMALL(1) .EQ. 1117925532 &
!                 .AND. SMALL(2) .EQ. -448790528) &
!             .OR.     (SMALL(2) .EQ. 1117925532 &
!                 .AND. SMALL(1) .EQ. -448790528)) THEN
!                IMACH(10) = 2
!                IMACH(14) = 53
!                IMACH(15) = -1021
!                IMACH(16) = 1024
!             ELSE IF ( SMALL(1) .EQ. -2065213935 &
!                .AND. SMALL(2) .EQ. 10752) THEN
!                IMACH(10) = 2
!                IMACH(14) = 56
!                IMACH(15) = -127
!                IMACH(16) = 127
!             ELSE IF ( SMALL(1) .EQ. 1267827943 &
!                .AND. SMALL(2) .EQ. 704643072) THEN
!                IMACH(10) = 16
!                IMACH(14) = 14
!                IMACH(15) = -64
!                IMACH(16) = 63
!             ELSE
!                WRITE(*,9010)
!                STOP 777
!                END IF
!             IMACH(11) = IMACH(14)
!             IMACH(12) = IMACH(15)
!             IMACH(13) = IMACH(16)
!          ELSE
!             RMACH = 1234567.
!             IF (SMALL(1) .EQ. 1234613304) THEN
!                IMACH(10) = 2
!                IMACH(11) = 24
!                IMACH(12) = -125
!                IMACH(13) = 128
!                IMACH(14) = 53
!                IMACH(15) = -1021
!                IMACH(16) = 1024
!                SC = 987
!             ELSE IF (SMALL(1) .EQ. -1271379306) THEN
!                IMACH(10) = 2
!                IMACH(11) = 24
!                IMACH(12) = -127
!                IMACH(13) = 127
!                IMACH(14) = 56
!                IMACH(15) = -127
!                IMACH(16) = 127
!                SC = 987
!             ELSE IF (SMALL(1) .EQ. 1175639687) THEN
!                IMACH(10) = 16
!                IMACH(11) = 6
!                IMACH(12) = -64
!                IMACH(13) = 63
!                IMACH(14) = 14
!                IMACH(15) = -64
!                IMACH(16) = 63
!                SC = 987
!             ELSE IF (SMALL(1) .EQ. 1251390520) THEN
!                IMACH(10) = 2
!                IMACH(11) = 24
!                IMACH(12) = -128
!                IMACH(13) = 127
!                IMACH(14) = 53
!                IMACH(15) = -1024
!                IMACH(16) = 1023
!             ELSE
!                DO 10 I3 = 1, 3
!                   J = SMALL(1) / 10000000
!                   K = SMALL(1) - 10000000*J
!                   IF (K .NE. T3E(I3)) GO TO 20
!                   SMALL(1) = J
!  10               CONTINUE
!                IMACH( 1) = 5
!                IMACH( 2) = 6
!                IMACH( 3) = 0
!                IMACH( 4) = 0
!                IMACH( 5) = 64
!                IMACH( 6) = 8
!                IMACH( 7) = 2
!                IMACH( 8) = 63
!                CALL I1MCR1(IMACH(9), K, 32767, 16777215, 16777215)
!                IMACH(10) = 2
!                IMACH(11) = 53
!                IMACH(12) = -1021
!                IMACH(13) = 1024
!                IMACH(14) = 53
!                IMACH(15) = -1021
!                IMACH(16) = 1024
!                GO TO 35
!  20            CALL I1MCR1(J, K, 16405, 9876536, 0)
!                IF (SMALL(1) .NE. J) THEN
!                   WRITE(*,9020)
!                   STOP 777
!                   END IF
!                IMACH(1) = 5
!                IMACH(2) = 6
!                IMACH(3) = 102
!                IMACH(4) = 6
!                IMACH(5) = 46
!                IMACH(6) = 8
!                IMACH(7) = 2
!                IMACH(8) = 45
!                CALL I1MCR1(IMACH(9), K, 0, 4194303, 16777215)
!                IMACH(10) = 2
!                IMACH(11) = 47
!                IMACH(12) = -8188
!                IMACH(13) = 8189
!                IMACH(14) = 94
!                IMACH(15) = -8141
!                IMACH(16) = 8189
!                GO TO 35
!                END IF
!             END IF
!          IMACH( 1) = 5
!          IMACH( 2) = 6
!          IMACH( 3) = 7
!          IMACH( 4) = 6
!          IMACH( 5) = 32
!          IMACH( 6) = 4
!          IMACH( 7) = 2
!          IMACH( 8) = 31
!          IMACH( 9) = 2147483647
!  35      SC = 987
!          END IF
!  9010 FORMAT(/' Adjust autodoubled I1MACH by uncommenting data'/ &
!       ' statements appropriate for your machine and setting'/ &
!       ' IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.')
!  9020 FORMAT(/' Adjust I1MACH by uncommenting data statements'/ &
!       ' appropriate for your machine.')
!       IF (I .LT. 1  .OR.  I .GT. 16) GO TO 40
!       I1MACH = IMACH(I)
!       RETURN
!  40   WRITE(*,*) 'I1MACH(I): I =',I,' is out of bounds.'
!       STOP
!       END
!       SUBROUTINE I1MCR1(A, A1, B, C, D)
!       INTEGER A, A1, B, C, D
!       A1 = 16777216*B + C
!       A = 16777216*A1 + D
!       END


!       program main
!          interface
!             integer function i1mach(i)
!                integer i
!             end function
!          end interface
!          print *, i1mach(4)
!       end program


double precision function d1mach(i)
integer i
integer small(2)
integer large(2)
integer right(2)
integer diver(2)
integer log10(2)
integer sc, cray1(38), j
common /d9mach/ cray1
save small, large, right, diver, log10,sc
double precision dmach(5)
equivalence (dmach(1),small(1))
equivalence (dmach(2),large(1))
equivalence (dmach(3),right(1))
equivalence (dmach(4),diver(1))
equivalence (dmach(5),log10(1))
end


program main
   interface
      double precision function d1mach(i)
         integer i
      end function
   end interface
   print *, d1mach(4)
end program