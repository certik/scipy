double precision function d1mach(i)
use iso_c_binding, only: c_loc, c_f_pointer
integer i
! integer small(2)
! integer large(2)
! integer right(2)
! integer diver(2)
! integer log10(2)
integer sc, cray1(38), j
common /d9mach/ cray1
! save small, large, right, diver, log10,
save sc
! double precision dmach(5)
! equivalence (dmach(1),small(1))
! equivalence (dmach(2),large(1))
! equivalence (dmach(3),right(1))
! equivalence (dmach(4),diver(1))
! equivalence (dmach(5),log10(1))

! implicit none
DOUBLE PRECISION, target :: DMACH(5)
INTEGER*4, pointer :: SMALL(:)
INTEGER*4, pointer :: LARGE(:)
INTEGER*4, pointer :: RIGHT(:)
INTEGER*4, pointer :: DIVER(:)
INTEGER*4, pointer :: LOG10(:)
call c_f_pointer(c_loc(dmach(1)), small, [2])
call c_f_pointer(c_loc(dmach(2)), large, [2])
call c_f_pointer(c_loc(dmach(3)), right, [2])
call c_f_pointer(c_loc(dmach(4)), diver, [2])
call c_f_pointer(c_loc(dmach(5)), log10, [2])
DATA SC/0/
if (sc .ne. 987) then
   dmach(1) = 1.d13
   if (small(1) .eq. 1117925532 .and. small(2) .eq. -448790528) then
      small(1) = 1048576
      small(2) = 0
      large(1) = 2146435071
      large(2) = -1
      right(1) = 1017118720
      right(2) = 0
      diver(1) = 1018167296
      diver(2) = 0
      log10(1) = 1070810131
      log10(2) = 1352628735
   else if ( small(2) .eq. 1117925532 .and. small(1) .eq. -448790528) then
      small(2) = 1048576
      small(1) = 0
      large(2) = 2146435071
      large(1) = -1
      right(2) = 1017118720
      right(1) = 0
      diver(2) = 1018167296
      diver(1) = 0
      log10(2) = 1070810131
      log10(1) = 1352628735
   else if ( small(1) .eq. -2065213935 .and. small(2) .eq. 10752) then
      small(1) = 128
      small(2) = 0
      large(1) = -32769
      large(2) = -1
      right(1) = 9344
      right(2) = 0
      diver(1) = 9472
      diver(2) = 0
      log10(1) = 546979738
      log10(2) = -805796613
   else if ( small(1) .eq. 1267827943 .and. small(2) .eq. 704643072) then
      small(1) = 1048576
      small(2) = 0
      large(1) = 2147483647
      large(2) = -1
      right(1) = 856686592
      right(2) = 0
      diver(1) = 873463808
      diver(2) = 0
      log10(1) = 1091781651
      log10(2) = 1352628735
   else if ( small(1) .eq. 1120022684 .and. small(2) .eq. -448790528) then
      small(1) = 1048576
      small(2) = 0
      large(1) = 2147483647
      large(2) = -1
      right(1) = 1019215872
      right(2) = 0
      diver(1) = 1020264448
      diver(2) = 0
      log10(1) = 1072907283
      log10(2) = 1352628735
   else if ( small(1) .eq. 815547074 .and. small(2) .eq. 58688) then
      small(1) = 16
      small(2) = 0
      large(1) = -32769
      large(2) = -1
      right(1) = 15552
      right(2) = 0
      diver(1) = 15568
      diver(2) = 0
      log10(1) = 1142112243
      log10(2) = 2046775455
   else
      dmach(2) = 1.d27 + 1
      dmach(3) = 1.d27
      large(2) = large(2) - right(2)
      if (large(2) .eq. 64 .and. small(2) .eq. 0) then
         ! Ideally lfortran should come here
         cray1(1) = 67291416
         do 10 j = 1, 20
            cray1(j+1) = cray1(j) + cray1(j)
10               continue
         cray1(22) = cray1(21) + 321322
         do 20 j = 22, 37
            cray1(j+1) = cray1(j) + cray1(j)
20               continue
         if (cray1(38) .eq. small(1)) then
            call i1mcry(small(1), j, 8285, 8388608, 0)
            small(2) = 0
            call i1mcry(large(1), j, 24574, 16777215, 16777215)
            call i1mcry(large(2), j, 0, 16777215, 16777214)
            call i1mcry(right(1), j, 16291, 8388608, 0)
            right(2) = 0
            call i1mcry(diver(1), j, 16292, 8388608, 0)
            diver(2) = 0
            call i1mcry(log10(1), j, 16383, 10100890, 8715215)
            call i1mcry(log10(2), j, 0, 16226447, 9001388)
         else
            write(*,9000)
            stop 779
            end if
      else
         write(*,9000)
         ! Problem is here, lfortran comes here and stops
         stop 779
         end if
      end if
   sc = 987
   end if
if (dmach(4) .ge. 1.0d0) stop 778
if (i .lt. 1 .or. i .gt. 5) then
   write(*,*) 'd1mach(i): i =',i,' is out of bounds.'
   stop
   end if
d1mach = dmach(i)
return
9000 FORMAT(/' Adjust D1MACH by uncommenting data statements appropriate for your machine.')
end
subroutine i1mcry(a, a1, b, c, d)
integer a, a1, b, c, d
a1 = 16777216*b + c
a = 16777216*a1 + d
end

program main
   interface
      double precision function d1mach(i)
         integer i
      end function
   end interface
   print *, d1mach(4)
end program