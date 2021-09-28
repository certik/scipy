      module X
      implicit none

      integer, parameter :: dp = kind(0.d0)

      contains

      real(dp) function alnorm(X, b)
      real(dp), intent(in) :: X
      logical, intent(in) :: b
      alnorm = 1.0
      end function

      real(8) function prho(n, is, ifault)
c
c        Algorithm AS 89   Appl. Statist. (1975) Vol.24, No. 3, P377.
c       
c        To evaluate the probability of obtaining a value greater than or
c        equal to is, where is=(n**3-n)*(1-r)/6, r=Spearman's rho and n
c        must be greater than 1
c
c     Auxiliary function required: ALNORM = algorithm AS66
c
      integer, intent(in) :: n, is
      integer, intent(out) :: ifault
      integer :: l(6)
      double precision b, x, y, z, u
      double precision, parameter :: zero = 0.0d0, one=1.0d0,
     $    two=2.0d0, six=6.0d0
      double precision, parameter :: c1=0.2274d0, c2 =0.2531d0,
     $  c3 =0.1745d0, c4 =0.0758d0, c5 =0.1033d0, c6 =0.3932d0,
     $  c7 =0.0879d0, c8 =0.0151d0, c9 =0.0072d0, c10=0.0831d0,
     $  c11=0.0131d0, c12=0.00046d0
      integer :: js, nfac, i, ifr, m, ise, n1, mt, nn
      integer :: goto
c
c        Test admissibility of arguments and initialize
c
      prho = one
      ifault = 1
      if (n .le. 1) return
      ifault = 0
      if (is .le. 0) return
      prho = zero
      if (is .gt. n * (n * n -1) / 3) return
      js = is
      if (js .ne. 2 * (js / 2)) js = js + 1
      if (n .gt. 6) goto=6
c
c        Exact evaluation of probability
c
      nfac = 1
      do 1 i = 1, n
        nfac = nfac * i
        l(i) = i
    1 continue
      prho = one / real(nfac,dp)
      if (js .eq. n * (n * n -1) / 3) return
      ifr = 0
      do 5 m = 1,nfac
        ise = 0
        do 2 i = 1, n
          ise = ise + (i - l(i)) ** 2
    2   continue
        if (js .le. ise) ifr = ifr + 1
        n1 = n
    3   mt = l(1)
        nn = n1 - 1
        do 4 i = 1, nn
          l(i) = l(i + 1)
    4   continue
        l(n1) = mt
        if (l(n1) .ne. n1 .or. n1 .eq. 2) goto=5
        n1 = n1 - 1
        if (m .ne. nfac) goto=3
    5 continue
      prho = real(ifr,dp) / real(nfac,dp)
      return
c
c        Evaluation by Edgeworth series expansion
c
    6 b = one / real(n,dp)
      x = (six * (real(js,dp) - one) * b / (one / (b * b) -one) -
     $  one) * sqrt(one / b - one)
      y = x * x
      u = x * b * (c1 + b * (c2 + c3 * b) + y * (-c4
     $  + b * (c5 + c6 * b) - y * b * (c7 + c8 * b
     $  - y * (c9 - c10 * b + y * b * (c11 - c12 * y)))))
c
c      Call to algorithm AS 66
c
      prho = u / exp(y / two) + alnorm(x, .true.)
      if (prho .lt. zero) prho = zero
      if (prho .gt. one) prho = one
      return
      end

      end module
