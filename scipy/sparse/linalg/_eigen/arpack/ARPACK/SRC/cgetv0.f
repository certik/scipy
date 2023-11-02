      subroutine cgetv0 
     &   (  n, resid, rnorm )

      include   'debug.h'
      include   'stat.h'

          call svout (logfil, 1, rnorm0, ndigit, 
     &                '_getv0: re-orthonalization ; rnorm0 is')
          call svout (logfil, 1, rnorm, ndigit, 
     &                '_getv0: re-orthonalization ; rnorm is')
         call svout (logfil, 1, rnorm, ndigit,
     &        '_getv0: B-norm of initial / restarted starting vector')
         call cvout (logfil, n, resid, ndigit,
     &        '_getv0: initial / restarted starting vector')
     
      return

      end
