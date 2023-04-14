      
      subroutine calculate_cholesky_randoms(ibeg_hyd_atm3,ndummy,
     &                                      mdummy,mithrd,mirep)
                                      
      use allocatable_arrays         
      USE MKL_VSL_TYPE
      USE MKL_VSL

      implicit real(a-h,o-z)        

C note that mdummy here identifies the HI system, i.e. the 3rd dimension
C of the diffusion and cholesky decompositions

C use sgemm - documentation for which remains a marvel of obfuscation

C sgemm : C=alpha*A*B+beta*C
C       : call sgemm('T','N',M,N,K,ALPHA,A,M,B,K,BETA,C,M)
C                                          ^   ^        ^
C note that we're setting the leading dimension terms (marked above) 
C to max_hyd_atms3 since the arrays are held in memory this way

      call sgemm(
     &    'N',                        ! keep A unchanged so we go across its rows
     &    'N',                        ! keep B unchanged so we go down its cols
     &    ndummy,                     ! #rows of A
     &    num_hyd_stp,                ! #cols of B
     &    num_hyd_atms3(mdummy),      ! #cols of A 
     &    1.0,                        ! scale AB by this
     &    s_a_r(ibeg_hyd_atm3,1,mdummy,mirep), ! use w/'N'; strt at row n col 1 wrk across
     &    max_hyd_atms3,              ! do entire col?
     &    h_u_r(1,1,mdummy,mirep),          ! use w/'N'; strt at row 1 col 1 wrk down
     &    max_hyd_atms3,              ! do entire col?
     &    0.0,                        ! multiply original C by this
     &    h_c_r(ibeg_hyd_atm3,1,mdummy,mirep),     ! strt at row n col 1 wrk across
     &    max_hyd_atms3)              ! do entire col?

      return
      end

