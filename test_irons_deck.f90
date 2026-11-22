! Test 1 for the Irons (1970) frontal solver: the paper's own specimen data
! (Appendix III), run through ZIPP with the card-input ELMAK exactly as
! printed in Appendix II.
!
!   Problem 1  NEWRHS=1, NTIREX=0 (output by elements), NRMAX=2, MAXTAP=10,
!              MAXELT=20. Element 2 is fragmented into two segments.
!              The paper states the assembled system:
!                 2x2                 = 1
!                       9x3 + 2x4     = 3
!                       2x3 + 6x4 + x5 = 5
!                              x4 + 3x5 = 7
!   Re-solution NEWRHS=2, element right-hand sides presented in reverse order
!              (continuation segment first). The paper states the RHS:
!                 (4,1), (7,7), (11,14), (4,8)
!   Problem 2  NEWRHS=1, NTIREX=1 (entire output), one element with a repeated
!              nickname (5,89,5); the paper states the assembled system:
!                 6x5 + 4x89 = 4
!                 4x5 + 4x89 = 2
!
! Each answer is compared with the exact solution (computed here in quad
! precision from the systems quoted in the paper), and the residual
! b - A x of the frontal answers is computed in quad precision.

program test_irons_deck
  use zcap
  implicit none
  real(16) :: A1(4,4), B1(4,3), X1(4,3), A2(2,2), B2(2,1), X2(2,1)
  real(8)  :: xf(4,3), xf2(2)
  integer  :: n1(4) = (/ 2, 3, 4, 5 /), n2(2) = (/ 5, 89 /)
  integer  :: i, j, nfail
  real(8)  :: err, res, worst_err, worst_res

  nfail = 0
  call reset_capture()
  do i = 1, 4
    j = reg(n1(i))
  end do
  j = reg(89)

  ! ---- the card deck of Appendix III ---------------------------------------
  open (40, file='irons_appendix3.dat', status='replace')
  write (40,'(a)') '1 0 2 10 20'              ! NEWRHS NTIREX NRMAX MAXTAP MAXELT
  write (40,'(a)') '3   2 4 5'                ! KUREL, LVABL
  write (40,'(a)') '3   3 4 5'
  write (40,'(a)') '1   3'
  write (40,'(a)') '0   0'                    ! pseudo-element terminates data
  write (40,'(a)') '-1  9  2. 0. 2. 0. 0. 1.  1. 2. 3.'   ! NGUIDE, LZ, EL
  write (40,'(a)') ' 1  4  3. 2. 4. 0.'
  write (40,'(a)') '-2  5  1. 2.  2. 3. 4.'
  write (40,'(a)') '-1  2  6.  1.'
  write (40,'(a)') '2'                        ! NEWRHS = 2: re-solution follows
  write (40,'(a)') '-1  2  4.  7.'            ! elements in reverse order
  write (40,'(a)') '-1  6  3. 4. 2.  0. 8. 3.'
  write (40,'(a)') ' 2  3  1. 6. 5.'          ! continuation segment first
  write (40,'(a)') '-1  3  4. 7. 2.'          ! 'first' segment, with BED
  write (40,'(a)') '-1'                       ! NEWRHS = -1: a second problem
  write (40,'(a)') '1 1 1 40 20'
  write (40,'(a)') '3   5 89 5'               ! repeated nickname
  write (40,'(a)') '0   0'
  write (40,'(a)') '-1  9  2. 2. 4. 1. 2. 2.  1. 2. 3.'
  write (40,'(a)') '0'                        ! NEWRHS = 0: no more problems
  close (40)

  NCARD = 41; NPRT = 42; NTP1 = 43; NTP2 = 44
  open (NCARD, file='irons_appendix3.dat', status='old')
  open (NPRT,  file='irons_appendix3.lis', status='replace')
  open (NTP1,  form='unformatted', status='scratch')
  open (NTP2,  form='unformatted', status='scratch')
  call zipp()
  close (NCARD); close (NPRT); close (NTP1); close (NTP2)

  print '(a)', ' Irons (1970) Appendix III specimen data, run through ZIPP'
  print '(a)', ' (full printer listing in irons_appendix3.lis)'
  if (jwhere_err /= 0) then
    print '(a,i0)', ' ZIPP stopped with JWHERE = ', jwhere_err
    stop 1
  end if
  if (npass /= 3) then
    print '(a,i0)', ' expected 3 back-substitution passes, got ', npass
    stop 1
  end if

  ! ---- exact systems quoted in the paper -----------------------------------
  A1 = 0
  A1(1,1) = 2
  A1(2,2) = 9; A1(2,3) = 2
  A1(3,2) = 2; A1(3,3) = 6; A1(3,4) = 1
  A1(4,3) = 1; A1(4,4) = 3
  B1(:,1) = (/ 1, 3, 5, 7 /)
  B1(:,2) = (/ 4, 7, 11, 4 /)
  B1(:,3) = (/ 1, 7, 14, 8 /)
  call qsolve(A1, B1, X1, 4, 3)
  A2 = reshape((/ 6, 4, 4, 4 /), (/ 2, 2 /)); B2(:,1) = (/ 4, 2 /)
  call qsolve(A2, B2, X2, 2, 1)

  do i = 1, 4
    xf(i,1) = ans(idxof(n1(i)), 1, 1)                ! problem 1
    xf(i,2) = ans(idxof(n1(i)), 1, 2)                ! re-solution, rhs 1
    xf(i,3) = ans(idxof(n1(i)), 2, 2)                ! re-solution, rhs 2
  end do
  do i = 1, 2
    xf2(i) = ans(idxof(n2(i)), 1, 3)                 ! problem 2
  end do

  worst_err = 0; worst_res = 0
  print '(/a)', ' Problem 1 (output by elements) and its re-solution with 2 RHS'
  print '(a)',  '   nickname   solution        exact         | re-solve RHS 1   exact' // &
                '          | re-solve RHS 2   exact'
  do i = 1, 4
    print '(3x,i5,3(3x,f14.10,1x,f14.10))', n1(i), (xf(i,j), real(X1(i,j),8), j=1,3)
  end do
  do j = 1, 3
    err = maxval(abs(xf(:,j) - real(X1(:,j),8))) / maxval(abs(real(X1(:,j),8)))
    res = maxval(abs(real(B1(:,j) - matmul(A1, real(xf(:,j),16)), 8))) / maxval(abs(real(B1(:,j),8)))
    print '(3x,a,i1,a,es9.2,a,es9.2)', 'RHS set ', j, ':  relative error =', err, &
          '   relative residual ||b-Ax||/||b|| =', res
    worst_err = max(worst_err, err); worst_res = max(worst_res, res)
  end do

  print '(/a)', ' Problem 2 (entire output, element nicknames 5,89,5)'
  print '(a)',  '   nickname   solution        exact'
  do i = 1, 2
    print '(3x,i5,3x,f14.10,1x,f14.10)', n2(i), xf2(i), real(X2(i,1),8)
  end do
  err = maxval(abs(xf2 - real(X2(:,1),8))) / maxval(abs(real(X2(:,1),8)))
  res = maxval(abs(real(B2(:,1) - matmul(A2, real(xf2,16)), 8))) / maxval(abs(real(B2(:,1),8)))
  print '(3x,a,es9.2,a,es9.2)', 'relative error =', err, &
        '   relative residual ||b-Ax||/||b|| =', res
  worst_err = max(worst_err, err); worst_res = max(worst_res, res)
  print '(3x,a,i0,a)', 'entire output vector has MAXNIC = ', maxnic_seen, ' entries' // &
        ' (unused nicknames give zero)'

  print '(/a,i0,a,es9.2)', ' Pivots taken: ', npivots, ';  largest Irons roundoff criterion CRIT = ', critmax
  if (worst_err > 1.d-14 .or. worst_res > 1.d-14) nfail = nfail + 1
  if (any(abs(xf2 - (/ 1.d0, -0.5d0 /)) > 1.d-14)) nfail = nfail + 1
  if (nfail == 0) then
    print '(/a)', ' ALL APPENDIX III CHECKS PASSED'
  else
    print '(/a)', ' APPENDIX III CHECK FAILED'
    stop 1
  end if

contains

  subroutine qsolve(A, B, X, n, nr)        ! Gauss with partial pivoting, quad
    integer, intent(in) :: n, nr
    real(16), intent(in) :: A(n,n), B(n,nr)
    real(16), intent(out) :: X(n,nr)
    real(16) :: W(n,n), Y(n,nr), f, tw(n), ty(nr)
    integer :: i, k, p
    W = A; Y = B
    do k = 1, n-1
      p = k - 1 + maxloc(abs(W(k:n,k)), 1)
      tw = W(k,:); W(k,:) = W(p,:); W(p,:) = tw
      ty = Y(k,:); Y(k,:) = Y(p,:); Y(p,:) = ty
      do i = k+1, n
        f = W(i,k)/W(k,k)
        W(i,k:n) = W(i,k:n) - f*W(k,k:n)
        Y(i,:) = Y(i,:) - f*Y(k,:)
      end do
    end do
    do i = n, 1, -1
      X(i,:) = (Y(i,:) - matmul(W(i,i+1:n), X(i+1:n,:))) / W(i,i)
    end do
  end subroutine qsolve

end program test_irons_deck
