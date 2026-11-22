! Shared test support for the Irons (1970) frontal solver.
!  - module zcap: the ZIPP blank COMMON and /ZIPIO/ (same layout as in
!    zipp.f), a nickname <-> index map, and storage for the answers that
!    ZIPP hands back through FRNOUT.
!  - subroutine FRNOUT: the hook called from ZIPP (see note 6 in zipp.f).
module zcap
  implicit none
  ! ---- ZIPP blank common (zipp.f) --------------------------------------
  real(8) :: ELPA(3000), ELCOR(3,20)
  integer :: LVABL(30), LDEST(30), MVABL(60)
  integer :: INITL, NTIREX, NEWRHS, NELEM, NELEMZ, KUREL, LPREQ, &
             LZ, NELZ, NBAXO, NBZ, KL, LDES, NSTRES
  common ELPA, LVABL, LDEST, MVABL, ELCOR, INITL, NTIREX, NEWRHS, &
         NELEM, NELEMZ, KUREL, LPREQ, LZ, NELZ, NBAXO, NBZ, KL, LDES, NSTRES
  integer :: NCARD, NPRT, NTP1, NTP2
  common /ZIPIO/ NCARD, NPRT, NTP1, NTP2

  ! ---- capture storage ---------------------------------------------------
  integer, parameter :: MXV = 600, MXR = 6, MXP = 10, MXNICK = 100000
  integer :: nvar = 0, nickof(MXV), idxof(MXNICK) = 0
  real(8) :: ans(MXV, MXR, MXP)            ! answers by variable, rhs, pass
  logical :: got(MXV, MXR, MXP)
  real(8) :: entire(3000, MXP)             ! raw ELPA for entire output
  integer :: npass = 0, nrhs_pass(MXP) = 0, kind_pass(MXP) = 0
  integer :: maxpa_seen = 0, maxnic_seen = 0, jwhere_err = 0
  integer :: ncritwarn = 0, npivots = 0
  real(8) :: critmax = 0.d0
  integer :: nic_critmax = 0

contains

  subroutine reset_capture()
    nvar = 0; idxof = 0; got = .false.; ans = 0.d0; entire = 0.d0
    npass = 0; nrhs_pass = 0; kind_pass = 0
    maxpa_seen = 0; maxnic_seen = 0; jwhere_err = 0
    ncritwarn = 0; npivots = 0; critmax = 0.d0; nic_critmax = 0
  end subroutine reset_capture

  integer function reg(nic)                ! register a nickname, return its index
    integer, intent(in) :: nic
    if (idxof(nic) == 0) then
      nvar = nvar + 1
      nickof(nvar) = nic
      idxof(nic) = nvar
    end if
    reg = idxof(nic)
  end function reg

end module zcap

! ----------------------------------------------------------------------------
subroutine frnout(kind, nrhs, maxnic, maxpa, jwhere, crit, pivot, nic)
  ! KIND = 1 output by elements: answers for element NELEM are in
  !          ELPA(KL + (J-1)*KUREL), nicknames in LVABL
  !        2 entire output: answer for nickname N, rhs J is ELPA(N+(J-1)*MAXNIC)
  !        3 fatal error JWHERE
  !        5 a pivot has been taken; CRIT is Irons' roundoff criterion
  use zcap
  implicit none
  integer, intent(in) :: kind, nrhs, maxnic, maxpa, jwhere, nic
  real(8), intent(in) :: crit, pivot
  integer :: k, j, i
  maxpa_seen = maxpa
  maxnic_seen = maxnic
  select case (kind)
  case (1)
    if (NELEM == NELEMZ) then                ! first element met in back-substitution
      npass = npass + 1
      nrhs_pass(npass) = nrhs
      kind_pass(npass) = 1
    end if
    do k = 1, KUREL
      i = reg(LVABL(k))
      do j = 1, nrhs
        ans(i, j, npass) = ELPA(k + (j-1)*KUREL)
        got(i, j, npass) = .true.
      end do
    end do
  case (2)
    npass = npass + 1
    nrhs_pass(npass) = nrhs
    kind_pass(npass) = 2
    entire(1:maxnic*nrhs, npass) = ELPA(1:maxnic*nrhs)
    do i = 1, nvar
      if (nickof(i) <= maxnic) then
        do j = 1, nrhs
          ans(i, j, npass) = ELPA(nickof(i) + (j-1)*maxnic)
          got(i, j, npass) = .true.
        end do
      end if
    end do
  case (3)
    jwhere_err = jwhere
  case (5)
    npivots = npivots + 1
    if (crit > critmax) then
      critmax = crit
      nic_critmax = nic
    end if
    if (crit > 1.d4 .or. pivot < 0.d0) ncritwarn = ncritwarn + 1
  end select
end subroutine frnout
