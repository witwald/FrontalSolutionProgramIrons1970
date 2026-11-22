! Test 2 for the Irons (1970) frontal solver ZIPP: finite element models.
!
! ELMAK here is a user element routine, as the paper intends ("in a real
! application it would be created in ELMAK"). In the pre-program it writes
! each element's stiffness and loads to tape 1, split into segments of at most
! MAXELT words. During back-substitution it uses the re-solution facility for
! iterative refinement, as the paper suggests: as soon as an element's answers
! appear it computes the element residual R_e = X_e - C_e x_e (accumulated in
! quad precision) and writes it to tape 1 in reversed segment order as the
! right-hand side for the next re-solution.
!
!  A  Hoit (1986) Figure 1 frame: 18 equations, element order and DOF
!     numbers of the figure; lateral and gravity load cases; entire output
!     and output by elements.
!  B  3-bay x 15-storey frame (180 equations) with springs, one spring having
!     a repeated nickname; 3 load cases.
!       B1  sequential nicknames, large buffer, whole elements
!       B2  scrambled nicknames up to 99991, 60-word buffer (many tape
!           blocks), MAXELT = 20 (every frame element fragmented)
!       B3  nearly rigid beams (ill-conditioned): Irons' CRIT warning and
!           two passes of iterative refinement by re-solution
!  C  Cantilever of 8-node quadrilaterals (the elements with midside nodes
!     that motivated the frontal method), supports by the 10**20 penalty
!     of Appendix IV 1(b) with reactions recovered as R = -10**20 x.
!     Checked against beam theory and by equilibrium.
!  D  Unsupported frame ("flying structure"): ZIPP must stop with its
!     JWHERE = 13 or 14 diagnostic.
!
! Accuracy measures for every solution (all accumulated in quad precision):
!   residual   max over load cases of ||b - A x|| / ||b||
!   error      ||x - x_ref|| / ||x_ref||, x_ref = dense solve in quad precision
!   correction ||dx|| / ||x|| from the next refinement pass (an error estimate
!              that needs no reference solution)
!   CRIT       Irons' diagonal-decay roundoff criterion (largest over pivots)

module femodel
  implicit none
  integer, parameter :: MXE = 300, MXK = 16, MXL = 3
  integer :: nel = 0, nload = 0, segmax = 1000
  integer :: kur(MXE), nik(MXK, MXE)
  real(8) :: ck(MXK, MXK, MXE), xr(MXK, MXL, MXE)
  integer :: nseg_written = 0, nres_written = 0
contains
  subroutine add_elem(k, nicks, c, r)
    integer, intent(in) :: k, nicks(k)
    real(8), intent(in) :: c(k,k), r(k,MXL)
    nel = nel + 1
    kur(nel) = k
    nik(1:k, nel) = nicks
    ck(:,:,nel) = 0.d0; ck(1:k,1:k,nel) = c
    xr(:,:,nel) = 0.d0; xr(1:k,:,nel) = r
  end subroutine add_elem
end module femodel

! ----------------------------------------------------------------------------
subroutine elmak
  use zcap
  use femodel
  implicit none
  integer :: e, k, i, l, j, n, ns, s, i0, i1, q, ix
  real(8) :: el(3000)
  real(16) :: xq(MXK), r
  e = NELEM
  k = KUREL
  if (INITL /= 0) then
    ! ---- pre-program: element stiffness (upper triangle by columns) + loads
    n = 0
    do j = 1, k
      do i = 1, j
        n = n + 1; el(n) = ck(i, j, e)
      end do
    end do
    do j = 1, NEWRHS
      do i = 1, k
        n = n + 1; el(n) = xr(i, j, e)
      end do
    end do
    ns = (n + segmax - 1)/segmax
    do s = 1, ns                                   ! segments in natural order
      i0 = (s-1)*segmax + 1; i1 = min(s*segmax, n)
      if (s == 1) then
        write (NTP1) KUREL, LPREQ, (LVABL(i), LDEST(i), i = 1, KUREL), i1-i0+1, el(i0:i1)
      else
        write (NTP1) i1-i0+1, el(i0:i1)
      end if
    end do
    nseg_written = nseg_written + ns
  else
    ! ---- back-substitution: residual right-hand sides for the re-solution
    if (NEWRHS <= 0) return
    n = 0
    do j = 1, NEWRHS
      do i = 1, k                                  ! accumulated solution so far
        ix = idxof(LVABL(i))
        xq(i) = 0
        do q = 1, npass
          xq(i) = xq(i) + real(ans(ix, j, q), 16)
        end do
      end do
      do i = 1, k
        r = real(xr(i, j, e), 16)
        do l = 1, k
          r = r - real(ck(i, l, e), 16)*xq(l)
        end do
        n = n + 1; el(n) = real(r, 8)
      end do
    end do
    ns = (n + segmax - 1)/segmax
    do s = ns, 1, -1                               ! reversed: BED segment last
      i0 = (s-1)*segmax + 1; i1 = min(s*segmax, n)
      if (s == 1) then
        write (NTP1) KUREL, LPREQ, (LVABL(i), LDEST(i), i = 1, KUREL), i1-i0+1, el(i0:i1)
      else
        write (NTP1) i1-i0+1, el(i0:i1)
      end if
    end do
    nres_written = nres_written + ns
  end if
end subroutine elmak

! ============================================================================

program test_irons_fe
  use zcap
  use femodel
  implicit none
  integer :: nfail = 0
  integer(8) :: iseed = 12345
  real(16), allocatable :: Aq(:,:), Bq(:,:), Xq(:,:)
  integer :: nq
  ! frame geometry shared by builders
  real(8) :: xn(200), yn(200)
  integer :: idn(200,3)
  real(8) :: xq_saved(MXV, MXR), last_len(MXE+1)
  integer :: semiband_seq = 0, semiband_scr = 0, semiband_last = 0

  print '(a)', ' Irons (1970) frontal solver ZIPP - finite element tests'
  print '(a)', ' (accuracy measures are accumulated in quad precision; listings in irons_fe_*.lis)'

  call test_a()
  call test_b()
  call test_c()
  call test_d()

  print *
  if (nfail == 0) then
    print '(a)', ' ALL FINITE ELEMENT CHECKS PASSED'
  else
    print '(a,i0,a)', ' ', nfail, ' CHECK(S) FAILED'
    stop 1
  end if

contains

  ! ==========================================================================
  ! A: Hoit (1986) Figure 1 frame
  ! ==========================================================================
  subroutine test_a()
    integer :: nd, lc
    real(8) :: rel
    print '(/a)', ' A. Two-bay, two-storey frame of Hoit (1986) Fig. 1 (18 equations, 10 elements)'
    call build_fig1()
    call run('a1', ntirex=1, nrmax=2, maxtap=2000, maxelt=100, nref=0, label='entire output')
    call report(0, 1.d-13, .true.)
    print '(a)', '     Displacements (entire output):          load case 1: lateral        ' // &
                 '        load case 2: gravity on beams'
    print '(a)', '     node      u            v          theta       |      u            v          theta'
    do nd = 4, 9
      print '(5x,a,1x,3es13.5,3x,a,3es13.5)', node_name(nd), &
        (ans(idxof(idn(nd,lc)), 1, 1), lc = 1, 3), '|', (ans(idxof(idn(nd,lc)), 2, 1), lc = 1, 3)
    end do
    call run('a2', ntirex=0, nrmax=2, maxtap=2000, maxelt=100, nref=0, label='output by elements')
    call report(0, 1.d-13, .true.)
    ! the two output modes must agree exactly
    rel = maxval(abs(ans(1:nvar,1:2,1) - xq_saved(1:nvar,1:2)))
    call check('A  entire output and output by elements give identical answers', rel == 0.d0)
  end subroutine test_a

  character(len=10) function node_name(nd)
    integer, intent(in) :: nd
    character(len=10), parameter :: nm(9) = (/ 'sup-left  ', 'sup-mid   ', 'sup-right ', &
      'mid-left  ', 'mid-mid   ', 'mid-right ', 'top-left  ', 'top-mid   ', 'top-right ' /)
    node_name = nm(nd)
  end function node_name

  subroutine build_fig1()
    integer :: e
    nel = 0; nload = 2
    ! nodes 1-3 supports; DOF numbers (u,v,theta) as in Hoit's Figure 1
    xn(1:9) = (/ 0.d0, 6.d0, 12.d0, 0.d0, 6.d0, 12.d0, 0.d0, 6.d0, 12.d0 /)
    yn(1:9) = (/ 0.d0, 0.d0, 0.d0, 4.d0, 4.d0, 4.d0, 8.d0, 8.d0, 8.d0 /)
    idn = 0
    idn(4,:) = (/ 16,17,18 /); idn(5,:) = (/ 7,8,9 /);    idn(6,:) = (/ 4,5,6 /)
    idn(7,:) = (/ 13,14,15 /); idn(8,:) = (/ 10,11,12 /); idn(9,:) = (/ 1,2,3 /)
    ! columns 1-6 then beams 7-10: the element numbering of the figure
    call frame_elem(1,4,1); call frame_elem(4,7,1); call frame_elem(2,5,1)
    call frame_elem(5,8,1); call frame_elem(3,6,1); call frame_elem(6,9,1)
    call frame_elem(4,5,2); call frame_elem(7,8,2); call frame_elem(5,6,2)
    call frame_elem(8,9,2)
    ! load case 1: lateral loads 10 at roof, 5 at first floor (left side)
    call nodal_load(13, 1, 10.d0)
    call nodal_load(16, 1, 5.d0)
    ! load case 2: uniform gravity load w = 2 on all beams (fixed-end forces)
    do e = 7, 10
      call beam_udl(e, 2, 2.d0)
    end do
  end subroutine build_fig1

  ! ==========================================================================
  ! B: 3-bay x 15-storey frame with springs
  ! ==========================================================================
  
  subroutine test_b()
    integer :: ib, k, e
    real(8) :: xsave(MXV, MXR), d
    integer :: maxpa1, ib_nick(MXV)
    print '(/a)', ' B. Three-bay, fifteen-storey frame with springs (180 equations, 108 elements, 3 load cases)'

    call build_grid(3, 15, 1.d0, .false.)
    call run('b1', ntirex=0, nrmax=3, maxtap=2400, maxelt=1000, nref=1, &
             label='sequential nicknames, whole elements')
    call report(1, 1.d-12, .true.)
    maxpa1 = maxpa_seen
    do k = 1, nvar
      ib_nick(k) = nickof(k)
    end do
    xsave(1:nvar,1:3) = ans(1:nvar,1:3,1)

    call build_grid(3, 15, 1.d0, .true.)
    call run('b2', ntirex=0, nrmax=3, maxtap=60, maxelt=20, nref=1, &
             label='scrambled nicknames, 60-word buffer, MAXELT=20')
    call report(1, 1.d-12, .true.)
    print '(5x,a,i0,a,i0,a)', 'elements written as ', nseg_written, ' tape segments (', nel, ' elements)'
    ! same physical answers whatever the nicknames, and the same front
    d = 0.d0
    do k = 1, nvar
      d = max(d, maxval(abs(ans(idxof(scramble(ib_nick(k))), 1:3, 1) - xsave(k, 1:3))))
    end do
    d = d / maxval(abs(xsave(1:nvar,1:3)))
    print '(5x,a,i0,a,i0,a,i0,a)', 'largest front (MAXPA) = ', maxpa_seen, ' for both numberings;' // &
          ' semi-bandwidth: ', semiband_seq, ' (sequential), ', semiband_scr, ' (scrambled)'
    print '(5x,a,es9.2)', 'answers, scrambled vs sequential nicknames: max relative difference ', d
    call check('B  front width independent of nickname order (Irons, note (c))', maxpa_seen == maxpa1)
    call check('B  same answers for scrambled nicknames and fragmented elements', d < 1.d-12)

    call build_grid(3, 15, 1.d4, .false.)
    call run('b3', ntirex=0, nrmax=3, maxtap=2400, maxelt=1000, nref=2, &
             label='nearly rigid beams (EA x 10**4), 2 refinement passes')
    call report(2, 1.d-12, .false.)
    call check('B3 Irons CRIT warns of possible roundoff trouble (CRIT > 1E4)', ncritwarn > 0)
    e = 0
  end subroutine test_b

  integer function scramble(k)                     ! one-to-one map onto 1..99991
    integer, intent(in) :: k
    scramble = int(mod(int(k,8)*48271_8, 99991_8)) + 1
  end function scramble

  subroutine build_grid(nbay, nsto, beamfac, scr)
    integer, intent(in) :: nbay, nsto
    real(8), intent(in) :: beamfac
    logical, intent(in) :: scr
    integer :: i, j, nd, eq, e1, e2, e
    real(8) :: c(2,2), r(2,MXL)
    iseed = 12345                                  ! same random loads every build
    nel = 0; nload = 3; eq = 0; idn = 0
    do j = 0, nsto
      do i = 0, nbay
        nd = j*(nbay+1) + i + 1
        xn(nd) = 5.d0*i + 0.3d0*j
        yn(nd) = 3.5d0*j
        if (j > 0) then
          idn(nd,:) = (/ eq+1, eq+2, eq+3 /); eq = eq + 3
        end if
      end do
    end do
    if (scr) then
      do nd = 1, (nbay+1)*(nsto+1)
        do i = 1, 3
          if (idn(nd,i) > 0) idn(nd,i) = scramble(idn(nd,i))
        end do
      end do
    end if
    do j = 1, nsto                                 ! storey by storey: the front
      do i = 0, nbay                               ! advances up the building
        call frame_elem((j-1)*(nbay+1)+i+1, j*(nbay+1)+i+1, 1)
        if (i == 0) call nodal_load(idn(j*(nbay+1)+1, 1), 1, 2.d0 + 0.2d0*j)  ! wind
      end do
      do i = 1, nbay
        call frame_elem(j*(nbay+1)+i, j*(nbay+1)+i+1, 2, beamfac)
        e = nel
        call beam_udl(e, 2, 3.d0)                  ! gravity
        call random_rhs(e, 3)                      ! case 3: random element loads
      end do
      if (j == 1) then                             ! long-range spring
        e1 = idn(nsto*(nbay+1)+1, 1); e2 = idn(1*(nbay+1)+nbay+1, 1)
        c = reshape((/ 4.d2, -1.5d2, -1.5d2, 3.d2 /), (/2,2/)); r = 0.d0
        r(:,3) = (/ 0.5d0, -0.25d0 /)
        call add_elem(2, (/ e1, e2 /), c, r)
      end if
      if (j == 8) then                             ! repeated nickname
        e1 = idn(8*(nbay+1)+2, 2)
        c = reshape((/ 7.d2, -2.d2, -2.d2, 5.d2 /), (/2,2/)); r = 0.d0
        r(:,3) = (/ 0.3d0, 0.2d0 /)
        call add_elem(2, (/ e1, e1 /), c, r)
      end if
      if (j == 6) then                             ! short rotational spring
        e1 = idn(5*(nbay+1)+3, 3); e2 = idn(6*(nbay+1)+3, 3)
        c = reshape((/ 9.d1, -3.d1, -3.d1, 6.d1 /), (/2,2/)); r = 0.d0
        call add_elem(2, (/ e1, e2 /), c, r)
      end if
    end do
  end subroutine build_grid

  ! ==========================================================================
  ! C: 8-node quadrilateral cantilever with penalty supports
  ! ==========================================================================
  
  subroutine test_c()
    integer, parameter :: NX = 10, NY = 2
    real(8), parameter :: L = 10.d0, H = 2.d0, E = 1000.d0, NU = 0.3d0, P = 1.d0
    real(8), parameter :: PEN = 1.d20
    integer :: nodeid(0:2*NX, 0:2*NY), nn, i, j, ex, ey, a, nodes(8), nk(16)
    real(8) :: xy(2,8), kel(16,16), rel(16,MXL), c1(1,1), r1(1,MXL)
    real(8) :: ei, gk, vtip, vtg, utip, sumr(2,2), appl(2,2), xs, vb, vfe, xcol
    real(8) :: rx, ry
    integer :: ns, k
    print '(/a)', ' C. Cantilever, 10 x 2 eight-node quadrilaterals, plane stress, penalty supports'
    nel = 0; nload = 2
    ! node numbering column by column
    nn = 0
    do i = 0, 2*NX
      do j = 0, 2*NY
        if (mod(i,2) == 1 .and. mod(j,2) == 1) then
          nodeid(i,j) = 0
        else
          nn = nn + 1; nodeid(i,j) = nn
        end if
      end do
    end do
    ! supports first: single-variable elements C = 10**20 (Appendix IV 1(b))
    c1 = PEN; r1 = 0.d0
    do j = 0, 2*NY
      call add_elem(1, (/ 2*nodeid(0,j)-1 /), c1, r1)
      call add_elem(1, (/ 2*nodeid(0,j)   /), c1, r1)
    end do
    ns = nel
    ! elements column by column
    do ex = 1, NX
      do ey = 1, NY
        nodes = (/ nodeid(2*ex-2,2*ey-2), nodeid(2*ex,2*ey-2), nodeid(2*ex,2*ey), &
                   nodeid(2*ex-2,2*ey),   nodeid(2*ex-1,2*ey-2), nodeid(2*ex,2*ey-1), &
                   nodeid(2*ex-1,2*ey),   nodeid(2*ex-2,2*ey-1) /)
        xy(1,:) = (/ ex-1.d0, ex*1.d0, ex*1.d0, ex-1.d0, ex-0.5d0, ex*1.d0, ex-0.5d0, ex-1.d0 /) * (L/NX)
        xy(2,:) = (/ ey-1.d0, ey-1.d0, ey*1.d0, ey*1.d0, ey-1.d0, ey-0.5d0, ey*1.d0, ey-0.5d0 /) * (H/NY)
        call quad8_k(xy, E, NU, kel)
        rel = 0.d0
        if (ex == NX) call tip_tractions(xy, P, H, rel)
        do a = 1, 8
          nk(2*a-1) = 2*nodes(a)-1; nk(2*a) = 2*nodes(a)
        end do
        call add_elem(16, nk, kel, rel)
      end do
    end do
    call run('c', ntirex=0, nrmax=2, maxtap=2400, maxelt=40, nref=1, &
             label='MAXELT=40: each element in 5 segments')
    call report(1, 1.d-12, .true.)
    print '(5x,a,i0,a,i0)', 'largest front (MAXPA) = ', maxpa_seen, ';  semi-bandwidth of the numbering = ', &
          semiband_seq
    ! beam theory: Timoshenko & Goodier end deflection with parabolic shear
    ei = E*H**3/12.d0
    vtg = P*L**3/(3.d0*ei) + (4.d0+5.d0*NU)*P*L/(2.d0*E*H)
    k = idxof(2*nodeid(2*NX,NY))
    vtip = -(ans(k,1,1) + ans(k,1,2))
    k = idxof(2*nodeid(2*NX,NY)-1)
    utip = ans(k,2,1) + ans(k,2,2)
    print '(5x,a)', 'Load case 1, tip shear P = 1 (parabolic), deflection of the beam axis:'
    print '(5x,a)', '     x       frontal FE       beam theory     ratio'
    gk = E/(2.d0*(1.d0+NU))
    do i = 2, 2*NX, 4
      xcol = i*(L/(2*NX))
      k = idxof(2*nodeid(i,NY))
      vfe = -(ans(k,1,1) + ans(k,1,2))
      vb = P*xcol**2*(3.d0*L - xcol)/(6.d0*ei) + (4.d0+5.d0*NU)*P*xcol/(2.d0*E*H)
      print '(5x,f7.2,2x,es14.6,2x,es14.6,2x,f7.4)', xcol, vfe, vb, vfe/vb
    end do
    xcol = L
    print '(5x,f7.2,2x,es14.6,2x,es14.6,2x,f7.4,a)', xcol, vtip, vtg, vtip/vtg, '  <- tip'
    print '(5x,a,es12.5,a,es12.5,a,f7.4)', 'Load case 2, tip tension P = 1: tip extension ', utip, &
          ';  PL/EA = ', P*L/(E*H), ';  ratio ', utip/(P*L/(E*H))
    call check('C  tip deflection within 2% of Timoshenko-Goodier', abs(vtip/vtg - 1.d0) < 0.02d0)
    call check('C  tip extension within 3% of PL/EA', abs(utip/(P*L/(E*H)) - 1.d0) < 0.03d0)
    ! reactions R = -10**20 x at supported variables; equilibrium with the loads
    sumr = 0.d0; appl = 0.d0
    do k = 1, ns
      a = nik(1,k)
      i = 2 - mod(a,2)                             ! 1 = x (odd nickname), 2 = y
      do j = 1, 2
        sumr(i,j) = sumr(i,j) - PEN*(ans(idxof(a),j,1) + ans(idxof(a),j,2))
      end do
    end do
    do k = ns+1, nel
      do a = 1, kur(k)
        i = 2 - mod(nik(a,k),2)
        appl(i,:) = appl(i,:) + xr(a,1:2,k)
      end do
    end do
    print '(5x,a)', 'Reactions recovered as R = -1E20 x (Appendix IV 1(b)) against the applied loads:'
    print '(5x,a)', '              sum Rx          sum Ry          applied Fx      applied Fy'
    do j = 1, 2
      print '(5x,a,i1,1x,4es16.8)', 'load case ', j, sumr(1,j), sumr(2,j), appl(1,j), appl(2,j)
    end do
    rx = maxval(abs(sumr + appl))
    print '(5x,a,es9.2)', 'equilibrium error max|sum R + sum F| = ', rx
    call check('C  reactions balance the applied loads (|R+F| < 1E-10)', rx < 1.d-10)
    ry = 0
  end subroutine test_c

  subroutine quad8_k(xy, E, NU, kel)                ! 8-node serendipity, 3x3 Gauss
    real(8), intent(in) :: xy(2,8), E, NU
    real(8), intent(out) :: kel(16,16)
    real(8), parameter :: g(3) = (/ -0.7745966692414834d0, 0.d0, 0.7745966692414834d0 /)
    real(8), parameter :: w(3) = (/ 5.d0/9.d0, 8.d0/9.d0, 5.d0/9.d0 /)
    real(8) :: D(3,3), dn(2,8), jac(2,2), ji(2,2), dj, bm(3,16), dxy(2,8)
    integer :: ig, jg, a
    D = 0.d0
    D(1,1) = 1.d0; D(1,2) = NU; D(2,1) = NU; D(2,2) = 1.d0; D(3,3) = (1.d0-NU)/2.d0
    D = D*E/(1.d0-NU**2)
    kel = 0.d0
    do ig = 1, 3
      do jg = 1, 3
        call q8_deriv(g(ig), g(jg), dn)
        jac = matmul(dn, transpose(xy))
        dj = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
        ji = reshape((/ jac(2,2), -jac(2,1), -jac(1,2), jac(1,1) /), (/2,2/)) / dj
        dxy = matmul(ji, dn)
        bm = 0.d0
        do a = 1, 8
          bm(1,2*a-1) = dxy(1,a); bm(2,2*a) = dxy(2,a)
          bm(3,2*a-1) = dxy(2,a); bm(3,2*a) = dxy(1,a)
        end do
        kel = kel + matmul(transpose(bm), matmul(D, bm)) * dj * w(ig)*w(jg)
      end do
    end do
  end subroutine quad8_k

  subroutine q8_deriv(xi, eta, dn)
    real(8), intent(in) :: xi, eta
    real(8), intent(out) :: dn(2,8)
    real(8), parameter :: xa(8) = (/ -1.d0, 1.d0, 1.d0, -1.d0, 0.d0, 1.d0, 0.d0, -1.d0 /)
    real(8), parameter :: ya(8) = (/ -1.d0, -1.d0, 1.d0, 1.d0, -1.d0, 0.d0, 1.d0, 0.d0 /)
    integer :: a
    do a = 1, 8
      if (a <= 4) then
        dn(1,a) = xa(a)*(1.d0+eta*ya(a))*(2.d0*xi*xa(a)+eta*ya(a))/4.d0
        dn(2,a) = ya(a)*(1.d0+xi*xa(a))*(xi*xa(a)+2.d0*eta*ya(a))/4.d0
      else if (xa(a) == 0.d0) then
        dn(1,a) = -xi*(1.d0+eta*ya(a))
        dn(2,a) = ya(a)*(1.d0-xi**2)/2.d0
      else
        dn(1,a) = xa(a)*(1.d0-eta**2)/2.d0
        dn(2,a) = -eta*(1.d0+xi*xa(a))
      end if
    end do
  end subroutine q8_deriv

  subroutine tip_tractions(xy, P, H, rel)          ! right edge = nodes 2, 6, 3
    real(8), intent(in) :: xy(2,8), P, H
    real(8), intent(inout) :: rel(16,MXL)
    real(8), parameter :: g(3) = (/ -0.7745966692414834d0, 0.d0, 0.7745966692414834d0 /)
    real(8), parameter :: w(3) = (/ 5.d0/9.d0, 8.d0/9.d0, 5.d0/9.d0 /)
    integer :: ig, a
    integer, parameter :: en(3) = (/ 2, 6, 3 /)
    real(8) :: et, nf(3), y, ds, yc, c, tau
    ds = (xy(2,3) - xy(2,2))/2.d0
    c = H/2.d0
    do ig = 1, 3
      et = g(ig)
      nf = (/ et*(et-1.d0)/2.d0, 1.d0-et**2, et*(et+1.d0)/2.d0 /)
      y = sum(nf*(/ xy(2,2), xy(2,6), xy(2,3) /))
      yc = y - c
      tau = 1.5d0*P/H*(1.d0 - (yc/c)**2)          ! parabolic shear, total P
      do a = 1, 3
        rel(2*en(a),   1) = rel(2*en(a),   1) - nf(a)*tau*ds*w(ig)   ! downward
        rel(2*en(a)-1, 2) = rel(2*en(a)-1, 2) + nf(a)*(P/H)*ds*w(ig) ! tension
      end do
    end do
  end subroutine tip_tractions

  ! ==========================================================================
  ! D: an unsupported structure must be diagnosed
  ! ==========================================================================
  
  subroutine test_d()
    integer :: e
    print '(/a)', ' D. Fig. 1 frame with its supports removed ("flying structure")'
    call build_fig1()
    ! give the bottom nodes free DOFs so nothing is earthed
    idn(1,:) = (/ 19,20,21 /); idn(2,:) = (/ 22,23,24 /); idn(3,:) = (/ 25,26,27 /)
    nel = 0
    call frame_elem(1,4,1); call frame_elem(4,7,1); call frame_elem(2,5,1)
    call frame_elem(5,8,1); call frame_elem(3,6,1); call frame_elem(6,9,1)
    call frame_elem(4,5,2); call frame_elem(7,8,2); call frame_elem(5,6,2)
    call frame_elem(8,9,2)
    call nodal_load(13, 1, 10.d0)
    call run('d', ntirex=1, nrmax=2, maxtap=2000, maxelt=100, nref=0, label='expect JWHERE 13 or 14', &
             expect_error=.true.)
    print '(5x,a,i0,a,es9.2)', 'ZIPP stopped with JWHERE = ', jwhere_err, ',  CRIT = ', critmax
    call check('D  singular stiffness diagnosed (JWHERE = 13 or 14)', jwhere_err == 13 .or. jwhere_err == 14)
    e = 0
  end subroutine test_d

  ! ==========================================================================
  ! Frame element helpers
  ! ==========================================================================
  
  subroutine frame_elem(n1, n2, grp, fac)
    integer, intent(in) :: n1, n2, grp
    real(8), intent(in), optional :: fac
    real(8) :: sk(6,6), Ea, AR, ZI, c(6,6), r(6,MXL)
    integer :: dof(6), keep(6), k, i
    Ea = 2.0d5
    if (grp == 1) then
      AR = 1.2d-2; ZI = 2.0d-4
    else
      AR = 8.0d-3; ZI = 3.5d-4
      if (present(fac)) AR = AR*fac
    end if
    call frame_k(xn(n1), yn(n1), xn(n2), yn(n2), Ea, AR, ZI, sk)
    dof = (/ idn(n1,1), idn(n1,2), idn(n1,3), idn(n2,1), idn(n2,2), idn(n2,3) /)
    k = 0
    do i = 1, 6                                     ! supported DOFs are deleted
      if (dof(i) > 0) then                          ! (Appendix IV, 1(a))
        k = k + 1; keep(k) = i
      end if
    end do
    c(1:k,1:k) = sk(keep(1:k), keep(1:k))
    r = 0.d0
    call add_elem(k, dof(keep(1:k)), c(1:k,1:k), r(1:k,:))
  end subroutine frame_elem

  subroutine nodal_load(nick, lc, f)                ! add to the latest element
    integer, intent(in) :: nick, lc                 ! that contains the variable
    real(8), intent(in) :: f
    integer :: e, a
    do e = nel, 1, -1
      do a = 1, kur(e)
        if (nik(a,e) == nick) then
          xr(a,lc,e) = xr(a,lc,e) + f
          return
        end if
      end do
    end do
  end subroutine nodal_load

  subroutine beam_udl(e, lc, w)                      ! horizontal beam, downward w
    integer, intent(in) :: e, lc
    real(8), intent(in) :: w
    real(8) :: L, fe(6)
    integer :: a
    L = beam_len(e)
    fe = (/ 0.d0, -w*L/2.d0, -w*L**2/12.d0, 0.d0, -w*L/2.d0, w*L**2/12.d0 /)
    do a = 1, 6
      xr(a,lc,e) = xr(a,lc,e) + fe(a)
    end do
  end subroutine beam_udl

  real(8) function beam_len(e)
    integer, intent(in) :: e
    beam_len = last_len(e)
  end function beam_len

  subroutine random_rhs(e, lc)
    integer, intent(in) :: e, lc
    integer :: a
    do a = 1, kur(e)
      xr(a,lc,e) = 4.d0*(rnd() - 0.5d0)
    end do
  end subroutine random_rhs

  subroutine frame_k(x1, y1, x2, y2, E, AR, ZI, SK)
    real(8), intent(in) :: x1, y1, x2, y2, E, AR, ZI
    real(8), intent(out) :: SK(6,6)
    real(8) :: L, c, s, KL(6,6), T(6,6)
    L = sqrt((x2-x1)**2 + (y2-y1)**2); c = (x2-x1)/L; s = (y2-y1)/L
    last_len(nel+1) = L
    KL = 0.d0
    KL(1,1) = E*AR/L;  KL(4,4) = KL(1,1); KL(1,4) = -KL(1,1); KL(4,1) = KL(1,4)
    KL(2,2) = 12*E*ZI/L**3; KL(5,5) = KL(2,2); KL(2,5) = -KL(2,2); KL(5,2) = KL(2,5)
    KL(2,3) = 6*E*ZI/L**2; KL(3,2) = KL(2,3); KL(2,6) = KL(2,3); KL(6,2) = KL(2,3)
    KL(3,5) = -KL(2,3); KL(5,3) = KL(3,5); KL(5,6) = -KL(2,3); KL(6,5) = KL(5,6)
    KL(3,3) = 4*E*ZI/L; KL(6,6) = KL(3,3); KL(3,6) = 2*E*ZI/L; KL(6,3) = KL(3,6)
    T = 0.d0
    T(1,1) = c; T(1,2) = s; T(2,1) = -s; T(2,2) = c; T(3,3) = 1.d0
    T(4:6,4:6) = T(1:3,1:3)
    SK = matmul(transpose(T), matmul(KL, T))
  end subroutine frame_k

  double precision function rnd()
    iseed = mod(iseed*16807_8, 2147483647_8)
    rnd = dble(iseed)/2147483647.d0
  end function rnd

  ! ==========================================================================
  ! Run ZIPP on the current model
  ! ==========================================================================
  
  subroutine run(tag, ntirex, nrmax, maxtap, maxelt, nref, label, expect_error)
    character(len=*), intent(in) :: tag, label
    integer, intent(in) :: ntirex, nrmax, maxtap, maxelt, nref
    logical, intent(in), optional :: expect_error
    integer :: e, a, k
    character(len=64) :: fdeck, flist
    call reset_capture()
    do e = 1, nel
      do a = 1, kur(e)
        k = reg(nik(a,e))
      end do
    end do
    segmax = maxelt
    nseg_written = 0; nres_written = 0
    semiband_last = semiband()
    if (tag(1:1) == 'b' .and. tag /= 'b2') semiband_seq = semiband_last
    if (tag == 'b2') semiband_scr = semiband_last
    if (tag == 'c') semiband_seq = semiband_last
    fdeck = 'irons_fe_' // tag // '.dat'
    flist = 'irons_fe_' // tag // '.lis'
    open (50, file=trim(fdeck), status='replace')
    write (50,'(5i8)') nload, ntirex, nrmax, maxtap, maxelt
    do e = 1, nel
      write (50,'(i4,16i7)') kur(e), nik(1:kur(e), e)
    end do
    write (50,'(a)') '0 0'
    do k = 1, nref
      write (50,'(i4)') nload                      ! residual re-solutions
    end do
    write (50,'(i4)') 0                            ! no more problems
    close (50)
    NCARD = 51; NPRT = 52; NTP1 = 53; NTP2 = 54
    open (NCARD, file=trim(fdeck), status='old')
    open (NPRT,  file=trim(flist), status='replace')
    open (NTP1,  form='unformatted', status='scratch')
    open (NTP2,  form='unformatted', status='scratch')
    call zipp()
    close (NCARD); close (NPRT); close (NTP1); close (NTP2)
    print '(3x,a,1x,a,a,a)', tag, '(', trim(label), ')'
    if (.not. present(expect_error)) then
      if (jwhere_err /= 0) then
        print '(5x,a,i0)', 'ZIPP stopped with JWHERE = ', jwhere_err
        nfail = nfail + 1
        return
      end if
      call dense_reference()
      xq_saved(1:nvar,1:nload) = ans(1:nvar,1:nload,1)
    end if
  end subroutine run

  integer function semiband()
    integer :: e
    semiband = 0
    do e = 1, nel
      semiband = max(semiband, maxval(nik(1:kur(e),e)) - minval(nik(1:kur(e),e)))
    end do
  end function semiband

  subroutine dense_reference()                     ! quad precision A, B, x_ref
    integer :: e, a, b, ia, ib, i, k, p
    real(16) :: f
    real(16), allocatable :: W(:,:), Y(:,:), tw(:), ty(:)
    nq = nvar
    if (allocated(Aq)) deallocate(Aq, Bq, Xq)
    allocate (Aq(nq,nq), Bq(nq,nload), Xq(nq,nload))
    Aq = 0; Bq = 0
    do e = 1, nel
      do a = 1, kur(e)
        ia = idxof(nik(a,e))
        Bq(ia,:) = Bq(ia,:) + real(xr(a,1:nload,e), 16)
        do b = 1, kur(e)
          ib = idxof(nik(b,e))
          Aq(ia,ib) = Aq(ia,ib) + real(ck(a,b,e), 16)
        end do
      end do
    end do
    allocate (W(nq,nq), Y(nq,nload), tw(nq), ty(nload))
    W = Aq; Y = Bq
    do k = 1, nq-1
      p = k - 1 + maxloc(abs(W(k:nq,k)), 1)
      if (p /= k) then
        tw = W(k,:); W(k,:) = W(p,:); W(p,:) = tw
        ty = Y(k,:); Y(k,:) = Y(p,:); Y(p,:) = ty
      end if
      do i = k+1, nq
        f = W(i,k)/W(k,k)
        if (f /= 0) then
          W(i,k:nq) = W(i,k:nq) - f*W(k,k:nq)
          Y(i,:) = Y(i,:) - f*Y(k,:)
        end if
      end do
    end do
    do i = nq, 1, -1
      Xq(i,:) = (Y(i,:) - matmul(W(i,i+1:nq), Xq(i+1:nq,:))) / W(i,i)
    end do
  end subroutine dense_reference

  ! Accuracy table for every pass: residual, error vs quad reference, correction
  
  subroutine report(nref, tol, require_tol)
    integer, intent(in) :: nref
    real(8), intent(in) :: tol
    logical, intent(in) :: require_tol
    integer :: p, j
    real(16) :: xs(nq, nload), r(nq, nload)
    real(8) :: res, err, cor, err1, errl
    print '(5x,a,i0,a,i0,a,i0,a,es9.2,a,i0)', 'equations ', nq, ', load cases ', nload, &
          ', largest front ', maxpa_seen, ', largest CRIT ', critmax, ', CRIT warnings ', ncritwarn
    print '(5x,a)', 'pass  solution           ||b-Ax||/||b||   ||x-x_ref||/||x_ref||   ||dx||/||x||'
    xs = 0
    do p = 1, npass
      do j = 1, nload
        xs(:,j) = xs(:,j) + real(ans(1:nq,j,p), 16)
      end do
      r = Bq - matmul(Aq, xs)
      res = 0; err = 0
      do j = 1, nload
        res = max(res, real(maxval(abs(r(:,j)))/maxval(abs(Bq(:,j))), 8))
        err = max(err, real(maxval(abs(xs(:,j) - Xq(:,j)))/maxval(abs(Xq(:,j))), 8))
      end do
      if (p < npass) then
        cor = 0
        do j = 1, nload
          cor = max(cor, maxval(abs(ans(1:nq,j,p+1))) / real(maxval(abs(xs(:,j))), 8))
        end do
        print '(5x,i3,2x,a,es14.2,es20.2,es20.2)', p, solname(p), res, err, cor
      else
        print '(5x,i3,2x,a,es14.2,es20.2,a)', p, solname(p), res, err, '        -'
      end if
      if (p == 1) err1 = err
      errl = err
    end do
    if (require_tol) then
      call check(trim(tagline(nref)) // ' accuracy against quad-precision reference', errl < tol)
    else
      call check('refinement by re-solution reduces the error at least 1000-fold', &
                 errl < 1.d-3*err1 .and. errl < tol)
    end if
  end subroutine report

  character(len=18) function solname(p)
    integer, intent(in) :: p
    if (p == 1) then
      solname = 'frontal solve     '
    else
      write (solname, '(a,i1,a)') '+ refinement ', p-1, '    '
    end if
  end function solname

  character(len=24) function tagline(nref)
    integer, intent(in) :: nref
    if (nref == 0) then
      tagline = '   solution'
    else
      tagline = '   solution + refinement'
    end if
  end function tagline

  subroutine check(label, ok)
    character(len=*), intent(in) :: label
    logical, intent(in) :: ok
    print '(3x,a,1x,a)', merge('ok  ', 'FAIL', ok), label
    if (.not. ok) nfail = nfail + 1
  end subroutine check

end program test_irons_fe
