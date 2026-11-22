Irons (1970) frontal solver - Visual Studio 2022 solution
=========================================================

Requires: Visual Studio 2022 with the Intel Fortran Compiler (ifx) and its
Visual Studio integration (Intel oneAPI HPC Toolkit / Fortran Essentials).
Platform: x64 (ifx has no 32-bit target).

Open FrontalSolutionProgramIrons1970.sln.

It holds two console projects that share src\:

  test_irons_deck  zcapture.f90, zipp.f, elmak_card.f, test_irons_deck.f90
                   The paper's Appendix III specimen data, card-input ELMAK.
  test_irons_fe    zcapture.f90, zipp.f, test_irons_fe.f90
                   Finite element tests with a user ELMAK (in the .f90 file).

The two projects must not share elmak_card.f: test_irons_fe supplies its own
ELMAK.

Build: Build > Build Solution (Debug or Release, x64).
Run:   right-click a project > Set as Startup Project, then Ctrl+F5 (Start
       Without Debugging) so the console window stays open.

Output files (decks irons_*.dat, ZIPP listings irons_*.lis) are written to
the working directory - by default the project folder when started from
Visual Studio.

Debug configuration: bounds checking, stack frame checks and traceback on.
ZIPP's error exit (JWHERE traps, exercised by test D) divides by zero on
purpose, as Irons' artifice to get a trace; the default floating-point
setting (/fpe:3) does not trap on it.
