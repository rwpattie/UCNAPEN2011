C  Common blocks are defined in the included file 'pmcomms.f'.

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C         PPPPP   EEEEEE  N    N  M    M    AA    IIII  N    N         C
C         P    P  E       NN   N  MM  MM   A  A    II   NN   N         C
C         P    P  E       N N  N  M MM M  A    A   II   N N  N         C
C         PPPPP   EEEE    N  N N  M    M  AAAAAA   II   N  N N         C
C         P       E       N   NN  M    M  A    A   II   N   NN         C
C         P       EEEEEE  N    N  M    M  A    A  IIII  N    N         C
C                                                                      C
C                                                   (version 2011).    C
C                                                                      C
C  This program performs Monte Carlo simulation of electron-photon     C
C  showers in material structures described with the constructive      C
C  quadric geometry package 'PENGEOM.F'.                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2011)                                     C
C  Copyright (c) 2001-2011                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  In the default mode, the program assumes that primary particles of a
C  given kind, KPARP, are emitted from a point or an extended source,
C  either with fixed energy SE0 or with a specified (histogram-like)
C  energy spectrum. The initial direction of the primary particles is
C  sampled uniformly within a circle of the unit sphere (conical beam),
C  or within a 'rectangular' window (pyramidal beam). Alternatively, the
C  program can read the initial state variables of 'primary' particles
C  from a pre-calculated phase-space file. This option may be used to
C  split the simulation of complex problems into several consecutive
C  stages. The program also admits radiation sources defined by the
C  user, although this requires some editing work.
C
C  >>>>>>>> NOTE: All energies and lengths are given in eV and cm,
C                 respectively.
C
C
C  ************  Structure of the input data file.
C
C  Each line in the input data file consists of a 6-character keyword
C  (columns 1-6) followed either by numerical data (in free format) or
C  by a character string, which start at the 8th column. Keywords are
C  explicitly used/verified by the program (which is case sensitive!).
C  Notice also that the order of the data lines is important. The
C  keyword '______' (6 blanks) indicates comment lines, these can be
C  placed anywhere in the input file. The program ignores any text
C  following the first blank after the last numerical datum, or after
C  the character string, in each line (thus, in the table given below,
C  the comments in square brackets are ignored by the program). Lines
C  with some keywords (e.g., 'SPECTR', 'IPSFN') can appear an arbitrary
C  number of times, limited only by the allocated amount of memory.
C
C  The program assigns default values to many input variables; lines
C  that declare default values may be eliminated from the input file.
C
C
C  The structure of the input file is the following,
C
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C  TITLE  Title of the job, up to 65 characters.
C         . (the dot prevents editors from removing trailing blanks)
C         >>>>>>>> Source definition.
C  SKPAR  KPARP    [Primary particles: 1=electron, 2=photon, 3=positron]
C         KPARP=0 activates a user-defined SOURCE model.
C  SENERG SE0              [Initial energy (monoenergetic sources only)]
C  SPECTR Ei,Pi                 [E bin: lower-end and total probability]
C  SGPOL  SP1,SP2,SP3          [Stokes parameters for polarised photons]
C  SPOSIT SX0,SY0,SZ0                        [Coordinates of the source]
C  SBOX   SSX,SSY,SSZ                            [Source box dimensions]
C  SBODY  KB                [Active source body; one line for each body]
C  SCONE  THETA,PHI,ALPHA                  [Conical beam; angles in deg]
C  SPYRAM THETAL,THETAU,PHIL,PHIU      [Rectangular beam; angles in deg]
C         .
C         >>>>>>>> Input phase-space file (psf).
C  IPSFN  psf-filename.ext         [Input psf name, up to 20 characters]
C  IPSPLI NSPLIT                                      [Splitting number]
C  WGTWIN WGMIN,WGMAX         [Weight window, RR & spl of psf particles]
C  EPMAX  EPMAX                 [Maximum energy of particles in the psf]
C         .
C         >>>>>>>> Material data and simulation parameters.
C                  Up to MAXMAT materials; 2 lines for each material.
C  MFNAME mat-filename.ext               [Material file, up to 20 chars]
C  MSIMPA EABS(1:3),C1,C2,WCC,WCR              [EABS(1:3),C1,C2,WCC,WCR]
C         .
C         >>>>>>>> Geometry and local simulation parameters.
C  GEOMFN geo-filename.ext               [Geometry file, up to 20 chars]
C  DSMAX  KB,DSMAX(KB)              [KB, maximum step length in body KB]
C  EABSB  KB,EABSB(1:3,KB)   [KB, local absorption energies, EABSB(1:3)]
C         .
C         >>>>>>>> Interaction forcing.
C  IFORCE KB,KPAR,ICOL,FORCER,WLOW,WHIG  [KB,KPAR,ICOL,FORCER,WLOW,WHIG]
C         .
C         >>>>>>>> Emerging particles. Energy and angular distributions.
C  NBE    EL,EU,NBE                      [Energy window and no. of bins]
C  NBANGL NBTH,NBPH          [Nos. of bins for the angles THETA and PHI]
C         .
C         >>>>>>>> Impact detectors (up to 25 different detectors).
C         IPSF=0; no psf is created.
C         IPSF=1; a psf is created (for only one detector).
C         IDCUT=0; tracking is discontinued at the detector entrance.
C         IDCUT=1; the detector does not affect the tracking.
C         IDCUT=2; the detector does not affect tracking, the energy
C                  distribution of particle fluence (averaged over the
C                  volume of the detector) is calculated.
C  IMPDET EL,EU,NBE,IPSF,IDCUT      [E-window, no. of bins, IPSF, IDCUT]
C  IDSPC  spc-impdet-##.dat               [Spectrum file name, 20 chars]
C  IDPSF  psf-impdet-##.dat            [Phase-space file name, 20 chars]
C  IDFLNC fln-impdet-##.dat       [Fluence spectrum file name, 20 chars]
C  IDBODY KB                       [Active body; one line for each body]
C  IDKPAR KPAR               [Kind of detected particles, one line each]
C         .
C         >>>>>>>> Energy-deposition detectors (up to 25).
C  ENDETC EL,EU,NBE                   [Energy window and number of bins]
C  EDSPC  spc-enddet-##.dat        [Output spectrum file name, 20 chars]
C  EDBODY KB                       [Active body; one line for each body]
C         .
C         >>>>>>>> Dose distribution.
C  GRIDX  XL,XU                 [X coordinates of the dose box vertices]
C  GRIDY  YL,YU                 [Y coordinates of the dose box vertices]
C  GRIDZ  ZL,ZU                 [Z coordinates of the dose box vertices]
C  GRIDBN NDBX,NDBY,NDBZ                               [Numbers of bins]
C         .
C         >>>>>>>> Job properties.
C  RESUME dump1.dmp               [Resume from this dump file, 20 chars]
C  DUMPTO dump2.dmp                  [Generate this dump file, 20 chars]
C  DUMPP  DUMPP                                 [Dumping period, in sec]
C         .
C  RSEED  ISEED1,ISEED2           [Seeds of the random-number generator]
C  NSIMSH DSHN                     [Desired number of simulated showers]
C  TIME   TIMEA                       [Allotted simulation time, in sec]
C         .
C  END                                  [Ends the reading of input data]
C  ....+....1....+....2....+....3....+....4....+....5....+....6....+....
C
C
C  The following listing describes the function of each of the keywords,
C  the accompanying data and their default values. For clarity, blanks
C  in keywords are indicated as '_'.
C
C  TITLE_ : Title of the job (up to 65 characters).
C             DEFAULT: none (the input file must start with this line)
C
C           The TITLE string is used to mark dump files. To prevent the
C           improper use of wrong resuming files, change the title each
C           time you modify basic parameters of your problem. The code
C           will then be able to identify the inconsistency and to print
C           an error message before stopping.
C
C  >>>>>>>> Source definition.
C
C  SKPAR_ : Kind of primary particle KPARP (1=electrons, 2=photons or
C           3=positrons).
C             DEFAULT: KPARP=1
C           If KPARP=0, the initial states of primary particles are
C           set by subroutine SOURCE, to be provided by the user. An
C           example of that subroutine, corresponding to a 60-Co source
C           (two gamma rays in each nuclear deexcitation), is included
C           at the end of the present source file.
C
C  SENERG : For a monoenergetic source, initial energy SE0 of primary
C           particles.
C             DEFAULT: SE0=1.0E6
C
C  SPECTR : For a source with continuous (stepwise constant) spectrum,
C           each 'SPECTR' line gives the lower end-point of an energy
C           bin of the source spectrum (Ei) and the associated relative
C           probability (Pi), integrated over the bin. Up to NSEM=1000
C           lines, in arbitrary order. The upper end of the spectrum is
C           defined by entering a line with Ei equal to the upper energy
C           end point and with a negative Pi value.
C             DEFAULT: none
C
C  SGPOL_ : This line activates the simulation of polarisation effects
C           in the scattering of photons (electrons and positrons are
C           assumed to be unpolarised). SP1, SP2, SP3 are the Stokes
C           parameters of primary photons, which define the degrees of
C           linear polarisation at 45 deg azimuth, of circular
C           polarisation, and of linear polarisation at zero azimuth,
C           respectively. It is assumed that secondary photons are
C           emitted with null polarisation (SP1=SP2=SP3=0).
C             DEFAULT: none
C
C  SPOSIT : Coordinates of the source centre.
C             DEFAULT: SX0=SY0=0.0, SZ0=0.0D0
C
C  SBOX__ : Extended source box. The source has uniform activity within
C           the volume of a right prism centred at the point (SX0,SY0,
C           SZ0) and whose sides have lengths SSX, SSY and SSZ.
C             DEFAULT: SSX=SSY=SSZ=0.0
C
C  In the case of a extended source, the active volume can be restricted
C  to that of a body or a set of bodies, which must be defined as parts
C  of the geometry. The activity of the source is assumed to be uniform
C  within the volume of the intersection of the active bodies and the
C  source box. Note that the initial coordinates of primary particles
C  are sampled by the rejection method; the sampling efficiency is equal
C  to the fraction of the source box volume that is occupied by active
C  bodies.
C
C  To define each active source body, add the following line:
C
C  SBODY_ : Active source body. One line for each body.
C             DEFAULT: None
C           The program stops if the source box has not been defined
C           previously
C
C  The initial direction of primary particles is sampled unifomrly
C  within a circle on the unit sphere (conical beam), or within a
C  'rectangular' window (pyramidal beam).
C
C  SCONE_ : Conical source beam. Polar and azimuthal angles of the
C           beam axis direction, THETA and PHI, and angular aperture,
C           ALPHA, in deg.
C             DEFAULTS: THETA=0.0, PHI=0.0, ALPHA=0.0
C
C           The case ALPHA=0 defines a monodirectional source, and ALPHA
C           =180 deg corresponds to an isotropic source.
C
C  SPYRAM : Pyramidal source beam. Limiting polar and azimuthal angles
C           of the source beam window, (THETAL,THETAU)x(PHIL,PHIU), in
C           deg.
C             DEFAULTS: THETAL=0.0, THETAU=0.0, PHIL=0.0, PHIU=0.0
C
C           The case THETAL=THETAU, PHIL=PHIU defines a monodirectional
C           source. To define an isotropic source, set THETAL=0, THETAU=
C           180, PHIL=0 and PHIU=360.
C
C           --> Notice that the default source is a pencil beam that
C           moves upwards along the Z-axis.
C
C  >>>>>>>> Input phase-space file.
C
C  The initial state variables of primary particles can be read directly
C  from a set of pre-calculated phase-space files (psf). When this
C  option is active, previous definitions about the source are ignored.
C  Photons from the psf's are assumed to be unpolarized.
C
C  IPSFN_ : Name of an input psf (up to 20 characters).
C             DEFAULT: none
C           Up to 100 psf's may be declared. They are read sequentially.
C
C  The input psf is in ASCII format. Each line defines the initial state
C  of a particle; it contains the following quantities in free format
C  (and in the order they are listed here):
C    -- KPAR, type of particle (1, electron; 2, photon; 3, positron).
C    -- E, energy.
C    -- X,Y,Z, position coordinates.
C    -- U,V,W, direction cosines.
C    -- WGHT, weight.
C    -- ILB(1),ILB(2),ILB(3),ILB(4), a set of indices that provide
C           information on how the particle was generated (see the file
C           'manual.txt').
C    -- NSHI, incremental shower number (difference between the shower
C           numbers of the present particle and the one preceding it
C           in the psf).
C  Phase-space files can be generated by running PENMAIN using an impact
C  detector with the flag IPSF=1 (see below).
C
C  Because of the limited size of the psf's, the results of analogue
C  simulations tend to be 'too noisy'. This can be partially corrected
C  by splitting the particles from the psf.
C
C  IPSPLI : Splitting number. Each particle in the psf's will be split
C           into NSPLIT equivalent particles, with weights equal to
C           WGHT/NSPLIT.
C             DEFAULT: NSPLIT=1 (no splitting)
C
C  --> WARNING: Notice that there is a 'latent' uncertainty in the psf,
C  which sets a limit to the accuracy that can be attained by using
C  large splitting numbers.
C
C  WGTWIN : Weight window, (WGMIN,WGMAX). Particles in the phase-space
C           file that have initial weights WGHT less than WGMIN will be
C           subjected to Russian roulette, and those with WGHT larger
C           than WGMAX will be split. Note that the weight window has
C           preference over the splitting option, i.e., a particle will
C           be split into NSPLIT or less particles only if the latter
C           have weights larger than WGMIN.
C             DEFAULTS: WGMIN=1.0E-35, WGMAX=1.0E35  (no action)
C
C  EPMAX_ : Maximum energy (in eV) of particles in the psf's.
C           EPMAX is the upper limit of the energy interval covered by
C           the simulation lookup tables. To minimise interpolation
C           errors, EPMAX should not be much larger than the maximum
C           energy actually occurring during the simulation.
C
C           When the initial state variables of particles are read from
C           a psf, this parameter is required to initialise PENELOPE and
C           is critical; the code crashes if it finds a particle that
C           has energy larger than EPMAX.
C             DEFAULT: EPMAX=1.0E9 (interpolation is not optimal)
C
C  >>>>>>>> Material data and simulation parameters.
C
C  Each material is defined by introducing the following _two_ lines;
C
C  MFNAME : Name of a PENELOPE input material data file (up to 20
C           characters). This file must be generated in advance by
C           running the program MATERIAL.
C             DEFAULT: none
C
C  MSIMPA : Set of simulation parameters for this material; absorption
C           energies, EABS(1:3,M), elastic scattering parameters, C1(M)
C           and C2(M), and cutoff energy losses for inelastic collisions
C           and bremsstrahlung emission, WCC(M) and WCR(M).
C             DEFAULTS: EABS(1,M)=EABS(3,M)=0.01*EPMAX,
C                       EABS(2,M)=0.001*EPMAX
C                       C1(M)=C2(M)=0.1, WCC=EABS(1,M), WCR=EABS(2,M)
C             EPMAX is the upper limit of the energy interval covered
C             by the simulation lookup tables.
C
C  Note that we must declare a separate material data file name and a
C  set of simulation parameters for each material. The label (material
C  number) asigned by PENELOPE to each material is determined by the
C  ordering of the material data files in the PENMAIN input file. That
C  is, the first, second, ... materials are assigned the labels 1, 2,
C  ... These labels are also used in the geometry definition file.
C
C  The original programs in the distribution package allow up to 10
C  materials. This number can be increased by changing the value of the
C  parameter MAXMAT in the original source files.
C
C  >>>>>>>> Geometry and local simulation parameters.
C
C  GEOMFN : PENGEOM geometry definition file name (a string of up to
C           20 characters).
C             DEFAULT: none.
C
C           --> The geometry definition file can be debugged/visualised
C           with the viewers GVIEW2D and GVIEW3D (operable only under
C           Windows).
C
C  DSMAX_ : Maximum step length DSMAX(KB) of electrons and positrons in
C           body KB. This parameter is important only for thin bodies;
C           it should be given a value of the order of one tenth of the
C           body thickness or less.
C             DEFAULT: DSMAX=1.0E20 (no step length control)
C
C  EABSB_ : Local absorption energies EABSB(KPAR,KB) of particles of
C           type KPAR in body KB. These values must be larger than
C           EABS(KPAR,M), where M is the material of body KB. When the
C           particle is moving within body KB, the absorption energy
C           EABS(KPAR,M) is temporarily set equal to EABSB(KPAR,KB).
C           Thus, the simulation of the particle history is discontinued
C           when the energy becomes less than EABSB(KPAR,KB). This
C           feature can be used, e.g., to reduce the simulation work in
C           regions of lesser interest.
C             DEFAULTS: EABSB(KPAR,KB)=EABS(KPAR,M)  (no action)
C
C  >>>>>>>> Interaction forcing.
C
C  IFORCE : Activates forcing of interactions of type ICOL of particles
C           KPAR in body KB. FORCER is the forcing factor, which must
C           be larger than unity. WLOW and WHIG are the lower and upper
C           limits of the weight window where interaction forcing is
C           applied.
C             DEFAULT: no interaction forcing
C
C           If the mean free path for real interactions of type ICOL is
C           MFP, the program will simulate interactions of this type
C           (real or forced) with an effective mean free path equal to
C           MFP/FORCER.
C
C           TRICK: a negative input value of FORCER, -FN, is  assumed to
C           mean that each particle should interact, on average and
C           approximately, +FN times in a path length equal to the range
C           of that kind of particle with E=EPMAX. This is very useful
C           to generate x-ray spectra from bulk samples.
C
C  The real effect of interaction forcing on the efficiency is not easy
C  to predict. Please, do tentative runs with different FORCER values
C  and check the efficiency gain (or loss!).
C
C  >>>>>>>> Energy and angular distributions of emerging particles.
C
C  NBE___ : Limits EL and EU of the interval where energy distributions
C           of emerging particles are tallied. Number of energy bins
C           (.LE. 1000).
C             DEFAULT: EL=0.0, EU=EPMAX, NBE=500
C
C  NBANGL : Numbers of bins for the polar angle THETA and the azimuthal
C           angle PHI, respectively, NBTH and NBPH (.LE. 180)
C             DEFAULT: NBTH=90, NBPH=1 (azimuthal average)
C
C           NOTE: In the output files, the terms 'upbound' and
C           'downbound' are used to denote particles that leave the
C           material system moving upwards (W>0) and downwards (W<0),
C           respectively.
C
C  >>>>>>>> Impact detectors.
C
C  Each impact detector consists of a set of active bodies, which must
C  have been defined as parts of the geometry. The output spectrum is
C  the energy distribution of particles that entered any of the active
C  bodies coming from a body that is not active (i.e. that is not part
C  of the detector). Notice that a detected particle can re-enter the
C  detector volume and, consequently, be 'counted' several times (except
C  when the flag IDCUT is set equal to 0, see below).
C
C  Active bodies cannot be void, because the geometry routines would not
C  stop particles at their limiting surfaces. In case you need to define
C  detectors outside the material system, fill them with an arbitrary
C  material of very small density to avoid perturbing the transport
C  process.
C
C  To define each impact detector, insert the following block of lines;
C
C  IMPDET : Starts the definition of a new detector. Up to 25 different
C           detectors can be considered.
C           EL and EU are the lower and upper limits of the energy
C             window covered by the impact detector.
C           NBE is the number of bins in the output energy spectrum of
C             the detector (.LE. 1000). If NBE is positive, energy bins
C             have uniform width, DE=(EU-EL)/NBE. When NBE is negative,
C             the bin width increases geometrically with the energy,
C             i.e., the energy bins have uniform width in a logarithmic
C             scale.
C
C           The integer flag IPSF serves to activate the creation of a
C           phase-space file (psf), which contains the state variables
C           of all particles that enter the detector. Use this option
C           with care, because psf's may grow very fast.
C           IPSF=0; no psf is created.
C           IPSF=1; the psf is created. Only one PSF can be created in
C             each simulation run.
C
C           The integer flag IDCUT allows discontinuing the tracking of
C           particles that enter the detector.
C           IDCUT=0; the simulation of a particle is discontinued when
C             it enters the detector (useful to stop the simulation of
C             particles recorded in a psf).
C           IDCUT=1; the presence of the detector does not affect the
C             tracking of particles.
C           IDCUT=2; the presence of the detector does not affect the
C             tracking of particles. The distribution of particle
C             fluence with respect to energy (integrated over the volume
C             of the detector) is tallied. The calculated distribution
C             has dimensions of length/energy.
C
C             DEFAULTS: None
C
C  IDPSF_ : Name of the output phase-space file (up to 20 characters).
C             DEFAULT: 'psf-impdet-##.dat'
C
C  IDSPC_ : Name of the output energy spectrum file (20 characters).
C             DEFAULT: 'spc-impdet-##.dat'
C
C  IDFLNC : Name of the output file with the energy distribution of
C           particle fluence (20 characters). This file is generated
C           only when IDCUT=2.
C             DEFAULT: 'fln-impdet-##.dat'
C
C  IDBODY : Active body of the detector. One line for each active body.
C             DEFAULT: None
C           --> Notice that a body cannot be part of more than one
C           impact detector.
C
C  IDKPAR : Kind of particle that is detected (1=electrons, 2=photons or
C           3=positrons). One line for each kind.
C             DEFAULT: All particles are detected
C
C  >>>>>>>> Energy-deposition detectors.
C
C  Each energy-deposition detector consists of a set of active bodies,
C  which must have been defined as parts of the geometry. The output
C  spectrum is the distribution of absorbed energy (per primary shower)
C  in the active bodies.
C
C           *** WARNING: The energy-deposition spectrum may be strongly
C           biased when interaction forcing is applied within the
C           detector bodies.
C
C  To define each energy-deposition detector insert the following block
C  of lines;
C
C  ENDETC : Starts the definition of a new energy-deposition detector.
C           Up to 25 different detectors can be considered.
C           EL and EU are the lower and upper limits of the energy
C             window covered by the detector.
C           NBE is the number of bins in the output energy spectrum of
C             the detector (.LE. 1000). If NBE is positive, energy bins
C             have uniform width, DE=(EU-EL)/NBE. When NBE is negative,
C             the bin width increases geometrically with the energy,
C             i.e., the energy bins have uniform width in a logarithmic
C             scale.
C
C  EDSPC_ : Name of the output spectrum file (up to 20 characters).
C             DEFAULT: 'spc-enddet-##.dat'
C
C  EDBODY : Active body KB of the detector. One line for each active
C           body.
C             DEFAULT: None
C           --> Notice that a body cannot be part of more than one
C           energy-deposition detector.
C
C  >>>>>>>> Dose map.
C
C  Optionally, the program can calculate the dose distribution inside a
C  parallelepiped (dose box) whose edges are parallel to the axes of the
C  laboratory frame. The dose box is defined by giving the coordinates
C  of its vertices. The dose is tallied using a uniform orthogonal grid
C  with NDBX, NDBY and NDBZ bins (= voxels) along the directions of
C  the respective coordinate axes. These numbers should be odd, to make
C  sure that the 'central' axes (i.e., lines that join the centres of
C  oposite faces of the box) go through the centres of a row of voxels.
C
C  GRIDX_ : X-coordinates of the vertices of the dose box.
C             DEFAULT: None
C  GRIDY_ : Y-coordinates of the vertices of the dose box.
C             DEFAULT: None
C  GRIDZ_ : Z-coordinates of the vertices of the dose box.
C             DEFAULT: None
C  GRIDBN : Numbers of bins NDBX, NDBY, and NDBZ in the X, Y and Z
C             directions, respectively. These values shold be odd and
C             .LE. 101.
C             DEFAULTS: NDBX=25, NDBY=25, NDBZ=25
C
C  --> The grid defined here to calculate the dose distribution can be
C  used to tally other 3D distributions (e.g. the space distribution of
C  inner-shell ionisation events, used in electron-probe microanalysis).
C  This, however, requires to edit the present source file.
C
C  >>>>>>>> Job properties.
C
C  RESUME : The program will read the dump file named `dump1.dmp' (up to
C           20 characters) and resume the simulation from the point
C           where it was left. Use this option very, _VERY_ carefully.
C           Make sure that the input data file is fully consistent with
C           the one used to generate the dump file.
C             DEFAULT: off
C
C  DUMPTO : Generate a dump file named 'dump2.dmp' (up to 20 characters)
C           after completing the simulation run. This allows the
C           simulation to be resumed later on to improve statistics.
C             DEFAULT: off
C
C           NOTE: If the file 'dump2.dmp' already exists, it is
C           overwritten.
C
C  DUMPP_ : When the DUMPTO option is activated, simulation results are
C           written in the output files every DUMPP seconds. This option
C           is useful to check the progress of long simulations. It also
C           allows the program to be run with a long execution time and
C           to be stopped when the required statistical uncertainty has
C           been reached.
C             DEFAULT: DUMPP=1.0E15
C
C  RSEED_ : Seeds of the random-number generator.
C             DEFAULT: ISEED1=1, ISEED2=1
C
C  NSIMSH : Desired number of simulated showers.
C             DEFAULT: DSHN=2.0E9
C
C  TIME__ : Allotted simulation time, in sec.
C             DEFAULT: TIMEA=2.0E9
C
C  END___ : Ends the reading of the input file. This line is needed only
C           when the TIME__ line is missing.
C
C  The program is aborted when an incorrect input datum is found. The
C  conflicting quantity usually appears in the last line of the output
C  file 'penmain.dat'. If the trouble is with arrays having dimensions
C  smaller than required, the program indicates how the problem can be
C  solved (this usually requires editing the source file, be careful).
C
C  The clock subroutine (TIMER) may have to be adapted to your specific
C  computer-compiler configuration; standard Fortran 77 does not provide
C  timing tools. However, the routines in module TIMER.F do work for
C  many Fortran compilers.
C
C  ************  Generating the executable PENMAIN and running it.
C
C  To generate the executable binary file PENMAIN.EXE, compile and link
C  the Fortran source files PENMAIN.F, PENELOPE.F, PENGEOM.F, PENVARED.F
C  and TIMER.F. For example, if you are using the GNU Fortran compiler
C  under Windows, place these five files in the same directory, open a
C  command window and from that directory enter the command
C    `gfortran -Wall -Os PENMAIN.F -o PENMAIN.EXE'
C  (The same, with file names in lowercase, should work under Linux).
C
C  To run PENMAIN, you have to generate an input data file, let's call
C  it PENMAIN.IN, and the corresponding geometry definition and material
C  data files. Place these files and the binary file PEMAIN.EXE in the
C  same directory and, from there, issue the command
C    `PENMAIN.EXE < PENMAIN.IN'
C
C  The calculated distributions are written in separate files, whose
C  names have the extension '.dat'. These files are in a format suited
C  for direct visualisation with GNUPLOT (version 4.2).
C
C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      INCLUDE 'pmcomms.f'
      INCLUDE 'ucnapenmain.h'
C
      CALL HLIMIT(NWPAWC) ! A FEW PAW FORUMS SUGGEST THIS FOR MAKING
      IQUEST(10) = 65000  ! BIGGER NTUPLES. KINDA WORKS BUT STILL 
                          ! LIMITED.
       CALL HROPEN(1,'EVENT','test.hbook','NQE',1024,IERR)
! c      if(NTYPE.NE.10)THEN
       CALL HBOOKN(34,'DECAYS',nentries,'//EVENT',1024,DTAGS) 
c      ELSE 
c        CALL HBOOKN(34,'DECAYS',ngentries,'//EVENT',1024,GETAGS)
c      ENDIF
C      ELSE IF(DTYPE.EQ.2)THEN
c        CALL HBOOKN(34,'DECAYS',nCDentries,'//EVENT',1024,CDTAGS)
c      ENDIF
      CALL HBOOK1(100,'INTIAL BETA DECAY SPECTRUM',300,0.,3000.,0.)
      CALL HBOOK1(110,'THROWN ASYMMETRY',100,0.,1000.,0.)
      CALL HBOOK1(120,'MEASURED ASYMMETRY (RAW)',100,0.,1000.,0.)
      CALL HBOOK1(130,'RAW MWPC SPECTRUM - EAST',100,0.,1000.,0.)
      CALL HBOOK1(140,'RAW MWCP SPECTRUM - WEST',100,0.,1000.,0.)
      CALL HBOOK1(150,'RAW PMT SPECTRUM - EAST',100,0.,1000.,0.)
      CALL HBOOK1(160,'RAW PMT SPECTRUM - WEST',100,0.,1000.,0.)
      CALL HBOOK1(170,'RAW MYLAR SPECTRUM - EAST',100,0.,1000.,0.)
      CALL HBOOK1(180,'RAW MYLAR SPECTRUM - WEST',100,0.,1000.,0.)
      CALL HBOOK1(190,'SPECTRUM HITS',100,0.,1000.,0.)
      CALL HBOOK1(195,'PHOTON ELECTRONS',36,0.,360.,0.)
      CALL HBOOK1(200,'BACKSCATTER ENERGY SPECTRUM',100,0.,1000.,0.)
      CALL HBOOK2(210,'BACKSCATTER ENERGY VS ANGLE SPECTRUM',100,0.,
     +                1000.,100,0.,1.,0.)
      CALL HBOOK1(220,'ANGLUAR DISTRIBUTION OF BACKSCATTERS',100,0.,
     +                1.,0.)

C
C  ****  Read input files and initialize the simulation packages.
C
      CALL PMRDR
      IF(JOBEND.NE.0) GO TO 103 ! The simulation was already completed.
C
C  ****  Simulation of a new shower and scoring.
C
  101 CONTINUE
      CALL SHOWER
      IF(JOBEND.NE.0) GO TO 102  ! The simulation is completed.
C
C  ****  End the simulation after the allotted time or after completing
C        DSHN showers.
C
      CALL TIMER(TSEC)
      IF(TSEC.LT.TSECA.AND.SHN.LT.DSHN) THEN
C  ****  Write partial results after each dumping period.
        IF(LDUMP) THEN
          IF(TSEC-TSECAD.GT.DUMPP) THEN
            TSECAD=TSEC
            TSIM=TSIM+CPUTIM()-CPUT0
            CALL PMWRT(-1)
            WRITE(6,1001) SHN
 1001       FORMAT(3X,'Number of simulated showers =',1P,E14.7)
            CPUT0=CPUTIM()
            GO TO 101
          ENDIF
        ENDIF
        GO TO 101
      ENDIF
C
  102 CONTINUE
      TSIM=TSIM+CPUTIM()-CPUT0
  103 CONTINUE
      CALL PMWRT(1)
      WRITE(6,1002) SHN
 1002 FORMAT(3X,'Number of simulated showers =',1P,E14.7,
     1  /3X,'*** END ***')
      CALL HROUT(0,ICYCLE,' ')
      CALL HPRINT(34)
      CALL HPRINT(150)
      CALL HPRINT(195)
      CALL HREND('EVENT')
     
      STOP
      END

C  *********************************************************************
C                       SUBROUTINE SHOWER
C  *********************************************************************
      SUBROUTINE SHOWER
C
C  Simulates a new shower and keeps score of relevant quantities.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      LOGICAL LINTF
      INCLUDE 'pmcomms.f'
      EXTERNAL RAND
C
C  ------------------------  Shower simulation starts here.
C
C  ************  Primary particle counters.
C
      
  101 CONTINUE
c     write(6,*) ISEED1,ISEED2
      DO I=1,3
        DPRIM(I)=0.0D0
        DO K=1,3
          DSEC(K,I)=0.0D0
        ENDDO
      ENDDO
C
      DO KB=1,NBODY
        DEBO(KB)=0.0D0  ! Energies deposited in the various bodies KB.
      ENDDO
C
      IF(NDIDEF.GT.0) THEN
        DO KD=1,NDIDEF
          DEDI(KD)=0.0D0  ! Energy collected by each detector.
        ENDDO
      ENDIF
      IEXIT=0
C
      CALL CLEANS          ! Cleans the secondary stack.
C
C  **********  Set the initial state of the primary particle.
C
  201 CONTINUE
      IF(KPARP.EQ.0) THEN
C  ****  User-defined source.
        CALL SOURCE
        SHN=SHN+1.0D0
        N=N+1
        IF(N.GT.2000000000) N=N-2000000000
        IF(LSPEC) THEN
          KE=E*RDSHE+1.0D0
          SHIST(KE)=SHIST(KE)+1.0D0
        ENDIF
c     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,0p,i3)')
c    1    MOD(N,100),KPAR,ILB(1),X,Y,Z,W,E,IBODY
      ELSE IF(LPSF) THEN
C  ****  Phase-space file.
        CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
        IF(KODEPS.EQ.-1) THEN
c       write(35,'(i2,1p,8e13.5,i3,i2,i2,i9,i9,2i3)')
c    1    KPAR,E,X,Y,Z,U,V,W,WGHT,ILB(1),ILB(2),ILB(3),ILB(4),
c    2    NSHI,ISEC,KODEPS
          CLOSE(IPSFI)
          NPSN=NPSN+1
          RLREAD=0.0D0
          IF(NPSN.GT.NPSF) THEN
            JOBEND=4
            RETURN
          ENDIF
          OPEN(IPSFI,FILE=PSFI(NPSN))
          WRITE(6,3040) SHN
          WRITE(6,1904) PSFI(NPSN)
          WRITE(26,1904) PSFI(NPSN)
          KODEPS=0
          CALL RDPSF(IPSFI,NSHI,ISEC,KODEPS)
c       write(35,'(i2,1p,8e13.5,i3,i2,i2,i9,i9,2i3)')
c    1    KPAR,E,X,Y,Z,U,V,W,WGHT,ILB(1),ILB(2),ILB(3),ILB(4),
c    2    NSHI,ISEC,KODEPS
          GO TO 101
        ENDIF
 3040   FORMAT(2X,'Number of simulated showers =',1P,E14.7)
 1904   FORMAT(/3X,'+ The phase-space file ',A20,' is opened.')
C
        RLREAD=RLREAD+1.0D0
        SHN=SHN+DBLE(NSHI)
        N=N+NSHI
        IF(N.GT.2000000000) N=N-2000000000
        ILB(5)=0
        IPOL=0  ! Particles in the phase-space file are unpolarized.
      ELSE
C  ****  External source.
        SHN=SHN+1.0D0
        N=N+1
        IF(N.GT.2000000000) N=N-2000000000
        KPAR=KPARP
        WGHT=1.0D0
C  ----  Initial position ...
        IF(LEXSRC) THEN
          IF(LEXBD) THEN
            NTRIAL=0
  301       CONTINUE
            X=SX0+(RAND(1.0D0)-0.5D0)*SSX
            Y=SY0+(RAND(2.0D0)-0.5D0)*SSY
            Z=SZ0+(RAND(3.0D0)-0.5D0)*SSZ
            CALL LOCATE
            NTRIAL=NTRIAL+1
            IF(NTRIAL.GT.200) THEN
              WRITE(26,'(3X,''WARNING: the sampling of initial '',
     1          ''positions may be very inefficient.'')')
              WRITE(6,'(3X,''WARNING: the sampling of initial '',
     1          ''positions may be very inefficient.'')')
            ENDIF
            IF(IXSBOD(IBODY).EQ.0) GO TO 301
          ELSE
            X=SX0+(RAND(1.0D0)-0.5D0)*SSX
            Y=SY0+(RAND(2.0D0)-0.5D0)*SSY
            Z=SZ0+(RAND(3.0D0)-0.5D0)*SSZ
          ENDIF
        ELSE
          X=SX0
          Y=SY0
          Z=SZ0
        ENDIF
C  ----  Initial direction ...
        IF(LSCONE) THEN
          CALL GCONE(U,V,W)  ! Conical beam.
        ELSE
          W=CTHL+RAND(4.0D0)*DCTH  ! Rectangular beam.
          UV=SQRT(1.0D0-W*W)
          PHI=PHIL+RAND(5.0D0)*DPHI
          U=UV*COS(PHI)
          V=UV*SIN(PHI)
        ENDIF
C  ----  Initial energy ...
        IF(LSPEC) THEN
          RN=RAND(6.0D0)*NSEB+1
          K=INT(RN) ! Continuous spectrum. E sampled by Walker's method.
          RNF=RN-K
          IF(RNF.GT.FS(K)) THEN
            KE=IAS(K)
          ELSE
            KE=K
          ENDIF
          E=ES(KE)+RAND(7.0D0)*(ES(KE+1)-ES(KE))
          SHIST(KE)=SHIST(KE)+1.0D0
        ELSE
          E=E0      ! Monoenergetic source.
          SHIST(1)=SHIST(1)+1.0D0
        ENDIF
        ILB(1)=1  ! Identifies primary particles.
        ILB(2)=0
        ILB(3)=0
        ILB(4)=0
        ILB(5)=0
        IF(KPAR.EQ.2) THEN
          IF(LGPOL) THEN
            IPOL=1  ! Photon polarisation.
            SP1=SP10
            SP2=SP20
            SP3=SP30
          ELSE
            IPOL=0
          ENDIF
        ELSE
          IPOL=0
        ENDIF
      ENDIF
      CALL START           ! Starts simulation in current medium.
C
C  ****  Check if the trajectory intersects the material system.
C
  302 CONTINUE
      CALL LOCATE
C     
      CALL TPEMF0(ULDV,ULDE,ULEM,DSMAX)
      IF(MAT.EQ.0) THEN
        IBODYL=IBODY
        CALL JUMP(DSMAX,DS)
        CALL TPEMF1(DS,DSEF,NCROSS)
        !CALL STEP(1.0D30,DSEF,NCROSS)
        IF(MAT.EQ.0) THEN  ! The particle does not enter the system.
          IF(W.GT.0) THEN
            IEXIT=1        ! Labels emerging upbound particles.
          ELSE
            IEXIT=2        ! Labels emerging downbound particles.
          ENDIF
          GO TO 104        ! Exit.
        ENDIF

C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Impact detectors.
        IDET=KDET(IBODY)
        IF(IDET.NE.0) THEN
          IF(KDET(IBODYL).NE.IDET.AND.KKDI(IDET,KPAR).EQ.1) THEN
C
            IF(IPSF(IDET).EQ.1) THEN
              NSHJ=SHN-RLAST
              CALL WRPSF(IPSFO,NSHJ,0)
              RWRITE=RWRITE+1.0D0
              RLAST=SHN
            ENDIF
C
            DEDI(IDET)=DEDI(IDET)+E*WGHT
            IF(LDILOG(IDET)) THEN
              IE=1.0D0+(LOG(E)-EDILL(IDET))*RBDIEL(IDET)
            ELSE
              IE=1.0D0+(E-EDIL(IDET))*RBDIE(IDET)
            ENDIF
            IF(IE.GT.0.AND.IE.LE.NDICH(IDET)) THEN
              IF(N.NE.LDIT(IDET,IE)) THEN
                DIT(IDET,IE)=DIT(IDET,IE)+DITP(IDET,IE)
                DIT2(IDET,IE)=DIT2(IDET,IE)+DITP(IDET,IE)**2
                DITP(IDET,IE)=WGHT
                LDIT(IDET,IE)=N
              ELSE
                DITP(IDET,IE)=DITP(IDET,IE)+WGHT
              ENDIF
              IF(N.NE.LDIP(IDET,IE,KPAR)) THEN
                DIP(IDET,IE,KPAR)=DIP(IDET,IE,KPAR)+DIPP(IDET,IE,KPAR)
                DIP2(IDET,IE,KPAR)=
     1              DIP2(IDET,IE,KPAR)+DIPP(IDET,IE,KPAR)**2
                DIPP(IDET,IE,KPAR)=WGHT
                LDIP(IDET,IE,KPAR)=N
              ELSE
                DIPP(IDET,IE,KPAR)=DIPP(IDET,IE,KPAR)+WGHT
              ENDIF
            ENDIF
            IF(IDCUT(IDET).EQ.0) THEN
              DEBO(IBODY)=DEBO(IBODY)+E*WGHT
              IEXIT=3
              GO TO 104
            ENDIF
          ENDIF
        ENDIF
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      ENDIF
C
      IF(E.LT.EABSB(KPAR,IBODY)) THEN  ! The energy is too low.
        DEP=E*WGHT
        DEBO(IBODY)=DEBO(IBODY)+DEP
        IF(LDOSEM) THEN  ! Particle inside the dose box.
          IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
     1       (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
     1       (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
            I1=1.0D0+(X-DXL(1))*RBDOSE(1)
            I2=1.0D0+(Y-DXL(2))*RBDOSE(2)
            I3=1.0D0+(Z-DXL(3))*RBDOSE(3)
            DOSP=DEP*RHOI(MAT)
            IF(N.NE.LDOSE(I1,I2,I3)) THEN
              DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
              DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
              DOSEP(I1,I2,I3)=DOSP
              LDOSE(I1,I2,I3)=N
            ELSE
              DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DOSP
            ENDIF
C
            IF(N.NE.LDDOSE(I3)) THEN
              DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
              DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
              DDOSEP(I3)=DOSP
              LDDOSE(I3)=N
            ELSE
              DDOSEP(I3)=DDOSEP(I3)+DOSP
            ENDIF
          ENDIF
        ENDIF
        IEXIT=3
        GO TO 104
      ENDIF
C  ---------------------------------------------------------------------
C  ------------------------  Track simulation begins here.
C

C  >>>>>>>>>>>>  Particle splitting and Russian roulette  >>>>>>>>>>>>>>
C  ...  only for particles read from phase-space files.
      IF(LPSF) THEN
C  ----  Russian roulette (weight window).
        IF(WGHT.LT.WGMIN) THEN
          PKILL=1.0D0-WGHT*RWGMIN
          CALL VKILL(PKILL)
          IF(E.LT.1.0D0) GO TO 201
        ENDIF
C  ----  Splitting.
        NSPL1=MAX(MIN(NSPLIT,INT(MIN(1.0D5,WGHT*RWGMIN))),1)
        IF(NSPL1.GT.1) THEN
          CALL VSPLIT(NSPL1)
C  ----  Energy is locally deposited in the material.
          DEP=(NSPL1-1)*E*WGHT
          DEBO(IBODY)=DEBO(IBODY)+DEP
          IF(LDOSEM) THEN  ! Particle inside the dose box.
            IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
     1         (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
     1         (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
              I1=1.0D0+(X-DXL(1))*RBDOSE(1)
              I2=1.0D0+(Y-DXL(2))*RBDOSE(2)
              I3=1.0D0+(Z-DXL(3))*RBDOSE(3)
              DOSP=DEP*RHOI(MAT)
              IF(N.NE.LDOSE(I1,I2,I3)) THEN
                DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
                DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
                DOSEP(I1,I2,I3)=DOSP
                LDOSE(I1,I2,I3)=N
              ELSE
                DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DOSP
              ENDIF
C
              IF(N.NE.LDDOSE(I3)) THEN
                DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
                DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
                DDOSEP(I3)=DOSP
                LDDOSE(I3)=N
              ELSE
                DDOSEP(I3)=DDOSEP(I3)+DOSP
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C  <<<<<<<<<<<<  Particle splitting and Russian roulette  <<<<<<<<<<<<<<

  102 CONTINUE
      EABS(KPAR,MAT)=EABSB(KPAR,IBODY)
C
C  ************  The particle energy is less than EABS.
C
      IF(E.LT.EABS(KPAR,MAT)) THEN  ! The particle is absorbed.
        DEP=E*WGHT
C  ----  Energy is locally deposited in the material.
        DEBO(IBODY)=DEBO(IBODY)+DEP
        IF(LDOSEM) THEN  ! Particle inside the dose box.
          IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
     1       (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
     1       (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
            I1=1.0D0+(X-DXL(1))*RBDOSE(1)
            I2=1.0D0+(Y-DXL(2))*RBDOSE(2)
            I3=1.0D0+(Z-DXL(3))*RBDOSE(3)
            DOSP=DEP*RHOI(MAT)
            IF(N.NE.LDOSE(I1,I2,I3)) THEN
              DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
              DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
              DOSEP(I1,I2,I3)=DOSP
              LDOSE(I1,I2,I3)=N
            ELSE
              DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DOSP
            ENDIF
C
            IF(N.NE.LDDOSE(I3)) THEN
              DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
              DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
              DDOSEP(I3)=DOSP
              LDDOSE(I3)=N
            ELSE
              DDOSEP(I3)=DDOSEP(I3)+DOSP
            ENDIF
          ENDIF
        ENDIF
        IEXIT=3                     ! Labels absorbed particles.
        GO TO 104                   ! Exit.
      ENDIF
C
C
  103 CONTINUE
      IBODYL=IBODY
c     write(6,'(''n,kpar,gen,x,y,z,w,e,ibody='',3i3,1p,5e11.3,0p,i3)')
c    1    MOD(N,100),KPAR,ILB(1),X,Y,Z,W,E,IBODY
C
      CALL TPEMF0(ULDV,ULDE,ULEM,DSMAX)
      IF(LFORCE(IBODY,KPAR).AND.((WGHT.GE.WLOW(IBODY,KPAR)).AND.
     1  (WGHT.LE.WHIG(IBODY,KPAR)))) THEN
        CALL JUMPF(DSMAX(IBODY),DS)  ! Interaction forcing.
        LINTF=.TRUE.
      ELSE
        CALL JUMP(DSMAX(IBODY),DS)   ! Analogue simulation.
        LINTF=.FALSE.
      ENDIF
      CALL TPEMF1(DS,DSEF,NCROSS)
      !CALL STEP(DS,DSEF,NCROSS)      ! Determines step end position.

C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Energy distribution of fluence.
      IDET=KDET(IBODYL)
      IF(IDET.NE.0) THEN
        IF(IDCUT(IDET).EQ.2) THEN
          IF(KKDI(IDET,KPAR).EQ.1) THEN
            IF(LDILOG(IDET)) THEN
              IE=1.0D0+(LOG(E)-EDILL(IDET))*RBDIEL(IDET)
            ELSE
              IE=1.0D0+(E-EDIL(IDET))*RBDIE(IDET)
            ENDIF
            IF(IE.GT.0.AND.IE.LE.NDICH(IDET)) THEN
              IF(N.NE.LFST(IDET,IE)) THEN
                FST(IDET,IE)=FST(IDET,IE)+FSTP(IDET,IE)
                FST2(IDET,IE)=FST2(IDET,IE)+FSTP(IDET,IE)**2
                FSTP(IDET,IE)=WGHT*DSEF
                LFST(IDET,IE)=N
              ELSE
                FSTP(IDET,IE)=FSTP(IDET,IE)+WGHT*DSEF
              ENDIF
              IF(N.NE.LFSP(IDET,IE,KPAR)) THEN
                FSP(IDET,IE,KPAR)=FSP(IDET,IE,KPAR)+FSPP(IDET,IE,KPAR)
                FSP2(IDET,IE,KPAR)=
     1              FSP2(IDET,IE,KPAR)+FSPP(IDET,IE,KPAR)**2
                FSPP(IDET,IE,KPAR)=WGHT*DSEF
                LFSP(IDET,IE,KPAR)=N
              ELSE
                FSPP(IDET,IE,KPAR)=FSPP(IDET,IE,KPAR)+WGHT*DSEF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C  ----  Check whether the particle is outside the enclosure.
      IF(MAT.EQ.0) THEN
        IF(W.GT.0) THEN
          IEXIT=1        ! Labels emerging upbound particles.
        ELSE
          IEXIT=2        ! Labels emerging downbound particles.
        ENDIF
        GO TO 104        ! Exit.
      ENDIF

C  ----  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ----  Impact detectors.
      IDET=KDET(IBODY)
      IF(IDET.NE.0) THEN
        IF(KDET(IBODYL).NE.IDET.AND.KKDI(IDET,KPAR).EQ.1) THEN
C
          IF(IPSF(IDET).EQ.1) THEN
            NSHJ=SHN-RLAST
            CALL WRPSF(IPSFO,NSHJ,0)
            RWRITE=RWRITE+1.0D0
            RLAST=SHN
          ENDIF
C
          DEDI(IDET)=DEDI(IDET)+E*WGHT
          IF(LDILOG(IDET)) THEN
            IE=1.0D0+(LOG(E)-EDILL(IDET))*RBDIEL(IDET)
          ELSE
            IE=1.0D0+(E-EDIL(IDET))*RBDIE(IDET)
          ENDIF
          IF(IE.GT.0.AND.IE.LE.NDICH(IDET)) THEN
            IF(N.NE.LDIT(IDET,IE)) THEN
              DIT(IDET,IE)=DIT(IDET,IE)+DITP(IDET,IE)
              DIT2(IDET,IE)=DIT2(IDET,IE)+DITP(IDET,IE)**2
              DITP(IDET,IE)=WGHT
              LDIT(IDET,IE)=N
            ELSE
              DITP(IDET,IE)=DITP(IDET,IE)+WGHT
            ENDIF
            IF(N.NE.LDIP(IDET,IE,KPAR)) THEN
              DIP(IDET,IE,KPAR)=DIP(IDET,IE,KPAR)+DIPP(IDET,IE,KPAR)
              DIP2(IDET,IE,KPAR)=
     1            DIP2(IDET,IE,KPAR)+DIPP(IDET,IE,KPAR)**2
              DIPP(IDET,IE,KPAR)=WGHT
              LDIP(IDET,IE,KPAR)=N
            ELSE
              DIPP(IDET,IE,KPAR)=DIPP(IDET,IE,KPAR)+WGHT
            ENDIF
          ENDIF
          IF(IDCUT(IDET).EQ.0) THEN
            DEBO(IBODY)=DEBO(IBODY)+E*WGHT
            IEXIT=3
            GO TO 104
          ENDIF
        ENDIF
      ENDIF
C  ----  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

C  ----  If the particle has crossed an interface, restart the track in
C        the new material.
      IF(NCROSS.GT.0) GO TO 102    ! The particle crossed an interface.
C  ----  Simulate next event.
      IF(LINTF) THEN
        CALL KNOCKF(DE,ICOL)       ! Interaction forcing is active.
      ELSE
        CALL KNOCK(DE,ICOL)        ! Analogue simulation.
      ENDIF
      DEP=DE*WGHT
C  ----  Energy is locally deposited in the material.
      DEBO(IBODY)=DEBO(IBODY)+DEP
!       IF(LDOSEM) THEN  ! The particle is inside the dose box.
!         IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
!      1     (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
!      1     (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
!           I1=1.0D0+(X-DXL(1))*RBDOSE(1)
!           I2=1.0D0+(Y-DXL(2))*RBDOSE(2)
!           I3=1.0D0+(Z-DXL(3))*RBDOSE(3)
!           DOSP=DEP*RHOI(MAT)
!           IF(N.NE.LDOSE(I1,I2,I3)) THEN
!             DOSE(I1,I2,I3)=DOSE(I1,I2,I3)+DOSEP(I1,I2,I3)
!             DOSE2(I1,I2,I3)=DOSE2(I1,I2,I3)+DOSEP(I1,I2,I3)**2
!             DOSEP(I1,I2,I3)=DOSP
!             LDOSE(I1,I2,I3)=N
!           ELSE
!             DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)+DOSP
!           ENDIF
! C
!           IF(N.NE.LDDOSE(I3)) THEN
!             DDOSE(I3)=DDOSE(I3)+DDOSEP(I3)
!             DDOSE2(I3)=DDOSE2(I3)+DDOSEP(I3)**2
!             DDOSEP(I3)=DOSP
!             LDDOSE(I3)=N
!           ELSE
!             DDOSEP(I3)=DDOSEP(I3)+DOSP
!           ENDIF
!         ENDIF
!       ENDIF
C
      IF(E.LT.EABS(KPAR,MAT)) THEN  ! The particle has been absorbed.
        IEXIT=3                     ! Labels absorbed particles.
        GO TO 104                   ! Exit.
      ENDIF
C
      GO TO 103
C  ------------------------  The simulation of the track ends here.
C  ---------------------------------------------------------------------
  104 CONTINUE
C
C  ************  Increment particle counters.
C
      IF(ILB(1).EQ.1) THEN
        DPRIM(IEXIT)=DPRIM(IEXIT)+WGHT
      ELSE
        DSEC(KPAR,IEXIT)=DSEC(KPAR,IEXIT)+WGHT
      ENDIF
C
C      IF(IEXIT.LT.3) THEN
C  ****  Energy distribution of emerging particles.
!         K=1.0D0+(E-EMIN)*RBSE
!         IF(K.GT.0.AND.K.LE.NBE) THEN
!           IF(N.NE.LPDE(KPAR,IEXIT,K)) THEN
!             PDE(KPAR,IEXIT,K)=PDE(KPAR,IEXIT,K)+PDEP(KPAR,IEXIT,K)
!             PDE2(KPAR,IEXIT,K)=
!      1        PDE2(KPAR,IEXIT,K)+PDEP(KPAR,IEXIT,K)**2
!             PDEP(KPAR,IEXIT,K)=WGHT
!             LPDE(KPAR,IEXIT,K)=N
!           ELSE
!             PDEP(KPAR,IEXIT,K)=PDEP(KPAR,IEXIT,K)+WGHT
!           ENDIF
!         ENDIF
! C  ****  Angular distribution of emerging particles.
!         THETA=ACOS(W)
!         KTH=1.0D0+THETA*RA2DE*RBSTH
!         IF(ABS(U).GT.1.0D-16) THEN  ! Azimuthal bin number corrected.
!            PHI=ATAN2(V,U)
!         ELSE IF(ABS(V).GT.1.0D-16) THEN
!            PHI=ATAN2(V,U)
!         ELSE
!            PHI=0.0D0
!         ENDIF
!         IF(PHI.LT.0.0D0) PHI=TWOPI+PHI
!         KPH=1.0D0+PHI*RA2DE*RBSPH
!         IF(N.NE.LPDA(KPAR,KTH,KPH)) THEN
!           PDA(KPAR,KTH,KPH)=PDA(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)
!           PDA2(KPAR,KTH,KPH)=PDA2(KPAR,KTH,KPH)+PDAP(KPAR,KTH,KPH)**2
!           PDAP(KPAR,KTH,KPH)=WGHT
!           LPDA(KPAR,KTH,KPH)=N
!         ELSE
!           PDAP(KPAR,KTH,KPH)=PDAP(KPAR,KTH,KPH)+WGHT
!         ENDIF
C      ENDIF
C
C  ************  Any secondary left?
C
  202 CONTINUE
      CALL SECPAR(LEFT)
      IF(LEFT.GT.0) THEN
c     write(6,'(/''new secondary'')')
c     write(6,'(''n,kpar,ilb(1:4),ibody,mat='',i6,4i4,i10,2i4)')
c    1    N,KPAR,ILB(1),ILB(2),ILB(3),ILB(4),IBODY,MAT
c     write(6,'(''  ener,x,y,z='',1p,5e14.6)') E,X,Y,Z
c     write(6,'(''  wght,u,v,w='',1p,5e14.6)') WGHT,U,V,W
        IF(ILB(1).EQ.-1) THEN  ! Primary particle from SOURCE.
          ILB(1)=1  ! Energy is not removed from the site.
          IF(LSPEC) THEN
            KE=E*RDSHE+1.0D0
            SHIST(KE)=SHIST(KE)+1.0D0
          ENDIF
          GO TO 302
        ENDIF
        IF(E.GT.EABSB(KPAR,IBODY)) THEN
          DEBO(IBODY)=DEBO(IBODY)-E*WGHT  ! Energy is removed.
!           IF(LDOSEM) THEN  ! Particle inside the dose box.
!             IF((X.GT.DXL(1).AND.X.LT.DXU(1)).AND.
!      1         (Y.GT.DXL(2).AND.Y.LT.DXU(2)).AND.
!      1         (Z.GT.DXL(3).AND.Z.LT.DXU(3))) THEN
!               I1=1.0D0+(X-DXL(1))*RBDOSE(1)
!               I2=1.0D0+(Y-DXL(2))*RBDOSE(2)
!               I3=1.0D0+(Z-DXL(3))*RBDOSE(3)
!               DOSP=E*WGHT*RHOI(MAT)
!               DOSEP(I1,I2,I3)=DOSEP(I1,I2,I3)-DOSP
!               DDOSEP(I3)=DDOSEP(I3)-DOSP
!             ENDIF
!           ENDIF
        ELSE
          GO TO 202
        ENDIF

C  >>>>>>>>>>>>>>>>>>>>>>>  Russian roulette  >>>>>>>>>>>>>>>>>>>>>>>>>>
C  ****  Russian roulette for photons moving downstream.
cr      IF(KPAR.EQ.2.AND.W.LT.0.0D0) THEN
cr        IF(WGHT.LT.0.1D0) THEN
cr          PKILL=0.75D0
cr          CALL VKILL(PKILL)
cr          IF(E.LT.EABSB(KPAR,IBODY)) GO TO 202
cr        ENDIF
cr      ENDIF
C  <<<<<<<<<<<<<<<<<<<<<<<  Russian roulette  <<<<<<<<<<<<<<<<<<<<<<<<<<

        GO TO 102
      ENDIF
C
      IF(LPSF) THEN
        IF(ISEC.EQ.1) GO TO 201
      ENDIF
C
C  ----  Energies deposited in different bodies and detectors.
C
C  ----  Tallying the spectra from energy-deposition detectors.
      IF(NDEDEF.GT.0) THEN
        DO KD=1,NDEDEF
          DEDE(KD)=0.0D0
        ENDDO
        DO KB=1,NBODY
          IDET=KBDE(KB)
          IF(IDET.NE.0) THEN
            DEDE(IDET)=DEDE(IDET)+DEBO(KB)
          ENDIF
        ENDDO
C
        DO KD=1,NDEDEF
          TDED(KD)=TDED(KD)+DEDE(KD)
          TDED2(KD)=TDED2(KD)+DEDE(KD)**2
          IF(DEDE(KD).GT.1.0D-5) THEN
            IF(LDELOG(KD)) THEN
              KE=1.0D0+(LOG(DEDE(KD))-EDELL(KD))*RBDEEL(KD)
            ELSE
              KE=1.0D0+(DEDE(KD)-EDEL(KD))*RBDEE(KD)
            ENDIF
            IF(KE.GT.0.AND.KE.LE.NDECH(KD)) THEN
              DET(KD,KE)=DET(KD,KE)+1.0D0
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      DO KB=1,NBODY
        TDEBO(KB)=TDEBO(KB)+DEBO(KB)
        TDEBO2(KB)=TDEBO2(KB)+DEBO(KB)**2
      ENDDO
C  --  Average energies 'collected' by impact detectors.
      IF(NDIDEF.GT.0) THEN
        DO KD=1,NDIDEF
          TDID(KD)=TDID(KD)+DEDI(KD)
          TDID2(KD)=TDID2(KD)+DEDI(KD)**2
        ENDDO
      ENDIF
C  --  Final state counters.
      DO I=1,3
        PRIM(I)=PRIM(I)+DPRIM(I)
        PRIM2(I)=PRIM2(I)+DPRIM(I)**2
        DO K=1,3
          SEC(K,I)=SEC(K,I)+DSEC(K,I)
          SEC2(K,I)=SEC2(K,I)+DSEC(K,I)**2
        ENDDO
      ENDDO
C
C  ------------------------  The simulation of the shower ends here.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RETURN
      END









