J/A+A/539/A80       SDSS DR7 superclusters. The catalogues   (Liivamagi+, 2012)
================================================================================
SDSS DR7 superclusters. The catalogues.
    Liivamagi L.J., Tempel E., Saar E.
   <Astron. Astrophys. 539, A80 (2012)>
   =2012A&A...539A..80L
================================================================================
ADC_Keywords: Clusters, galaxy
Keywords: large-scale structure of the Universe - galaxies: clusters: general

Abstract:
    We have constructed a set of supercluster catalogues for the galaxies
    from the SDSS survey main and luminous red galaxy (LRG) flux-limited
    samples. To delineate superclusters, we calculated luminosity density
    fields using the B_3_-spline kernel of the radius of 8Mpc/h for the
    main sample and 16Mpc/h for the LRG sample and define regions with
    densities over a selected threshold as superclusters, while utilising
    almost the whole volume of both samples. We created two types of
    catalogues, one with an adaptive local threshold and a set of
    catalogues with different global thresholds. We describe the
    supercluster catalogues and their general properties. Using smoothed
    bootstrap, we find uncertainty estimates for the density field and use
    these to attribute confidence levels to the catalogue objects. We have
    also created a test catalogue for the galaxies from the Millennium
    simulation to compare the simulated and observed superclusters and to
    clarify the methods we use. We find that the superclusters are
    well-defined systems, and the properties of the superclusters of the
    main and LRG samples are similar. We also show that with adaptive
    local thresholds we get a sample of superclusters, the properties of
    which do not depend on their distance from the observer. The
    Millennium galaxy catalogue superclusters are similar to those
    observed.

Description:
    Excerpt of the SDSS DR7 supercluster catalogue containing four tables:
    two for the main and two for LRG samples. For both samples, one table
    contains the fixed threshold superclusters and the adaptively assigned
    threshold superclusters. Full dataset for all density thresholds is
    accessible at http://meteo.physic.ut.ee/~juhan/super/ .

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
main50.dat    274      982   Main sample D=5.0 mdens threshold catalog
mainadap.dat  284     1313   Main sample adaptive threshold catalog
lrg44.dat     264     3761   Luminous Red Galaxy (LRG) sample D=4.4 mdens
                             threshold catalog
lrgadap.dat   274     2701   Luminous Red Galaxy (LRG) sample adaptive
                             threshold catalog
--------------------------------------------------------------------------------

See also:
   J/A+A/514/A102 : SDSS DR7 groups of galaxies (Tago+, 2010)

Byte-by-byte Description of file: main50.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1- 12  A12   ---       UID       Name (NNN+NNN+NNNN)
  16- 20  I5    ---       ID        [1/1440] Short unique ID
  27- 30  I4    ---       Ngal      Number of galaxies
  38- 40  I3    ---       Ngr       Number of groups (G1)
  48- 52  I5    Mpc3      Vol       Supercluster volume (in (Mpc/h)^3^ units)
  53- 54  A2    ---       ---       [.0]
  61- 68  F8.3 10+10Lsun  dlum      Supercluster luminosity from density field
  74- 82  F9.4 10+10Lsun  slum      Sum of supercluster galaxy luminosities (G3)
  88- 96  F9.4 10+10Lsun  wlum      Sum of weighted supercluster galaxy
                                     luminosities
 105-110  F6.3  ---       Peak      Value of highest supercluster peak
                                     density (in mdens units) (G4)
     118  I1    ---       Nsps      Number of the density peaks in supercluster
 124-130  F7.3  deg       RA1deg    Right ascension (J2000) of highest density
                                     peak
 137-142  F6.3  deg       DE1deg    Declination (J2000) of highest density peak
 148-154  F7.3  Mpc       Dist1     Distance to highest density peak (Mpc/h)
 160-166  F7.3  deg       RA2deg    Right ascension (J2000) of mass centre
 173-178  F6.3  deg       DE2deg    Declination (J2000) of mass centre
 184-190  F7.3  Mpc       Dist2     Distance to mass centre (Mpc/h)
 195-202  F8.4  Mpc       Diam      Supercluster diameter (Mpc/h)
 209-214  I6    ---       mgID      Marker galaxy ID (G5)
 220-226  F7.3  deg       RA3deg    Right ascension (J2000) of marker galaxy
 233-238  F6.3  deg       DE3deg    Declination (J2000) of marker galaxy
 243-250  F8.6  ---       mgz       Redshift of marker galaxy
 256-262  F7.4  ---       Conf      Supercluster confidence estimate
     268  I1    ---       Edge      [0/1] 1=Mask boundary contact
     274  I1    ---       Comp      [0/5] Number of splitting objects
--------------------------------------------------------------------------------

Byte-by-byte Description of file: mainadap.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1- 12  A12   ---       UID       Name (NNN+NNN+NNNN)
  16- 20  I5    ---       ID        [1/4060] Short unique ID
  27- 29  F3.1  ---       D         Supercluster adaptive threshold
                                     (in mdens units) (G4)
      30  A1    ---       ---       [0]
  37- 40  I4    ---       Ngal      Number of galaxies
  48- 50  I3    ---       Ngr       Number of groups (G1)
  58- 62  I5    Mpc3      Vol       Supercluster volume (in (Mpc/h)^3^ units)
  63- 64  A2    ---       ---       [.0]
  71- 78  F8.3 10+10Lsun  dlum      Supercluster luminosity from density field
  84- 92  F9.4 10+10Lsun  slum      Sum of supercluster galaxy luminosities (G3)
  98-106  F9.4 10+10Lsun  wlum      Sum of weighted supercluster galaxy
                                     luminosities
 115-120  F6.3  ---       Peak      Value of highest supercluster peak density
                                     (in mdens units) (G4)
     128  I1    ---       Nsps      Number of the density peaks in supercluster
 134-140  F7.3  deg       RA1deg    Right ascension (J2000) of highest density
                                     peak
 147-152  F6.3  deg       DE1deg    Declination (J2000) of highest density peak
 158-164  F7.3  Mpc       Dist1     Distance to highest density peak (in Mpc/h)
 170-176  F7.3  deg       RA2deg    Right ascension (J2000) of mass centre
 183-188  F6.3  deg       DE2deg    Declination (J2000) of mass centre
 194-200  F7.3  Mpc       Dist2     Distance to mass centre (in Mpc/h)
 205-212  F8.4  Mpc       Diam      Supercluster diameter (in Mpc/h)
 219-224  I6    ---       mgID      Marker galaxy ID (in Tago et al 2010)
 230-236  F7.3  deg       RA3deg    Right ascension (J2000) of marker galaxy
 243-248  F6.3  deg       DE3deg    Declination (J2000) of marker galaxy
 253-260  F8.6  ---       mgz       Redshift of marker galaxy
 266-272  F7.4  ---       Conf      Supercluster confidence estimate
     278  I1    ---       Edge      [0/1] 1=Mask boundary contact
     284  I1    ---       Comp      [0/6] Number of splitting objects
--------------------------------------------------------------------------------

Byte-by-byte Description of file: lrg44.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1- 12  A12   ---       UID       Name (NNN+NNN+NNNN)
  16- 20  I5    ---       ID        [1/6115] Short unique ID
  27- 30  I4    ---       Ngal      Number of galaxies
  37- 42  I6    Mpc3      Vol       Supercluster volume (in (Mpc/h)^3^ units)
  43- 44  A2    ---       ---       [.0]
  51- 58  F8.3  10+10Lsun dlum      Supercluster luminosity from density field
  64- 72  F9.4  10+10Lsun slum      Sum of supercluster galaxy luminosities (G3)
  78- 86  F9.4  10+10Lsun wlum      Sum of weighted supercluster galaxy
                                     luminosities
  95-100  F6.3  ---       Peak      Value of highest supercluster peak density
                                     (in mdens unit) (G4)
 107-108  I2    ---       Nsps      Number of the density peaks in supercluster
 114-120  F7.3  deg       RA1deg    Right ascension (J2000) of highest density
                                     peak
 127-132  F6.3  deg       DE1deg    Declination (J2000) of highest density peak
 137-144  F8.3  Mpc       Dist1     Distance to highest density peak (in Mpc/h)
 150-156  F7.3  deg       RA2deg    Right ascension (J2000) of mass centre
 163-168  F6.3  deg       DE2deg    Declination (J2000) of mass centre
 173-180  F8.3  Mpc       Dist2     Distance to mass centre (in Mpc/h)
 185-192  F8.4  Mpc       Diam      Supercluster diameter (in Mpc/h)
 199-204  I6    ---       mgID      Marker galaxy ID
 210-216  F7.3  deg       RA3deg    Right ascension (J2000) of marker galaxy
 223-228  F6.3  deg       DE3deg    Declination (J2000) of marker galaxy
 233-240  F8.6  ---       mgz       Redshift of marker galaxy
 246-252  F7.4  ---       Conf      Supercluster confidence estimate
     258  I1    ---       Edge      [0/1] 1=Mask boundary contact
     264  I1    ---       Comp      [0/4] Number of splitting objects
--------------------------------------------------------------------------------

Byte-by-byte Description of file: lrgadap.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1- 12  A12   ---       UID       Name (NNN+NNN+NNNN)
  16- 20  I5    ---       ID        [1/10865] Short unique ID
  27- 29  F3.1  ---       D         Supercluster adaptive threshold
                                     (in mdens units) (G4)
      31  A1    ---       ---       [0]
  37- 40  I4    ---       Ngal      Number of galaxies
  47- 52  I6    Mpc3      Vol       Supercluster volume (in (Mpc/h)^3^ units)
  53- 54  A2    ---       ---       [.0]
  61- 68  F8.3  10+10Lsun dlum      Supercluster luminosity from density field
  74- 82  F9.4  10+10Lsun slum      Sum of supercluster galaxy luminosities (G3)
  88- 96  F9.4  10+10Lsun wlum      Sum of weighted supercluster galaxy
                                     luminosities
 105-110  F6.3  ---       Peak      Value of highest supercluster peak density
                                     (in mdens units) (G4)
 117-118  I2    ---       Nsps      Number of the density peaks in supercluster
 124-130  F7.3  deg       RA1deg    Right ascension (J2000) of highest density
                                     peak
 137-142  F6.3  deg       DE1deg    Declination (J2000) of highest density peak
 147-154  F8.3  Mpc       Dist1     Distance to highest density peak
 160-166  F7.3  deg       RA2deg    Right ascension (J2000) of mass centre
 173-178  F6.3  deg       DE2deg    Declination (J2000) of mass centre
 183-190  F8.3  Mpc       Dist2     Distance to mass centre
 195-202  F8.4  Mpc       Diam      Supercluster diameter
 209-214  I6    ---       mgID      Marker galaxy ID
 220-226  F7.3  deg       RA3deg    Right ascension (J2000) of marker galaxy
 233-238  F6.3  deg       DE3deg    Declination (J2000) of marker galaxy
 243-250  F8.6  ---       mgz       Redshift of marker galaxy
 256-262  F7.4  ---       Conf      Supercluster confidence estimate
     268  I1    ---       Edge      [0/1] Mask boundary contact
     274  I1    ---       Comp      [0/7] Number of splitting objects
--------------------------------------------------------------------------------

Global notes:
Note (G1): SDSS main sample group data is taken from Tago et al.
     (2010, Cat. J/A+A/514/A102)
Note (G3): all luminosities are in SDSS r-band
Note (G4): mean density value "mdens" is
   * for the main sample: mdens=1.526e-2*10^10^h^-1^Lsun(Mpc/h)^-3
   * for the LRG sample:  mdens=8.148e-4*10^10^h^-1^Lsun(Mpc/h)^-3
Note (G5): Galaxy IDs are from Tago et al. (2010, Cat. J/A+A/514/A102)
--------------------------------------------------------------------------------

Acknowledgements:
    Lauri Juhan Liivamagi, wire(at)ut.ee

================================================================================
(End)    Lauri Juhan Liivamagi [Tartu Obs.], Patricia Vannier [CDS]  03-Jan-2012
