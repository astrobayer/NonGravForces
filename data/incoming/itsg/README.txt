
*******************************************************************
         WELCOME TO THE IfG FTP PRODUCT SERVER
         https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/
*******************************************************************

          Institute of Geodesy - Satellite orbit products

                   10. August 2022
                   revised on
                   16. October 2025

                   Graz University of Technology
                   Institute of Geodesy

                   Torsten Mayer-Guerr
                   mayer-guerr@tugraz.at
                   +43 316 873 6359

Adress:            Steyrergasse 30/III
                   8010 Graz, Austria

Introduction:
-------------

On this ftp server you can find following satellite orbit products produced at the Institute of Geodesy:

- Kinematic orbit
- Covariance
- Reduced dynamic orbit
- Attitude
- Non-conservative forces
- Neutral density

Access is granted without any registration and free of charge. All provided products can be used for scientific
research or any other application. Nevertheless, if you publish any results based on our products please
add an acknowledgment.

The satellite orbit products were processed with the GROOPS (Gravity Recovery Object Oriented Programming System)
Software.

GROOPS Open source version available at https://github.com/groops-devs/groops

Description of the processing strategy can be found in: Suesser-Rechberger, B., Krauss, S., Strasser, S., & Mayer-Guerr, T. (2022).
                                                        Improved precise kinematic LEO orbits based on the raw observation approach.
                                                        Advances in Space Research, 69(10), 3559-3570. https://doi.org/10.1016/j.asr.2022.03.014


                                                        Krauss, S., Behzadpour, S., Temmer, M., & Lhotka, C. (2020).
                                                        Exploring thermospheric variations triggered by severe geomagnetic storm
                                                        on 26 August 2018 using GRACE Follow-On data. Journal of Geophysical Research:
                                                        Space Physics, 125, e2019JA027731. https://doi.org/10.1029/2019JA027731

Data sampling:
--------------
- Kinematic orbit and covariance data sampling corresponds to the measurement epochs and is not equidistant.
  The epochs are given in GPS time and corrected for receiver clock errors.
  POSITIONS and COVARIANCES ARE NOT INTERPOLATED TO INTEGER SECONDS NOR A CONSTANT SAMPLING!!
- Reduced dynamic orbit, attitude and non-conservative forces data have a 10 s sampling. For easier interpolation
  both midnight epochs are included.
- Neutral density data have a sampling of 60s for release 1.0. TanDEM-X and TerraSAR-X have a 30s sampling. For release 2.0, CHAMP
  data have a sampling of 10s. For release 1.5, GRACE data have a sampling of 5s. GRACE-FO data also uses a sampling of 5s.


File naming convention:
------------------------
The filenames consist of three parts:satellite name, satellite orbit product and date.

e.g.: CHAMP_kinematicOrbit_2007-01-01.txt.gz
      CHAMP_kinematicOrbitCovariance_2007-01-01.txt.gz
      CHAMP_reducedDynamicOrbit_2007-01-01.txt.gz
      CHAMP_attitude_2007-01-01.txt.gz
      CHAMP_nonConservativeForces_2007-01-01.txt.gz
      CHAMP_neutralDensity_ACC_2007-01-01.txt.gz


File format description:
------------------------
All files are compressed.
GROOPS format is used which is a simple ASCII format.

- kinematicOrbit file: Contains the kinematic orbit position coordinates frame per epoch.
                       Position coordinates are given in the celestial reference frame.

  column    | content                             | unit
------------|-------------------------------------|---------------------------
     1      | GPS time                            | MJD
     2      | x-component of position             | m
     3      | y-component of position             | m
     4      | z-component of position             | m
 ............................................................................

- kinematicOrbitCovariance file: Contains the elements of the 3x3 covariance matrix per epoch in the following way: xx yy zz xy xz yz.
                                 Elements are given in the celestial reference frame.

  column    | content         | unit
------------|-----------------|------
     1      | GPS time        | MJD
     2      | xx-value        | [m^2]
     3      | yy-value        | [m^2]
     4      | zz-value        | [m^2]
     5      | xy-value        | [m^2]
     6      | xz-value        | [m^2]
     7      | yz-value        | [m^2]
 ............................................................................

- Reduced dynamic orbit file: Contains the position coordinates and the velocity components of the reduced dynamic orbit per epoch.
                              Position coordinates and velocities are given in the celestial reference frame.

  column    | content                    | unit
------------|----------------------------|------
     1      | GPS time                   | MJD
     2      | x-component of position    | m
     3      | y-component of position    | m
     4      | z-component of position    | m
     5      | x-component of velocity    | m/s
     6      | y-component of velocity    | m/s
     7      | z-component of velocity    | m/s
 ............................................................................

- attitude file: Contains the rotation from the satellite reference frame to the celestial reference frame in terms of quaternions per epoch.

  column    | content                    | unit
------------|----------------------------|------
     1      | GPS time                   | MJD
     2      | scalar part                | -
     3      | qx-value                   | -
     4      | qy-value                   | -
     5      | qz-value                   | -
 ............................................................................

- nonConservativeForces file: Contains the accelerations due to non-conservative forces per epoch.
                              Accelerations are given in the satellite reference frame.

 column     | content                            | unit
------------|------------------------------------|------
     1      | GPS time                           | MJD
     2      | x-component of the acceleration    | m/s^2
     3      | y-component of the acceleration    | m/s^2
     4      | z-component of the acceleration    | m/s^2
 ............................................................................

- neutralDensity file release 1.0 (CHAMP, GRACE): Contains the estimated neutral density per epoch. Abbreviation ACC in the file name indicates
                                   that the density estimate is based on on-board accelerometer measurements. The identifier KIN
                                   points out that the solution is based on kinematic orbit information.

  column    | content                             | unit
------------|-------------------------------------|------
     1      | GPS time                            | MJD
     2      | neutral density                     | kg/m^3
 ............................................................................

- neutralDensity files release 1.5 / 2.0 (CHAMP, GRACE, GRACE-FO, TanDEM-X, TerraSAR-X):
                                   Contains the estimated neutral density per epoch and the satellite's position (GRS80) and local solar time (LST).
                                   CHAMP, GRACE, GRACE-FO density estimates are derived from accelerometer measurements.
                                   TanDEM-X and TerraSAR-X density estimates are derived from kinematic/POD observations.

  column    | content                             | unit
------------|-------------------------------------|------
     1      | GPS time                            | MJD
     2      | Longitude                           | degree
     3      | Latitude                            | degree
     4      | Altitude                            | m
     5      | local solar time                    | hour (decimal)
     6      | neutral density                     | kg/m^3
 ............................................................................
