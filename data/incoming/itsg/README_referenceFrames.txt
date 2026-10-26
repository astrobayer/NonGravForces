Reference frames
----------------
We estimate our own GNSS constellation (GSP/Galileo)
based on the IGS groundstation network (about 400 stations per day)
using the raw observation approach on daily basis.
(orbit/clocks/biases, station coordinates, Earth Rotation Parameters, Troposphere, Ionopshere).

The IGS core station coordinates are aligned to ITRF2020/IGS20 (or updated versions)
via strict No-Net-Rotation and No-Net-Translation constraint
(CoN, Center of Network).

The GNSS and LEO orbits are estimated in a space fixed (quasi-inertial) frame
in the Center of Mass (CoM).

Earth Rotation parameters are co-estimated in GNSS network adjustment.

The instantaneous transformation between CoN and CoM is realized by
+ the annual and semiannual ITRF2020 model
+ center of mass correction of ocean tidal loading (FES2022)

All times are GPS times.

............................................................................

File formats:

operational/geocenter/ (since 2025-01-01):
-----------------------
The geocenter model  must be added to CoM coordinates in TRF to get coordinates in CoN.

Data columns of gecenterCoM2CoN.txt.gz
- Time:   GPS [MJD]
- data0:  x in TRF [m]
- data1:  y in TRF [m]
- data2:  z in TRF [m]

............................................................................

operational/earthRotation/ (since 2025-01-01):
---------------------------
Estimated Earth Rotation including all models and corrections
following IERS2010 conventions
(pmsdnut2, utlibr, eraXy06, HF-EOP model: desai_jgrb51665-sup-0002-ds01)

Data columns of eop.txt.gz
- Time:   GPS [MJD]
- data0:  xp [rad]
- data1:  xp [rad]
- data2:  sp [rad]
- data3:  UT1-UTC [sconds]
- data4:  length of day (LOD) [seconds]
- data5:  X [rad]
- data6:  Y [rad]
- data7:  S [rad]

Data columns of rotationTRF2CRF.txt.gz
- Time:   GPS [MJD]
- data0:  xx  [-]
- data1:  xy  [-]
- data2:  xz  [-]
- data3:  yx  [-]
- data4:  yy  [-]
- data5:  yz  [-]
- data6:  zx  [-]
- data7:  zy  [-]
- data8:  zz  [-]
