# fetchData probe report, 2026-10-26 18:10 UTC

Written by `tools/fetchData.py` on a CI runner (hosted: GitHub; classic: self-hosted). Staging only: see the script's docstring. Each line is a request and what came back.

## Dst (O-9)

- `https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2.text` -> 200, 18021 bytes, sha256 `be4aca7bdf030a3eb6ed584fb6960913e3be909d5627b2a69fe93a4765962d7f` -- OMNI2 format description
  - committed `omni2/omni2.text` (18021 bytes)
- `https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_2007.dat` -> 200, 2873280 bytes, sha256 `2bdae448301e6de3863e78643178628e3bcd39dee2aea111fe523d1e393217e7`
  - window 2007-12-25..2008-01-07 from 2007: 168 hourly lines
  - committed `omni2/omni2_2007_2007-12-25_2008-01-07.dat` (55104 bytes)
- `https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_2008.dat` -> 200, 2881152 bytes, sha256 `3115c8f3c1d28eaca6a4ada84ca2f36865d36d2b2db1c4f8a4c4ba667ac1e274`
  - window 2007-12-25..2008-01-07 from 2008: 168 hourly lines
  - committed `omni2/omni2_2008_2007-12-25_2008-01-07.dat` (55104 bytes)
- `https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_2019.dat` -> 200, 2873280 bytes, sha256 `8f4e57f313687b4fb86865dbb2df06adf414c7e496b99893b86013d91744b48e`
  - window 2019-12-25..2020-01-07 from 2019: 168 hourly lines
  - committed `omni2/omni2_2019_2019-12-25_2020-01-07.dat` (55104 bytes)
- `https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_2020.dat` -> 200, 2881152 bytes, sha256 `11c2de98f1fd1ebe0156ca36c439fd4063d28c95e3706ec8752e7f59eea41700`
  - window 2019-12-25..2020-01-07 from 2020: 168 hourly lines
  - committed `omni2/omni2_2020_2019-12-25_2020-01-07.dat` (55104 bytes)
- `https://wdc.kugi.kyoto-u.ac.jp/dst_final/200712/dst0712.for.request` -> 200, 3752 bytes, sha256 `7b3ff06cf0cc4560497e8ca9030d9b67f6f118e396ff40657daef86b0b7fc432`
  - committed `kyoto/dst0712.for.request` (3752 bytes)
- `https://wdc.kugi.kyoto-u.ac.jp/dst_final/200801/dst0801.for.request` -> 200, 3752 bytes, sha256 `a09c56d390e89bad6d7dd89a3db7858503e52649f6087fd57aee594b8707b23b`
  - committed `kyoto/dst0801.for.request` (3752 bytes)
- `https://wdc.kugi.kyoto-u.ac.jp/dst_final/201912/dst1912.for.request` -> 200, 3752 bytes, sha256 `8813b818ffaa2a04521be4ed0f225a8164e490c520b4d4c408dc1c0211f8e99a`
  - committed `kyoto/dst1912.for.request` (3752 bytes)
- `https://wdc.kugi.kyoto-u.ac.jp/dst_final/202001/dst2001.for.request` -> 200, 3752 bytes, sha256 `f7545be24ea25aecb985bb5b73f760cd620d2897bbe041bcd6145761c45da483`
  - committed `kyoto/dst2001.for.request` (3752 bytes)

## F30, the 30 cm solar radio flux (O-9)

- `https://solar.nro.nao.ac.jp/norp/data/daily/` -> failed -- URLError: <urlopen error timed out>
- `https://solar.nro.nao.ac.jp/norp/data/daily/TYKW-NoRP_dailyflux.txt` -> URLError: <urlopen error timed out>
listing of `ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/`: `radio_flux_absolute_observation.txt`, `radio_flux_adjusted_observation.txt`
- `ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_absolute_observation.txt` -> , 5260151 bytes, sha256 `50c8df9ab77ebe0be6387a0d596e21956e71f5016458d4252790d8d6514cf46f`
  - committed `cls/radio_flux_absolute_observation.txt.head.txt` (8375 bytes)
  - `radio_flux_absolute_observation.txt` daily windows: 426 dated lines
  - committed `cls/radio_flux_absolute_observation_2007-09_2008-03_2019-09_2020-03.txt` (83447 bytes)
- `ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_adjusted_observation.txt` -> , 5260147 bytes, sha256 `6423391a3bd0b9fb073170c0a896ff9b6c22f1db0e5ed6c554e7a3a358f5e422`
  - committed `cls/radio_flux_adjusted_observation.txt.head.txt` (8371 bytes)
  - `radio_flux_adjusted_observation.txt` daily windows: 426 dated lines
  - committed `cls/radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt` (83443 bytes)

## ITSG operational products: attitude, non-conservative forces, orbits (N-12, O-19)

- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/README.txt` -> 200, 8448 bytes, sha256 `09a1a4f320edd8125557d67a6ec0c43a657a8fbfb4aa6ba44f6c85f862c02d05`
  - committed `itsg/README.txt` (8448 bytes)
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/README_referenceFrames.txt` -> 200, 2000 bytes, sha256 `290f86f7248d6d8122fc3d4e7bb801a203104ba5a27800ef7dadc69200a0ccb9`
  - committed `itsg/README_referenceFrames.txt` (2000 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-1/attitude/2008/`: 366 entries; for 2008-01-01: `GRACE-1_attitude_2008-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-1/attitude/2008/GRACE-1_attitude_2008-01-01.txt.gz` -> 200, 438416 bytes, sha256 `b6689a4f2c0fc1f28a5db9eba82dd28a2f4a8f6bf1a77d53b2b38f264446b386`
  - `GRACE-1_attitude_2008-01-01.txt`: 1123673 bytes, sha256 `e51dfdc9a0fe89acaa3e95efac18bf34ad39e8a33d092f65cebf4cee9c612bdc`
  - committed `itsg/GRACE-1/attitude/GRACE-1_attitude_2008-01-01.txt.head.txt` (7363 bytes)
  - `GRACE-1_attitude_2008-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACE-1/attitude/GRACE-1_attitude_2008-01-01.txt` (78342 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-2/attitude/2008/`: 366 entries; for 2008-01-01: `GRACE-2_attitude_2008-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-2/attitude/2008/GRACE-2_attitude_2008-01-01.txt.gz` -> 200, 437714 bytes, sha256 `11f8319ce529933c3945fe75ad2027f24a5e77e6b4677c27b65904af01dcb1fe`
  - `GRACE-2_attitude_2008-01-01.txt`: 1123673 bytes, sha256 `f1c8bc021f8100fbd43370e84cce112f780280163dd841098831301e60316b58`
  - committed `itsg/GRACE-2/attitude/GRACE-2_attitude_2008-01-01.txt.head.txt` (7363 bytes)
  - `GRACE-2_attitude_2008-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACE-2/attitude/GRACE-2_attitude_2008-01-01.txt` (78342 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/attitude/2020/`: 359 entries; for 2020-01-01: `GRACEFO-1_attitude_2020-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/attitude/2020/GRACEFO-1_attitude_2020-01-01.txt.gz` -> 200, 439031 bytes, sha256 `6c845b358f91175aad28f87f70ef4605e0b0abe253b015fc88fc24c5dcd5ba72`
  - `GRACEFO-1_attitude_2020-01-01.txt`: 1123673 bytes, sha256 `1fa543419822011f00397f35de1d7a3dcbe6594d135cc672fda00a59955975ea`
  - committed `itsg/GRACEFO-1/attitude/GRACEFO-1_attitude_2020-01-01.txt.head.txt` (7363 bytes)
  - `GRACEFO-1_attitude_2020-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/attitude/GRACEFO-1_attitude_2020-01-01.txt` (78342 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/attitude/2020/`: 359 entries; for 2020-01-02: `GRACEFO-1_attitude_2020-01-02.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/attitude/2020/GRACEFO-1_attitude_2020-01-02.txt.gz` -> 200, 438820 bytes, sha256 `c7614304a1b30ef15c7ce2660f9c3c6a7690c0bbf3ac91cfc914b06b803d2f89`
  - `GRACEFO-1_attitude_2020-01-02.txt`: 1123673 bytes, sha256 `0ed571132c4d94aa2be6dd34d10f4df5f3ef9bf22fe5893e7f79d7470d135ba5`
  - committed `itsg/GRACEFO-1/attitude/GRACEFO-1_attitude_2020-01-02.txt.head.txt` (7363 bytes)
  - `GRACEFO-1_attitude_2020-01-02.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/attitude/GRACEFO-1_attitude_2020-01-02.txt` (78342 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/attitude/2020/`: 359 entries; for 2020-01-03: `GRACEFO-1_attitude_2020-01-03.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/attitude/2020/GRACEFO-1_attitude_2020-01-03.txt.gz` -> 200, 438903 bytes, sha256 `8776843d45c936a32cd8fc461a63e7e8371cc8e6c785ef400add4cc7a7582d4e`
  - `GRACEFO-1_attitude_2020-01-03.txt`: 1123673 bytes, sha256 `7ac48f84a06893f28cd3b289be3d5e05a7ee7b8c36dc7d883ed2298236c939d9`
  - committed `itsg/GRACEFO-1/attitude/GRACEFO-1_attitude_2020-01-03.txt.head.txt` (7363 bytes)
  - `GRACEFO-1_attitude_2020-01-03.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/attitude/GRACEFO-1_attitude_2020-01-03.txt` (78342 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-2/attitude/2020/`: 359 entries; for 2020-01-01: `GRACEFO-2_attitude_2020-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-2/attitude/2020/GRACEFO-2_attitude_2020-01-01.txt.gz` -> 200, 437144 bytes, sha256 `523f83d00bff30dc8193ed06e876cb24032e14da4b8d32c00913a3480ce2e1dd`
  - `GRACEFO-2_attitude_2020-01-01.txt`: 1123673 bytes, sha256 `71ffca04a3cdae3c4637ba0e5c8e7b4814fd902f750bd5f3e3687e3550c69c48`
  - committed `itsg/GRACEFO-2/attitude/GRACEFO-2_attitude_2020-01-01.txt.head.txt` (7363 bytes)
  - `GRACEFO-2_attitude_2020-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-2/attitude/GRACEFO-2_attitude_2020-01-01.txt` (78342 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-1/nonConservativeForces/2008/`: 366 entries; for 2008-01-01: `GRACE-1_nonConservativeForces_2008-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-1/nonConservativeForces/2008/GRACE-1_nonConservativeForces_2008-01-01.txt.gz` -> 200, 341164 bytes, sha256 `4e720c8b730435c94f17a334e8ec291e414673f557f2fd6b2acdcb89e295b985`
  - `GRACE-1_nonConservativeForces_2008-01-01.txt`: 898958 bytes, sha256 `5f3331c48d9c9a35c29bf3a68f9407d9ad93dcaef1695e661d8cad3391c6e27a`
  - committed `itsg/GRACE-1/nonConservativeForces/GRACE-1_nonConservativeForces_2008-01-01.txt.head.txt` (5910 bytes)
  - `GRACE-1_nonConservativeForces_2008-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACE-1/nonConservativeForces/GRACE-1_nonConservativeForces_2008-01-01.txt` (62693 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-2/nonConservativeForces/2008/`: 366 entries; for 2008-01-01: `GRACE-2_nonConservativeForces_2008-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-2/nonConservativeForces/2008/GRACE-2_nonConservativeForces_2008-01-01.txt.gz` -> 200, 340654 bytes, sha256 `8326d5563d5031f725f5bb858084df9b8b4266dc6f7372a5c9db79994c83849b`
  - `GRACE-2_nonConservativeForces_2008-01-01.txt`: 898958 bytes, sha256 `eb639781cb728094eb963c68ed7bc05cfb238717e6e62778eae0cb16e0f17610`
  - committed `itsg/GRACE-2/nonConservativeForces/GRACE-2_nonConservativeForces_2008-01-01.txt.head.txt` (5910 bytes)
  - `GRACE-2_nonConservativeForces_2008-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACE-2/nonConservativeForces/GRACE-2_nonConservativeForces_2008-01-01.txt` (62693 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/nonConservativeForces/2020/`: 359 entries; for 2020-01-01: `GRACEFO-1_nonConservativeForces_2020-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/nonConservativeForces/2020/GRACEFO-1_nonConservativeForces_2020-01-01.txt.gz` -> 200, 340331 bytes, sha256 `57023e3f6383e843a325b597b81c5e24624397fc4e9b135b481f90721338034b`
  - `GRACEFO-1_nonConservativeForces_2020-01-01.txt`: 898958 bytes, sha256 `a8b5ad30c98014c5859163f307aa489f493a15fe82f2cf490ee2e9e1e9482cc9`
  - committed `itsg/GRACEFO-1/nonConservativeForces/GRACEFO-1_nonConservativeForces_2020-01-01.txt.head.txt` (5910 bytes)
  - `GRACEFO-1_nonConservativeForces_2020-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/nonConservativeForces/GRACEFO-1_nonConservativeForces_2020-01-01.txt` (62693 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/nonConservativeForces/2020/`: 359 entries; for 2020-01-02: `GRACEFO-1_nonConservativeForces_2020-01-02.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/nonConservativeForces/2020/GRACEFO-1_nonConservativeForces_2020-01-02.txt.gz` -> 200, 340232 bytes, sha256 `5528e19a783b924add176a849d2904217437d60c7e7f517efff4ddda0b2b4914`
  - `GRACEFO-1_nonConservativeForces_2020-01-02.txt`: 898958 bytes, sha256 `25c984db98269d855518e2176837012c087755f8cb976ab60a333ca803a3db3d`
  - committed `itsg/GRACEFO-1/nonConservativeForces/GRACEFO-1_nonConservativeForces_2020-01-02.txt.head.txt` (5910 bytes)
  - `GRACEFO-1_nonConservativeForces_2020-01-02.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/nonConservativeForces/GRACEFO-1_nonConservativeForces_2020-01-02.txt` (62693 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/nonConservativeForces/2020/`: 359 entries; for 2020-01-03: `GRACEFO-1_nonConservativeForces_2020-01-03.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/nonConservativeForces/2020/GRACEFO-1_nonConservativeForces_2020-01-03.txt.gz` -> 200, 340753 bytes, sha256 `3a8316726b8fc71de60492affe4af63e4789a15e4f451e279e59b4595a4430c9`
  - `GRACEFO-1_nonConservativeForces_2020-01-03.txt`: 898958 bytes, sha256 `afef190515daa08664bb2b2e3de60a1bbe5e01b03762b2b82dc3a86275855092`
  - committed `itsg/GRACEFO-1/nonConservativeForces/GRACEFO-1_nonConservativeForces_2020-01-03.txt.head.txt` (5910 bytes)
  - `GRACEFO-1_nonConservativeForces_2020-01-03.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/nonConservativeForces/GRACEFO-1_nonConservativeForces_2020-01-03.txt` (62693 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-2/nonConservativeForces/2020/`: 359 entries; for 2020-01-01: `GRACEFO-2_nonConservativeForces_2020-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-2/nonConservativeForces/2020/GRACEFO-2_nonConservativeForces_2020-01-01.txt.gz` -> 200, 339753 bytes, sha256 `3d3490e1ed629867a6481cd8c80ae574c3a0dd6f36f5625ecd4a76dba9b95723`
  - `GRACEFO-2_nonConservativeForces_2020-01-01.txt`: 898958 bytes, sha256 `e6dbe3507f09c15ccc37413ccda09fbb5b70e289c3a7c1b47510a85766197c58`
  - committed `itsg/GRACEFO-2/nonConservativeForces/GRACEFO-2_nonConservativeForces_2020-01-01.txt.head.txt` (5910 bytes)
  - `GRACEFO-2_nonConservativeForces_2020-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-2/nonConservativeForces/GRACEFO-2_nonConservativeForces_2020-01-01.txt` (62693 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-1/reducedDynamicOrbit/2008/`: 366 entries; for 2008-01-01: `GRACE-1_reducedDynamicOrbit_2008-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-1/reducedDynamicOrbit/2008/GRACE-1_reducedDynamicOrbit_2008-01-01.txt.gz` -> 200, 661333 bytes, sha256 `40813ac1aeb1ad2177e1b19e48ced675d272b240ba88247a34289d8e4ca4801b`
  - `GRACE-1_reducedDynamicOrbit_2008-01-01.txt`: 2247258 bytes, sha256 `c8a77a8a17f28e7f9adaf428c0ba2ce90153358aadfc4b75090b6c23d9e53a2d`
  - committed `itsg/GRACE-1/reducedDynamicOrbit/GRACE-1_reducedDynamicOrbit_2008-01-01.txt.head.txt` (14638 bytes)
  - `GRACE-1_reducedDynamicOrbit_2008-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACE-1/reducedDynamicOrbit/GRACE-1_reducedDynamicOrbit_2008-01-01.txt` (156597 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-2/reducedDynamicOrbit/2008/`: 366 entries; for 2008-01-01: `GRACE-2_reducedDynamicOrbit_2008-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACE-2/reducedDynamicOrbit/2008/GRACE-2_reducedDynamicOrbit_2008-01-01.txt.gz` -> 200, 661454 bytes, sha256 `e3d9d7463f5e7ecfa34acae3bc09f64619151d00971242883d2b64d448a51aa6`
  - `GRACE-2_reducedDynamicOrbit_2008-01-01.txt`: 2247258 bytes, sha256 `dfa1eba16e4dcd326754784dbe14b5ac69f8743a63b5745f4d2de97128620424`
  - committed `itsg/GRACE-2/reducedDynamicOrbit/GRACE-2_reducedDynamicOrbit_2008-01-01.txt.head.txt` (14638 bytes)
  - `GRACE-2_reducedDynamicOrbit_2008-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACE-2/reducedDynamicOrbit/GRACE-2_reducedDynamicOrbit_2008-01-01.txt` (156597 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/reducedDynamicOrbit/2020/`: 359 entries; for 2020-01-01: `GRACEFO-1_reducedDynamicOrbit_2020-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/reducedDynamicOrbit/2020/GRACEFO-1_reducedDynamicOrbit_2020-01-01.txt.gz` -> 200, 659073 bytes, sha256 `32608c82721ac7ca8e97291687475a25f880738e7009c3319fe916e425868f7b`
  - `GRACEFO-1_reducedDynamicOrbit_2020-01-01.txt`: 2247258 bytes, sha256 `b5252ce26776ab4e45c5cbb0eda4be06fc1e34196a72056ae818ed93c9d495e5`
  - committed `itsg/GRACEFO-1/reducedDynamicOrbit/GRACEFO-1_reducedDynamicOrbit_2020-01-01.txt.head.txt` (14638 bytes)
  - `GRACEFO-1_reducedDynamicOrbit_2020-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/reducedDynamicOrbit/GRACEFO-1_reducedDynamicOrbit_2020-01-01.txt` (156597 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/reducedDynamicOrbit/2020/`: 359 entries; for 2020-01-02: `GRACEFO-1_reducedDynamicOrbit_2020-01-02.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/reducedDynamicOrbit/2020/GRACEFO-1_reducedDynamicOrbit_2020-01-02.txt.gz` -> 200, 658882 bytes, sha256 `24d9065d9895388a41a4dc68978504d6b2e64e7bb8a826578822181594a1daf1`
  - `GRACEFO-1_reducedDynamicOrbit_2020-01-02.txt`: 2247258 bytes, sha256 `d20e91e00d24734fd7bbd0a23d9391098436b06d2dac3c81967540ef8e36e2a9`
  - committed `itsg/GRACEFO-1/reducedDynamicOrbit/GRACEFO-1_reducedDynamicOrbit_2020-01-02.txt.head.txt` (14638 bytes)
  - `GRACEFO-1_reducedDynamicOrbit_2020-01-02.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/reducedDynamicOrbit/GRACEFO-1_reducedDynamicOrbit_2020-01-02.txt` (156597 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/reducedDynamicOrbit/2020/`: 359 entries; for 2020-01-03: `GRACEFO-1_reducedDynamicOrbit_2020-01-03.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-1/reducedDynamicOrbit/2020/GRACEFO-1_reducedDynamicOrbit_2020-01-03.txt.gz` -> 200, 658883 bytes, sha256 `dc46ecf50d9841629a12a9ea5c47e00dbf6a68e86bb2448e5ec193260072524c`
  - `GRACEFO-1_reducedDynamicOrbit_2020-01-03.txt`: 2247258 bytes, sha256 `b24e43b4c2d864e315d0e98ac08cd32e7d42e85caba7de71a6872824ddf14883`
  - committed `itsg/GRACEFO-1/reducedDynamicOrbit/GRACEFO-1_reducedDynamicOrbit_2020-01-03.txt.head.txt` (14638 bytes)
  - `GRACEFO-1_reducedDynamicOrbit_2020-01-03.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-1/reducedDynamicOrbit/GRACEFO-1_reducedDynamicOrbit_2020-01-03.txt` (156597 bytes)
`https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-2/reducedDynamicOrbit/2020/`: 359 entries; for 2020-01-01: `GRACEFO-2_reducedDynamicOrbit_2020-01-01.txt.gz`
- `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/GRACEFO-2/reducedDynamicOrbit/2020/GRACEFO-2_reducedDynamicOrbit_2020-01-01.txt.gz` -> 200, 659035 bytes, sha256 `e90356d8869d8662dbaf8832b3ef57bd7139a4a69aeaf04e7c20999f449426d4`
  - `GRACEFO-2_reducedDynamicOrbit_2020-01-01.txt`: 2247258 bytes, sha256 `9845cd43a9b2c2d3dbde82d95c4907037bf2159937972a15391ffcd41705297d`
  - committed `itsg/GRACEFO-2/reducedDynamicOrbit/GRACEFO-2_reducedDynamicOrbit_2020-01-01.txt.head.txt` (14638 bytes)
  - `GRACEFO-2_reducedDynamicOrbit_2020-01-01.txt`: 600 of 8641 epochs kept; count rewritten: count line after the rule
  - committed `itsg/GRACEFO-2/reducedDynamicOrbit/GRACEFO-2_reducedDynamicOrbit_2020-01-01.txt` (156597 bytes)
