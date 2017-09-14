!in_gr_common

include dtb/elements.dtb
include dtb/logk_eq36_hkf.dtb
include dtb/kinetics_usgs.dtb

SYSTEM
  TdgC 25
  SI  INERT    SIO2(AQ)
  NA  INERT    NA+
  K   INERT    K+
  AL  INERT    AL+3
  H   MOBILE   H+    4.30
  CL  BALANCE  CL-
END SYSTEM

DYNAMIC.ROCK
  QTZ   species QUARTZ-A   kinetics QUARTZ       mode D radius 1.0E-6 mole 4.
  ALB   species ALBITE-LOW kinetics ALBITE       mode D radius 1.0E-6 mole 1.
  KFS   species MICROCLINE kinetics K-FELDSPAR   mode D radius 1.0E-6 mole 1.
  MUS   species MUSCOVITE  kinetics MUSCOVITE    mode D radius 1.0E-9 mole 1.E-3
  KAOL  species KAOLINITE  kinetics KAOLINITE    mode D radius 1.0E-9 mole 1.E-3
  PARAG species PARAGONITE kinetics PARAGONITE   mode D radius 1.0E-9 mole 1.E-3
  GIBBS species GIBBSITE   kinetics GIBBSITE     mode D radius 1.0E-9 mole 1.E-3
END DYNAMIC.ROCK

  QTZ   species QUARTZ-A   kinetics QUARTZ       mode D radius 1.0E-6 volume 27.50
  ALB   species ALBITE-LOW kinetics ALBITE       mode D radius 1.0E-6 volume 35.75
  KFS   species MICROCLINE kinetics K-FELDSPAR   mode D radius 1.0E-6 volume 35.75
  MUS   species MUSCOVITE  kinetics MUSCOVITE    mode D radius 1.0E-9 volume 0.001
  KAOL  species KAOLINITE  kinetics KAOLINITE    mode D radius 1.0E-9 volume 0.001
  PARAG species PARAGONITE kinetics PARAGONITE   mode D radius 1.0E-9 volume 0.001
  GIBBS species GIBBSITE   kinetics GIBBSITE     mode D radius 1.0E-9 volume 0.001

DYNAMIC
  POROSITY 0.20
  TUNIT    YEAR
END DYNAMIC

SOLVENT
  MODEL DAV_1
END SOLVENT
