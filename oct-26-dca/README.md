## Input files
- MC = LHC25d1/kDpmjetSingleA, run 544123, (/alice/sim/2025/LHC25d1/kDpmjetSingleA/0/544123/AOD/*/AO2D.root)
- data = Derived data from train 291050, run 544123, (/alice/cern.ch/user/a/alihyperloop/jobs/0088/hy_884121/AOD/*/AO2D.root)

## Configurations
- MC = `conf-spectra-mc.json`. It uses the singleGapProducer as service task, not 100% sure that the config is correct, but it runs fine.
- Data = `conf-spectra-data.json`. Config related only to my task, derived data taken from train reported on the ([Analysis Note](https://alice-notes.web.cern.ch/system/files/notes/analysis/1626/2025-03-15-AN_spectra_UPC-4.pdf))

## Fit to DCA spectra
Done using the macro `fitDCASpectra.c`. It reads TH2 that contains DCAxy vs pT, differently for TPC and TOF PID.