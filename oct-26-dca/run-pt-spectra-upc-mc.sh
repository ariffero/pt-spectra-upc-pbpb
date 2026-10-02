o2-analysis-ud-sgcand-producer -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json \
| o2-analysis-ft0-corrected-table -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json \
| o2-analysis-pid-tof-merge -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json \
| o2-analysis-tracks-extra-v002-converter -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json \
| o2-analysis-event-selection-service -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json \
| o2-analysis-trackselection -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json \
| o2-analysis-propagationservice -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json \
| o2-analysis-pid-tpc-service -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json \
| o2-analysis-ud-pt-spectra-inclusive-upc -b --shm-segment-size 4000000000 --configuration json://conf-spectra-mc.json --aod-file @file_list_mc.txt

mv AnalysisResults.root AnalysisResults-spectra-mc.root