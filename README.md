# eDITH_RiverDNA

R code supporting "Modelling eDNA transport in river networks reveals highly resolved spatio-temporal patterns of freshwater biodiversity", by Luca Carraro, Rosetta C. Blackman, Florian Altermatt

## Content

- `RUN_MODEL.r`: extract river network, process landscape data, read eDNA data, and run eDITH model for each genus and season. Write output in the `calibration` folder.
- `MAIN.r`: read output from the `calibration` folder, analyze output and produce figures.
- `calibration`: folder containing eDNA raw data and output from the eDITH model.
- `support`: folder containing catchment data and support R functions.
