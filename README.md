# thylacine

> **Status:** developed in the open as a working record. The analyses support the peer-reviewed paper below; later changes are not peer-reviewed unless stated.

<b>R code for analysis of the Tasmanian Thylacine Sighting Records Database (TTSRD)</b>

- Analyses can be driven through the method script, written in R (it requires loading the functions script).
- R code for running all analyses are in /full_analysis, and just the spatial modelling in /spatial_analysis
- The reduced CSV form of the database (loaded in the R scripts), along with GIS layers for the spatial analysis, are in /data
- The revised database (v1.1.1, 1223 records, after some museum-record errors and duplicates were removed) is in /ttsrd_database as CSV files and as `TTSRD_1.1.1_public.xlsx`, which omits the `Observer` field

## Source documents and observer identities

The scanned attachments for each record (original sighting descriptions, letters and agency forms) and the full spreadsheets with observer names (v1.0.0 with all 1237 records, v1.1.0 and v1.1.1) are no longer kept in this repository. Many of the documents carry personal details of the people who reported sightings, such as names, signatures, addresses and phone numbers, so they are held in a restricted-access Zenodo record, [doi:10.5281/zenodo.23876288](https://doi.org/10.5281/zenodo.23876288), and released on request for non-commercial research. They were removed from this repository's history on 2 November 2026; the analysis code and results are unchanged.

<i>Associated peer-reviewed paper:</i>  
Brook, B.W., Sleightholme, S.R., Campbell, C.R., Jarić, I. & Buettel, J.C. (2023). Resolving when (and where) the Thylacine went extinct. <i>Science of The Total Environment</i>, 877, 162878. [https://doi.org/10.1016/j.scitotenv.2023.162878](https://doi.org/10.1016/j.scitotenv.2023.162878)

<i>Archived code and database:</i> release v1.0.0 of this repository is archived at [doi:10.5281/zenodo.23935032](https://doi.org/10.5281/zenodo.23935032). Please cite the paper and this archive; `CITATION.cff` holds both.
