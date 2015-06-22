
# MASSSPECHISTORY

This is a suite of scripts that will generate a visual QC analysis of a Thermo mass-spectrometer, thus providing a visual history of the mass-spectrometer.

The scripts perform:

1. Batch processing of mzML/raw files for PSM analysis using the Morpheus search engine.
2. Excalibur script to print report for iRT peptides
3. Parsers to extract diagnostic variables from Morpheus search results and Excalibur reports.
4. Construction of a visual webapp to display the data. This can be put in a static directory in a web-server for public consumption, or viewed directly on a local machine.

Morpheus was chosen because of its ease of installation across platforms, it's open-source license and the robustness of its output files.


## Morpheus wrapper

## Windows batch files

## MSMS parsing

simple spectra, PSM peptide, protein counts

## top peptide analysis

## Excalibur diagnostic iRT peptide analysis


## Make charts

Javascript/HTML

nvd3 graphing library
converting summary.txt/reports.txt into data
loading into HTML

## automation
machine copy
windows server because of RAW
linux web/server
shared network drive
can also use winscp/ftp



