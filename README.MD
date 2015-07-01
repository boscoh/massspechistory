
# MASSSPECHISTORY

This is a suite of scripts that will generate a visual QC analysis of a Thermo mass-spectrometer, thus providing a visual history of the mass-spectrometer's performance.

The scripts perform:

1. Batch analysis of mzML/raw files for peptide-spectrum matches using the Morpheus search engine.
2. Excalibur script to generate reports for MS peak detecton for iRT peptides
3. Parsers to extract diagnostic variables from Morpheus search results and Excalibur reports.
4. Construction of a javascript/HTML5 webapp to display the data. 
5. Automate the process and report on a web-server and email.


## E. Coli and Hela digest

We use E. Coli and Hela digest to evaluate the performance of the mass-spectrometer. We have chosen Morpheus as it is easily installable across both Windows/Mac/Linux platforms. Morpheus is open-source and the output files it provides a convenient mode that works directly with Thermo RAW files.

We have written a Python wrapper `morpheus.py` around Morpheus. On Windows, Morpheus works with the default.Net module. On Mac and Linux, Morpheus requires the RAW files to be converted to mzML files, and Morpheus requires the installation of Mono.

`morpheus.py` provides some documentation for the command-line version of Morpheus, as well as a clean params dictionary api for the function. It adds some extra functionality such as creating the output directory, and switching between versions depending on the extension of the input file. `morpheus.py` also includes a simple batch processing function.

## Analysis of MS/MS searches of digests

### Default MSMS parameters

Morpheus files are quite simple to parse. They are essentially .tsv files that is read as a list of dictionaries in Python. As well are interested in the overall performance of a digest analysis, we only make use of the `summary.txt` file that summarizes the number of PSMs, peptides, and protein groups identified. Further analysis of individual PSM will be deferred to the iRT peptides below.

### Common peptides
Even though the PSM fluctuates, much of the PSM is due to a lot of poor quality matches that is ultimately not particularly informative. A better way of evaluating the peptides identified in the PSM is to consider only commonly identified peptides. 

To identify the 300 most commonly identified peptides, we make use of the  many samples of the Hela and E. Coli digests that we have run in the past. For each peptide, we can define an acceptable range for fraction of ion products identified, and intensity scores. 

We can then re-evaluate each MS/MS search of the digest in terms of how they recapitulate the 300 most common peptides. 

We believe this is the best indication of the core peptide identification capacity of the instrument.

### Mass calibration

### Percentage of unique peptides to PSMs

## Excalibur diagnostic iRT peptide analysis
For detailed analysis of MS1 spectra, we used commercial iRT peptides to analyze retention times and peak intensities. We created an excalibur script that matched the MS1 spectra of peaks associated with the iRT peptides and generated a report. This is read into the reporting system.

## Charts display
For our analysis, RAW files and iRT peptide reports were saved at every time point. These were read into our reporting system as a function of time.

We generated the charts as interactive javascript/HTML5 graphs on a web-page. This is saved both internally and copied to a public website. The web-app is mobile friendly and renders well on mobile devices. We have used the nvd3, canvas.js and google charts library.

## Reporting Outliers
In order to streamline the analysis process, we wrote a script to identify outliers with the intention of sending email alerts. The email alerts allow the facility manager to be selectively informed of quality failure. 

To analyze the data we ran a two pass averaging. The first pass finds the average and standard deviation of the past data set. Only values that were above the 95% one-tail limit were kept for the second round. The average and standard deviation of this second set was used to define the acceptable lower limit.

An outlier was flagged if a particular value was below this lower limit.

For the iRT peptides, we have loosened the requirement to flag an outlier only if more than 3 irT peptide retention times are outliers.

## Automation

The installation on our systems exhibits a mix of Windows and Linux servers, connected with a Samba network drive. We expect other systems to be similarly arranged.

We have included sample scripts, more for convenience than anything, to demonstrate the automation proces:

- in Windows, using robocopy.exe utility to copy files internally
- batch files with the Windows Task Scheduler
- batch file to batch process RAW files with Morpheus
- batch files to mount and unmount Samba network drives
- using WinSCP to transfer files
- crontabs for Linux
- shell scripts to copy files

