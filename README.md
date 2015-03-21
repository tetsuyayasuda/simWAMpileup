# simWAMpileup

## Overview
simWAMpileup is a Monte-Carlo simulation code based on C++, which is reproduce the on-board analog electronics of Suzaku/WAM to correct pile-up effect.
To read and write the Flexible Imaging Transport System (FITS) data format files, the SLLIB and the SFITSIO libraries [http://www.ir.isas.jaxa.jp/~cyamauch/sli/index.html] are utilized.
Inputs of this software are a background spectrum (FITS format), a light curve (ASCII format), a response matrix (FITS format), and a model spectrum (FITS format), and output is a model spectrum affected by the pile-up effect (FITS format).

## Requirement
- C++ compiler
- Cmake
- SLI Libraries/SFITSIO [http://www.ir.isas.jaxa.jp/~cyamauch/sli/index.html]
- ROOT [https://root.cern.ch/drupal/]

## Test Environment
- CentOS 6.6
- ROOT ver5.99/05
- SLLIB 1.4.2
- SFITSIO 1.4.2

## Usage
### Compile
$ ls  
 source/  
$ mkdir build run  
$ cd build  
$ cmake ../source  
$ make  
$ cd ../run  
$ ln -s ../build/simWAMpileup  


### Preparation  
You should prepare a source spectrum on the FITS format.

[example]

$ xspec  
xspec> model powerlaw  
xspec> fakeit none  
(then you should choice [Use counting statistics] no)  

### Run  
./readResponse [resp:FITS] [bgd:FITS] [src:FITS] [incidentCountRate:int] [lightcurve.qdp:ASCII] [overFlowRate:int]

[inputs]
- resp: a FITS file of response matrix 
- bgd:  a FITS file of background spectrum
- src:  a FITS file of source spectrum before the pileup effect affected
- incidentCountRate: int value of total counts of the source plus the background components in the WAM entire energy band for 1 second
- lightcurve: a ASCII file of space-separated table, 2 columns (time and counts) and n rows.
- overflowRate: int value of over flow. default value is 700.

[output]
- Dead.fits : estimated spectrum affected by the pileup effect


## Licence
You can use this software under the MIC license.

## Author
- Tetsuya Yasuda ( yasuda AT heal.phy.saitama-u.ac.jp )
- Saitama University 
