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
./readResponse [resp] [bgd] [src] [incidentCountRate] [option:lightcurve.qdp] [option:overFlowRate]


## Install

## Licence
You can use this software under the MIC license.

## Author
- Tetsuya Yasuda ( yasuda AT heal.phy.saitama-u.ac.jp )
- Saitama University 
