# Laser Analysis software

This software is to analyse the waveform output of the laser. It plots the gain as a function of x and y, and the minumum voltage as a function of x and y.

## Getting Started

Clone the repo and set the paths for the waveform and TCTAnalyse libraries.

### Prerequisites

TCTAnalyse Version 2.0 available here http://www.particulars.si/TCTAnalyse/TCTAnalyse-Downloads.html

ROOT available here https://root.cern.ch/downloading-root

### Installing

Unzip the files and execute the type of analysis you would like

```
root plot.C
```

## Deployment

For the analysis files:

Set the path in line 3 to point at your PSTCT .sl or .dll file

Set the path in line 5 to point at your waveform.rtct file


