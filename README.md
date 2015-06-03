# Dijets
Investigation of dijet asymmetries

This repo is subordinate to my [Analysis](https://github.com/kpedro88/Analysis) repo. It reimplements a few classes from that repo and uses many others directly.
It is necessary to check out both repos in order to compile the macros in this repo.

## Skimming

[input\_selection.txt](input/input\_selection.txt) defines the basic dijet selection and lists the available QCD samples.

To run interactively, applying the "dijet" selection to the "QCD" sample and writing the output tree to a folder "tree_dijet":
```
root -b -q -l 'KSkimDriver.C+("input/input_selection.txt","QCD","dijet"," root://cmseos.fnal.gov//store/user/pedrok/SUSY2015/Phys14_QCD_Pt-binned_PUPPI","tree")'
```

To recompile the driver without running:
```
root -b -q -l 'KSkimDriver.C++()'
```

## Histogram generation

The dijet asymmetry distributions are binned in several variables. This leads to a large number of histograms, so it is best to create them all at once and store them in ROOT files.
The binning for these histograms is defined in [input\_bins.txt](input/input\_bins.txt).
```
root -b -q -l KHistDriver.C+
```

## Fitting

To fit all the histograms saved in the previous step, extrapolate the fit parameters to α = 0, and look at the extrapolated parameter trends vs. p<sub>T</sub>:
```
root -b -l -q 'dijet_comp.C+(1)'
```
This saves all the plots of the fits, extrapolations, and trends in the "plots" folder, and then makes a zip of that folder for easy download.