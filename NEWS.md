
#Todo:
- Make sure that a related unknown is positioned last and change warning message to
 "Since fst>0, we highly recommend that only the last unknown is specified as a related".
- A warning message is given when specifying relationships when fst>0:	
- Include a separate function getUpperLR for providing the theoretical maximum LR obtainable.

NEEDED: Print defined relationships for the unknowns into the report

#v1.1.0: (22.11.07)
- Added possibility to export DC profile to Deconvolution panel.
- Added progress-bar with a expected upper time limit for MLE optimizing (similar as EuroForMix).
- GUI windows provided to show different results (param, logLik, LR per marker).
- Model validation is no longer excecuted when creating report (instead it is stored after each user calculation).
- Now only "normal points" are shown in the model validation figure.
- Including "Related reference" to number of typed alleles (may affect testing).

- Bug fixes (thanks to Sandor Furedi for finding these issues):
- Fixed bug when creating report: getReportText-L77: Wrong environment variable "Freqfile" used (should be "optFreqfile"). 
- Fixed bug causing Hp to be also be run for deconvolution when no references are conditioned on.

#v1.0.1: (21.06.22)
 - Fixed bug when creating report: updating getReportText file.
 - Included test "test_logLik1pRepsExtendedError": Expect an error only if a sample contains more markers than specified in a kit (must exist in frequency data).
 - The number of threads is set to zero by default (meaning using all available).
 - OPENMP is optional: Makevars.win and *.cpp scripts were modified.
 - Added license file: GNU Lesser General Public License v3.0

#v1.0.0:
 - Release version (used in publication)
 - Fixed GUI issue when number of related + references exceeds the number of contibutors.

#v0.7:
 - Fixed GUI issue when number of condtional references exceeds the number of contibutors.
 - All references forwarded to model specification are counted as typed individuals (same as from EFM v3.4.x).

#v0.6:
 - Changing package name from EFMreps to EFMrep
 - Updated calcloglik_cumprob:L242 to change lower integration term ('AT-1' instead of 'AT') , fixes issue of observation on AT .
 - Fixed bug in efm2:L698 (getRepGrpID helpfunction); the evidence name when extracted from GUItable was splitted if it contains white space (this is now avoided).
 - Fixed bug in contLikMLE2:L301, where FW stutter param was not properly proposed
 - In fitgammamodel2: Make pre-fitting more robust for highly degraded samples (zero-obs exchanged with smallest obs/2).
 - Found bug in prepareC2:L246 causing a negative stutter from-index for marker dropout.
 - prepareData2 function now throws an error if one of the frequencies become negative (may happen for the Q-allele).
 - Fixed small bug in GUI when configurating evidence specific models: Store only global settings to file (Global settings after Evid settings).
 
#v0.5:
 - Include Rmath instead of using BH::Boost

#v0.4: Implements relatedness
 - Including tests
 - Create vertical scroll for replicates at data import (issue for screens with low resolution).

#v0.3: Small changes:
 - Fixed bug when showing no references.
 - Fixed bug when showing PCA and all replicates is missing the marker.

#v0.2: User-improvements:
- Global settings and selected population freq file is stored and remembered.
- All calculations/results are stored throughout the session (before using REFRESH) and can be highligheted at any time (useful for model comparison).

#v0.1: Initiated development version of EFMreps