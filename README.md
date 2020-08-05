Requirements
============

  - IDL version >6.4(?) under UNIX (Windows not tested)
  - IDL Astronomy User's Library -- http://idlastro.gsfc.nasa.gov
  - Coyote Library by D. Fanning -- http://www.dfanning.com
  - MPFIT package by C. B. Markwardt -- http://purl.com/net/mpfit

  - time has to be in UT in the raw files otherwise include the offset in the sorting routine
  - the sampling rate HAS TO BE 1 second!



Installation
============

  - choose a main directory, e.g., $HOME$/Lightmeter/ (this example is used in the following)
  - copy the routines (*.pro) to $HOME$/Lightmeter/
  - copy the folder 'fix' to $HOME$/Lightmeter/
  - create folder 'data_original/' in $HOME$/Lightmeter/
  - create folder 'data_original/LOCATION/', the location is the name of the observing site, e.g., Paranal

  - open 'SETUP.txt' in $HOME$/Lightmeter/fix/ and modify the main path with, e.g., $HOME$/Lightmeter/

  - open 'LOCATION.txt' in $HOME$/Lightmeter/fix/ and enter the complete location properties using the given format - DO NOT CHANGE THE FORMAT
    the location has to be identical with the name given to the folder in 'data_original/LOCATION/'!


Running the Pipeline
====================

  - copy the raw data to $HOME$/Lightmeter/data_original/LOCATION/

  - open a console a got the the main directory, e.g., $HOME$/Lightmeter/
  - compile and run 'lightmeter'
  - choose one of the offered sites by typing the associated number
  - type 's' for sorting the raw data, this can take a while depending on the number of nights
  - call 'lightmeter' again and type 'r' for the reduction process and wait...
  - seveal folders cotaining different output will created



Modification History
====================

version 2.0 -- under construction
  changed radius of moon slightly in lunar_model.pro (Astronomical Almanac 2010, K7)
  changed AU in selen_coord_sun.pro (Astronomical Almanac 2010, K6)
  taking dynamical time offset into account for ephemeris computation (+66s, in SETUP file)
  possible automatic detection of time offsets - not activated in this version
  new computation of sun rise and offset times and twilight times
  new requirements for reduction, i.e., data have to be available for all times, not only during night
  taking temperature dependence and non-linearity of Lightmeter into account
  bug fix in lunar_model.pro (disk reflectance)
  sun_mod.pro removed, using only sunpos/sunpos_mod.pro
  solar and lunar model used to calibrate Lightmeter, i.e. getting physical units out of counts
  computation of LOSSAM-like data for day time
  producing of a global result file per astronomical night
  global file including parameter limits for irradiance fit
  using global file for starting values for the fit
  ...and more...more or less a rewrite of everything


version 1.7
  using MPFIT package as fit routine

version 1.6
  lunar_model: computation of moon altitude using eq2hor.pro which includes correction for refraction, aberration, nutation, precession...
  rename of axis labels
  locations included in the file names of the plots 

version 1.5
  handles non-linearity of lightmeter by considering only certain range of counts (32k-100k)
  NO correction of counts above the highest threshold is applied, i.e. its only for the purpose to derive a more or less trustful estimate of the conversion factor
  first and last 2 hours of night removed in order to get rid of zodiacal light (which is not visible during the whole year with same intensity)

version 1
  user friendly, fully automatized
  uses global setup files
  possibility to investigate non-linearity of lightmeter by binned model fitting

beta_version
  first version without setup files and super user friendly behavior
