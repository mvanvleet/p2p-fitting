<Date>
Mon Sep  1 11:24:30 CDT 2014

<Purpose>

<Relevant Notebook Page(s)>
MVV-I, p. 51-52

<Method>
    
<Directory Contents>
main_dispersion:
    Code to compute dispersion coefficients for a molecular system (computed
    previously using camcasp) which can employ constraints on the frequency
    dependent polarizabilities used for certain atomtypes. Original code
    written by jmcdaniel. 
    Source: jmcdaniel; source files located in source_dispersion or in
        /home/jmcdaniel/apps/camcasp-5.5/fit_p2p/source 
    Usage:
        $ main_dispersion ifile path_relative_to/home/mvanvleet/  dispersion_constraint_file
    An example constraint file is provided in this directory. Each hard
    constraint is listed by giving its frequency dependent polarizabilities.
main_drude:
    Similar code to as above, but instead computes drude oscillator
    coefficients for a molecular system.
    Source: jmcdaniel; source files located in source_drude or in
        /home/jmcdaniel/apps/camcasp-5.5/fit_p2p/source_fit_static_drude
    Usage:
        $ main_drude ifile path_relative_to/home/mvanvleet/  drude_constraint_file
    An example constraint file is provided in this directory. Unlike the
    dispersion constraint file, the drude constraint file contains both hard
    constraints and initial guesses. A guess (or hard constraint) is required for
    each atom type. Guesses are marked by positive polarizability values, while
    hard constraints are denoted by negative polarizability values.

<System Requirements>
Various custom basis sets can be found in ~/basis_sets, and are necessary to
run all of these jobs.

<Cautions/Warnings/Bugs>

<Results/Conclusions>

