This file gives a basic info on how to proceed with the simulations for cochlear implant development.

CMCD, ClassD and ClassDE amplifier LTSpice files are for testing purposes to provide with an estimate of the efficiency and performance.

The LTSpice files corresponding to ClassE amplifier consist of different variations.
Without coils, With a resonant driver, With MOS switch, with coils but without rectifier and a complete ClassE amplifier using GaNFET (ClassE_v2.asc/v3.asc).

ClassE amplifier has been thoroughly tested. MATLAB files are written to theoretically find the appropriate L,C values for a given ClassE design (two matlab files: i) without rectifier and ii) with rectifier)

************************
ClassE_v3.asc :
--> This LTSPice file imports a txt file named "ClassE_Rect_Param" which includes all the parameters in the circuit. This allows for a complete control over the design by editing the parameters file.
--> This parameters file can in turn be written using the MATLAB file "ClassE_rect.m"
--> The MATLAB file "Feedback_analysis.m" can be used to extract the efficiencies using the corrected equations of M. Baker et. al, which can be used to tweak the passive components values in the Class E amplifier design.

