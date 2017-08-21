****************************************************************************************************
README.txt
Creation Date: 2017-08-21
Authors: Georges Al Makdessi, John Kildea, Robert Maglieri, Logan Montgomery

This file explains:
1) The file structure of the project, including inputs & outputs
2) Instructions for executing the program.

Please refer to dependencies.txt for dependencies
****************************************************************************************************

----------------------------------------------------------------------------------------------------
1) File structure
----------------------------------------------------------------------------------------------------
- Progam code: main.cpp

- Settings file: mlem.cfg
    - Settings file for the algorithm and NNS specs:
        mlem_cutoff                 max # of iterations of MLEM algorithm

        nns_normalization           Vendor-provided normalization factor for the NNS used (varies
                                    with NNS)

        mlem_max_error              The maximum allowed error in reconstructed measurement values
                                    from the true measurment values (in nC). Calculated for each 
                                    MLEM iteration for each measured value (8). Once all values are 
                                    within this max error, MLEM algorithm will stop.

        f_factor                    Vendor-provided calibration factor to convert measured current
                                    (i.e. measured charge divided by duration) to counts per second
                                    (cps). Value in units of [fA / cps]. Value was verified 
                                    experimentally.

- Files read by the program are included in input/
    - The contents of these files should remain constant for most situations
        energy_bins.csv             Energy bins corresponding to other input files [MeV].

        icrp_conversions.csv        Factors to convert neutron flux [n cm^-2 s^-1] to neutron
                                    ambient dose equivalent rate [pSv s^-1]. Factors in units of
                                    [pSv cm^2] for each bin in energy_bins.csv. Values were linearly
                                    interpolated from tabulated values in ICRP 74 Table A.42 (1st
                                    column).

        input_spectrum.csv          The "guess" spectrum input into the MLEM unfolding algorithm.
                                    Flux values in units of [n cm^-2 s^-1] for each bin in
                                    energy_bins.csv. Currently use a step function (high value at
                                    low energies and low value at high energies, step at thermal).

        nns_response.csv            Vendor-provided NNS relative response as a function of energy
                                    for each level of moderation. Line 1: response for bare probe,
                                    Line 2: response for 1 moderator, etc. Units [cm^2].

- Files written to by the program are located in output/
        output_spectrum.csv         The unfolded spectrum generated from executing this program is
                                    appended to this file. Energy bins are also included. The new
                                    spectrum is appended as a column (i.e. a new element is added
                                    to each row). Format this way for direct import into spreadsheet
                                    program

        output_dose.csv             The measured total neutron ambient dose equivalent rate and its
                                    uncertainty are appended to this file. The new values are
                                    inserted as a new row.

        REPORT_*.txt                Report file generated each time the program is executed. The
                                    report contains all pertinent inputs & outputs of the program.

        FIGURE_*.png                Figure file generated each time the program is executed. The
                                    figure displays the spectrum and its uncertainty plotted as a
                                    function of energy.

----------------------------------------------------------------------------------------------------
2) Execution instructions
----------------------------------------------------------------------------------------------------
- Input measured data:
    - Create file input/measurements.csv
    - Refer to input/template_measurements.csv for structure of the file
        - Line 1: A string representing the irradiation conditions
        - Line 2: The measurement duration in seconds
        - Lines 3 through 10: Measured charge values in nC (decreasing number of moderators)

- Adjust settings (if desired/necessary):
    - edit mlem.cfg or input/*.csv

- Compile the program:
    $ make

- Execute the program:
    $ ./main.exe