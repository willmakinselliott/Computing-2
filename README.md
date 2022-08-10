# Computing-2
Tasks for coursework/overview of what code does
1. Read in the data file and produce a plot of absorption versus magnetic field for a frequency of 20 GHz, with appropriate labelling of the axes. You will need to annotate this plot and then include it in your report. If your code handles more than one frequency, you should only plot the data for 20 GHz.
2. Find the location in magnetic field of the peak absorption for each frequency in your data file. You should annotate your plot to show the located peak for the 20 GHz data.
3. Determine the width of the peak you found in the previous task where width is defined in the Lorentzian function. Annotate your plot of the 20 GHz data with this value too.
Both the peak position and peak width should be quoted with uncertainities.
The remaining tasks are only possible if you select data with multiple frequencies.
4. When multiple frequencies are selected in the data,
  A. Find the peak position and width for each frequency
  B. Make a plot (with error bars) of frequency versus peak position
  C. Fit the frequency vs peak position with the Kittel equation and obtain values of
    a. the gyromagnetic ratio, 
    b. hence Lande g-factor,
    c. the anisotropy field,
    d. saturation Magnetisation
  D. Make a plot, with error bars of peak width, versus frequency and fit this with a straight line according to equation at the end of secgtion 2.3 to obtain a value for the Gilbert Damping constant .
5. Your code must return the results it has found, along with an estimate of the standard error for each, in a dictionary. The details of the keys and what the corresponding values should be are given in the template python file (see section 4).
6. Write a report on your code. Your report should have two sections:
A. A results section that tabulates all the parameters found by your code from your data set, with
errors. You may manually adjust the precision you quote your answers to, so that errors are quoted to 1 s.f. and the result is quoted to a precision that is consistent with the error (as you would do for labs). This section should also include all plots produced by your code.
B. A flow chart that describes how your code works. It should describe the functions you have written and the overall logic of your program. It is not necessary to document every single line,
but the flow chart should be sufficient to allow someone else to write a program that would work
the same way as your code.
