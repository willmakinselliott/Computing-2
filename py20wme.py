 # Your imports go here
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import peak_widths
import scipy.integrate
from scipy.constants import Planck as h
from scipy.constants import hbar
# Your functions go here
def ProcessData(filename):
    """
    This function opens the file called, splits each line into a list of its components then removes the lines containing items
    other than just data so the data can be analysed.
   
    Firstly, the noise level is determined and then an average of this value is removed from all absorption values.
   
    Secondly, the 20GHz data is located, the peak width and magnetic field at peak absorption are then located using a curve fit with the
    absorption_fit function inside this one. The uncertainty in these values are also calcualted.
   
    The next steps involve data for all frequencies. The peak position and width are found in the same way as for 20 GHz but for all frequencies.
    frequncy vs peak position is then fitted with the kittel_eqn_freq function inside this one. From this the gyromagnetic ratio, lande g factor,
    ainsotropy field and saturation magnetisation. uncertainties for these values are also given.
   
    Finally peak width vs frequncy is fitted to the peak_width_fit function iside of this one. From that fit the intrinsic line width, gilbert damping
    parameter and their uncertainties are calculated.
   
    The function returns a dictionary called results containting peak position at 20GHz in Am^-1, peak width at 20GHz in Am^-1, gyromagentic ratio
    in rads^-1 T^-1, lande g factor, ainsotropy field in Am^-1, saturation magnetisation in Am^-1, intrinsic line width in Am^-1,
    gilbert daming parameter and all associated uncertainties.  
    """
    # Your code goes here
    file = open(filename) #opens the file called and splits it into a list whose elements are a list containing every value/word from each line as a string
    list_of_lines = []
    for line in file:
        line_strip = line.strip()
        line_split = line_strip.split()
        list_of_lines.append(line_split)

    index_END = list_of_lines.index(['&END']) #finds the index of &END in the file(where the metadata ends) so the meta data can be fremoved from the file when data is analysed

    frequencies_raw = list_of_lines[index_END+1][3:] #creates a list of the frewunceies in the data by reading the first line of the data and stripping GHz from the value
    frequencies_string = [i.strip('GHz') for i in frequencies_raw]
    frequencies = []
    for i in range(len(frequencies_string)):
        frequencies.append(float(frequencies_string[i]))
   
    freq_values_in_Hz = [] #this converts the values of frequency from GHz to Hz as in equations value in Hz are used
    for i in range(len(frequencies)):
        freq_values_in_Hz.append(frequencies[i]*10**9)
       
    only_data_list_of_lines = list_of_lines[index_END+2:] #this removes the frequnency row from the data, so it only contains the mag field and absorbtion values

    mag_field_in_data = [] #this takes the first value from each element in only_data_list_of_lines to create a list of mag field values
    for i in range(len(only_data_list_of_lines)):
        mag_field_in_data.append(only_data_list_of_lines[i][0])

    float_mag_field_in_data = [] #this converts all the values of mag field from string to float so the values can be used in plots and calcualtions
    for i in range(len(mag_field_in_data)):
        float_mag_field_in_data.append(float(mag_field_in_data[i]))

    list_of_data_without_mag_field_values = [] #this creates a list of lits for just absorbtion values by removing the first element of each list(the mag field vaules)
    for i in range(len(only_data_list_of_lines)):
        list_of_data_without_mag_field_values.append(only_data_list_of_lines[i][1:])
    transposed_data = np.transpose(list_of_data_without_mag_field_values) #gives list of lists of absobtion values for each freq rather than absorption vlaues for each mag field value

    float_transposed_data = [] #as the abosrbtion values in the lists are all strings they must be converted to floats so they can be used in calculations
    for i in range(len(transposed_data)):
        individual_lists_float_transposed_data = []
        for j in range(len(transposed_data[i])):
            individual_lists_float_transposed_data.append(float(transposed_data[i][j]))
        float_transposed_data.append(individual_lists_float_transposed_data)

    noise_finder = [] #this takes the transposed data and makes its values all floats the same as the for loop above however this is just so that the noise level can be calcualted
    for i in range(len(transposed_data)):
        individual_noise_finder = []
        for j in range(len(transposed_data[i])):
            individual_noise_finder.append(float(transposed_data[i][j]))
        noise_finder.append(individual_noise_finder)

    abs_first_half_frequencies = noise_finder[:int(len(noise_finder)/2)] #this takes the absorbtion values for the first half of the frequencies
    abs_second_half_frequencies = noise_finder[int(len(noise_finder)/2):] #this takes the absorbtion values for the remaining frequnecies

    last_100_values = [] #this creates a list containing the last 100 values of absorbtion for the first half of the frequnecies
    for i in range(len(abs_first_half_frequencies)):
        last_100_values.append(abs_first_half_frequencies[i][-100:])

    first_100_values = [] #this creates a list containing the first 100 values of absorbtion for the second half of the frequencies
    for i in range(len(abs_second_half_frequencies)):
        first_100_values.append(abs_second_half_frequencies[i][:100])

    average_noise_1st_half_freq = [] #creates a list of the average value of the last 100 values for each of the first half of frequnecies
    for i in range(len(last_100_values)):
        average_noise_1st_half_freq.append(np.mean(last_100_values[i]))

    average_noise_2nd_half_freq = [] #creates a list of the average value of the first 100 values for each of the second half of frequnecies
    for i in range(len(first_100_values)):
        average_noise_2nd_half_freq.append(np.mean(first_100_values[i]))

    average_noise = average_noise_1st_half_freq + average_noise_2nd_half_freq #this combines the 2 lists above

    float_transposed_data_minus_average_noise = [] #this removes the average noise level for each frequncies from every absorbtion value at that frequnency
    for i in range(len(float_transposed_data)):
        individual_float_transposed_data_minus_average_noise = []
        for j in range(len(float_transposed_data[i])):
            individual_float_transposed_data_minus_average_noise.append(float_transposed_data[i][j]-average_noise[i])
        float_transposed_data_minus_average_noise.append(individual_float_transposed_data_minus_average_noise)

    absorption_for_20GHz = float_transposed_data_minus_average_noise[frequencies.index(20.0)] #this is a list of the absorbtion values at 20GHz with the average noise level removed
   
    #P0 values for 20 GHz absorption fit(lorentzian function)
    #magnetic field peak predictions, value of mag field at peak absorption(finds highest value of absorbtion and its corresponding value of mag field)
    max_absorption_20GHz = np.amax(absorption_for_20GHz)
    max_absorption_20GHz_index = int(np.where(max_absorption_20GHz == absorption_for_20GHz)[0])
    mag_field_at_peak_20GHz = float_mag_field_in_data[max_absorption_20GHz_index]
    #peak width prediction, fwhm(halves peak mag field value then finds the mag field values closest to the halved value on each side of the peak the difference between)
    half_peak_20GHz = max_absorption_20GHz/2
    first_half_split_20GHz = absorption_for_20GHz[:absorption_for_20GHz.index(max_absorption_20GHz)]
    second_half_split_20GHz = absorption_for_20GHz[absorption_for_20GHz.index(max_absorption_20GHz)+1:]
    _1st_half_index_of_absorption_closest_to_half_peak = absorption_for_20GHz.index(min(first_half_split_20GHz, key=lambda x:abs(x-half_peak_20GHz)))
    _2nd_half_index_of_absorption_closest_to_half_peak = absorption_for_20GHz.index(min(second_half_split_20GHz, key=lambda x:abs(x-half_peak_20GHz)))
    mag_field_at_given_index_first_20GHz = float_mag_field_in_data[_1st_half_index_of_absorption_closest_to_half_peak]
    mag_field_at_given_index_second_20GHz = float_mag_field_in_data[_2nd_half_index_of_absorption_closest_to_half_peak]
    peak_width_20GHz = mag_field_at_given_index_second_20GHz - mag_field_at_given_index_first_20GHz
    #I0 prediction, integrate under curve
    integral_20GHz = scipy.integrate.trapz(absorption_for_20GHz, float_mag_field_in_data)

    def absorption_fit(float_mag_field_in_data, peak_width, I_0, mag_field_at_peak_absorption):
        '''
        This function returns values of absortion from a lorentzian function given by, I = (I0/2*pi)*(delH/((H-H0)**2+(delH/2)**2)
        where:
        I = absorption values in arb units
        I0 = area under curve in a square units
        delH = peak width in Am^-1
        H = magnetic field in Am^-1
        H0 = magnetic field at peak in Am^-1
        '''
       
        first_term = I_0/(2*np.pi)
        second_term = peak_width/((float_mag_field_in_data - mag_field_at_peak_absorption)**2 + (peak_width/2)**2)
        return first_term * second_term
    #curve fitting the absorbtion for 20GHz, values of peak width, I0 and peak mag field contaied within GHz20opt and there uncertainties in GHz20cov
    GHz20opt,GHz20cov = curve_fit(absorption_fit, float_mag_field_in_data, absorption_for_20GHz, p0=(peak_width_20GHz, integral_20GHz, mag_field_at_peak_20GHz), maxfev = 50000)
    absorption_fit_errors = np.sqrt(np.diag(GHz20cov))
   
    #P0 values for all frequencies (lortetzian function)
    #magnetic field peak predictions, value of mag field at peak absorption(finds highest value of absorbtion and its corresponding value of mag field)
    max_values_of_absorption = []
    for i in range(len(float_transposed_data_minus_average_noise)):
        max_values_of_absorption.append(np.amax(float_transposed_data_minus_average_noise[i]))

    index_position_of_max_absorption = []
    for i in range(len(max_values_of_absorption)):
        index_position_of_max_absorption.append(np.where(max_values_of_absorption[i] == float_transposed_data_minus_average_noise))

    refined_index_of_max_absorption = []
    for i in range(len(index_position_of_max_absorption)):
        refined_index_of_max_absorption.append(int(index_position_of_max_absorption[i][1]))

    corresponding_values_of_mag_field = []
    for i in range(len(max_values_of_absorption)):
        corresponding_values_of_mag_field.append(float_mag_field_in_data[refined_index_of_max_absorption[i]])

    #peak width prediction, fwhm(halves peak mag field value then finds the mag field values closest to the halved value on each side of the peak the difference between)
    half_peak = [] #finds half peak mag field value for each frequency
    for i in range(len(max_values_of_absorption)):
        half_peak.append(max_values_of_absorption[i]/2)

    first_half_split = [] #this makes a list of the absorbtion values up to the peak value for each frequency
    for i in range(len(float_transposed_data_minus_average_noise)):
        first_half_split.append(float_transposed_data_minus_average_noise[i][:float_transposed_data_minus_average_noise[i].index(max_values_of_absorption[i])])

    second_half_split = [] #this makes a list of the absorbtion values after the peak value for each frequency
    for i in range(len(float_transposed_data_minus_average_noise)):
        second_half_split.append(float_transposed_data_minus_average_noise[i][float_transposed_data_minus_average_noise[i].index(max_values_of_absorption[i])+1:])

    index_of_closest_value_to_half_peak_first_half = [] #this makes a list of the index of the value closest to the half peak value for the first split of the data
    for i in range(len(float_transposed_data_minus_average_noise)):
        index_of_closest_value_to_half_peak_first_half.append(float_transposed_data_minus_average_noise[i].index(min(first_half_split[i], key=lambda x:abs(x-half_peak[i]))))

    index_of_closest_value_to_half_peak_second_half = [] #this makes a list of the index of the value closest to the half peak value for the second split of the data
    for i in range(len(float_transposed_data_minus_average_noise)):
        index_of_closest_value_to_half_peak_second_half.append(float_transposed_data_minus_average_noise[i].index(min(second_half_split[i], key=lambda x:abs(x-half_peak[i]))))

    mag_fields_at_given_indexes_first = [] #this makes a list of the mag field value at the index of the absorbtion value closest to the half peak value for the first half split
    for i in range(len(index_of_closest_value_to_half_peak_first_half)):
        mag_fields_at_given_indexes_first.append(float_mag_field_in_data[index_of_closest_value_to_half_peak_first_half[i]])

    mag_fields_at_given_indexes_second = [] #this makes a list of the mag field value at the index of the absorbtion value closest to the half peak value for the second half split
    for i in range(len(index_of_closest_value_to_half_peak_second_half)):
        mag_fields_at_given_indexes_second.append(float_mag_field_in_data[index_of_closest_value_to_half_peak_second_half[i]])

    peak_width_prediction = [] #this subtracts the mag field values from the first half split from the second haf split to give a fwhm prediction
    for i in range(len(mag_fields_at_given_indexes_first)):
        peak_width_prediction.append(mag_fields_at_given_indexes_second[i]-mag_fields_at_given_indexes_first[i])

    #I0 prediction, integrate under curve
    integral_absorption_vs_mag_field = []
    for i in range(len(float_transposed_data_minus_average_noise)):
        integral_absorption_vs_mag_field.append(scipy.integrate.trapz(float_transposed_data_minus_average_noise[i], float_mag_field_in_data))
       
    #curve fitting for the absorbtion fit(lorentzian function)
    list_of_peak_positions = []
    peak_position_uncertainties = []
    list_of_widths = []
    width_uncertainties = []
    for i in range(len(float_transposed_data_minus_average_noise)):
        lorentzopt,lorentzcov = curve_fit(absorption_fit, float_mag_field_in_data,float_transposed_data_minus_average_noise[i],p0=(peak_width_prediction[i],integral_absorption_vs_mag_field[i],corresponding_values_of_mag_field[i]), maxfev = 50000)
        list_of_peak_positions.append(lorentzopt[2])
        peak_position_uncertainties.append(np.sqrt(np.diag(lorentzcov))[2])
        list_of_widths.append(lorentzopt[0])
        width_uncertainties.append(np.sqrt(np.diag(lorentzcov))[0])

    def kittel_eqn_freq(list_of_peak_positions, gyromagnetic_ratio, anisotropy_field, saturation_magnetisation):
        '''
        This function returns values of frequency given by the kittel equation, v(H) = (u0*gamma/2*pi)/sqrt((H+HK)+(H+HK+MS))
        where:
        v = frequency in Hz
        u0 = 4 pi * 10 ** -7
        gamma = gyromagnetic ratio in rads^-1 T^-1
        H = magnetic field vlaues at peaks in Am^-1
        HK = anisotropy field in Am^-1
        MS = saturation magnetisation in Am^-1
        '''
        pre_factor = ((4*np.pi)*(10**-7)*gyromagnetic_ratio)/(2*np.pi)
        root_term = (list_of_peak_positions + anisotropy_field)*(list_of_peak_positions + anisotropy_field + saturation_magnetisation)
        kittel_freq = pre_factor * np.sqrt(root_term)    
        return kittel_freq    
   
    #curve fitting for kittel equation, values of gyromangetic ratio, ainsotropy field and saturation magnetisation conatined within kittelopt and there uncertainties in kittelcov
    kittelopt,kittelcov = curve_fit(kittel_eqn_freq, list_of_peak_positions, freq_values_in_Hz, p0=((1.75*9.27400968*10**-24)/(hbar),20000,60000),maxfev = 50000)
    kittel_eqn_errors = np.sqrt(np.diag(kittelcov))

    lande_g_factor = (kittelopt[0]*(hbar))/(9.27400968*10**-24) #calculation of lande g factor from the gyromagnetic ratio
    lande_g_factor_error = (kittel_eqn_errors[0]*hbar)/(9.27400968*10**-24) #error in lande g factor given by error in gyromagnetic ratio error multiplided by the same factors when the value is calcualted
   
    def peak_width_fit(freq_values_in_Hz, intrinsic_width, gilbert_damping):
        '''
        This function returns values of peak width in A/m given by the equation, dehH = delH_0+(alpha*4*pi/u_0*sqrt(3))*(v/gamma)
        where:
        delH = peak width in Am^-1
        delH_0 = intrinsic peak width in Am^-1
        alpha = gilbert damping parameter
        u_0 = 4*pi*10**-7
        v = frequency in Hz
        gamma = gyromagnetic ratio in rads^-1 T^-1
        '''
        constant = (4*np.pi)/(np.sqrt(3)*4*np.pi*(10**-7))
        factor = (gilbert_damping*freq_values_in_Hz)/kittelopt[0]
        return intrinsic_width + constant*factor
   
    #curve fitting for the width fit, values of intrisic peak width and gilbert damping parameter are contained within widthopt and there uncertainties in widthcov
    widthopt,widthcov = curve_fit(peak_width_fit, freq_values_in_Hz, list_of_widths, p0=(3000, 0.05))
    peak_width_errors = np.sqrt(np.diag(widthcov))
    
    plot_values_absorption = absorption_fit(float_mag_field_in_data, GHz20opt[0], GHz20opt[1], GHz20opt[2])
    plottable_values_of_freq = kittel_eqn_freq(list_of_peak_positions, kittelopt[0], kittelopt[1], kittelopt[2])
    plottable_widths = []
    for i in range(len(freq_values_in_Hz)):
        plottable_widths.append(peak_width_fit(freq_values_in_Hz[i], widthopt[0], widthopt[1]))
   
    plt.plot(float_mag_field_in_data, float_transposed_data[72], "r.",label='data')
    plt.plot(float_mag_field_in_data, plot_values_absorption,label='fit')
    plt.xlabel("magnetic field (A/M)")
    plt.ylabel("absorption (arb units)")
    plt.title("py20wme: absorption vs magnetic field for 20GHz")
    plt.annotate("Io = 797000 ± 800",(800000,40))
    plt.annotate("Ho = 623854 ± 6 A/m",(800000,37))
    plt.annotate("delH = 11570 ± 20 A/m",(800000,34))
    plt.annotate('',xy=(0.64*10**6,44),xytext=(0.8*10**6,41),arrowprops=dict(arrowstyle='->'),)
    plt.legend(loc="upper left")
    plt.show()
   
    plt.errorbar(list_of_peak_positions, freq_values_in_Hz, yerr=None, xerr=peak_position_uncertainties, fmt='b.',label='data',elinewidth=1,capsize=3)
    plt.plot(list_of_peak_positions, plottable_values_of_freq,'r', label="fit")
    plt.xlabel("magnetic field (A/M)")
    plt.ylabel("frequency (Hz)")
    plt.title("py20wme: frequency vs peak position")
    plt.annotate("g = 1.655653 + 0.000004",(0.0,30000000000))
    plt.annotate("Hk = 30657 ± 8 A/m",(0.0,27500000000))
    plt.annotate("Ms = 66210 ± 20 A/m",(0.0,25000000000))
    plt.legend(loc="lower right")
    plt.show()
   
    plt.errorbar(freq_values_in_Hz, list_of_widths, yerr=width_uncertainties, xerr=None, fmt='b.',label='data',elinewidth=1,capsize=3)
    plt.plot(freq_values_in_Hz, plottable_widths,'r',label='fit')
    plt.xlabel("frequency (Hz)")
    plt.ylabel("peak width (A/m)")
    plt.title("py20wme: peak width vs frequency")
    plt.annotate("deltaH = 3495 ± 4 A/m",(2500000000,16000))
    plt.annotate("alpha = 0.010193 ± 0.000005 ",(2500000000,15000))
    plt.legend(loc="lower right")
    plt.show()
   
    results={"20GHz_peak":GHz20opt[2], #peak position at 20GHz (A m^-1)
             "20GHz_peak_error":absorption_fit_errors[2], # uncertainity in above (A m^-1)
             "20GHz_width":GHz20opt[0], #Delta H for 20 GHz (A m^-1)
             "20GHz_width_error":absorption_fit_errors[0],  # uncertainity in above (A m^-1)
             "gamma": kittelopt[0], #your gamma value (rad s^-1 T^-1)
             "gamma_error": kittel_eqn_errors[0],  # uncertainity in above (rad s^-1 T^-1)
             "g": lande_g_factor, #Your Lande g factor (dimensionless number)
             "g_error": lande_g_factor_error,  # uncertainity in above (dimensionless number)
             "Hk": kittelopt[1], #Your value for the anisotropy field (A m^-1)
             "Hk_error": kittel_eqn_errors[1],  # uncertainity in above (A m^-1)
             "Ms": kittelopt[2], #Your value for M_s (A m^-1)
             "Ms_error": kittel_eqn_errors[2],  # uncertainity in above (A m^-1)
             "DeltaH": widthopt[0], #Intrinsic line width (A m^-1)
             "DeltaH_error": peak_width_errors[0],  # uncertainity in above (A m^-1)
             "alpha": widthopt[1], #Gilbert Damping parameter (dimensionless number)
             "alpha_error": peak_width_errors[1] } # uncertainity in above (dimensionless number)

    return results
#ProcessData("assessment_data_py20wme.dat")