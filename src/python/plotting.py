import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fit_log_line(x, y):
    log_y = np.log(y)
    coeffs = np.polyfit(x, log_y, 1, full=False)  # Ensure it's an array
    return coeffs.tolist()  # Convert to list to prevent unpacking issues


def plotAmp(filtered_data, plot=True, shear=False, transverse_axis="y"):
    row = filtered_data.iloc[0]  # Fix: Get first row correctly
    
    # Convert from strings (if necessary) to lists
    def ensure_array(val):
        return np.array(eval(val)) if isinstance(val, str) else np.array(val)
    
    x_parra = ensure_array(row["initial_distance_from_oscillation_output_x_fft"])
    y = ensure_array(row["amplitude_vector_x"])
    
    coeffs = fit_log_line(x_parra, y)
    yIntercept_amp_x, slope_amp_x = coeffs
    y_x = y

    if shear:
        x_parra = ensure_array(row["initial_distance_from_oscillation_output_y_fft"])
        y = ensure_array(row["amplitude_vector_y"])
        coeffs = fit_log_line(x_parra, y)
        yIntercept_amp_x, slope_amp_x = coeffs
        y_x = y

    if transverse_axis == "y":
        x_perp = ensure_array(row["initial_distance_from_oscillation_output_y_fft"])
        y = ensure_array(row["amplitude_vector_y"])
    else:
        x_perp = ensure_array(row["initial_distance_from_oscillation_output_z_fft"])
        y = ensure_array(row["amplitude_vector_z"])

    coeffs = fit_log_line(x_perp, y)
    yIntercept_amp_y, slope_amp_y = coeffs

    if shear:
        x_perp = ensure_array(row["initial_distance_from_oscillation_output_x_fft"])
        y = ensure_array(row["amplitude_vector_x"])
        coeffs = fit_log_line(x_perp, y)
        yIntercept_amp_y, slope_amp_y = coeffs

    if plot:
        plt.figure()
        plt.scatter(x_parra, y, label=r"$A_{||}(x)$")
        plt.scatter(x_perp, y, label=r"$A_{\perp}(x)$")

        plt.plot(x_parra, np.exp(yIntercept_amp_x + slope_amp_x * x_parra), linestyle="--", color="blue")
        plt.plot(x_perp, np.exp(yIntercept_amp_y + slope_amp_y * x_perp), linestyle="--", color="red")

        plt.yscale("log")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$A(x)$")
        plt.legend()
        plt.grid(True)
        plt.show()

    return np.exp(yIntercept_amp_y) / np.exp(yIntercept_amp_x)
