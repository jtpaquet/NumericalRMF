import numpy as np
from scipy.optimize import curve_fit
from scipy.fft import fft, fftfreq
from scipy.signal import butter, filtfilt
from scipy.integrate import simpson

def gaussian(x, A, mu, sigma, H): 
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2)) + H


def sig_figs(m, dm):
    # float: m = measurement
    # float: dm = error
    if np.isfinite(m) and np.isfinite(dm) and dm>0:
        exponent = -int(np.floor(np.log10(abs(dm))))
        dm=np.round(dm,exponent)
        exponent = -int(np.floor(np.log10(abs(dm))))
        m=np.round(m,exponent)
        if exponent <= 0:
            m=int(m)
            dm=int(dm)
    return m, dm


### PLOT PARAMETERS

from cycler import cycler

plot_options = {
    "axes.grid" : False,
    'axes.prop_cycle' : cycler('color', 'krbgcmy'),
    'axes.unicode_minus' : False, # Otherwise, bug with negative number in axes
    'errorbar.capsize' : 3,
    'figure.dpi':100,
    'figure.figsize':(12,9),
    'font.family' : 'serif',
    'font.sans-serif': 'Computer Modern Sans Serif',
    'font.serif': 'cmr10',
    'font.size': 20,
    'image.cmap' : 'jet',
    'lines.color' : 'k',
    'lines.linewidth' : 0.7,
    'lines.markerfacecolor': 'w',
    'lines.markersize': 10,
    'markers.fillstyle': 'none',
    'mathtext.fontset':'cm',
    'scatter.edgecolors' : 'face',
    'scatter.marker' : 'o',
    'text.usetex' : False,
    'xtick.direction' : 'in',
    'xtick.major.size' : 4.0, 
    'xtick.minor.visible' : True,
    'xtick.top': True,
    'ytick.direction' : 'in',
    'ytick.major.size' : 4.0, 
    'ytick.minor.visible' : True,
    'ytick.right': True,
    }



def sine(t, a, w, phi, dc):
    return a*np.sin(w*t+phi)+dc

def fit_sine(t,x,pad=1):
    xf, mag, phase = get_fft(t,x, pad=pad)
    popt, pcov = np.array([0, np.nan, 0, 0]), np.ones((4,4))
    a0  = abs(x.max() - x.min())/2
    a0 = np.sqrt(2)*np.sqrt(np.mean((x-x.mean())**2))
    # a0 = mag.max()
    f0 = xf[mag.argmax()]
    w0 = 2*np.pi*f0
    phi0 = phase[mag.argmax()]
    dc0 = x.mean()
    p0 = (a0, w0, phi0, dc0)
    # Fit
    popt, pcov = curve_fit(sine, t, x, p0=p0)
    a, w, phi, dc = popt
    if a < 0:
        a = abs(a)
        phi = phi - np.pi
    phi = phi % (2*np.pi)
    if phi > np.pi:
        phi = phi - 2*np.pi
    popt = [a, w, phi, dc]
    return popt, pcov


def get_fft(x,y,pad=1):
    y_ = np.tile(y, 1)
    y_ = np.pad(y, (len(y)*pad), 'mean') # Pads mean to have smaller bins on fft
    y_ = y_ - y_.mean()
    N = len(y_)
    dt = (x[2]-x[0])/2
    dt = np.diff(x).mean()
    yf = fft(y_)
    xf = fftfreq(N, dt)
    xf, yf = xf[:N//2], yf[:N//2]
    magnitude = 2 * np.abs(yf) / N
    phase = np.angle(yf, deg=True)
    return xf, magnitude, phase

# Low-pass filter design (Butterworth)
def lowpass(data, cutoff, fs, order=4):
    nyq = 0.5 * fs
    norm_cutoff = cutoff / nyq
    b, a = butter(order, norm_cutoff, btype='low')
    return filtfilt(b, a, data)

def get_amp_phi(t,x,f0=0):
    dt = np.diff(t).mean()
    fs = int(1/dt)
    if f0 == 0:
        xf, mag, phase = get_fft(t,x)
        f0 = xf[mag.argmax()]
        # f0 = find_peak(xf[np.argmax(mag)-20:np.argmax(mag)+20], mag[np.argmax(mag)-20:np.argmax(mag)+20])
        # f0 = zero_crossing(t,x)
        # print(f"f0={f0/1e3:.2f} kHz")

    # Reference signals for demodulation
    ref_cos = np.cos(2 * np.pi * f0 * t)
    ref_sin = np.sin(2 * np.pi * f0 * t)

    # Multiply signal by reference (mixing)
    I_raw = x * ref_cos
    Q_raw = x * ref_sin

    # Apply low-pass filter
    I = lowpass(I_raw, cutoff=f0/5, fs=fs)
    Q = lowpass(Q_raw, cutoff=f0/5, fs=fs)

    # Estimate amplitude and phase
    amplitude_est = 2*np.sqrt(I**2 + Q**2)
    phase_est = np.arctan2(Q, I)

    # Average values
    amplitude_avg = np.mean(amplitude_est)
    phase_avg = np.mean(phase_est)

    # Print results
    # print(f"f0={f0/1e3:.2f} kHz", f"Estimated Amplitude: {amplitude_avg:.3f}", f"Estimated Phase: {np.degrees(phase_avg):.2f} degrees")

    # import matplotlib.pyplot as plt
    # plt.plot(t, x, label='Signal')
    # plt.show()

    # plt.plot(t, I_raw, label='I_raw')
    # plt.plot(t, Q_raw, label='I_raw')
    # plt.show()

    # plt.plot(t, I, label='I_raw')
    # plt.plot(t, Q, label='I_raw')
    # plt.show()

    # plt.plot(t, amplitude_est, label='I_raw')
    # plt.plot(t, phase_est, label='I_raw')
    # plt.show()
 
    return f0, amplitude_avg, phase_avg


def fit_sine_lock_in(t,x):
    popt, pcov = np.array([0, np.nan, 0, 0]), np.ones((4,4))
    f0, a0, phi0 = get_amp_phi(t,x)
    w0 = 2*np.pi*f0
    dc0 = x.mean()
    p0 = (a0, w0, phi0, dc0)
    # Fit
    popt, pcov = curve_fit(sine, t, x, p0=p0)
    a, w, phi, dc = popt
    if a < 0:
        a = abs(a)
        phi = phi - np.pi
    phi = phi % (2*np.pi)
    if phi > np.pi:
        phi = phi - 2*np.pi
    popt = [a, w, phi, dc]
    return popt, pcov


def find_peak(f,x):
    q1 = simpson(y=f*x,x=f)
    q2 = simpson(y=x,x=f)
    return q1/q2


def find_integral(f,x):
    q1 = simpson(y=f,x=f)
    delta_f = f.max() - f.min()
    return q1/delta_f


def fit_sine_rms(t,x):
    xf, mag, phase = get_fft(t,x,pad=1)
    a0 = np.sqrt(2)*np.sqrt(np.mean((x-x.mean())**2))
    dc0 = x.mean()
    f0 = xf[mag.argmax()]
    f0 = find_freq(t,x)
    w0 = 2*np.pi*f0
    phi0 = np.radians(phase[mag.argmax()])
    x_ = x - x.mean()
    x_ = x_/a0
    p0 = (w0, phi0)
    # Fit
    try:
        popt, pcov = curve_fit(lambda t, w, phi: sine(t, 1, w, phi, 0), t, (x-dc0)/a0, p0=(w0, phi0))
    except RuntimeError:
        try:
            popt, pcov = curve_fit(lambda t, w, phi: sine(t, 1, w, phi, 0), t, (x-dc0)/a0, p0=(w0, phi0+np.pi/2))
        except RuntimeError:
            try:
                popt, pcov = curve_fit(lambda t, w, phi: sine(t, 1, w, phi, 0), t, (x-dc0)/a0, p0=(w0, phi0-np.pi/2))
            except RuntimeError:
                print("Optimal parameters not found")
                popt = p0
    w, phi = popt
    try:
        popt, pcov = curve_fit(lambda t, a: sine(t, a, w, phi, 0), t, (x-dc0), p0=(a0))
        a_fit = popt[0]
    except RuntimeError:
        print("Optimal parameters not found")
        a_fit = a0
    popt = [a0, w, phi, dc0]
    pcov = [0,0,0,0]
    return popt, pcov


def zero_crossing(t,x):
    a0 = np.sqrt(2)*np.sqrt(np.mean((x-x.mean())**2))
    Dt = t.max() - t.min()
    zc_x = np.copy(x-x.mean())/a0
    zc_mask = np.abs(zc_x) < 1/10
    zc_x[zc_mask] = 0
    zc_x = np.sign(zc_x)
    zc_x[zc_x == 0] = 1
    f0 = np.abs(np.diff(zc_x)).sum()/(4*Dt)
    if f0 < 1e3:
        import matplotlib.pyplot as plt
        plt.plot(t*1e3, np.copy(x-x.mean())/a0)
        plt.plot(t*1e3, zc_x)
        plt.xlim([0, 10/f0])
        plt.show()
    return f0

def find_freq(t,x):
    f0, a0, phi0 = get_amp_phi(t,x)
    
    dt = np.diff(t).mean()
    fs = int(1/dt)

    freq = []
    amp = []
    phi = []
    for f_ in np.linspace(f0*0.98, f0*1.02, 101):
        _, amplitude_avg, phase_avg = get_amp_phi(t,x,f0=f_)
        freq.append(f_)
        amp.append(amplitude_avg)
        # phi.append(phase_avg)

    f0 = find_peak(np.array(freq), np.array(amp))
    # print(f"f0={f0/1e3:.3f} kHz")

    return f0


def new_sine_fit(t, x):
    f0 = find_freq(t, x)
    # Initial amplitude and phase guess
    amp, phi = rms_amp_and_mean_phase(x)
    amp, phi = np.sqrt(2) * np.std(x - np.mean(x)), 0
    
    # Ensure initial phase corresponds to positive peak
    x0 = x - np.mean(x)
    if np.mean(x0 * np.sin(2*np.pi*f0*t + phi)) < 0:
        phi += np.pi

    phi = (phi + np.pi) % (2*np.pi) - np.pi  # wrap to [-π, π]
    
    p0 = [amp, 2*np.pi*f0, phi, np.mean(x)]
    if (t.max() - t.min())/f0 > 100:
        mask = (t.min() <= t) & (t <= t.min()+100/f0)
    else:
        mask = t.min() <= t
    popt, pcov = curve_fit(sine, t[mask], x[mask], p0=p0)
    amp, w, phi, dc = popt

    # Force positive amplitude
    if amp < 0:
        amp = -amp
        phi += np.pi

    # Final phase wrap
    phi = (phi + np.pi) % (2*np.pi) - np.pi

    return [amp, w, phi, dc]

from scipy.signal import hilbert

def rms_amp_and_mean_phase(x):
    x0 = x - np.mean(x)
    z = hilbert(x0)

    amp = np.sqrt(2) * np.std(x0)
    phi = np.angle(np.mean(z / np.abs(z)))  # circular mean phase

    return amp, phi


def ellipse_from_sine_params(
    bdot_x_amp, bdot_x_phi,
    bdot_z_amp, bdot_z_phi
):
    # amplitudes
    a1 = bdot_x_amp
    a2 = bdot_z_amp

    # phase difference
    dphi = bdot_x_phi - bdot_z_phi

    # precompute
    c = np.cos(dphi)

    # eigenvalues (squared semi-axes)
    term1 = 0.5 * (a1**2 + a2**2)
    term2 = 0.5 * np.sqrt((a1**2 - a2**2)**2 + 4 * a1**2 * a2**2 * c**2)

    lambda1 = term1 + term2  # major^2
    lambda2 = term1 - term2  # minor^2

    # semi-axes
    a = np.sqrt(lambda1)  # semi-major
    b = np.sqrt(lambda2)  # semi-minor

    # eccentricity
    # guard against tiny numerical negatives
    if lambda1 > 0:
        e = np.sqrt(max(0.0, 1.0 - lambda2 / lambda1))
    else:
        e = 0.0

    # tilt angle (radians)
    theta = 0.5 * np.arctan2(
        2 * a1 * a2 * c,
        (a1**2 - a2**2)
    )

    return {
        "lambda1": lambda1,
        "lambda2": lambda2,
        "semi_major": a,
        "semi_minor": b,
        "eccentricity": e,
        "theta_rad": theta,
        "theta_deg": np.degrees(theta),
    }


def ellipse_from_timeseries(x, z):
    x = np.asarray(x)
    z = np.asarray(z)

    # remove mean (important!)
    x = x - np.mean(x)
    z = z - np.mean(z)

    # covariance elements
    Cxx = np.mean(x * x)
    Czz = np.mean(z * z)
    Cxz = np.mean(x * z)

    # covariance matrix
    C = np.array([[Cxx, Cxz],
                  [Cxz, Czz]])

    # eigen decomposition
    vals, vecs = np.linalg.eigh(C)

    # sort (largest first)
    order = np.argsort(vals)[::-1]
    lambda1, lambda2 = vals[order]
    v1 = vecs[:, order[0]]

    # semi-axes
    a = np.sqrt(2 * lambda1)
    b = np.sqrt(2 * lambda2)

    # eccentricity
    e = np.sqrt(max(0.0, 1 - lambda2 / lambda1)) if lambda1 > 0 else 0.0

    # angle of major axis
    theta = np.arctan2(v1[1], v1[0])

    # Phase shift
    rho = Cxz / np.sqrt(Cxx * Czz)   # correlation coefficient
    rho = np.clip(rho, -1, 1)
    phi = np.arccos(rho)

    return {
        "lambda1": lambda1,
        "lambda2": lambda2,
        "semi_major": a,
        "semi_minor": b,
        "eccentricity": e,
        "theta_rad": theta,
        "theta_deg": np.degrees(theta),
        "phi_deg": np.degrees(phi),
    }

def ellipse_fit_from_param(ellipse_params):
    t = np.linspace(0, 2*np.pi, 500)

    a = ellipse_params["semi_major"]
    b = ellipse_params["semi_minor"]
    theta = ellipse_params["theta_rad"]

    xc = 0
    zc = 0

    # ellipse in its principal-axis frame
    x_ell = a * np.cos(t)
    z_ell = b * np.sin(t)

    # rotate
    R = np.array([[np.cos(theta), -np.sin(theta)],
                [np.sin(theta),  np.cos(theta)]])

    xy_rot = R @ np.vstack([x_ell, z_ell])

    # translate
    x_fit = xy_rot[0] + xc
    z_fit = xy_rot[1] + zc

    return x_fit, z_fit