import matplotlib.pyplot as plt
import numpy as np
import pickle

from matplotlib.widgets import Button, Slider
from src.config.setup import *
from src.navigation.data_processing import NavDataArrayIndices as IDX, generate_trial_curve


USE_TEST_DATA = False


with open(DATA_PATH + SAVED_DATA_FILE, "rb") as file:
    saved_nav_data = pickle.load(file)
if USE_TEST_DATA:
    with open(DATA_PATH + TEST_DATA_FILE, "rb") as file:
        test_nav_data = pickle.load(file)
else:
    test_nav_data = generate_trial_curve(saved_nav_data)

# The parametrized function to be plotted
def f(freq, amplitude, frequency):
    time = saved_nav_data[:, IDX.t]
    return freq + amplitude + (time - np.min(time)) * frequency

t = saved_nav_data[:, IDX.f]

# Define initial parameters
init_amplitude = 0
init_frequency = 0

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = ax.plot(saved_nav_data[:, IDX.t], f(t, init_amplitude, init_frequency), lw=0, marker=".", label="Measured")
line0, = ax.plot(test_nav_data[:, IDX.t], test_nav_data[:, IDX.f], lw=0, marker=".", label="Simulated")
ax.set_xlabel('Time [s]')

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

plt.legend()

# Make a horizontal slider to control the frequency.
axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
freq_slider = Slider(
    ax=axfreq,
    label='Drift [Hz/s]',
    valmin=-5,
    valmax=5,
    valinit=init_frequency,
)

# Make a vertically oriented slider to control the amplitude
axamp = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
amp_slider = Slider(
    ax=axamp,
    label="Offset [Hz]",
    valmin=-20e3,
    valmax=20e3,
    valinit=init_amplitude,
    orientation="vertical"
)


# The function to be called anytime a slider's value changes
def update(val):
    line.set_ydata(f(t, amp_slider.val, freq_slider.val))
    fig.canvas.draw_idle()


# register the update function with each slider
freq_slider.on_changed(update)
amp_slider.on_changed(update)

# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.
resetax = fig.add_axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')


def reset(event):
    freq_slider.reset()
    amp_slider.reset()
button.on_clicked(reset)

plt.show()