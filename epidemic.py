
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population, N.
N = 250
# Initial number of infected and recovered individuals, I0 and R0.
B10,B20,B30, I0 = 5, 0,0,0
# Everyone else, S0, is susceptible to infection initially.
D0 = N - B10 - I0 - B20-B30
# Contact rate, alfa, and mean recovery rate, beta, (in 1/days).
alfa, beta = 0.7, 0.3
# A grid of time points (in days)
t = np.linspace(0, 100, 100)

# The SIR model differential equations.
def deriv(y, t, N, alfa, beta):
    D, B1,B2,B3, I = y
    dDdt = -alfa * D * B1 / N
    dB1dt = alfa * D * B1 / N - beta * B2
    dB2dt = beta*B1 - beta*B2
    dB3dt = beta*B2 - beta*B3
    dIdt = beta * B3
    return dDdt, dB1dt,dB2dt,dB3dt, dIdt

# Initial conditions vector
y0 = D0, B10,B20,B30, I0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, alfa, beta))
D, B1,B2,B3, I = ret.T
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, D, 'b', alpha=0.5, lw=2, label='Dovzetni')
ax.plot(t, B1, 'r', alpha=0.5, lw=2, label='Inkubacijska doba')
ax.plot(t, B2, 'orange', alpha=0.5, lw=2, label='Kužnost')
ax.plot(t, B3, 'y', alpha=0.5, lw=2, label='Karantena')
ax.plot(t, I, 'g', alpha=0.5, lw=2, label='Imuni')
ax.set_xlabel('Število dni')
ax.set_ylabel('Število ljudi')
ax.set_ylim(0,250)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=1, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(True)
plt.show()
    