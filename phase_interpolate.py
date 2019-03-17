import numpy as np
import matplotlib.pyplot as plt
from math import pi as PI
from math import e
from inverter_model import*

STEPS = 16
OVERSAMPLE = 16
freq = 1.0

UIS = 16
SAMPLES_PER_UI = 500
n_samples = UIS*SAMPLES_PER_UI
time = np.arange(n_samples)/float(SAMPLES_PER_UI)

VMIN= 0
VMAX = 1.8

def make_exp_clock(delay, t_ui, uis, samples, _min, _max, initial_sol_uis=10):
    n_per_ui = samples/float(uis)
    t_step = t_ui/n_per_ui
    tau = t_ui/5.0
    wf = np.zeros(samples + int(initial_sol_uis*n_per_ui))
    last_index = 0
    vi = _min
    vf = _max

    while True:
        for n in range(int(n_per_ui)):
            if n + last_index < samples + int(initial_sol_uis*n_per_ui):
                wf[last_index+n] = vi + (vf-vi)*(1-e**(-(n*t_step)/tau))
            else:
                break
        if n + last_index >= samples + int(initial_sol_uis*n_per_ui):
            break
        vi = wf[last_index+n]
        vf = _max if vf < (_min+_max)/2.0 else _min
        last_index += int(n_per_ui)
    wf_new = np.zeros(samples)
    wf_new = wf[int(initial_sol_uis*n_per_ui):]
    wf = wf_new
    if delay:
        start_index = int(n_per_ui*delay)
        wf_new = np.zeros(samples)
        wf_new[start_index:] = wf[:-start_index]
        wf_new[:start_index] = wf[samples-start_index:]
        wf = wf_new
    return wf

nfet_params = {
    "W" : 1.0,
    "L" : 1.0,
    "K" : 1.0,
    "VTO" : 0.5,
    "LAMBDA" : 0.1,
}

pfet_params = {
    "W" : 1.0,
    "L" : 1.0,
    "K" : 1.0,
    "VTO" : -0.5,
    "LAMBDA" : 0.1,
}

def interpolate(I, Q, rI, rQ, vo_i, C, t_step, vdd, vss):
    vo = vo_i
    i = 0.0
    interpolated = np.zeros(len(I))
    for n, sample in enumerate(I):
        vo += (i/C)*t_step
        i_I = inverter_iv(I[n], vo, vdd, vss, nfet_params, pfet_params)
        i_Q = inverter_iv(Q[n], vo, vdd, vss, nfet_params, pfet_params)
        i = (rI*i_I + rQ*i_Q)
        if vo >= vdd and i > 0:
            i = 0.0
            vo = vdd
        elif vo <= vss and i < 0:
            i = 0.0
            vo = vss
        interpolated[n] = vo
    return interpolated

C = 0.1

I = make_exp_clock(0.0, 1.0, UIS, 2*n_samples, VMIN, VMAX)
Q = make_exp_clock(0.5, 1.0, UIS, 2*n_samples, VMIN, VMAX)

#interp0 = interpolate(I, Q, 1.0, 0.0, (VMAX-VMIN)/2.0,  0.1, 0.001, VMAX, VMIN)
#interp1 = interpolate(I, Q, 0.707, 0.707, (VMAX-VMIN)/2.0,  0.1, 0.001, VMAX, VMIN)
#interp2 = interpolate(I, Q, 0.0, 1.0, (VMAX-VMIN)/2.0,  0.1, 0.001, VMAX, VMIN)
#plt.plot(interp0)
#plt.plot(interp1)
#plt.plot(interp2)
#plt.axhline(pfet_params["VTO"]+VMAX)
#plt.axhline(nfet_params["VTO"]+VMIN)
#plt.plot(I)
#plt.plot(Q)
#plt.show()

#plt.figure(1)
crossings = []
for n in range(STEPS*OVERSAMPLE):
    #phase_n = I*(1-(n/float(STEPS*OVERSAMPLE))) + Q*(n/float(STEPS*OVERSAMPLE))
    rI = np.cos(PI*n/(2.0*STEPS*OVERSAMPLE))
    rQ = np.sin(PI*n/(2.0*STEPS*OVERSAMPLE))
    phase_n = interpolate(I, Q, rI, rQ, (VMAX-VMIN)/2.0, C, 0.001, VMAX, VMIN)[n_samples:]
    #plt.plot(time, phase_n-(VMAX-VMIN)/2.0)
    crossings.append(np.where(np.diff(np.sign(phase_n-(VMAX-VMIN)/2.0)))[0][0])

crossings = np.array(crossings)/float(SAMPLES_PER_UI)
offset = crossings[0]

linearized = np.interp(np.arange(STEPS*OVERSAMPLE)/float(2*STEPS*OVERSAMPLE), crossings-offset, np.arange(STEPS*OVERSAMPLE))

reduced_linearized = [v/float(OVERSAMPLE) for (n,v) in enumerate(linearized) if n%OVERSAMPLE==0]

#plt.plot(crossings)
#plt.show()

plt.figure(1)
#plt.plot(crossings)
lin_crossings = []
for n in reduced_linearized:
    #phase_n = I*(1-(n/float(STEPS))) + Q*(n/float(STEPS))
    rI = np.cos(PI*n/(2.0*STEPS))
    rQ = np.sin(PI*n/(2.0*STEPS))
    phase_n = interpolate(I, Q, rI, rQ, (VMAX-VMIN)/2.0, C, 0.001, VMAX, VMIN)[n_samples:]
    plt.plot(time, phase_n)
    lin_crossings.append(np.where(np.diff(np.sign(phase_n-(VMAX-VMIN)/2.0)))[0][0])
lin_crossings = np.array(lin_crossings)/float(SAMPLES_PER_UI)
plt.figure(2)
#plt.plot(crossings)
plt.plot(lin_crossings)
plt.plot(np.arange(STEPS), np.arange(STEPS)/float(2*STEPS))

print("\n* Calculated ratio for currents")
print("\tn_step\tphase\trel\tDNL\tI\tQ")
for n, val in enumerate(reduced_linearized):
    #I = 1 - (val/float(STEPS))
    #Q = val/float(STEPS)
    I = np.cos(PI*n/(2.0*STEPS))
    Q = np.sin(PI*n/(2.0*STEPS))
    print("\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f"%(n, lin_crossings[n], lin_crossings[n]-lin_crossings[0],(lin_crossings[n]-lin_crossings[0])-(n/float(2*STEPS)), I, Q))


plt.show()
