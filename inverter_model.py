import numpy as np
import matplotlib.pyplot as plt
from golden_section import *

def dict_comb(*dicts):
    comb_dict = {}
    for _dict in dicts:
        for k, v in _dict.items():
            comb_dict[k] = v
    return comb_dict

def nmos_iv_level1(vg, vs, vd, VTO, K, W, L, LAMBDA):
    veff = vg-vs-VTO
    vds = vd-vs
    if veff <= 0.0: # cutoff
        return 0.0
    elif vds < veff: # linear mode
        return K*(W/L)*(veff*vds-0.5*vds**2)
    else:   # saturation
        # CJM Eq 1.84 (pg 27) says to use 1.0 + LAMBDA*(vds-vedd)
        # Razavi Eq 2.29 (pg 23) says 1.0 + LAMBDA*vds, WITH discontinuity,
        #   says advanced models should be used to remove discontinuity
        return 0.5*K*(W/L)*veff**2*(1.0+LAMBDA*(vds-veff))

def pmos_iv_level1(vg, vs, vd, VTO, K, W, L, LAMBDA):
    veff = vg-vs-VTO
    vds = vd-vs
    if veff >= 0.0: # cutoff
        return 0.0
    elif vds > veff: # linear mode
        return -K*(W/L)*(veff*vds-0.5*vds**2)
    else:   # saturation
        return -0.5*K*(W/L)*veff**2*(1.0-LAMBDA*(vds-veff))

def inverter_iv(vi, vo, vdd, vss, nfet_params, pfet_params, nmos=nmos_iv_level1, pmos=pmos_iv_level1, *args, **kwargs):
    #solve op point
    vo = (vdd-vss)/2.0
    vo_last = 0.0
    i = 0.0
    i_n = nmos(vg=vi, vs=vss, vd = vo, **nfet_params)
    i_p = pmos(vg=vi, vs=vdd, vd = vo, **pfet_params)
    # params = dict_comb(pmos_op_point, pfet_params)
    # vd_pmos = gss_find_arg(pmos_iv_level1, arg="vd", target=-i, params=params, _min=vss, _max=vdd, conv_tol=conv_tol/10.0)
    #vo_last = vo
    #vo = vdd-vd_pmos
    #print i_n, i_p
    return -i_n-i_p

v_nmos_iv_level1 = np.vectorize(nmos_iv_level1, otypes=[float])
v_pmos_iv_level1 = np.vectorize(pmos_iv_level1, otypes=[float])
v_inverter_iv = np.vectorize(inverter_iv, otypes=[float])

nfet_params = {
    "W" : 1.0,
    "L" : 1.0,
    "K" : 1.0,
    "VTO" : 1.0,
    "LAMBDA" : 0.01,
}

pfet_params = {
    "W" : 1.0,
    "L" : 1.0,
    "K" : 1.0,
    "VTO" : -1.0,
    "LAMBDA" : 0.01,
}

#VDD = 5.0
#VSS = 0.0

#vg = np.linspace(VSS, VDD, 1000)
#i_n = v_nmos_iv_level1(vg=vg, vs=VSS, vd=VDD, **nfet_params)
#i_p = v_pmos_iv_level1(vg=vg, vs=VDD, vd=VSS, **pfet_params)

#plt.figure(0)
#plt.plot(vg, i_n)
#plt.plot(vg, i_p)

#plt.figure(1)
#vd = np.linspace(VSS, VDD, 1000)
#i_n = v_nmos_iv_level1(vg=2.5, vs=VSS, vd=vd, **nfet_params)
#i_p = v_pmos_iv_level1(vg=2.5, vs=VDD, vd=vd, **pfet_params)

#vi = np.linspace(VSS, VDD, 10)
#vo = v_inverter_iv(vi=vi, vo=2.5, vdd=VDD, vss=VSS, nfet_params=nfet_params, pfet_params=pfet_params, nmos=nmos_iv_level1, pmos=pmos_iv_level1)


#plt.plot(vi, vo)
#plt.show()

