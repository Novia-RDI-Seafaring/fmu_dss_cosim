# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 18:54:20 2025

@author: emrkay
"""
import numpy as np

# --- your f(Zc) from complex impedances ---
def f_of_Zc(Zc, Zk, Zm):
    Zc_abs, Zk_abs, Zm_abs = np.abs(Zc), np.abs(Zk), np.abs(Zm)
    dck = np.angle(Zc) - np.angle(Zk)
    dcm = np.angle(Zc) - np.angle(Zm)

    # numerical safety: clamp tiny negatives before sqrt
    def s(x): return np.sqrt(np.maximum(x, 0.0))

    num_k = s(Zc_abs**2 - 2*Zc_abs*Zk_abs*np.cos(dck) + Zk_abs**2)
    num_m = s(Zc_abs**2 - 2*Zc_abs*Zm_abs*np.cos(dcm) + Zm_abs**2)
    den_k = s(Zc_abs**2 + 2*Zc_abs*Zk_abs*np.cos(dck) + Zk_abs**2)
    den_m = s(Zc_abs**2 + 2*Zc_abs*Zm_abs*np.cos(dcm) + Zm_abs**2)

    return s((num_k * num_m) / (den_k * den_m))

# --- helpers ---
def polar(mag, deg):                  # mag∠deg° -> complex
    return mag * np.exp(1j * np.deg2rad(deg))

def f_Zc(Zc, Zk, Zm):
    """
    Compute f(Zc, θc) from the given formula,
    where Zc, Zk, Zm are complex impedances.

    Parameters
    ----------
    Zc, Zk, Zm : complex
        Complex impedances

    Returns
    -------
    f : float
    """

    # Magnitudes
    Zc_abs = np.abs(Zc)
    Zk_abs = np.abs(Zk)
    Zm_abs = np.abs(Zm)

    # Angles
    theta_c = np.angle(Zc)
    theta_k = np.angle(Zk)
    theta_m = np.angle(Zm)

    dck = theta_c - theta_k
    dcm = theta_c - theta_m

    # Numerator terms
    num_k = np.sqrt(Zc_abs**2 - 2*Zc_abs*Zk_abs*np.cos(dck) + Zk_abs**2)
    num_m = np.sqrt(Zc_abs**2 - 2*Zc_abs*Zm_abs*np.cos(dcm) + Zm_abs**2)

    # Denominator terms
    den_k = np.sqrt(Zc_abs**2 + 2*Zc_abs*Zk_abs*np.cos(dck) + Zk_abs**2)
    den_m = np.sqrt(Zc_abs**2 + 2*Zc_abs*Zm_abs*np.cos(dcm) + Zm_abs**2)

    # Final expression
    return np.sqrt((num_k * num_m) / (den_k * den_m))



def list_like(dct_keys, substrings=()):
    ks = [str(k) for k in dct_keys]
    if not substrings:
        return ks
    ll = []
    for k in ks:
        low = k.lower()
        if any(s.lower() in low for s in substrings):
            ll.append(k)
    return ll

def ensure_vars_exist(model, required):
    """Check required variable names exist; if not, print a helpful list and raise."""
    names = set(model.get_model_variables().keys())
    missing = [v for v in required if v not in names]
    if missing:
        print("\n[!] Missing variables in FMU:", model.get_name())
        print("    Missing:", missing)
        print("    Available (filtered on pcc/inj/slack/delta):")
        print("   ", list_like(names, ('pcc', 'inj', 'slack', 'delta')))
        raise RuntimeError(f"FMU {model.get_name()} is missing variables: {missing}")

def complex_from(model, re_name, im_name):
    return complex(float(model.get(re_name)), float(model.get(im_name)))

def set_complex(model, re_name, im_name, value: complex):
    model.set(re_name, float(np.real(value)))
    model.set(im_name, float(np.imag(value)))

def niAE(y, y_ref, DT, EPS):
    """Normalized Integral of Absolute Error (0..1), where 1 is perfect."""
    y = np.asarray(y, dtype=float)
    y_ref = np.asarray(y_ref, dtype=float)
    denom = np.trapz(np.abs(y_ref), dx=DT)
    if denom < EPS:
        return float('nan')
    numer = np.trapz(np.abs(y - y_ref), dx=DT)
    return 1.0 - (numer / denom)




def fmu_is_cs(m):  # CS FMUs have do_step
    return hasattr(m, "do_step")

def init_cs(m, start=0.0):
    # Recover from any prior error
    try: m.reset()
    except Exception: pass
    m.setup_experiment(start_time=start)
    m.enter_initialization_mode()
    # set default inputs here
    for nm in ("Iinj_re","Iinj_im"):
        if nm in m.get_model_variables(): m.set(nm, 0.0)
    m.exit_initialization_mode()
    # time reference for our local stepping
    return start

def init_me(m, start=0.0):
    try: m.reset()
    except Exception: pass
    m.setup_experiment(start_time=start)
    m.initialize()  # transitions to step-capable state
    # set default inputs after initialize for many ME FMUs
    for nm in ("Iinj_re","Iinj_im"):
        if nm in m.get_model_variables(): m.set(nm, 0.0)
    # return time and current states (may be empty)
    try: x = m.get_continuous_states()
    except Exception: x = np.array([])
    return start, x

def step_cs(m, t, dt):
    m.do_step(current_t=t, step_size=dt)
    return t+dt

def step_me(m, t, dt, x):
    # very light explicit Euler; adequate for our tiny probe steps
    try:
        m.time = t
    except Exception:
        pass
    if x.size:
        m.set_continuous_states(x)
        dx = m.get_derivatives()
        x  = x + dt*dx
        m.set_continuous_states(x)
    t = t + dt
    try: m.completed_integrator_step()
    except Exception: pass
    try: m.event_update()
    except Exception: pass
    return t, x

def pick_Zc(Zk, Zm):
    def lam_mag(Zc):
        Gk = (Zc - Zk)/(Zc + Zk)
        Gm = (Zc - Zm)/(Zc + Zm)
        return abs(Gk*Gm)**0.5
    cands = {
        "Zk": Zk,
        "Zm": Zm,
        "mid": 0.5*(Zk+Zm),
    }
    best = min(cands.items(), key=lambda kv: lam_mag(kv[1]))
    return best[0], best[1], {name: lam_mag(z) for name,z in cands.items()}