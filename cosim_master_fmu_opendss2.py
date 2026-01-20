import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import opendssdirect as dss
from pyfmi import load_fmu
from helper_functs import dict_to_complex

# ---------------------------------------------------------------------------
# User-configurable paths & bases
# ---------------------------------------------------------------------------

FMU_PATH   = "Two_Areas_TS_FMU.fmu"
DSS_MASTER = "34Bus/ieee34Mod2_forcosim2.dss"   # adjust path if needed

SBASE_MVA = 100.0     # common 3φ base power
VBASE_TS  = 400.0     # TS PCC base voltage (kV, line-to-line)
VBASE_DS  = 400.0     # DS PCC base voltage (kV, line-to-line)

T_FINAL = 10.0        # simulation horizon (s)
DT      = 0.007  # macro step size (s)

tfault1 = .4 
tfault2 = .5



# ---------------------------------------------------------------------------
# DSS helper: everything returned in *TS* per-unit
# ---------------------------------------------------------------------------
def DSS_pu(Vs_ts_pu: complex) -> tuple[complex, complex]:
    """
    Given TS-side sending voltage Vs (pu on TS base),
    drive OpenDSS Vsource.Source and return:
        Im_pu_ts : interface current (pu on TS base, phase current)
        Vm_pu_ts : DS PCC voltage (pu on TS base, phase voltage)
    """
    # 1) Convert Vs_ts_pu to DS pu for OpenDSS
    V_LL_TS_kV = Vs_ts_pu * VBASE_TS
    Vs_ds_pu   = V_LL_TS_kV / VBASE_DS

    mag_ds  = abs(Vs_ds_pu)
    ang_deg = np.degrees(np.angle(Vs_ts_pu))

    # 👉 FIX: edit the actual source element, Vsource.Source
    dss.Text.Command(f"Edit Vsource.Source pu={mag_ds} angle={ang_deg}")

    dss.Solution.SolveNoControl()

    # PCC bus is called PCC in the .dss file
    dss.Circuit.SetActiveBus("PCC")

    # 2) Voltage at PCC (phase A, LN, volts)
    Vm_vec   = dss.Bus.Voltages()
    V_phase  = complex(Vm_vec[0], Vm_vec[1])

    Vbase_LN_TS_V = (VBASE_TS * 1e3) / np.sqrt(3.0)
    Vm_pu_ts      = V_phase / Vbase_LN_TS_V

    # 3) Interface current → TS pu
    P_kW, Q_kvar = dss.Circuit.TotalPower()
    S_W = (P_kW + 1j * Q_kvar) * 1e3

    I_phase_A = np.conj(S_W / (3.0 * V_phase))
    Ibase_phase_TS = (SBASE_MVA * 1e6) / (3.0 * Vbase_LN_TS_V)
    Im_pu_ts       = I_phase_A / Ibase_phase_TS

    return Im_pu_ts, Vm_pu_ts



# ---------------------------------------------------------------------------
# Fictitious line + DS fixed-point solve (paper Sec. 5)
# ---------------------------------------------------------------------------

def solveLineOpenDSS(Ehm, Vs0_pu, Zc_pu, Zk_pu, eps=1e-4, max_iter=20):
    """
    Fixed-point inner loop for DS + fictitious line.

    All quantities are TS per-unit:
      Ehm    : incoming history term from TS-FMU (E-type, voltage-like, terminal m)
      Vs0_pu : initial guess for sending voltage at TS terminal (k side)
      Zc_pu  : characteristic impedance of fictitious line
      Zk_pu  : TS Thevenin impedance at terminal k

    Returns
    -------
    Vs_pu   : converged sending voltage (pu)
    Im_pu   : DS interface current (pu)
    Vm_pu   : DS PCC voltage (pu)
    it      : number of iterations
    """
    Vs    = Vs0_pu
    Im_pu = 0j
    Vm_pu = 0j

    for it in range(1, max_iter + 1):
        # DS solve for current Vs

        Im_pu_raw, Vm_pu = DSS_pu(Vs)
        Im_pu = -Im_pu_raw   # Flip sign to match Bergeron sign convention
        # Bergeron relation at DS terminal (m):

        # => expected Vm from Ehm & Im:
        Vm_berg = Ehm - Zc_pu * Im_pu
        # print('it:'+str(it)+'  '+str(abs(Vm_berg - Vm_pu)))
        # print('Vm_berg:'+str(Vm_berg))
        # print('Vm_pu:'+str(Vm_pu))
        # Convergence on DS voltage
        if abs(Vm_berg - Vm_pu) < eps:
            break

        # Update sending voltage using TS Thevenin:
        # Vs = Vm_berg - Zk * Im
        Vs = Vm_berg - Zk_pu * Im_pu

    return Vs, Im_pu, Vm_pu, it


def solveEhk(Vs, Ehm, t, dt, Zc_pu, Zk_pu, Vs_tm1, Vs_tm2,
             eps=1e-3, max_iter=20):
    """
    Predictor–corrector for fictitious line + DS (corresponds to solveEhk in Fig.4)

    Inputs
    ------
    Vs      : current sending voltage (pu, TS base)
    Ehm     : incoming history term from TS-FMU (E-type) at terminal m
    t, dt   : time and macro-step size
    Zc_pu   : characteristic impedance (pu)
    Zk_pu   : TS Thevenin impedance (pu)
    Vs_tm1  : Vs at previous macro step
    Vs_tm2  : Vs at step before previous

    Returns
    -------
    Vs_new  : updated sending voltage (pu)
    Ehk     : outgoing history term (E-type) at terminal k (for FMU input)
    Im_pu   : DS interface current (pu)
    Vm_pu   : DS PCC voltage (pu)
    it      : inner fixed-point iterations
    """
    # 2-step predictor for Vs
    if t > 2 * dt and Vs_tm2 != 0:
        Vs_est = (Vs_tm1 ** 2) / Vs_tm2
    else:
        Vs_est = Vs_tm1

    Vs_new, Im_pu, Vm_pu, it = solveLineOpenDSS(
        Ehm=Ehm,
        Vs0_pu=Vs_est,
        Zc_pu=Zc_pu,
        Zk_pu=Zk_pu,
        eps=eps,
        max_iter=max_iter,
    )

    # E-type history at TS terminal k: Ehk = Vk + Zc * Ik
    # Here we approximate Vk ≈ Vs_new (TS terminal voltage)
    #Ehk = Vs_new + Zc_pu * Im_pu
    Ehk = Vs_new - Zc_pu * Im_pu
    return Vs_new, Ehk, Im_pu, Vm_pu, it


# ---------------------------------------------------------------------------
# Main co-simulation
# ---------------------------------------------------------------------------

    # ----- Load coupling parameters from JSON -----
params_path = Path("coupling_params.json")
if not params_path.exists():
    raise FileNotFoundError("coupling_params.json not found.")

data   = json.load(params_path.open())
Zc_pu  = dict_to_complex(data["Zc_pu"])

Zk_pu  = dict_to_complex(data["Zk_pu"])
Vk0    = dict_to_complex(data["Vk0"])
Ik0    = dict_to_complex(data["Ik0"])
Vs_init = dict_to_complex(data["Vs_pu"])

# Yc for optional H→E conversion
if Zc_pu == 0:
    raise ValueError("Zc_pu is zero in JSON; cannot run co-simulation.")
Yc = 1.0 / Zc_pu
#Yc_pu = 1.0 / Zc_pu   # passed as an argument
Hk0 = Yc * Vk0 + Ik0
Hm0 = Yc * Vk0 - Ik0

Yc_re = Yc.real
Yc_im = Yc.imag
# Interpret Hk0/Hm0 either as H (current-like) or E (voltage-like)

Ehm0 = Hk0 / Yc   # TS->DS history (E at m)
Ehk0 = Hm0 / Yc   # DS->TS history (E at k)


# ----- Initialize DSS -----

dss.Basic.ClearAll()
dss.Text.Command(f'Redirect "{DSS_MASTER}"')
dss.Text.Command("CalcVoltageBases")

# 👉 disable regulator/cap control iterations for co-sim
dss.Text.Command("Set controlmode=OFF")
dss.Text.Command("Set maxcontroliter=1")

# Do one initial snap solve
dss.Solution.SolveNoControl()
# ----- Load FMU -----
fmu_TS = load_fmu(FMU_PATH, kind="cs")
# Optional: see FMU logs if something goes wrong
# fmu_TS.set_debug_logging(True)

# ----- Initialize FMU -----
fmu_TS.reset()
fmu_TS.setup_experiment(start_time=0.0, stop_time=T_FINAL)
fmu_TS.enter_initialization_mode()
#set fault times
fmu_TS.set('pwFault.t1',float(tfault1))
fmu_TS.set('pwFault.t2',float(tfault2))
# Set fictitious line admittance (from Zc)
fmu_TS.set("Yc_re", float(Yc_re))
fmu_TS.set("Yc_im", float(Yc_im))
# FMU expects current-type history H at its input

Hin0 = Hm0          # DS→TS history as current (what FMU wants)

fmu_TS.set("hin_re", float(Hin0.real))
fmu_TS.set("hin_im", float(Hin0.imag))

fmu_TS.exit_initialization_mode()
# ----- Initial conditions for line state -----
Vs     = Vs_init
Vs_tm1 = Vs_init
Vs_tm2 = Vs_init

# For the external fictitious line:
Ehm = Ehm0   # FMU→DS history (read from FMU's hout at t=0 in the paper)
Ehk = Ehk0   # DS→FMU history (we just used this as initial hin)

# ----- Storage for plots -----
results = {
    "t":   [],
    "Vts": [],
    "Vds": [],
    "Vss": [],
    "g1_ang":[],
    "g4_ang":[],
    "g3_ang":[]
}

t = 0.0

# -------------------------------------------------------------------
# Time stepping loop (implements Fig. 4 pseudo-code)
# -------------------------------------------------------------------
eps = 1e-3
max_iter=20
while t < T_FINAL:
    # 1) Independent solutions
    status_fmu = fmu_TS.do_step(current_t=t, step_size=DT)
    if status_fmu != 0:
        print(f"FMU error (status={status_fmu}) at t={t}")
        break

    Vs_new, Ehk, Im_pu, Vm_pu, it = solveEhk(
        Vs=Vs,
        Ehm=Ehm,
        t=t,
        dt=DT,
        Zc_pu=Zc_pu,
        Zk_pu=Zk_pu,
        Vs_tm1=Vs_tm1,
        Vs_tm2=Vs_tm2,
        eps = eps,
        max_iter=max_iter
    )

    # Update Vs history
    Vs_tm2 = Vs_tm1
    Vs_tm1 = Vs
    Vs     = Vs_new
    #print('Vs:'+str(Vs))
    # # 2) Get outgoing history from FMU (hout → Ehm for next step)

    Hk_out = fmu_TS.get("hout_re") + 1j * fmu_TS.get("hout_im")  # current-like
    

    Ehm = Hk_out / Yc    # convert H → E

    
    H_in = Yc * Ehk   # convert E → H
    # 3) Feed new DS→TS history to FMU
    fmu_TS.set("hin_re", float(H_in.real))
    fmu_TS.set("hin_im", float(H_in.imag))
    # 4) Read TS PCC voltage for plotting
    Vk_re = fmu_TS.get("Vpcc_re")
    Vk_im = fmu_TS.get("Vpcc_im")
    Vts   = Vk_re + 1j * Vk_im

    # 5) Store TS and DS PCC voltages (TS base pu)
    results["t"].append(t)
    results["Vts"].append(abs(Vts))
    results["Vds"].append(abs(Vm_pu))
    
    results["g1_ang"].append(float(fmu_TS.get('g1.g1.delta')))
    results["g3_ang"].append(float(fmu_TS.get('g3.g3.delta')))
    t += DT

# -------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------
plt.figure(figsize=(10, 5))
plt.plot(results["t"], results["Vts"], label="FMU Voltage (TS pu)")
plt.plot(results["t"], results["Vds"], label="DS Voltage mapped to TS pu")
plt.title("Voltage Coupling")
plt.xlabel("Time [s]")
plt.ylabel("Voltage [pu on TS base]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

plt.figure(figsize=(6, 6))
plt.plot(results["t"], np.array(results["g1_ang"])-np.array(results["g3_ang"]), label="")

plt.title("relative angle btw G1 and G3")
plt.xlabel("Time [s]")
plt.ylabel("Angle [rad]")
plt.grid(True)
#plt.legend()
plt.tight_layout()
plt.show()
print("✔ Co-simulation complete.")


# if __name__ == "__main__":
#     main()
