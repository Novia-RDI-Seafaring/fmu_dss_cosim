# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 12:59:41 2025

@author: emrkay
"""

import numpy as np
from pyfmi import load_fmu
import opendssdirect as dss
from Zs import compute_Zs_pu
from Zc import compute_Zc_pu
from Zm import compute_Zm_pu_from_dss
import json
from pathlib import Path
from Zk_pu_from_fmu import compute_Zk_pu_from_fmu
import os

Vbase_TS_kV = 400.0   # PCC base in the FMU
Vbase_DS_kV = 400.0    # PCC base in OpenDSS

fmu_path = "Two_Areas_TS_FMU.fmu"
dss_master = "34Bus/ieee34Mod2_forcosim2.dss"
# ------------------------------
# Helper functions
# ------------------------------

def complex_from_mag_angle(mag, ang_deg):
    """mag ∠ ang_deg -> complex."""
    return mag * np.exp(1j * np.deg2rad(ang_deg))

def mag_angle_from_complex(z):
    """complex -> (mag, ang_deg)."""
    return abs(z), np.rad2deg(np.angle(z))


def get_dss_slack_power_pu(S_base_MVA, slack_vsource_name="PCC"):
    """
    Return net complex power *absorbed by the DS* at the interface (Vsource),
    in per unit. Positive = load seen by the TS.
    """
    # Make Vsource.<name> the active element
    dss.Vsources.Name(slack_vsource_name)
    el_name = f"Vsource.{slack_vsource_name}"
    dss.Circuit.SetActiveElement(el_name)

    # Powers() -> [P1, Q1, P2, Q2, P3, Q3, ...] in kW/kvar
    pq = dss.CktElement.Powers()
    P_kW = pq[0] + pq[2] + pq[4]
    Q_kvar = pq[1] + pq[3] + pq[5]

    S_MVA_src = (P_kW + 1j*Q_kvar) / 1000.0   # source injection
    S_MVA_load = -S_MVA_src                   # power absorbed by DS

    return S_MVA_load / S_base_MVA

def set_dss_slack_voltage_pu(V_pu, Vbase_kV, slack_vsource_name=None):
    """
    Set OpenDSS slack/source voltage magnitude and angle (in pu).
    Assumes 1-phase or positive-sequence representation.
    """
    mag, ang_deg = mag_angle_from_complex(V_pu)

    if slack_vsource_name is not None:
        # Activate that Vsource explicitly (if you use multiple Vsources)
        dss.Vsources.Name(slack_vsource_name)
    else:
        dss.Vsources.First()  # assume the first one is the slack

    # In OpenDSS, Vsources.kV is the phase-to-neutral kV
    
    dss.Vsources.PU(mag)
    dss.Vsources.AngleDeg(ang_deg)

def get_dss_bus_voltage_pu(bus_name, Vbase_kV):
    """
    Returns bus positive-sequence voltage in pu as complex.
    For 1-phase or positive-sequence, we just take the first phase.
    """
    dss.Circuit.SetActiveBus(bus_name)
    vmag_angle = dss.Bus.VMagAngle()  # [V1_mag, V1_ang, V2_mag, V2_ang, ...]
    V1_mag = vmag_angle[0]
    V1_ang = vmag_angle[1]
    V1 = complex_from_mag_angle(V1_mag, V1_ang)
    Vbase_V = (Vbase_kV / np.sqrt(3)) * 1000.0
    
    #Vbase_V = Vbase_kV * 1000.0   # line-to-neutral base
    return V1 / Vbase_V

# ------------------------------
# TS–DS initialization mode
# ------------------------------

def initialize_ts_ds_coupling(
        fmu_path,
        dss_master_file,
        #interface_bus_ts,
        interface_bus_ds,
        Zs_pu,
        #Zc_pu,
        S_base_MVA,
        #Vbase_kV,
        tol=1e-3,
        max_iter=50,
        slack_vsource_name=None,
        verbose=True
    ):
    """
    Implements Section 5.1 initialization:
    - Iterative power flow between TS (FMU) and DS (OpenDSS)
    - Computes initial historical currents Hk(0), Hm(0)

    Returns:
        {
          "Vs_pu": Vs_pu,
          "Vpcc_ts_pu": Vpcc_ts_pu,
          "Vpcc_ds_pu": Vpcc_ds_pu,
          "Ik_pu": Ik_pu,
          "Yc_pu": Yc_pu,
          "Hk0": Hk0,
          "Hm0": Hm0,
          "fmu": fmu   # initialized FMU object (at the converged point)
        }
    """
    deltaVpccs=[];Vpcc_ts_pus=[];Vpcc_ds_pus=[]
    # --- Load FMU (transmission system) ---
    fmu = load_fmu(fmu_path)

    # --- Prepare OpenDSS (distribution system) ---
    dss.Basic.ClearAll()
    dss.Text.Command(f"Compile {dss_master_file}")

    # Initial guess Vs = 1∠0 pu
    Vs_pu = 1.0 + 0.0j

    # Characteristic admittance of fictitious line in pu

    # Iterative TS–DS power flow
    for it in range(max_iter):
        
        # 1) Set DS slack voltage to Vs (pu)

        # Vs_pu is on TS base. Convert to DS base before sending to OpenDSS.
        Vs_pu_DS = Vs_pu #* (Vbase_TS_kV / Vbase_DS_kV)
        set_dss_slack_voltage_pu(Vs_pu_DS, Vbase_DS_kV, slack_vsource_name="source")

        # fmu.exit_initialization_mode()
        # 2) DS power flow -> SDS
        dss.Solution.Solve()

        SDS_pu = get_dss_slack_power_pu(S_base_MVA, slack_vsource_name="source")
        if verbose:
            print(f"  SDS_pu(load at PCC) = {SDS_pu.real:.3f} + j{SDS_pu.imag:.3f} pu")

        
        fmu.reset()
        fmu.setup_experiment(start_time=0.0)
        fmu.enter_initialization_mode()
        fmu.set("Sds_re", SDS_pu.real)
        fmu.set("Sds_im", SDS_pu.imag)

    
        try:
            fmu.exit_initialization_mode()
        except Exception as e:
            print(f"FMU failed at iteration {it} with SDS_pu={SDS_pu}: {e}")
            try:
                for line in fmu.get_log():
                    print(line)
            except Exception:
                pass
            raise

        # 4) Read TS interface voltage and current (phasors in pu)
        # >>> ADAPT VARIABLE NAMES <<<
        Vpcc_ts_re = fmu.get("Vpcc_re")   # pu
        Vpcc_ts_im = fmu.get("Vpcc_im")
        Ipcc_ts_re = fmu.get("Im_re")   # pu, current leaving TS into line
        Ipcc_ts_im = fmu.get("Im_im")

        Vpcc_ts_pu = Vpcc_ts_re + 1j * Vpcc_ts_im
        Ik_pu = Ipcc_ts_re + 1j * Ipcc_ts_im  # current at terminal k (TS side)
        #Ik_pu = np.conjugate(SDS_pu/Vpcc_ts_pu)
        print(str(np.conjugate(SDS_pu/Vpcc_ts_pu)) + '   '+str(Ik_pu))
        # 5) DS interface voltage (bus on DS side) in pu

        Vpcc_ds_pu_DS = get_dss_bus_voltage_pu(interface_bus_ds, Vbase_DS_kV)

        # Convert DS pu → TS pu so it is comparable to Vpcc_ts_pu:
        Vpcc_ds_pu = Vpcc_ds_pu_DS #* (Vbase_DS_kV / Vbase_TS_kV)

        # 6) Update source voltage Vs = Vpcc_ts - Zs * Im
        # Here Im is current at DS side, with opposite sign of Ik in steady state.
        # For the iterative coupling of Section 5.1 they use:

        # but you can also write in terms of Ik with consistent sign convention.
        Im_pu = -Ik_pu   # by convention Im = -Ik at steady state
        Vs_new = Vpcc_ts_pu - Zs_pu * Im_pu
        print("Im_pu:"+str(Im_pu))
        # 7) Convergence check: |V_PCC_OM - V_PCC_OpenDSS| < eps
        mismatch = abs(Vpcc_ts_pu - Vpcc_ds_pu)

        if verbose:
            print("Iter " +str(it)+": ΔV_pcc = "+ str(mismatch) + " pu, Vs = "+str(abs(Vs_new))+" pu")
        deltaVpccs.append(mismatch)
        Vpcc_ts_pus.append(Vpcc_ts_pu)
        Vpcc_ds_pus.append(Vpcc_ds_pu)
        Vs_pu = Vs_new
        print("Vs_pu:"+str(Vs_pu))
        if mismatch < tol:
            if verbose:
                print(f"Converged in {it+1} iterations.")
            break
    else:
        print(f"TS–DS initialization did not converge in {max_iter} iterations.")



    # At this point, FMU is initialized at t=0 with consistent TS–DS conditions
    return {
        "Vs_pu": Vs_pu,
        "Vpcc_ts_pu": Vpcc_ts_pu,
        "Vpcc_ds_pu": Vpcc_ds_pu,
        "Ik_pu": Ik_pu,
        "deltaVpccs": deltaVpccs,   # add this
        "Vpcc_ts_pus":Vpcc_ts_pus,
        "Vpcc_ds_pus":Vpcc_ds_pus,
        # "Hk0": Hk0,
        # "Hm0": Hm0,
        "fmu": fmu,
    }





Zk_pudict = compute_Zk_pu_from_fmu(
    fmu_path="Two_Areas_TS_FMU.fmu",
    S_test_pu=0.1 + 0j,   # e.g. 0.1 pu real power load
    hin_pu=0+0j           # usually zero during PF-style init
)
Zk_pu = Zk_pudict[0]



cwd = os.getcwd()
S_base_MVA = 100.0
#Vbase_kV = 400.0 #/ np.sqrt(3)  # example: 69-kV line-to-line -> phase-to-neutral
Zs_pu_TS, Zs_ohm, bus_kV = compute_Zs_pu(
    dss_master_file=dss_master,
    interface_bus="pcc",
    S_base_MVA=S_base_MVA,
    Vbase_target_kV=Vbase_TS_kV,   # put Zs on TS base
)

print("Zs_ohm =", Zs_ohm)
print("Zs_pu on TS base =", Zs_pu_TS)

os.chdir(cwd)
# Then in your initialization:
Zs_pu = Zs_pu_TS

res = initialize_ts_ds_coupling(
    fmu_path=fmu_path,
    dss_master_file=dss_master,
    #interface_bus_ts="bus7_ts",      # adapt to your FMU / naming
    interface_bus_ds="pcc",       # DS interface bus name in OpenDSS
    Zs_pu=Zs_pu,
    S_base_MVA=S_base_MVA,

    tol=2e-4,

    max_iter=1000,

    slack_vsource_name=None,
    
    verbose=True
)

from matplotlib import pyplot as plt
# Now you can enter your co-simulation integration loop,
# using res["Hk0"], res["Hm0"] as initial historical terms.
mismatch = res["deltaVpccs"]

plt.figure()

Vpcc_ts_mag = [abs(v) for v in res["Vpcc_ts_pus"]]
Vpcc_ds_mag = [abs(v) for v in res["Vpcc_ds_pus"]]

#plt.plot(mismatch, marker='.', label="|V_TS - V_DS| (pu)")
plt.plot(Vpcc_ts_mag, marker='.', label="|V_TS| (pu, TS base)")
plt.plot(Vpcc_ds_mag, marker='.', label="|V_DS| (pu, TS base)")
plt.legend()

plt.xlabel("Iteration")
plt.ylabel("Voltage mismatch |V_TS − V_DS| (pu)")
plt.title("TS–DS Initialization Mismatch per Iteration")
plt.grid(True, which="both")
plt.show()


# --- Zm from DS (OpenDSS) ---
Zm_pu, Zm_ohm, Vbase_DS_kV = compute_Zm_pu_from_dss(
    dss_master_file=cwd+'\\'+dss_master,
    interface_bus="pcc",
    S_base_MVA=S_base_MVA,
)



Zc_pu, lam = compute_Zc_pu(Zk_pu, Zm_pu)

print("Chosen Zc_pu:", Zc_pu)
print("|lambda|:", lam)

# after convergence

Vk0 = res["Vpcc_ts_pus"][-1][0]
Ik0 = res["Ik_pu"][0]




def complex_to_dict(z: complex) -> dict:
    return {"re": z.real, "im": z.imag}

data = {
    "Zc_pu": complex_to_dict(Zc_pu),
    "Vk0":   complex_to_dict(Vk0),
    "Ik0":   complex_to_dict(Ik0),
    "Zk_pu": complex_to_dict(Zk_pu),
    "Vs_pu": complex_to_dict(res["Vs_pu"][0])
}
os.chdir(cwd)
out_path = Path("coupling_params.json")
with out_path.open("w") as f:
    json.dump(data, f, indent=2)


