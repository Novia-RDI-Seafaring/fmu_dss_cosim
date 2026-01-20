# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 16:32:46 2025

@author: emrkay
"""

import opendssdirect as dss
import numpy as np

def compute_Zs_pu(
    dss_master_file: str,
    interface_bus: str,
    S_base_MVA: float,
    Vbase_target_kV: float = None,
):
    """
    Compute Zs (Thevenin/source impedance between OpenDSS slack and PCC)
    from the OpenDSS model using the positive-sequence short-circuit
    impedance at the PCC bus.

    Parameters
    ----------
    dss_master_file : str
        Path to the OpenDSS master file (e.g. 'Master.dss').
    interface_bus : str
        Name of the DS interface bus in OpenDSS (e.g. 'pcc').
    S_base_MVA : float
        System MVA base used for per-unit (same base as TS FMU).
    Vbase_target_kV : float, optional
        Voltage base (kV) on which you want Zs in per-unit.
        If None, returns Zs in per-unit on the PCC bus base kV.

    Returns
    -------
    Zs_pu : complex
        Source/Thevenin impedance in per-unit on the chosen voltage base.
    Zs_ohm : complex
        Source/Thevenin impedance in ohms.
    bus_kV_base : float
        Base kV of the interface bus (line-to-line).
    """

    # Clear any existing circuit and compile
    dss.Basic.ClearAll()
    dss.Text.Command(f"Compile {dss_master_file}")

    # Run a fault study or refresh Zsc so Zsc1 is available
    # Either of the following strategies is OK:
    # 1) FaultStudy:
    dss.Text.Command("CalcVoltageBases")
    dss.Text.Command("Solve mode=faultstudy")
    # or:
    # 2) Use ZscRefresh on the active bus after a normal solve.
    # dss.Solution.Solve()

    # Set the active bus to the interface node
    dss.Circuit.SetActiveBus(interface_bus)

    # Make sure Zsc is computed for this bus
    dss.Bus.ZscRefresh()

    # Positive-sequence short-circuit impedance looking into the bus (ohms)
    Zs_ohm1 = dss.Bus.Zsc1()

    # Base voltage at this bus (kV, line-to-line)
    bus_kV_base = dss.Bus.kVBase()

    # Impedance base (ohms) on the bus base
    # Z_base = V_base^2 / S_base (with V in kV, S in MVA)
    Z_base_bus_ohm = (bus_kV_base ** 2) / S_base_MVA
    Zs_ohm = Zs_ohm1[0]+1j*Zs_ohm1[1]
    # Per-unit Zs on the bus base
    Zs_pu_bus = Zs_ohm / Z_base_bus_ohm

    if Vbase_target_kV is None or np.isclose(Vbase_target_kV, bus_kV_base):
        # Return per-unit on the bus base
        return Zs_pu_bus, Zs_ohm, bus_kV_base

    # Convert to a different voltage base (same S_base)
    # Z_pu,new = Z_pu,old * (V_old^2 / V_new^2)
    Zs_pu_target = Zs_pu_bus * (bus_kV_base ** 2) / (Vbase_target_kV ** 2)

    return Zs_pu_target, Zs_ohm, bus_kV_base
