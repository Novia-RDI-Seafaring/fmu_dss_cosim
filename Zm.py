# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 18:51:09 2025

@author: emrkay
"""

import opendssdirect as dss
from typing import Tuple

def compute_Zm_pu_from_dss(
    dss_master_file: str,
    interface_bus: str,
    S_base_MVA: float,
) -> Tuple[complex, complex, float]:
    """
    Compute Zm (DS Thevenin impedance) in per-unit from OpenDSS.

    Parameters
    ----------
    dss_master_file : str
        Path to the OpenDSS master file.
    interface_bus : str
        Name of the interface/PCC bus in OpenDSS.
    S_base_MVA : float
        Common per-unit MVA base.

    Returns
    -------
    Zm_pu : complex
        DS Thevenin impedance in per-unit (on the bus base).
    Zm_ohm : complex
        DS Thevenin impedance in ohms.
    bus_kV_base : float
        Base kV (line-to-line) of the interface bus in OpenDSS.
    """
    dss.Basic.ClearAll()
    dss.Text.Command(f"Compile {dss_master_file}")

    # Set voltage bases and run a fault study so that Zsc is available
    dss.Text.Command("CalcVoltageBases")
    dss.Text.Command("Solve mode=faultstudy")

    dss.Circuit.SetActiveBus(interface_bus)

    # Ensure Zsc matrices are refreshed
    dss.Bus.ZscRefresh()

    # Positive-sequence short-circuit impedance at the bus (ohms)
    Zm_ohm1 = dss.Bus.Zsc1()

    # Bus base voltage in kV (line-to-line)
    bus_kV_base = dss.Bus.kVBase()

    # Z_base = V_base^2 / S_base (kV^2 / MVA => ohms)
    Z_base_ohm = (bus_kV_base ** 2) / S_base_MVA
    Zm_ohm = Zm_ohm1[0]+1j*Zm_ohm1[1]
    Zm_pu = Zm_ohm / Z_base_ohm
    return Zm_pu, Zm_ohm, bus_kV_base
