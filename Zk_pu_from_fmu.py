# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 13:05:08 2025

@author: emrkay
"""

from pyfmi import load_fmu
import numpy as np

def compute_Zk_pu_from_fmu(
    fmu_path: str,
    S_test_pu: complex = 0.1 + 0j,
    hin_pu: complex = 0+0j,
) -> complex:
    """
    Estimate Zk_pu (Thevenin impedance of the TS at the PCC) from a TS FMU.

    Assumes the FMU has:
        - Real input  Sds_re, Sds_im  (DS complex power at PCC, pu)
        - Real output Vpcc_re, Vpcc_im (PCC voltage, pu, TS base)
        - Real output Im_re, Im_im     (PCC current TS->DS, pu)

    Method:
        1) Run FMU initialization with Sds = 0   -> (V1, I1)
        2) Run FMU initialization with Sds = S_test_pu -> (V2, I2)
        3) Use Thevenin relation V = Vs_th - Zk * I:

           V1 = Vs_th - Zk * I1
           V2 = Vs_th - Zk * I2

           => Zk = (V1 - V2) / (I2 - I1)

    Parameters
    ----------
    fmu_path : str
        Path to the FMU file (e.g. "Two_Areas_TS_FMU.fmu").
    S_test_pu : complex, optional
        Test power (pu) to draw at the PCC in the second run.
        Should be small enough to stay near the operating point (e.g. 0.05–0.2 pu).
    hin_pu : complex, optional
        Value for the fictitious-line history input hin during initialization
        (usually 0 for pure PF-style initialization).

    Returns
    -------
    Zk_pu : complex
        Estimated Thevenin impedance in per-unit on the TS base.
    """

    fmu = load_fmu(fmu_path, kind="cs")

    def run_steady(Sds_pu: complex):
        """Run one steady-state init and return (V_pu, I_pu)."""
        fmu.reset()
        fmu.setup_experiment(start_time=0.0)
        fmu.enter_initialization_mode()

        # Set DS power and history (if present)
        fmu.set("Sds_re", float(Sds_pu.real))
        fmu.set("Sds_im", float(Sds_pu.imag))

        # If your FMU has hin_re/hin_im, set them; otherwise remove these lines
        fmu.set("hin_re", float(hin_pu.real))
        fmu.set("hin_im", float(hin_pu.imag))

        fmu.exit_initialization_mode()

        # Read PCC voltage and current (pu)
        V = fmu.get("Vpcc_re") + 1j * fmu.get("Vpcc_im")
        I = fmu.get("Im_re")   + 1j * fmu.get("Im_im")
        return V, I

    # 1) No-load case (approx open circuit)
    V1, I1 = run_steady(0+0j)

    # 2) Small test load
    V2, I2 = run_steady(S_test_pu)

    # Avoid division by zero / degenerate case
    dI = I2 - I1
    if abs(dI) < 1e-8:
        raise RuntimeError(
            f"Test current change too small (|I2 - I1|={abs(dI):.3e}). "
            "Try a larger |S_test_pu|."
        )

    # Zk = (V1 - V2) / (I2 - I1)
    Zk_pu = (V1 - V2) / dI
    return Zk_pu
