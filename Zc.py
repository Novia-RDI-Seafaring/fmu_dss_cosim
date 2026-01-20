# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 18:36:21 2025

@author: emrkay
"""

import cmath
import math
from typing import Tuple

def compute_Zc_pu(
    Zk_pu: complex,
    Zm_pu: complex,
    mag_min: float = 0.5,
    mag_max: float = 1.5,
    n_mag: int = 40,
    ang_span_deg: float = 30.0,
    n_ang: int = 61,
) -> Tuple[complex, float]:
    """
    Compute a suitable characteristic impedance Zc (in pu) for the fictitious line,
    given transmission-side and distribution-side Thevenin impedances Zk, Zm (in pu).

    The function numerically searches around Zk in polar coordinates and chooses
    the Zc that minimizes |lambda|, where:
        Gamma_k = (Zk - Zc) / (Zk + Zc)
        Gamma_m = (Zm - Zc) / (Zm + Zc)
        lambda^2 = Gamma_k * Gamma_m
    and we want |lambda| < 1 for stability, ideally as small as possible.

    Parameters
    ----------
    Zk_pu : complex
        Thevenin impedance of the transmission system at the interface (per-unit).
    Zm_pu : complex
        Thevenin impedance of the distribution system at the interface (per-unit).
    mag_min, mag_max : float
        Relative search bounds for |Zc| around |Zk|.
        Zc magnitude is swept from mag_min*|Zk| to mag_max*|Zk|.
    n_mag : int
        Number of magnitude points in the search grid.
    ang_span_deg : float
        Half-width of the angle search window around angle(Zk), in degrees.
        Total span is [-ang_span_deg, +ang_span_deg] around angle(Zk).
    n_ang : int
        Number of angle points in the search grid.

    Returns
    -------
    Zc_best : complex
        Selected Zc in per-unit.
    lambda_best : float
        Corresponding |lambda| value.
    """

    # Handle degenerate case
    if Zk_pu == 0:
        # As a desperate fallback, just return something small
        return 0.01 + 0j, 0.0

    # Polar representation of Zk
    Zk_mag = abs(Zk_pu)
    Zk_ang = cmath.phase(Zk_pu)  # radians

    # Magnitude and angle grids
    mag_vals = [Zk_mag * (mag_min + (mag_max - mag_min) * i / max(n_mag - 1, 1))
                for i in range(n_mag)]
    ang_min = Zk_ang - math.radians(ang_span_deg)
    ang_max = Zk_ang + math.radians(ang_span_deg)
    ang_vals = [ang_min + (ang_max - ang_min) * j / max(n_ang - 1, 1)
                for j in range(n_ang)]

    def lambda_mag(Zc: complex) -> float:
        # Avoid singularities where Zk + Zc or Zm + Zc ~ 0
        if abs(Zk_pu + Zc) < 1e-9 or abs(Zm_pu + Zc) < 1e-9:
            return float("inf")

        Gamma_k = (Zk_pu - Zc) / (Zk_pu + Zc)
        Gamma_m = (Zm_pu - Zc) / (Zm_pu + Zc)

        # lambda^2 = Gamma_k * Gamma_m  → lambda = sqrt(...)
        lam2 = Gamma_k * Gamma_m
        lam = cmath.sqrt(lam2)
        return abs(lam)

    Zc_best = Zk_pu
    lambda_best = lambda_mag(Zk_pu)

    # Grid search
    for mag in mag_vals:
        for ang in ang_vals:
            Zc_candidate = cmath.rect(mag, ang)
            lam_abs = lambda_mag(Zc_candidate)

            # Only accept candidates with |lambda| < 1
            if lam_abs < 1.0 and lam_abs < lambda_best:
                lambda_best = lam_abs
                Zc_best = Zc_candidate

    # If nothing better found (all |lambda| >= 1), just return Zk itself
    # (in practice this is usually okay).
    return Zc_best, lambda_best
