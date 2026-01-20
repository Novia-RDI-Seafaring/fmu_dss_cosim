# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 14:43:48 2025

@author: emrkay
"""


# ---------------------------------------------------------------------------
# Helper for JSON <-> complex
# ---------------------------------------------------------------------------

def dict_to_complex(d):
    return complex(d["re"], d["im"])
