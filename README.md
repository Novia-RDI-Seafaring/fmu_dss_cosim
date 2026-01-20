----------------
Unofficial implementation of the paper 

**"Co-simulation applied to power systems with high penetration of distributed energy resources" (Chagas & Tomim, 2022)**
----------------
**Requirements:** matplotlib, numpy, PyFMI, tkinter, json, pathlib, opendssdirect

**How to use:**
- Run Zc_GUI.py to start the GUI for studying the effects of Zc selection.

- Run fmu_dss_cosim_initialization.py for computing the coupling parameters for the FMU-OpenDSS co-simulation.

- Run cosim_master_fmu_opendss2.py to co-simulate FMU and OpenDSS.

**Note:** OpenModelica source files (.mo files) are provided for reference. Install OpenModelica and OpenIPSL if you want to export FMUS by compiling them. 
