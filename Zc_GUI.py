import tkinter as tk
from tkinter import ttk
import opendssdirect as dss
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.gridspec import GridSpec
from pyfmi import load_fmu

from plot_helper2 import Zc_eigenplotter   # uses ax=... variant

#FMU_DS = True
# ---------------------- system parameters ----------------------
TS_FMU = "IEEE_14_Buses_FMU.fmu"
#if FMU_DS:
DS_FMU = "IEEE_9_Buses_FMU.fmu"
# else:
#     DS_DSS = "IEEE_9_Buses.dss"

    # # Path to your DSS master file
    # dss_file = os.path.abspath(DS_DSS)

    # # 1) Compile the circuit
    # dss.Text.Command(f"compile [{dss_file}]")

fault_ts_t1 = 0.2
fault_ts_t2 = 0.3
fault_ds_t1 = 0.5
fault_ds_t2 = 0.6
# Thevenin impedances (pu)

Zk_mag = 0.151028 #0.200242 #0.0203  0.151028 ∠ -99.60
Zm_mag =  0.057539 #∠  #0.01465
Zk_ang = 80.40#-99.60#180#90.0
Zm_ang = 74.63#-105.37#90.0
r_max = 2*max(Zk_mag,Zm_mag)
Zk = Zk_mag * np.exp(1j * np.deg2rad(Zk_ang))
Zm = Zm_mag * np.exp(1j * np.deg2rad(Zm_ang))

Zc_init_mag = (Zk_mag) #/ 2
Zc_init_ang = (Zk_ang) #/ 2
Zc_init = Zc_init_mag * np.exp(1j * np.deg2rad(Zc_init_ang))

# timestep slider settings (4× range, 4× length)
H_DEFAULT = 1e-3
H_MIN = 1e-6
H_MAX = 5e-3
H_RES = 1e-5

# t_end slider settings
T_END_DEFAULT = 0.05
T_END_MIN = 0.01
T_END_MAX = 2.0
T_END_RES = 0.01


def get_VI(m):
    vr = m.get('V_re'); vi = m.get('V_im')
    ir = m.get('I_re'); ii = m.get('I_im')
    return complex(vr, vi), complex(ir, ii)

def DSS(Vs):
    pu_mag = abs(Vs) #/ Vbase
    angle_deg = np.degrees(np.angle(Vs))
    dss.Text.Command(f"Edit Vsource.PCC pu={pu_mag} angle={angle_deg}")
    dss.Solution.Solve()
    dss.Circuit.SetActiveBus("B1")
    Vm_vec  = dss.Bus.Voltages()

    # convert Vm_vec (real/imag pairs) to complex 3-phase; here simplify to pos-seq:
    Vm = complex(Vm_vec[0], Vm_vec[1])

    P_DSS, Q_DSS = dss.Circuit.TotalPower()
    Im = np.conj((P_DSS*1e3 + 1j*Q_DSS*1e3) / Vm)   # S = V * I*
    return Im,Vm

def solveLineOpenDSS(Ehm,Vs0,Zc,Zk,eps=0.1):
    Vs = Vs0
    h=0
    non_converged = True
    while non_converged:
        Imdss,Vmdss = DSS(Vs)
        Vmberg = Ehm+Zc*Imdss
        non_converged= (abs(Vmberg-Vmdss)>=eps)
        Vs = Vmberg-Zk*Imdss
        h=h+1
        print('solveLineOpenDSS:iteration:'+str(h))
        print(Vmberg)
        print(Vmdss)
    
    return Vs


def solveEhk(Vs,Ehm,t,h,Zc,Zk):

    if t>2*h:
        Vs_est = Vs*((t-h)**2)/(Vs*(t-2*h))
    else:
        Vs_est = Vs*(t-h)
    
    Vs = solveLineOpenDSS(Ehm,Vs_est,Zc,Zk)
    
    Im,Vm = DSS(Vs)
    Ehk = Vm + Zc*Im
    return Vs, Ehk
# ---------------------- co-simulation ----------------------
def run_cosim(Zc: complex, h: float, t_end: float):
    Yc = 1.0 / Zc
    # Yc_re = float(np.real(Yc))
    # Yc_im = float(np.imag(Yc))

    print(f"\nRunning co-sim with h = {h:.2e} s, "
          f"t_end = {t_end:.3f} s, "
          f"Zc = {abs(Zc):.4f} ∠ {np.rad2deg(np.angle(Zc)):.2f} deg")

    # load FMUs
    fmu_ts = load_fmu(TS_FMU, kind='CS')

    # instantiate & initialize
    fmu_ts.instantiate()

    fmu_ts.setup_experiment()
    

    fmu_ts.enter_initialization_mode()
    fmu_ts.set('pwFault2.t1',fault_ts_t1)
    fmu_ts.set('pwFault2.t2',fault_ts_t2)
    # set coupling admittance in both FMUs before exiting initialization
    # for fmu in (fmu_ts,):
    #     fmu.set('Yc_re', Yc_re)
    #     fmu.set('Yc_im', Yc_im)
    fmu_ts.time = 0

    # exit initialization so that FMUs reach their steady-state (load-flow) point
    fmu_ts.exit_initialization_mode()
    
    
    #if FMU_DS:
    fmu_ds = load_fmu(DS_FMU, kind='CS')
    fmu_ds.instantiate()
    fmu_ds.setup_experiment()
    fmu_ds.enter_initialization_mode()
    fmu_ds.set('pwFault2.t1',fault_ds_t1)
    fmu_ds.set('pwFault2.t2',fault_ds_t2)
    fmu_ds.time = 0
    fmu_ds.exit_initialization_mode()

    # time vector and storage
    ts = [0.0]

    # read initial steady-state V and I from both FMUs at t = 0
    V_k, I_k = get_VI(fmu_ts)


    # else: 
    #     V_m, I_m = get_VI_from_opendss(dss_ds)
    # Paper assumption: Vm = Vk and Im = -Ik (used only for Hin calculation)
    V_m = V_k
    I_m = -I_k

    Vts = [V_k]
    Vds = [V_m]

    # initial historical currents (Hin) as in the paper:
    # H_to_ds is computed from TS side quantities (V_k, I_k)
    # H_to_ts is computed from DS side quantities (V_m, I_m)
    #if FMU_DS:
    H_to_ds = Yc * V_k + I_k
    H_to_ts = Yc * V_m + I_m
    #else:
    Eh_to_ds = Zc * I_k + V_k*230000
        #Eh_to_ts = Zc * I_m + V_m 
    

    fmu_ts.set('Hin_re', float(np.real(H_to_ts)))
    fmu_ts.set('Hin_im', float(np.imag(H_to_ts)))
    #if FMU_DS:
    fmu_ds.set('Hin_re', float(np.real(H_to_ds)))
    fmu_ds.set('Hin_im', float(np.imag(H_to_ds)))
    # else:
    #     #INITIALIZATION MODE
    #     Vs = 1.0#V_k  # good starting guess
        
    
    tt = 0.0
    while tt < t_end - 1e-12:
        fmu_ts.do_step(current_t=tt, step_size=h, new_step=True)
        #if FMU_DS:
        fmu_ds.do_step(current_t=tt, step_size=h, new_step=True)
        # else:
        #     Vs, Ehk = solveEhk(Vs,Eh_to_ds,tt,h,Zc,Zk)

            


        # read updated PCC voltages and currents
        V_k, I_k = get_VI(fmu_ts)
        #if FMU_DS:
        V_m, I_m = get_VI(fmu_ds)

            
        Vts.append(V_k)
        Vds.append(V_m)

        # update historical currents for the next step
        H_to_ds = Yc * V_k + I_k
        #Eh_to_ds = -Zc * I_k + V_k
        #if FMU_DS:
        H_to_ts = Yc * V_m + I_m
        #else:
        #H_to_ts = (V_k-Ehk)/Zc

        fmu_ts.set('Hin_re', float(np.real(H_to_ts)))
        fmu_ts.set('Hin_im', float(np.imag(H_to_ts)))
        #if FMU_DS:
        fmu_ds.set('Hin_re', float(np.real(H_to_ds)))
        fmu_ds.set('Hin_im', float(np.imag(H_to_ds)))


        tt = round(tt + h, 12)
        ts.append(tt)

    fmu_ts.terminate()
    #if FMU_DS:
    fmu_ds.terminate()

    Vts_mag = np.abs(np.array(Vts))
    Vds_mag = np.abs(np.array(Vds))
    return ts, Vts_mag, Vds_mag


# ---------------------- Tk GUI ----------------------
def main():
    root = tk.Tk()
    root.title("Co-simulation GUI (Zc_eigenplotter + 2 sliders)")

    state = {
        "Zc": Zc_init,
        "h": H_DEFAULT,
        "t_end": T_END_DEFAULT,
    }

    main_frame = ttk.Frame(root)
    main_frame.pack(fill=tk.BOTH, expand=True)

    # one figure, two subplots
    fig = plt.Figure(figsize=(7, 5))
    gs = GridSpec(2, 1, height_ratios=[3, 1], figure=fig)
    
    ax_polar = fig.add_subplot(gs[0, 0], projection='polar')
    ax_time  = fig.add_subplot(gs[1, 0])
    
    # draw eigenplot INTO this axis (no new window)
    Zc_eigenplotter(Zc_init, Zk, Zm, r_max=r_max, ax=ax_polar)
    
    # marker for currently selected Zc
    theta0 = np.angle(Zc_init)
    r0 = np.abs(Zc_init)
    (marker,) = ax_polar.plot(theta0, r0, "ro", markersize=8)

    # time subplot
    ax_time.set_xlabel("Time [s]")
    ax_time.set_ylabel("Voltage magnitude [pu]")
    ax_time.set_title("PCC Voltage")
    ax_time.grid(True)

    (line_ts,) = ax_time.plot([], [], color='tab:blue',  label="|V| TS (PCC)")
    (line_ds,) = ax_time.plot([], [], color='tab:orange', label="|V| DS (PCC)")
    ax_time.legend()

    canvas = FigureCanvasTkAgg(fig, master=main_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # -------- two sliders stacked vertically --------
    slider_frame = ttk.Frame(main_frame, padding=5)
    slider_frame.pack(side=tk.BOTTOM, fill=tk.X)

    # --- t_end slider row (TOP) ---
    tend_row = ttk.Frame(slider_frame)
    tend_row.pack(side=tk.TOP, fill=tk.X, pady=2)

    ttk.Label(tend_row, text="Total time t_end [s]:").pack(side=tk.LEFT)

    t_end_var = tk.DoubleVar(value=T_END_DEFAULT)
    t_end_slider = tk.Scale(
        tend_row,
        variable=t_end_var,
        from_=T_END_MIN,
        to=T_END_MAX,
        resolution=T_END_RES,
        orient=tk.HORIZONTAL,
        showvalue=False,
        length=400,
    )
    t_end_slider.pack(side=tk.LEFT, padx=5)

    t_end_label = ttk.Label(tend_row, text=f"{T_END_DEFAULT:.2f}")
    t_end_label.pack(side=tk.LEFT, padx=5)

    # --- h slider row (BOTTOM) ---
    h_row = ttk.Frame(slider_frame)
    h_row.pack(side=tk.TOP, fill=tk.X, pady=2)

    ttk.Label(h_row, text="Timestep h [s]:").pack(side=tk.LEFT)

    h_var = tk.DoubleVar(value=H_DEFAULT)
    h_slider = tk.Scale(
        h_row,
        variable=h_var,
        from_=H_MIN,
        to=H_MAX,
        resolution=H_RES,
        orient=tk.HORIZONTAL,
        showvalue=False,
        length=400,
    )
    h_slider.pack(side=tk.LEFT, padx=5)

    h_label = ttk.Label(h_row, text=f"{H_DEFAULT:.1e}")
    h_label.pack(side=tk.LEFT, padx=5)
    # --------- update logic ---------
    def update_time_plot():
        Zc = state["Zc"]
        h = state["h"]
        t_end = state["t_end"]
        t, Vts_mag, Vds_mag = run_cosim(Zc, h, t_end)
        line_ts.set_data(t, Vts_mag)
        line_ds.set_data(t, Vds_mag)
        ax_time.relim()
        ax_time.autoscale_view()
        canvas.draw()

    # initial run
    update_time_plot()

    # --------- click handler for polar plot ---------
    def on_polar_click(event):
        if event.inaxes != ax_polar:
            return
        if event.xdata is None or event.ydata is None:
            return
        theta_click = event.xdata
        r_click = event.ydata
        state["Zc"] = r_click * np.exp(1j * theta_click)
        marker.set_data([theta_click], [r_click])
        canvas.draw()
        update_time_plot()

    canvas.mpl_connect("button_press_event", on_polar_click)

    # slider callbacks
    def on_h_move(val):
        h_label.config(text=f"{float(val):.1e}")

    def on_h_release(event):
        state["h"] = h_var.get()
        update_time_plot()

    def on_tend_move(val):
        t_end_label.config(text=f"{float(val):.2f}")

    def on_tend_release(event):
        state["t_end"] = t_end_var.get()
        update_time_plot()

    h_slider.config(command=on_h_move)
    h_slider.bind("<ButtonRelease-1>", on_h_release)

    t_end_slider.config(command=on_tend_move)
    t_end_slider.bind("<ButtonRelease-1>", on_tend_release)

    root.mainloop()


if __name__ == "__main__":
    main()
