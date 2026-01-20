import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge
from helpers import f_of_Zc


def Zc_eigenplotter(Zc, Zk, Zm, r_max=0.5, ax=None):
    """
    Plot f(Zc, theta_c) in the Zc-plane.

    If ax is None  -> create a new figure & polar axis (stand-alone use).
    If ax is given -> draw into that axis (for embedding in Tkinter, etc.).
    """

    # ----- grid over Zc -----
    r = np.linspace(0.0, r_max, 300)
    th = np.linspace(0.0, 2*np.pi, 1080, endpoint=False)
    TH, R = np.meshgrid(th, r)
    Zc_grid = R * np.exp(1j * TH)

    F = f_of_Zc(Zc_grid, Zk, Zm)

    vmin = 0.0
    vmax = 1.0
    F_clip = np.clip(F, vmin, vmax)
    F_plot = F_clip
    cb_label = r'$f(Z_c,\theta_c)$'
    levels = np.linspace(vmin, vmax, 16)

    # ----- figure / axis handling -----
    created_fig = False
    if ax is None:
        fig = plt.figure(figsize=(7.5, 7.5))
        ax = plt.subplot(111, projection='polar')
        created_fig = True
    else:
        fig = ax.figure

    # polar layout
    ax.set_theta_zero_location('E')
    ax.set_theta_direction(1)

    # filled contours
    cf = ax.contourf(TH, R, F_plot, levels=levels, extend='both')

    # colorbar attached to this figure
    cb = fig.colorbar(cf, ax=ax, pad=0.1)
    cb.set_label(cb_label)

    # threshold contour
    thr = 0.1
    ax.contour(TH, R, F, levels=[thr], colors='k', linewidths=1.8)

    # axis limits / ticks
    ax.set_rmax(r_max)
    ax.set_rticks([0.02, 0.05, 0.08, 0.10])
    ax.set_title(r'Polar contour of $f(Z_c,\theta_c)$', pad=20)

    # optional wedge
    add_wedge = False
    if add_wedge:
        w = Wedge(center=(0, 0), r=r_max,
                  theta1=300, theta2=360, width=r_max,
                  transform=ax.transData._b, alpha=0.08)
        ax.add_patch(w)

    # overlay markers
    ax.plot(np.angle(Zk), np.abs(Zk), '*', ms=7, label='$Z_k$', color='blue')
    ax.plot(np.angle(Zm), np.abs(Zm), '*', ms=7, label='$Z_m$', color='violet')
    ax.plot(np.angle(Zc), np.abs(Zc), '*', ms=7, label='$Z_c$', color='red')

    # legend OUTSIDE the polar disk so it doesn't overlap
    ax.legend(
        loc='center right',
        bbox_to_anchor=(-0.15, 0.5),  # move it to the left of the plot
        framealpha=0.9
    )

    if created_fig:
        fig.tight_layout()
        plt.show()

    return fig, ax
