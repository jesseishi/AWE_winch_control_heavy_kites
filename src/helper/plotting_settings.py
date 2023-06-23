# Linewidth in latex is 15.75 cm, which is 15.75/2.54 inch. I make it a bit larger which
# makes the text a bit smaller, but the figure a little bigger.
FULLSIZE = (15.75 / 2.54 + 3, 4)
PARTSIZE = (0.7 * (15.75 / 2.54 + 3), 4)

# TODO: Can do more customisation for how the plots look like if needed.
# From awe_workshop.
# # Configure the default plotting settings
# size=13
# params = {'legend.fontsize': 'large',
#           'figure.figsize': (8,5),
#           'axes.labelsize': size,
#           'axes.titlesize': size,
#           'xtick.labelsize': size*0.85,
#           'ytick.labelsize': size*0.85,
#           'axes.titlepad': 25}
# plt.rcParams.update(params)

# Seaborn sets the x and y axis labels to the column names, but we want them to be a bit
# more human-readable.
var_to_label = {
    "vw_mps": "wind speed [m/s]",
    "Lt_m": "tether length [m]",
    "phi_deg": "azimuth [deg]",
    "beta_deg": "elevation [deg]",
    "chi_deg": "course angle [deg]",
    "Ftk_N_star": "tether force kite [N]",
    "f_star": "reel-out factor [-]",
    "Ftg_N_star": "tether force [N]",
    "Ftg_N": "tether force [N]",
    "P_W_star": "power [W]",
    "P_W": "power [W]",
    "E_eff": "effective lift-to-drag ratio [-]",
    "vr_mps_star": "reel-out speed [m/s]",
    "vr_mps": "reel-out speed [m/s]",
    "f": "reel-out factor [-]",
    "": "",
}


def set_labels(ax):
    # If ax is an iterable, loop over all the individual axis.
    if hasattr(ax, "__iter__"):
        for ax_i in ax:
            set_labels(ax_i)

    else:
        ax.set_xlabel(var_to_label[ax.get_xlabel()])
        ax.set_ylabel(var_to_label[ax.get_ylabel()])

        # Also the legend.
        leg = ax.get_legend()
        if leg:
            if leg.get_visible():
                try:
                    leg.set_title(var_to_label[leg.get_title().get_text()])

                # When hue and style are set, Seaborn doesn't set the title but just puts
                # multiple texts in the legend. With some digging into the source code we
                # can still fit it.
                except Exception:
                    pass
                finally:
                    for T in leg.texts:
                        if T._text in var_to_label.keys():
                            T._text = var_to_label[T._text]
