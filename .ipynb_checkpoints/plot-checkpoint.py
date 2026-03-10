from libraries import *

plt.rcParams['font.family']  = 'DeJavu Serif'
plt.rcParams['font.serif']   = ['Times New Roman']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

df = pd.read_csv("dataset/balloon_sim_data.csv")

R  = 6371  # Earth's radius in km

params1 = {
    # Optical
    "wavelength": 1550e-9,     # [m]
    "W0": 0.15,                # Initial beam waist [m] (15 cm)
    
    # Apertures
    "rx_aperture_down": 0.75,  # Ground station [m]
    "rx_aperture_up": 0.40,    # HAP [m]
    "obs_ratio": 0.3,

    # Atmosphere & turbulence
    "Cn0": 1e-17,         # [m^(-2/3)] "Cn0": 9.6e-14,            # [m^(-2/3)]
    "u_rms": 1,           # [m/s]

    # Pointing & tracking
    "pointing_error": 1e-7,    # 1 μrad (paper)
    "tracking_efficiency": 0.85,

    # Detection
    "detector_eff": 0.85,

    # Netsquid
    "init_time": 1
}

params2 = {
    # Optical
    "wavelength": 1550e-9,     # [m]
    "W0": 0.15,                # Initial beam waist [m] (15 cm)
    
    # Apertures
    "rx_aperture_down": 0.75,  # Ground station [m]
    "rx_aperture_up": 0.40,    # HAP [m]
    "obs_ratio": 0.3,

    # Atmosphere & turbulence
    "Cn0": 1e-13,            # [m^(-2/3)] "Cn0": 9.6e-14,            # [m^(-2/3)]
    "u_rms": 10,             # [m/s]

    # Pointing & tracking
    "pointing_error": 1e-7,    # 1 μrad (paper)
    "tracking_efficiency": 0.85,

    # Detection
    "detector_eff": 0.85,

    # Netsquid
    "init_time": 1
}

def clip_eta(eta, eps=1e-16):
    """
    Enforce physical transmissivity bounds.
    """
    return min(eta, 1.0 - eps)


def plob_skr_kbps(eta, ts):
    """
    PLOB SKR in kbps (direct formula, no compute_skr).
    """
    eta = clip_eta(eta)
    return (
        -ts.ratesources
        * ts.sourceeff
        * math.log1p(-eta) / math.log(2)
    )

def geometric_eta(distance_km, ts):
    """
    Ideal geometric (diffraction-only) transmissivity.
    Assumes circular apertures and Airy divergence.
    """

    D_tx = ts.params["rx_aperture_up"]     # TX aperture (HAP)
    D_rx = ts.params["rx_aperture_down"]   # RX aperture (GS)
    lam  = ts.params["wavelength"]

    theta = 1.22 * lam / D_tx              # diffraction half-angle
    beam_diameter = theta * distance_km * 1e3

    eta_geo = (D_rx / (D_tx + beam_diameter))**2

    print(f"eta_geo: {eta_geo}")
    return clip_eta(eta_geo)

# def plot_skr_distance(n=6, d_min=25, d_max=250):
#     """
#     Plot skr over different LoS distances
#     df : pandas.DataFrame with ["Time_s", "Longitude_deg", "Latitude_deg", "Altitude_m"]
#     ts : transmittance simulation module/object with .theoretical_eff, .simulated_eff, .compute_skr
#     """

#     distances = [d for d in range(d_min, d_max)]
#     altitude  = d_min
    
#     skr_theory1     = []
#     # skr_sim1         = []
#     # skr_sim_plob1    = []
#     eta_theory1     = []
#     # eta_sim1         = []

#     skr_theory2     = []
#     # skr_sim2         = []
#     # skr_sim_plob2    = []
#     eta_theory2     = []
#     # eta_sim2         = []

#     eta_fiber1      = []
#     eta_fiber2      = []
#     skr_fiber1      = []
#     skr_fiber2      = []

#     eta_ideal       = []

#     skr_ideal  = []
#     # Compute SKRs
#     for idx_d, d in enumerate(distances):
#         #eta_t, eta_s, skr_t, skr_s, skr_pt, skr_ps = compute_point(d, 15, "downlink", 10)
#         eta_t1, skr_t1, skr_pt1 = compute_point_wsim(d, altitude, "downlink", 10, params1)
#         eta_t2, skr_t2, skr_pt2 = compute_point_wsim(d, altitude, "downlink", 10, params2)
        
#         #eta_t = ts.channel_theory("downlink", 0, 15, d, n)
#         #eta_s = ts.channel_simulation("downlink", 0, 15, d, n)

#         eta_theory1.append(eta_t1 * 100)
#         #eta_sim1.append(eta_s1 * 100)
#         skr_theory1.append(skr_t1 / 1000)
#         #skr_sim1.append(skr_s1 / 1000)

#         eta_theory2.append(eta_t2 * 100)
#         #eta_sim2.append(eta_s2 * 100)
#         skr_theory2.append(skr_t2 / 1000)
#         #skr_sim2.append(skr_s2 / 1000)
        
#         #skr_sim_plob.append(-ts.ratesources * ts.sourceeff * math.log2(1 - eta_s) / 1000)

#         fiber_alpha1 = 0.2  # dB/km
#         fiber_alpha2 = 0.35 # dB/km
#         fiber_loss1  = fiber_alpha1 * d + 3 + 0.15 * 2 + 0.5 * 2
#         fiber_loss2  = fiber_alpha2 * d + 3 + 0.3 * 2 + 0.75 * 2
#         fiber_eta1    = 10 ** (-fiber_loss1 / 10)
#         fiber_eta2    = 10 ** (-fiber_loss2 / 10)

#         eta_fiber1.append(fiber_eta1 * 100)
#         eta_fiber2.append(fiber_eta2 * 100)

#         skr_fiber1.append(ts.compute_skr(fiber_eta1) / 1000) #-ts.ratesources * ts.sourceeff * math.log1p(-fiber_eta) / math.log(2) / 1000)
#         skr_fiber2.append(ts.compute_skr(fiber_eta2) / 1000) #-ts.ratesources * ts.sourceeff * math.log1p(-fiber_eta) / math.log(2) / 1000)

#         #####################
#         # Elevation angle
#         R_RX = ts.params["rx_aperture_down"]
#         R_TX = ts.params["rx_aperture_up"]
#         W    = ts.params["wavelength"]
#         THETA = 1.22 * W / R_TX

#         # Losses
#         L_geo = 10 * math.log10((R_TX + d * 1000 * THETA)**2 / R_RX**2)
#         L_t = L_geo

#         eta_geo = 10 ** (-L_t / 10)

#         eta_ideal.append(eta_geo * 100)
#         skr_ideal.append(ts.compute_skr(eta_geo) / 1000) #-ts.ratesources * ts.sourceeff * math.log1p(-ETA) / math.log(2) / 1000)

#         # # Apertures (diameters in meters)
#         # D_TX = ts.params["rx_aperture_down"]   # 0.75 m
#         # D_RX = ts.params["rx_aperture_up"]     # 0.40 m
        
#         # # Wavelength
#         # lam = ts.params["wavelength"]
        
#         # # Diffraction divergence (half-angle)
#         # theta = 1.22 * lam / D_TX
        
#         # # Beam radius at receiver
#         # w = theta * d * 1000  # meters

#         # print(f"eta_geo: {10 * math.log10((D_TX + d * 1000 * theta)**2 / D_RX**2)}")
        
#         # # Geometric coupling efficiency
#         # eta_geo = min((D_RX / (2 * w))**2, 1.0)
        
#         # eta_ideal.append(eta_geo * 100)
#         # skr_ideal.append(ts.compute_skr(eta_geo) / 1000)
#         #####################

#         print(f"{idx_d} -> d: {d}, eta_t1: {eta_t1}, eta_t2: {eta_t2}, eta_i: {eta_geo}, skr_t1: {skr_theory1[-1]}, skr_t2: {skr_theory2[-1]}, skr_f1: {skr_fiber1[-1]}, skr_f2: {skr_fiber2[-1]}, skr_ip: {skr_ideal[-1]}")

#     # Plot
#     # plt.figure(figsize=(3,3))
#     # plt.plot(distances, eta_theory, label="HAP", color="blue", marker="o", markevery=20)
#     # #plt.plot(distances, eta_sim, label="Simulation", color="orange", marker="^", markevery=20)

#     # plt.plot(distances, eta_fiber1, linestyle="--", color="red")
#     # plt.plot(distances, eta_fiber2, linestyle="--", color="red")
#     # plt.plot(distances, eta_ideal, label="HAP (Geo)", linestyle="--", color="green")

#     # plt.fill_between(
#     #     distances,
#     #     eta_fiber1,
#     #     eta_fiber2,
#     #     color="red",
#     #     alpha=0.2,
#     #     label="Fiber"
#     # )

#     # plt.xlabel("Distance (km)", fontsize=12)
#     # plt.ylabel("Channel Efficiency (%)", fontsize=12)
#     # plt.grid(True)
#     # plt.legend(fontsize=12)
#     # plt.yscale("log")
#     # plt.savefig("trans_distance.svg", format="svg", dpi=300, bbox_inches="tight")
#     # plt.show()
    
#     # Plot
#     plt.figure(figsize=(3,3))
#     plt.plot(distances, skr_theory1, color="blue")
#     plt.plot(distances, skr_theory2, color="blue")
#     # plt.plot(distances, skr_sim, label="Simulation (DW)", color="orange", marker="^", markevery=20)
#     #plt.plot(distances, skr_theory_plob, label="HAP (PLOB)", color="green", marker="h", markevery=20)
#     # plt.plot(distances, skr_sim_plob, label="Simulation (PLOB)", color="red", marker="*", markevery=20)
#     plt.plot(distances, skr_fiber1, linestyle="--", color="red")
#     plt.plot(distances, skr_fiber2, linestyle="--", color="red")
#     plt.plot(distances, skr_ideal, label="HAP (Geo)", linestyle="--", color="green")

#     plt.fill_between(
#         distances,
#         skr_theory1,
#         skr_theory2,
#         color="blue",
#         alpha=0.2,
#         label="HAP"
#     )
    
#     plt.fill_between(
#         distances,
#         skr_fiber1,
#         skr_fiber2,
#         color="red",
#         alpha=0.2,
#         label="Fiber"
#     )

#     plt.xlabel("Distance (km)", fontsize=12)
#     plt.ylabel("Maximum SKR (Kbps)", fontsize=12)
#     plt.grid(True)
#     plt.legend(fontsize=12, borderpad=0.2, handletextpad=0.1)
#     plt.yscale("log")
#     plt.savefig("skr_distance.svg", format="svg", dpi=300, bbox_inches="tight")
#     plt.show()

def plot_skr_distance(d_min=25, d_max=250):
    """
    Plot SKR vs distance with:
      - DW bounds (HAP, Fiber)
      - PLOB bounds (HAP, Fiber)
      - Ideal geometric-loss-only HAP curve
    """

    distances_km = list(range(d_min, d_max))
    hap_altitude_km = d_min

    # =========================
    # Storage
    # =========================

    # --- DW bounds ---
    hap_dw_low, hap_dw_high = [], []
    fiber_dw_low, fiber_dw_high = [], []

    # --- PLOB bounds ---
    hap_plob_low, hap_plob_high = [], []
    fiber_plob_low, fiber_plob_high = [], []

    # --- Ideal geometric ---
    hap_geo_dw = []
    hap_geo_plob = []

    # =========================
    # Loop
    # =========================

    for d in distances_km:

        # ---------- HAP (DW, two parameter sets) ----------
        eta_hap_1, skr_dw_1, _ = compute_point_wsim(
            d, hap_altitude_km, "downlink", 10, params1
        )
        eta_hap_2, skr_dw_2, _ = compute_point_wsim(
            d, hap_altitude_km, "downlink", 10, params2
        )

        eta_hap_1 = clip_eta(eta_hap_1)
        eta_hap_2 = clip_eta(eta_hap_2)

        hap_dw_low.append(min(skr_dw_1, skr_dw_2))
        hap_dw_high.append(max(skr_dw_1, skr_dw_2))

        hap_plob_low.append(
            min(
                plob_skr_kbps(eta_hap_1, ts),
                plob_skr_kbps(eta_hap_2, ts)
            )
        )
        hap_plob_high.append(
            max(
                plob_skr_kbps(eta_hap_1, ts),
                plob_skr_kbps(eta_hap_2, ts)
            )
        )

        # ---------- Fiber ----------
        alpha_low  = 0.16  # dB/km
        alpha_high = 0.2   # dB/km

        loss_low  = alpha_low  * math.sqrt(d**2 - 25**2) + 3 + 0.15 * 2 + 0.50 * 2
        loss_high = alpha_high * math.sqrt(d**2 - 25**2) + 3 + 0.30 * 2 + 0.75 * 2

        eta_fiber_low  = clip_eta(10 ** (-loss_low / 10))
        eta_fiber_high = clip_eta(10 ** (-loss_high / 10))

        fiber_dw_low.append(ts.compute_skr(eta_fiber_low))
        fiber_dw_high.append(ts.compute_skr(eta_fiber_high))

        fiber_plob_low.append(plob_skr_kbps(eta_fiber_low, ts))
        fiber_plob_high.append(plob_skr_kbps(eta_fiber_high, ts))

        # ---------- Ideal geometric (HAP) ----------
        eta_geo = geometric_eta(d, ts)

        hap_geo_dw.append(ts.compute_skr(eta_geo))

    # =========================
    # Plot 1: DW bounds
    # =========================

    plt.figure(figsize=(3.2, 3.2))

    # plt.plot(
    #     distances_km,
    #     hap_dw_low,
    #     linestyle="-",
    #     color="blue",
    # )

    # plt.plot(
    #     distances_km,
    #     hap_dw_high,
    #     linestyle="-",
    #     color="blue",
    # )

    plt.plot(
        distances_km,
        hap_geo_dw,
        linestyle="-",
        color="black",
        marker="^", 
        markevery=25,
        label="HAP-Geo"
    )

    plt.fill_between(
        distances_km,
        fiber_dw_low,
        fiber_dw_high,
        facecolor="lightcoral",     # important
        edgecolor="black",
        hatch="o",
        alpha=0.8,
        label="Fiber"
    )
    plt.fill_between(
        distances_km,
        hap_dw_low,
        hap_dw_high,
        facecolor="lightblue",     # important
        edgecolor="black",
        hatch=".",
        alpha=0.8,
        label="HAP"
    )

    # plt.plot(
    #     distances_km,
    #     fiber_dw_low,
    #     linestyle="-",
    #     color="red",
    # )

    # plt.plot(
    #     distances_km,
    #     fiber_dw_high,
    #     linestyle="-",
    #     color="red",
    # )

    plt.yscale("log")
    plt.xlabel("LoS Distance (km)")
    plt.ylabel("SKR DW Bound (bps)")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig("skr_distance_dw.pdf", format='pdf', transparent=True, bbox_inches='tight', pad_inches=0.01)
    plt.savefig("skr_distance_dw.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()

    # =========================
    # Plot 2: PLOB bounds
    # =========================

    plt.figure(figsize=(3.2, 3.2))

    # plt.plot(
    #     distances_km,
    #     hap_plob_low,
    #     linestyle="-",
    #     color="blue",
    # )

    # plt.plot(
    #     distances_km,
    #     hap_plob_high,
    #     linestyle="-",
    #     color="blue",
    # )
    
    

    plt.plot(
        distances_km,
        hap_geo_plob,
        linestyle="-",
        color="black",
        marker="^", 
        markevery=25,
        label="HAP-Geo"
    )
    plt.fill_between(
        distances_km,
        fiber_plob_low,
        fiber_plob_high,
        facecolor="lightcoral",     # important
        edgecolor="black",
        hatch="o",
        alpha=0.8,
        label="Fiber"
    )
    plt.fill_between(
        distances_km,
        hap_plob_low,
        hap_plob_high,
        facecolor="lightblue",     # important
        edgecolor="black",
        hatch=".",
        alpha=0.8,
        label="HAP"
    )
    

    # plt.plot(
    #     distances_km,
    #     fiber_plob_low,
    #     linestyle="-",
    #     color="red",
    # )

    # plt.plot(
    #     distances_km,
    #     fiber_plob_high,
    #     linestyle="-",
    #     color="red",
    # )

    plt.yscale("log")
    plt.xlabel("LoS Distance (km)")
    plt.ylabel("SKR PLOB Bound (bps)")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig("skr_distance_plob.pdf", format='pdf', transparent=True, bbox_inches='tight', pad_inches=0.01)
    plt.savefig("skr_distance_plob.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()

def plot_connectivity_graph_3d(
    gnodes,
    hnodes,
    links,
    t=0,                 # time index for node positions
    cyl_radius=0.65,      # cylinder radius (lon/lat units)
    cyl_alpha=0.05       # transparency
):
    """
    3D connectivity graph with GS, HAP trajectories, links,
    and a transparent vertical cylinder.
    """

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, projection="3d")

    # ==========================
    # Plot Ground Stations (GS)
    # ==========================
    for gs in gnodes:
        ax.scatter(
            gs.lg, gs.la, 0,
            color="skyblue", marker="^", s=120, zorder=5
        )
        if hasattr(gs, "tag"):
            ax.text(gs.lg - 0.05, gs.la - 0.05, 0, gs.tag, fontsize=14)

    # ==========================
    # Plot HAP initial positions
    # ==========================
    for h in hnodes:
        ax.scatter(
            h.lg[t], h.la[t], h.H[t],
            color="red", s=10, zorder=6
        )
        if hasattr(h, "tag"):
            ax.text(h.lg[t] - 0.05, h.la[t] - 0.05, h.H[t], h.tag, fontsize=14)

    # ==========================
    # Plot links (no duplicates)
    # ==========================
    plotted_edges = set()
    for l in links:
        edge_key = frozenset([l.n1, l.n2])
        if edge_key in plotted_edges:
            continue
        plotted_edges.add(edge_key)

        def node_xyz(n):
            if isinstance(n.lg, list):
                return n.lg[t], n.la[t], n.H[t]
            else:
                return n.lg, n.la, 0

        x1, y1, z1 = node_xyz(l.n1)
        x2, y2, z2 = node_xyz(l.n2)

        ax.plot(
            [x1, x2], [y1, y2], [z1, z2],
            color="gray", linestyle="--", linewidth=0.6, alpha=0.6
        )

    # ==========================
    # Plot HAP trajectories
    # ==========================
    for h in hnodes:
        ax.plot(
            h.lg, h.la, h.H,
            color="red", linewidth=1, alpha=0.8
        )

    # ==========================
    # Transparent cylinder
    # ==========================
    # Cylinder center (mean lon/lat)
    all_lons = [gs.lg for gs in gnodes] + [h.lg[t] for h in hnodes]
    all_lats = [gs.la for gs in gnodes] + [h.la[t] for h in hnodes]
    z_max = max(max(h.H) for h in hnodes)

    cx = np.mean(all_lons)
    cy = np.mean(all_lats)

    theta = np.linspace(0, 2*np.pi, 60)
    z = np.linspace(0, z_max, 40)
    theta, z = np.meshgrid(theta, z)

    x_cyl = cx + cyl_radius * np.cos(theta)
    y_cyl = cy + cyl_radius * np.sin(theta)

    ax.plot_surface(
        x_cyl, y_cyl, z,
        color="blue",
        alpha=cyl_alpha,
        linewidth=0,
        shade=False
    )

    # ==========================
    # Axis labels and limits
    # ==========================
    ax.set_xlabel("Longitude", fontsize=14)
    ax.set_ylabel("Latitude", fontsize=14)
    ax.set_zlabel("Altitude (km)", fontsize=14)

    ax.set_xlim(min(all_lons) - cyl_radius, max(all_lons) + cyl_radius)
    ax.set_ylim(min(all_lats) - cyl_radius, max(all_lats) + cyl_radius)
    ax.set_zlim(0, z_max * 1.05)

    # Longitude / Latitude ticks every 0.5 deg
    x_ticks = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 4)
    y_ticks = np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 4)
    
    # Altitude ticks every ~25%
    z_ticks = np.linspace(0, z_max, 4)
    
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)

    ax.tick_params(labelsize=10)

    # ==========================
    # View angle (important)
    # ==========================
    ax.view_init(elev=15, azim=135)

    plt.tight_layout()
    plt.savefig("hap_qkd_trajectory_3d.pdf", format="pdf", bbox_inches="tight", pad_inches=0.01, transparent="True")
    plt.show()


def plot_skr_eta_stratotegic_real(n=6, d_min_t=0, d_max_t=86400):
    """
    Plot theoretical and simulated SKR for a balloon trajectory over time.
    df : pandas.DataFrame with ["Time_s", "Longitude_deg", "Latitude_deg", "Altitude_m"]
    ts : transmittance simulation module/object with .theoretical_eff, .simulated_eff, .compute_skr
    """
    df_filtered = (
        df[df["Time_s"] % 300 == 0]
        .drop_duplicates(subset=["Time_s"], keep="first")
        .reset_index(drop=True)
    )

    times = df_filtered["Time_s"].values / 3600
    mask = (times >= d_min_t) & (times <= d_max_t)
    times = times[mask]

    results = []  # list of tuples (d, eta_t, eta_s, skr_t, skr_s)

    la_rad_g = math.radians(49)
    lg_rad_g = math.radians(279)
    x_g = R * math.cos(la_rad_g) * math.cos(lg_rad_g)
    y_g = R * math.cos(la_rad_g) * math.sin(lg_rad_g)

    distances = []
    altitudes = []
    for _, row in df_filtered[mask].iterrows():
        lon = row["Longitude_deg"]
        lat = row["Latitude_deg"]
        alt = row["Altitude_m"] / 1000  # km

        # Balloon coordinates
        la_rad_h = math.radians(lat)
        lg_rad_h = math.radians(lon)
        x_h = R * math.cos(la_rad_h) * math.cos(lg_rad_h)
        y_h = R * math.cos(la_rad_h) * math.sin(lg_rad_h)

        # Horizontal distance
        d_los_hor = math.sqrt((x_h - x_g) ** 2 + (y_h - y_g) ** 2)

        # Elevation angle
        alpha = math.atan(alt / d_los_hor) if d_los_hor > 0 else math.pi / 2

        # LOS distance
        d_los = alt / math.sin(alpha)

        altitudes.append(alt)
        distances.append(d_los)

    skr_theory      = []
    #skr_sim         = []
    skr_theory_plob = []
    #skr_sim_plob    = []
    eta_theory      = []
    #eta_sim         = []
    if os.path.exists("skr_strato_data.pkl"):
        with open("skr_strato_data.pkl", "rb") as f:
            results = pickle.load(f)
        print(f"[INFO] Loaded checkpoint")

        # Iterate and access each entry
        for idx, entry in enumerate(results):
            d     = entry["d"]
            eta_t = entry["eta_t"]
            #eta_s = entry["eta_s"]
            skr_t = entry["skr_t"]
            #skr_s = entry["skr_s"]

            print(
                f"{idx} -> d: {d}, eta_t: {eta_t},"
                f"skr_t: {skr_t}"
            )

            #distances.append(d)
            eta_theory.append(eta_t)
            #eta_sim.append(eta_s)
            skr_theory.append(skr_t)
            #skr_sim.append(skr_s)

            skr_theory_plob.append(-ts.ratesources * ts.sourceeff * math.log2(1 - eta_t) / 1e6)
            #skr_sim_plob.append(-ts.ratesources * ts.sourceeff * math.log2(1 - eta_s) / 1000)
    else:
        # Compute SKRs
        print(f"{len(altitudes)}, {len(distances)}")
        for idx_d, d in enumerate(distances):
            eta_t, skr_t, skr_pt = compute_point(d, altitudes[idx_d], "downlink", 10)
            
            #eta_t = ts.channel_theory("downlink", 0, altitudes[idx_d], d, n)
            #eta_s = ts.channel_simulation("downlink", 0, altitudes[idx_d], d, n)

            eta_theory.append(eta_t * 100)
            #eta_sim.append(eta_s * 100)
            skr_theory.append(skr_t / 1e6)
            #skr_sim.append(skr_s / 1000)
            
            print(f"{idx_d} -> d: {d}, eta_t: {eta_t}, skr_t: {skr_theory[-1]}")

            skr_theory_plob.append(-ts.ratesources * ts.sourceeff * math.log2(1 - eta_t) / 1e6)
            #skr_sim_plob.append(-ts.ratesources * ts.sourceeff * math.log2(1 - eta_s) / 1000)

            results.append({
                "d": d,
                "eta_t": eta_t * 100,
                #"eta_s": eta_s,
                "skr_t": skr_t / 1e6,
                #"skr_s": skr_s,
            })
    
        with open("skr_strato_data.pkl", "wb") as f:
            pickle.dump(results, f)

    # Plot
    # plt.figure(figsize=(7,3))
    # plt.plot(times, eta_theory, label="Theory", color="blue", marker="o", markevery=20)
    # plt.plot(times, eta_sim, label="Simulation", color="orange", marker="^", markevery=20)

    # plt.xlabel("Time (s)", fontsize=12)
    # plt.ylabel("Transmittance (%)", fontsize=12)
    # plt.grid(True)
    # plt.legend(fontsize=12)
    # plt.savefig("trans_strato.svg", format="svg", dpi=300, bbox_inches="tight")
    # plt.show()
    
    # Plot
    plt.figure(figsize=(3,2))
    plt.plot(times, skr_theory, label="DW", color="blue", marker="o", markevery=40)
    #plt.plot(times, skr_sim, label="Simulation (DW)", color="orange", marker="^", markevery=20)
    plt.plot(times, skr_theory_plob, label="PLOB", color="green", marker="h", markevery=40)
    #plt.plot(times, skr_sim_plob, label="Simulation (PLOB)", color="red", marker="*", markevery=20)

    plt.xlabel("Time (hours)", fontsize=12)
    plt.ylabel("SKR (Mbps)", fontsize=12)
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.savefig("skr_strato.pdf", format='pdf', transparent=True, bbox_inches='tight', pad_inches=0.01)
    plt.show()

    # Plot
    plt.figure(figsize=(3,2))
    plt.plot(times, distances, color="purple", marker="o", markevery=40)

    plt.xlabel("Time (hours)", fontsize=12)
    plt.ylabel("Distance (km)", fontsize=12)
    plt.grid(True)
    plt.savefig("distance_strato.pdf", format='pdf', transparent=True, bbox_inches='tight', pad_inches=0.01)
    plt.show()

def plot_solution(solution):
    """
    Plot solution variables as time series (2D and 3D depending on dimension).
    Handles both numpy array style (solution[var] = np.array) 
    and dict style (solution[var][(i,t)] = value).
    """

    for var, data in solution.items():
        print(f"Plotting {var}...")

        # --- Case 1: Data is numpy array ---
        if isinstance(data, np.ndarray):
            if data.ndim == 1:   # shape (t,)
                plt.figure(figsize=(10,5))
                plt.plot(range(len(data)), data, label=f"{var}")
                plt.xlabel("Time")
                plt.ylabel(var)
                plt.title(f"{var} vs Time")
                plt.legend()
                plt.show()

            elif data.ndim == 2:   # shape (i, t)
                plt.figure(figsize=(10,5))
                for i in range(data.shape[0]):
                    plt.plot(range(data.shape[1]), data[i], label=f"{var}[{i}]")
                plt.xlabel("Time")
                plt.ylabel(var)
                plt.title(f"{var} vs Time")
                plt.legend()
                plt.show()

            elif data.ndim == 3:   # shape (i, j, t)
                fig = plt.figure(figsize=(8,6))
                ax = fig.add_subplot(111, projection="3d")
                for i in range(data.shape[0]):
                    for j in range(data.shape[1]):
                        ax.plot(range(data.shape[2]), [i]*data.shape[2], data[i,j], label=f"{var}[{i},{j}]")
                ax.set_xlabel("Time")
                ax.set_ylabel("i")
                ax.set_zlabel(var)
                plt.title(f"{var} vs Time (3D)")
                plt.show()

        # --- Case 2: Data is dict {(i,...,t): value} ---
        elif isinstance(data, dict):
            # Extract keys
            
            keys = list(data.keys())
            if not keys:
                continue

            # Detect dimension from tuple length
            key_len = len(keys[0])

            if key_len == 1:  # (t)
                times = sorted(k[0] for k in keys)
                vals = [data[(t,)] for t in times]
                plt.figure(figsize=(10,5))
                plt.plot(times, vals, label=var)
                plt.xlabel("Time")
                plt.ylabel(var)
                plt.title(f"{var} vs Time")
                plt.legend()
                plt.show()

            elif key_len == 2:  # (i,t)
                grouped = {}
                for (i,t), val in data.items():
                    grouped.setdefault(i, {})[t] = val
                plt.figure(figsize=(10,5))
                for i, tvals in grouped.items():
                    times = sorted(tvals.keys())
                    vals = [tvals[t] for t in times]
                    plt.plot(times, vals, label=f"{var}[{i}]")
                plt.xlabel("Time")
                plt.ylabel(var)
                plt.title(f"{var} vs Time")
                plt.legend()
                plt.show()

            elif key_len == 3:  # (i,j,t)
                grouped = {}
                for (i,j,t), val in data.items():
                    grouped.setdefault((i,j), {})[t] = val
                fig = plt.figure(figsize=(8,6))
                ax = fig.add_subplot(111, projection="3d")
                for (i,j), tvals in grouped.items():
                    times = sorted(tvals.keys())
                    vals = [tvals[t] for t in times]
                    ax.plot(times, [i]*len(times), vals, label=f"{var}[{i},{j}]")
                ax.set_xlabel("Time")
                ax.set_ylabel("i")
                ax.set_zlabel(var)
                plt.title(f"{var} vs Time (3D)")
                plt.show()

def print_solution(solution):
    if solution is None:
        print("No solution to display.")
        return
    
    for name, values in solution.items():
        print(f"\n-- {name} --")
        
        # Case 1: scalar (float, int, etc.)
        if np.isscalar(values):
            if values > 1e-6:
                print(f"{name} = {values}")
        
        # Case 2: numpy array
        elif isinstance(values, np.ndarray):
            if values.ndim == 3:
                for i in range(values.shape[0]):
                    for j in range(values.shape[1]):
                        for k in range(values.shape[2]):
                            if values[i, j, k] > 1e-6:
                                print(f"{name}[{i},{j},{k}] = {values[i, j, k]}")
            elif values.ndim == 2:
                for i in range(values.shape[0]):
                    for j in range(values.shape[1]):
                        if values[i, j] > 1e-6:
                            print(f"{name}[{i},{j}] = {values[i, j]}")
            elif values.ndim == 1:
                for i in range(len(values)):
                    if values[i] > 1e-6:
                        print(f"{name}[{i}] = {values[i]}")
        
        # Case 3: dict (recursively handle)
        elif isinstance(values, dict):
            for k, v in values.items():
                if np.isscalar(v):
                    if v > 1e-6:
                        print(f"{name}[{k}] = {v}")
                elif isinstance(v, np.ndarray):
                    print_solution({f"{name}[{k}]": v})  # recursive call
                else:
                    print(f"{name}[{k}] = {v}")
        
        # Case 4: catch-all
        else:
            print(f"{name} = {values}")

# --- Fog case ---
def plot_loss_fog():
    distance_range = range(15, 225)  # km
    alt = 15  # km
    
    # Light, Moderate, Heavy fog (km visibility)
    V_fog_values = [1, 0.5, 0.2]
    labels = ["Light Fog (1 km)", "Moderate Fog (0.5 km)", "Heavy Fog (0.2 km)"]

    fig, ax = plt.subplots(figsize=(8, 6))

    for V_fog, label in zip(V_fog_values, labels):
        d_values, L_values = [], []

        for d in distance_range:
            alpha = math.asin(alt / d)

            # Geometric & misalignment loss
            L_geo = 20 * max(math.log10((R_TX + d * 1000 * THETA) / R_RX), 0)
            L_ma  = 0.01 * d

            # Cloud layer penetration
            H_C = 0.5
            R_C = H_C / math.sin(alpha)

            # Fog extinction coefficient
            if V_fog > 50:
                U = 1.6
            elif 6 < V_fog <= 50:
                U = 1.3
            else:
                U = 0.585 * V_fog ** (1 / 3)

            L_fog = (3.91 / V_fog) * ((LAMBDA / (550e-9)) ** (-U)) * R_C

            # Total loss
            L_t = L_geo + L_ma + L_fog

            d_values.append(d)
            L_values.append(L_t)

        ax.plot(d_values, L_values, label=label)

    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Total Loss $L_t$ (dB)")
    ax.set_title("Total Channel Loss vs Distance under Fog")
    ax.legend()
    ax.grid(True)
    plt.show()


# --- Rain case ---
def plot_loss_rain():
    distance_range = range(15, 225)  # km
    alt = 15  # km
    
    # Light, Moderate, Heavy rain (mm/hr equivalent visibility index)
    V_rain_snow_values = [10, 4, 2]
    labels = ["Light Rain (10)", "Moderate Rain (4)", "Heavy Rain (2)"]

    fig, ax = plt.subplots(figsize=(8, 6))

    for V_rain_snow, label in zip(V_rain_snow_values, labels):
        d_values, L_values = [], []

        for d in distance_range:
            alpha = math.asin(alt / d)

            # Geometric & misalignment loss
            L_geo = 20 * max(math.log10((R_TX + d * 1000 * THETA) / R_RX), 0)
            L_ma  = 0.01 * d

            # Water layer penetration
            H_W = 5
            R_W = H_W / math.sin(alpha)

            L_rain = (2.8 / V_rain_snow) * R_W

            # Total loss
            L_t = L_geo + L_ma + L_rain

            d_values.append(d)
            L_values.append(L_t)

        ax.plot(d_values, L_values, label=label)

    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Total Loss $L_t$ (dB)")
    ax.set_title("Total Channel Loss vs Distance under Rain")
    ax.legend()
    ax.grid(True)
    plt.show()

def plot_loss_distance(p_all, fog, rain, snow):
    # Weather cases: (rain, fog, snow)
    # if p_all:
    #     cases = {
    #         "No Weather": (0, 0, 0),
    #         "Rain Only":  (1, 0, 0),
    #         "Fog Only":   (0, 1, 0),
    #         # "Snow Only":  (0, 0, 1),
    #     }
    # elif fog:
    #     cases = {"Fog Only": (0, 1, 0)}
    # elif rain:
    #     cases = {"Rain Only": (1, 0, 0)}
    # elif snow:
    #     cases = {"Snow Only": (0, 0, 1)}
        
    # Prepare plot with 3 subplots
    fig     = plt.figure(figsize=(8, 6))
    ax_loss = fig.add_subplot(311)
    ax_skr  = fig.add_subplot(312)

    distance_range = range(15, 225)  # in km

    alt = 15

    L_values = []
    d_values = []
    K_values = []
    K_fvalues = []
    for d in distance_range:
        # Elevation angle
        alpha = math.asin(alt / d)

        THETA = 1.22 * 1550e-9 / 0.4

        # Losses
        L_geo = 20 * max(math.log10((0.4 + d * 1000 * THETA) / 0.3), 0)
        L_ma  = 0.01 * d

        H_W = 5
        H_C = 0.5
        R_W = H_W / math.sin(alpha)
        R_C = H_C / math.sin(alpha)

        # if fog:
        #     U = 1.6 if V_fog > 50 else (1.3 if 6 < V_fog <= 50 else 0.585 * V_fog ** (1 / 3))
        # elif rain or snow:
        #     U = 1.6 if V_rain_snow > 50 else (1.3 if 6 < V_rain_snow <= 50 else 0.585 * V_rain_snow ** (1 / 3))
        # else:
            # U = 1.6

        # L_fog  = (3.91 / V_fog)       * ((LAMBDA / 550 * 1e9) ** (-U)) * R_C
        # L_snw  = (58   / V_rain_snow) * R_W
        # L_rain = (2.8  / V_rain_snow) * R_W

        # # Total loss
        # L_t = L_geo + L_ma + L_fog * fog + L_snw * snow + L_rain * rain
        
        L_t = L_geo + L_ma

        ETA = 10 ** (-L_t / 10)

        ETA_fiber = 10 ** (-0.02*d)

        # Key rate
        K_link  = -80e6 * math.log2(1 - ETA)
        K_flink = -80e6 * math.log2(1 - ETA_fiber)

        d_values.append(d)
        L_values.append(L_t)
        K_values.append(K_link)
        K_fvalues.append(K_flink)

    # Plot Loss vs Distance
    ax_loss.plot(d_values, L_values, label="")

    # Loss plot formatting
    ax_loss.set_xlabel("Distance (km)")
    ax_loss.set_ylabel("Total Loss L_t (dB)")
    ax_loss.legend()
    ax_loss.grid(True)

    # Plot Loss vs Distance
    ax_skr.plot(d_values, K_values,   label="HAP")
    ax_skr.plot(d_values, K_fvalues, label="Fiber")

    # Loss plot formatting
    ax_skr.set_xlabel("Distance (km)")
    ax_skr.set_ylabel("SKR")
    ax_skr.legend()
    ax_skr.grid(True)

    plt.yscale("log")
    plt.tight_layout()
    plt.show()


def plot_key_rate_stratotegic(p_all, fog, rain, snow):
    # Keep only rows where Time_s is an integer
    df_filtered = df[df["Time_s"] % 1 == 0].reset_index(drop=True)

    # Extract required columns
    result = df_filtered[["Time_s", "Longitude_deg", "Latitude_deg", "Altitude_m"]]

    # Static GS coordinates
    la_rad_g = math.radians(49)
    lg_rad_g = math.radians(279)
    x_g = R * math.cos(la_rad_g) * math.cos(lg_rad_g)
    y_g = R * math.cos(la_rad_g) * math.sin(lg_rad_g)

    # Weather cases: (rain, fog, snow)
    if p_all:
        cases = {
            "No Weather": (0, 0, 0),
            "Rain Only":  (1, 0, 0),
            "Fog Only":   (0, 1, 0),
            "Snow Only":  (0, 0, 1),
        }
    elif fog:
        cases = {"Fog Only": (0, 1, 0)}
    elif rain:
        cases = {"Rain Only": (1, 0, 0)}
    elif snow:
        cases = {"Snow Only": (0, 0, 1)}
        
    # Prepare plot with 3 subplots
    fig = plt.figure(figsize=(14, 16))
    ax3d   = fig.add_subplot(311, projection="3d")
    ax2d   = fig.add_subplot(312)
    ax_loss = fig.add_subplot(313)

    time_range = range(0, 86400)  # simulate full day (adjust if too slow)

    # loop through each weather case
    for label, (rain, fog, snow) in cases.items():
        d_values, t_values, K_values, L_values = [], [], [], []

        for t in time_range:
            row = result.loc[result["Time_s"] == t]
            if row.empty:
                continue

            lon = row["Longitude_deg"].values[0]
            lat = row["Latitude_deg"].values[0]
            alt = row["Altitude_m"].values[0] / 1000  # km

            # Balloon coordinates
            la_rad_h = math.radians(lat)
            lg_rad_h = math.radians(lon)
            x_h = R * math.cos(la_rad_h) * math.cos(lg_rad_h)
            y_h = R * math.cos(la_rad_h) * math.sin(lg_rad_h)

            # Horizontal distance
            d_los_hor = math.sqrt((x_h - x_g) ** 2 + (y_h - y_g) ** 2)

            # Elevation angle
            alpha = math.atan(alt / d_los_hor) if d_los_hor > 0 else math.pi / 2

            # LOS distance
            d_los = alt / math.sin(alpha)

            THETA = 1.22 * 1550e-9 / 0.4

            # Losses
            L_geo = 20 * max(math.log10((0.4 + d_los * 1000 * THETA) / 0.3), 0)
            L_ma  = 0.01 * d_los

            H_W = 5
            H_C = 0.5
            R_W = H_W / math.sin(alpha)
            R_C = H_C / math.sin(alpha)

            if fog:
                U = 1.6 if V_fog > 50 else (1.3 if 6 < V_fog <= 50 else 0.585 * V_fog ** (1 / 3))
            elif rain or snow:
                U = 1.6 if V_rain_snow > 50 else (1.3 if 6 < V_rain_snow <= 50 else 0.585 * V_rain_snow ** (1 / 3))
            else:
                U = 1.6

            L_fog  = (3.91 / V_fog)       * ((LAMBDA / 550 * 1e9) ** (-U)) * R_C
            L_snw  = (58   / V_rain_snow) * R_W
            L_rain = (2.8  / V_rain_snow) * R_W

            # Total loss
            L_t = L_geo + L_ma + L_fog * fog + L_snw * snow + L_rain * rain

            ETA = 10 ** (-L_t / 10)

            # Key rate
            K_link = -B * math.log2(1 - ETA)

            d_values.append(d_los)
            t_values.append(t)
            K_values.append(K_link)
            L_values.append(L_t)

        # Plot in 3D
        ax3d.plot(d_values, t_values, K_values, label=label)

        # Plot Key Rate vs Time
        ax2d.plot(t_values, K_values, label=label)

        # Plot Loss vs Time
        ax_loss.plot(t_values, L_values, label=label)

    # 3D plot formatting
    ax3d.set_xlabel("LoS Distance (km)")
    ax3d.set_ylabel("Time (s)")
    ax3d.set_zlabel("Max Key Rate (bps)")
    ax3d.set_title("Max Key Rate vs LoS Distance over Time")
    ax3d.legend()

    # 2D Key Rate plot formatting
    ax2d.set_xlabel("Time (s)")
    ax2d.set_ylabel("Key Rate (bps)")
    ax2d.set_title("Key Rate vs Time")
    ax2d.legend()
    ax2d.grid(True)

    # Loss plot formatting
    ax_loss.set_xlabel("Time (s)")
    ax_loss.set_ylabel("Total Loss L_t (dB)")
    ax_loss.set_title("Total Channel Loss vs Time")
    ax_loss.legend()
    ax_loss.grid(True)

    plt.tight_layout()
    plt.show()

def plot_transmittance_stratotegic_real(n=5, d_min_t=0, d_max_t=800): #86400
    """
    Plot theoretical and simulated transmittance for a balloon trajectory over time.
    df : pandas.DataFrame with ["Time_s", "Longitude_deg", "Latitude_deg", "Altitude_m"]
    ts : transmittance simulation module/object with .theoretical_eff and .simulated_eff
    """
    # Filter integer timestamps
    df_filtered = df[df["Time_s"] % 100 == 0].reset_index(drop=True)

    # Time range (0 to 86400 sec default)
    times = df_filtered["Time_s"].values
    mask = (times >= d_min_t) & (times <= d_max_t)
    times = times[mask]

    # Distances for each timestamp (LoS distance in km)
    distances = []
    for _, row in df_filtered[mask].iterrows():
        lat = math.radians(row["Latitude_deg"])
        lon = math.radians(row["Longitude_deg"])
        alt = row["Altitude_m"] / 1000  # km
        # Simplify: use alt as proxy distance (or replace with real LoS calc if needed)
        distances.append(alt)  
    
    # Compute efficiencies
    eta_theory = [ts.channel_theory("downlink", 0, 15, d, n) for d in distances]
    eta_sim    = [ts.channel_simulation("downlink", 0, 15, d, n) for d in distances]

    # Plot
    plt.figure(figsize=(12,6))
    plt.plot(times, eta_theory, label="Theoretical Transmittance", color="blue")
    plt.plot(times, eta_sim, label="Simulated Transmittance", color="red", linestyle="--")

    plt.xlabel("Time (s)")
    plt.ylabel("Transmittance (η)")
    plt.title("Transmittance vs Time (Full-day Balloon Trajectory)")
    plt.grid(True)
    plt.legend()
    plt.show()

# def plot_skr(dir, n, d_list, h_list):
#     """
#     Plot theoretical and simulated SKR for a balloon trajectory over time.
#     """
#     times = range(len(d_list))

#     # Compute SKRs
#     skr_theory = []
#     skr_sim    = []
#     skr_plob_theory = []
#     skr_plob_sim    = []
#     spinner = ['|', '/', '-', '\\']
#     for idx_d, d in enumerate(d_list):
#         eta_t = ts.channel_theory(direction=dir, gs_alt=0, balloon_alt=h_list[idx_d], distance=d, n_correction=n)
#         eta_s = ts.channel_simulation(direction=dir, gs_alt=0, balloon_alt=h_list[idx_d], distance=d, n_correction=n)

#         # eta_t = ts.theoretical_eff(distance=d, h_balloons=h_list[idx_d], n=n)
#         # eta_s = ts.simulated_eff(distance=d, h_balloons=h_list[idx_d], n=n)
        
#         skr_theory.append(ts.compute_skr(eta_t))
#         skr_sim.append(ts.compute_skr(eta_s))
#         skr_plob_theory.append(-ts.ratesources * ts.sourceeff * math.log2(1 - eta_t))
#         skr_plob_sim.append(-ts.ratesources * ts.sourceeff * math.log2(1 - eta_s))

#         # print(f"skr_plob_theory: {skr_plob_theory[idx_d]}, skr_plob_sim: {skr_plob_sim[idx_d]}")
#         sys.stdout.write("\rProcessing... " + spinner[idx_d % len(spinner)])
#         sys.stdout.flush()

#     # Plot
#     plt.figure(figsize=(12,6))
#     plt.plot(times, skr_theory, label="Theoretical SKR", color="green")
#     plt.plot(times, skr_sim, label="Simulated SKR", color="orange", linestyle="--")
#     plt.plot(times, skr_plob_theory, label="Theoretical SKR Upper Bound", color="blue", linestyle="-.")
#     plt.plot(times, skr_plob_sim, label="Simulated SKR Upper Bound", color="red", linestyle=":")

#     plt.xlabel("Time (s)")
#     plt.ylabel("SKR (bps)")
#     plt.title("Secret Key Rate vs Time (Full-day Balloon Trajectory)")
#     plt.grid(True)
#     plt.legend()
#     plt.show()2D

def plot_connectivity_graph(gnodes, hnodes, links):
    """
    Plot connectivity graph of Ground Stations (GS) and HAPs with trajectories
    using pure Matplotlib (no NetworkX drawing) to ensure proper longitude/latitude ticks.
    
    Parameters:
    -----------
    gnodes : list
        List of GS objects with attributes 'la' (latitude), 'lg' (longitude), and optional 'tag'.
    hnodes : list
        List of HAP objects with attributes 'la', 'lg', and optional 'tag'.
        If HAP has a trajectory, use full lists for 'la' and 'lg'.
    links : list
        List of link objects with attributes 'n1' and 'n2' (nodes from gnodes/hnodes).
    """
    plt.figure(figsize=(6, 4))
    
    # --- Plot GS nodes ---
    for gs_node in gnodes:
        plt.scatter(gs_node.lg, gs_node.la, color='skyblue', s=80, zorder=5, marker='^')
        # Optional: label the GS
        if hasattr(gs_node, 'tag'):
            plt.text(gs_node.lg + 0.04, gs_node.la + 0.04, gs_node.tag, fontsize=9)
    
    # --- Plot HAP nodes (initial position) ---
    for hap_node in hnodes:
        plt.scatter(hap_node.lg[0], hap_node.la[0], color='orange', s=5, zorder=5)
        if hasattr(hap_node, 'tag'):
            plt.text(hap_node.lg[0] - 0.4, hap_node.la[0] - 0.2, hap_node.tag, fontsize=9)
    
    # --- Plot edges without duplicates ---
    plotted_edges = set()
    for l in links:
        # Use frozenset to make the edge unordered (A-B same as B-A)
        edge_key = frozenset([l.n1, l.n2])
        if edge_key in plotted_edges:
            continue  # already plotted
        plotted_edges.add(edge_key)
    
        # Determine coordinates for nodes
        x = [l.n1.lg[0] if isinstance(l.n1.lg, list) else l.n1.lg,
             l.n2.lg[0] if isinstance(l.n2.lg, list) else l.n2.lg]
        y = [l.n1.la[0] if isinstance(l.n1.la, list) else l.n1.la,
             l.n2.la[0] if isinstance(l.n2.la, list) else l.n2.la]
        
        # Decide line style
        plt.plot(x, y, color='grey', linestyle='--', alpha=0.6, linewidth=0.5)
    
    # --- Plot HAP trajectories ---
    for hap_node in hnodes:
        plt.plot(hap_node.lg, hap_node.la, color='orange', linewidth=2, alpha=0.8)
    
    
    # --- Axis labels and limits ---
    all_lons = [gs.lg for gs in gnodes] + [hap.lg[0] for hap in hnodes]
    all_lats = [gs.la for gs in gnodes] + [hap.la[0] for hap in hnodes]
    plt.xlabel("Longitude", fontsize=13)
    plt.ylabel("Latitude", fontsize=13)
    plt.xlim(min(all_lons) - 0.2, max(all_lons) + 1.5)
    plt.ylim(min(all_lats) - 0.2, max(all_lats) + 0.3)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    
    # --- Legend ---
    custom_handles = [
        Line2D([], [], marker='^', color='skyblue', linestyle='None', markersize=6, label='GS'),
        Line2D([], [], marker='o', color='orange', linestyle='None', markersize=6, label='HAP')
    ]
    plt.legend(handles=custom_handles, loc='upper right', frameon=True, fontsize=13)
    
    plt.grid(True, alpha=0.3)

    # ==============================
    # Zoomed-in inset for one HAP
    # ==============================
    hap_zoom = hnodes[3]   # choose which HAP to zoom

    # Create inset axis
    ax = plt.gca()
    axins = inset_axes(
        ax,
        width="35%",   # relative size
        height="35%",
        loc="lower right",
        borderpad=1.2
    )

    # Plot trajectory inside inset
    axins.plot(hap_zoom.lg, hap_zoom.la,
               color='orange', linewidth=2)

    # Optional: mark start point
    axins.scatter(hap_zoom.lg[0], hap_zoom.la[0],
                  color='orange', s=2, zorder=5)

    # Set zoom window (tight bounds)
    margin = 0.05
    axins.set_xlim(min(hap_zoom.lg) - margin, max(hap_zoom.lg) + margin)
    axins.set_ylim(min(hap_zoom.la) - margin, max(hap_zoom.la) + margin)

    # Clean inset appearance
    axins.set_xticks([])
    axins.set_yticks([])
    axins.grid(True, alpha=0.3)

    # Draw rectangle on main plot to show zoomed region
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
    plt.savefig("hap_qkd_network.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()


def animate_hap_trajectories(times, lons_list, lats_list, hap_names):
    """
    Animate HAP trajectories over time using Plotly.

    Parameters
    ----------
    times : list or array
        List of time steps.
    lons_list : list of lists
        Each sublist is the longitude trajectory of a HAP over time.
    lats_list : list of lists
        Each sublist is the latitude trajectory of a HAP over time.
    hap_names : list of str
        Names/IDs of the HAPs (must match number of lons/lats sublists).
    """
    # Build DataFrame
    data = {"time": [], "lat": [], "lon": [], "hap": []}
    for hap_idx, hap_name in enumerate(hap_names):
        for t_idx, t in enumerate(times):
            data["time"].append(t)
            data["lat"].append(lats_list[hap_idx][t_idx])
            data["lon"].append(lons_list[hap_idx][t_idx])
            data["hap"].append(hap_name)
    
    df = pd.DataFrame(data)
    
    # Plot animated scatter
    fig = px.scatter(df, x="lat", y="lon", animation_frame="time", color="hap",
                     range_x=[min(map(min, lats_list)) - 0.1, max(map(max, lats_list)) + 0.1],
                     range_y=[min(map(min, lons_list)) - 0.1, max(map(max, lons_list)) + 0.1])
    fig.update_layout(xaxis_title="Latitude", yaxis_title="Longitude")
    fig.show()

























def compute_point(d, h, dir, n):
    """
    Worker function for one (distance, height).
    """
    eta_t = ts.channel_theory(direction=dir, gs_alt=0, hap_alt=h, distance=d, n_correction=n, params_in=ts.params)
    #eta_s = ts.channel_simulation(direction=dir, gs_alt=0, hap_alt=h, distance=d, n_correction=n)
    # print(f"dir: {dir}, h: {h}, d: {d}, n: {n}")
    # print(f"eta_t: {eta_t}, eta_s: {eta_s}")

    skr_t   = ts.compute_skr(eta_t)
    #skr_s   = ts.compute_skr(eta_s)
    skr_pt  = -ts.ratesources * ts.sourceeff * math.log2(1 - eta_t)
    #skr_ps  = -ts.ratesources * ts.sourceeff * math.log2(1 - eta_s)
    #return eta_t, eta_s, skr_t, skr_s, skr_pt, skr_ps
    return eta_t, skr_t, skr_pt

def compute_point_wsim(d, h, dir, n, params_in):
    """
    Worker function for one (distance, height).
    """
    eta_t = ts.channel_theory(direction=dir, gs_alt=0, hap_alt=h, distance=d, n_correction=n, params_in=params_in)

    skr_t   = ts.compute_skr(eta_t)
    skr_pt  = -ts.ratesources * ts.sourceeff * math.log2(1 - eta_t)
    return eta_t, skr_t, skr_pt

def plot_skr(dir, n, d_list, h_list, max_workers=24):
    """
    Parallelized SKR plotter across CPU cores with progress bar.
    """
    times = range(len(d_list))

    # Pack args
    tasks = [(d, h_list[idx], dir, n) for idx, d in enumerate(d_list)]
    skr_theory, skr_sim, skr_plob_theory, skr_plob_sim = [], [], [], []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_point, t) for t in tasks]
        for f in tqdm(as_completed(futures), total=len(futures), desc="Computing SKR"):
            skr_t, skr_s, skr_pt, skr_ps = f.result()
            skr_theory.append(skr_t * 1e-3)
            skr_sim.append(skr_s * 1e-3)
            skr_plob_theory.append(skr_pt * 1e-3)
            skr_plob_sim.append(skr_ps * 1e-3)

    print("\nAll points computed.")

    # Plot
    plt.figure(figsize=(12,6))
    plt.plot(times, skr_theory, label="SKR DW Bound (Theoretical)", color="green")
    plt.plot(times, skr_sim, label="SKR DW Bound (Simulation)", color="orange", linestyle="--")
    plt.plot(times, skr_plob_theory, label="SKR PLOB Bound (Theoretical)", color="blue", linestyle="-.")
    plt.plot(times, skr_plob_sim, label="SKR PLOB Bound (Simulation)", color="red", linestyle=":")

    plt.xlabel("Time (s)")
    plt.ylabel("SKR (kbps)")
    plt.grid(True)
    plt.legend()
    plt.savefig("skr.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()

# def plot_skr(dir, n, d_list, h_list, start_time=0, end_time=None, 
#               max_workers=24, result_file="skr_results.pkl"):
#     """
#     Parallelized SKR plotter across CPU cores with progress bar.
#     Supports incremental runs with (start_time, end_time) and persistent results.
#     """
#     if end_time is None:
#         end_time = len(d_list)

#     # Select only the slice for this run
#     d_list = d_list[start_time:end_time]
#     h_list = h_list[start_time:end_time]
#     times = list(range(start_time, end_time))

#     # Pack args
#     tasks = [(d, h_list[idx], dir, n) for idx, d in enumerate(d_list)]
#     new_theory, new_sim, new_plob_theory, new_plob_sim = [], [], [], []

#     print(f"Running interval {start_time} → {end_time} ({len(tasks)} points)")

#     # Parallel compute
#     with ProcessPoolExecutor(max_workers=max_workers) as executor:
#         futures = [executor.submit(compute_point, t) for t in tasks]
#         for f in tqdm(as_completed(futures), total=len(futures), desc="Computing SKR"):
#             skr_t, skr_s, skr_pt, skr_ps = f.result()
#             new_theory.append(skr_t * 1e-3)
#             new_sim.append(skr_s * 1e-3)
#             new_plob_theory.append(skr_pt * 1e-3)
#             new_plob_sim.append(skr_ps * 1e-3)

#     # print("len(times):", len(times))
#     # print("len(new_theory):", len(new_theory))
#     # print("len(new_sim):", len(new_sim))
#     # print("len(new_plob_theory):", len(new_plob_theory))

#     print("\nInterval computation done. Saving results...")

#     # === Load previous results if exist ===
#     if os.path.exists(result_file):
#         with open(result_file, "rb") as f:
#             data = pickle.load(f)
#         skr_theory = data["skr_theory"]
#         skr_sim = data["skr_sim"]
#         skr_plob_theory = data["skr_plob_theory"]
#         skr_plob_sim = data["skr_plob_sim"]
#         old_times = data["times"]
#     else:
#         skr_theory, skr_sim, skr_plob_theory, skr_plob_sim, old_times = [], [], [], [], []

#     # === Align results ===
#     # Ensure we’re extending only if the new times don’t overlap
#     existing_points = set(old_times)
#     for i, t in enumerate(times):
#         if t not in existing_points:
#             skr_theory.append(new_theory[i])
#             skr_sim.append(new_sim[i])
#             skr_plob_theory.append(new_plob_theory[i])
#             skr_plob_sim.append(new_plob_sim[i])
#             old_times.append(t)

#     all_times = old_times

#     # === Save updated results ===
#     with open(result_file, "wb") as f:
#         pickle.dump({
#             "skr_theory": skr_theory,
#             "skr_sim": skr_sim,
#             "skr_plob_theory": skr_plob_theory,
#             "skr_plob_sim": skr_plob_sim,
#             "times": all_times
#         }, f)

#     print(f"Results saved to {result_file} ({len(all_times)} total points).")

#     # === Plot cumulative results ===
#     plt.figure(figsize=(12,6))
#     plt.plot(all_times, skr_theory, label="SKR DW Bound (Theoretical)", color="green")
#     plt.plot(all_times, skr_sim, label="SKR DW Bound (Simulation)", color="orange", linestyle="--")
#     plt.plot(all_times, skr_plob_theory, label="SKR PLOB Bound (Theoretical)", color="blue", linestyle="-.")
#     plt.plot(all_times, skr_plob_sim, label="SKR PLOB Bound (Simulation)", color="red", linestyle=":")

#     plt.xlabel("Time (s)")
#     plt.ylabel("SKR (kbps)")
#     plt.grid(True)
#     plt.legend()
#     plt.savefig("skr.svg", format="svg", dpi=300, bbox_inches="tight")
#     plt.show()

def plot_skr(result_file="skr_results.pkl", outlier_factor=3.0):
    """
    Load and plot existing SKR data from result_file.
    Automatically filters outliers where a point jumps > outlier_factor × previous point.
    """

    if not os.path.exists(result_file):
        print(f"❌ No result file found at: {result_file}")
        return

    # === Load results ===
    with open(result_file, "rb") as f:
        data = pickle.load(f)

    skr_theory = data.get("skr_theory", [])
    skr_sim = data.get("skr_sim", [])
    skr_plob_theory = data.get("skr_plob_theory", [])
    skr_plob_sim = data.get("skr_plob_sim", [])
    times = data.get("times", [])

    if not times:
        print("⚠️ No data found in the file.")
        return

    print(f"Loaded {len(times)} data points from {result_file}")

    # === Define outlier filtering helper ===
    def remove_outliers(values):
        if not values:
            return values
        cleaned = [values[0]]
        for i in range(1, len(values)):
            if values[i] > outlier_factor * cleaned[-1]:
                # Outlier detected – replace with previous value
                cleaned.append(cleaned[-1])
            else:
                cleaned.append(values[i])
        return cleaned

    # === Apply outlier removal ===
    skr_theory = remove_outliers(skr_theory)
    skr_sim = remove_outliers(skr_sim)
    skr_plob_theory = remove_outliers(skr_plob_theory)
    skr_plob_sim = remove_outliers(skr_plob_sim)

    # === Convert time scale to hours ===
    times_hours = [t / 3600.0 for t in times]

    # === Plot cleaned results ===
    plt.figure(figsize=(8,3))
    plt.plot(times_hours, skr_sim, label="DW (Simulation)", color="orange", linestyle="--", marker="o", markersize=8, markevery=8640)
    plt.plot(times_hours, skr_plob_sim, label="PLOB (Simulation)", color="red", linestyle=":", marker="x", markersize=8, markevery=8640)
    plt.plot(times_hours, skr_theory, label="DW (Theory)", color="green", marker="s", markersize=8, markevery=8640)
    plt.plot(times_hours, skr_plob_theory, label="PLOB (Theory)", color="blue", linestyle="-.", marker="^", markersize=8, markevery=8640)

    plt.xlabel("Time (hours)")
    plt.ylabel("Maximum SKR (kbps)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("skr_cleaned.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()

    print("✅ Plot complete (outliers replaced if detected).")