from libraries import *

plt.rcParams['font.family']  = 'DeJavu Serif'
plt.rcParams['font.serif']   = ['Times New Roman']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype']  = 42

df = pd.read_csv("dataset/balloon_sim_data.csv")

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
    if p_all:
        cases = {
            "No Weather": (0, 0, 0),
            "Rain Only":  (1, 0, 0),
            "Fog Only":   (0, 1, 0),
            # "Snow Only":  (0, 0, 1),
        }
    elif fog:
        cases = {"Fog Only": (0, 1, 0)}
    elif rain:
        cases = {"Rain Only": (1, 0, 0)}
    elif snow:
        cases = {"Snow Only": (0, 0, 1)}
        
    # Prepare plot with 3 subplots
    fig     = plt.figure(figsize=(8, 6))
    ax_loss = fig.add_subplot(111)

    distance_range = range(15, 225)  # in km

    alt = 15

    # loop through each weather case
    for label, (rain, fog, snow) in cases.items():
        d_values, L_values = [], []

        for d in distance_range:
            # Elevation angle
            alpha = math.asin(alt / d)

            # Losses
            L_geo = 20 * max(math.log10((R_TX + d * 1000 * THETA) / R_RX), 0)
            L_ma  = 0.01 * d

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

            d_values.append(d)
            L_values.append(L_t)

        # Plot Loss vs Distance
        ax_loss.plot(d_values, L_values, label=label)

    # Loss plot formatting
    ax_loss.set_xlabel("Distance (km)")
    ax_loss.set_ylabel("Total Loss L_t (dB)")
    ax_loss.set_title("Total Channel Loss vs Distance")
    ax_loss.legend()
    ax_loss.grid(True)

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

            # Losses
            L_geo = 20 * max(math.log10((R_TX + d_los * 1000 * THETA) / R_RX), 0)
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

def plot_transmittance_stratotegic_real(n=5, d_min_t=0, d_max_t=8): #86400
    """
    Plot theoretical and simulated transmittance for a balloon trajectory over time.
    df : pandas.DataFrame with ["Time_s", "Longitude_deg", "Latitude_deg", "Altitude_m"]
    ts : transmittance simulation module/object with .theoretical_eff and .simulated_eff
    """
    # Filter integer timestamps
    df_filtered = df[df["Time_s"] % 1 == 0].reset_index(drop=True)

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
    eta_theory = [ts.theoretical_eff(distance=d, h_balloons=15, n=n) for d in distances]
    eta_sim    = [ts.simulated_eff(distance=d, h_balloons=15, n=n) for d in distances]

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

def plot_skr_stratotegic_real(n=5, d_min_t=0, d_max_t=8):
    """
    Plot theoretical and simulated SKR for a balloon trajectory over time.
    df : pandas.DataFrame with ["Time_s", "Longitude_deg", "Latitude_deg", "Altitude_m"]
    ts : transmittance simulation module/object with .theoretical_eff, .simulated_eff, .compute_skr
    """
    df_filtered = df[df["Time_s"] % 1 == 0].reset_index(drop=True)

    times = df_filtered["Time_s"].values
    mask = (times >= d_min_t) & (times <= d_max_t)
    times = times[mask]

    distances = []
    for _, row in df_filtered[mask].iterrows():
        alt = row["Altitude_m"] / 1000  # km
        distances.append(alt)

    # Compute SKRs
    skr_theory = []
    skr_sim    = []
    for d in distances:
        eta_t = ts.theoretical_eff(distance=d, h_balloons=15, n=n)
        eta_s = ts.simulated_eff(distance=d, h_balloons=15, n=n)
        skr_theory.append(ts.compute_skr(eta_t))
        skr_sim.append(ts.compute_skr(eta_s))

    # Plot
    plt.figure(figsize=(12,6))
    plt.plot(times, skr_theory, label="Theoretical SKR", color="green")
    plt.plot(times, skr_sim, label="Simulated SKR", color="orange", linestyle="--")

    plt.xlabel("Time (s)")
    plt.ylabel("SKR (bps)")
    plt.title("Secret Key Rate vs Time (Full-day Balloon Trajectory)")
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
    plt.ylim(min(all_lats) - 0.2, max(all_lats) + 0.2)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    
    # --- Legend ---
    custom_handles = [
        Line2D([], [], marker='^', color='skyblue', linestyle='None', markersize=6, label='GS'),
        Line2D([], [], marker='o', color='orange', linestyle='None', markersize=6, label='HAP')
    ]
    plt.legend(handles=custom_handles, loc='best', frameon=True, fontsize=13)
    
    plt.grid(True, alpha=0.3)

    # ==============================
    # Zoomed-in inset for one HAP
    # ==============================
    hap_zoom = hnodes[0]   # choose which HAP to zoom

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

























def compute_point(args):
    """
    Worker function for one (distance, height).
    """
    d, h, dir, n = args
    eta_t = ts.channel_theory(direction=dir, gs_alt=0, balloon_alt=h, distance=d, n_correction=n)
    eta_s = ts.channel_simulation(direction=dir, gs_alt=0, balloon_alt=h, distance=d, n_correction=n)
    # print(f"dir: {dir}, h: {h}, d: {d}, n: {n}")
    # print(f"eta_t: {eta_t}, eta_s: {eta_s}")

    skr_t   = ts.compute_skr(eta_t)
    skr_s   = ts.compute_skr(eta_s)
    skr_pt  = -ts.ratesources * ts.sourceeff * math.log2(1 - eta_t)
    skr_ps  = -ts.ratesources * ts.sourceeff * math.log2(1 - eta_s)
    return skr_t, skr_s, skr_pt, skr_ps

# def plot_skr(dir, n, d_list, h_list, max_workers=24):
#     """
#     Parallelized SKR plotter across CPU cores with progress bar.
#     """
#     times = range(len(d_list))

#     # Pack args
#     tasks = [(d, h_list[idx], dir, n) for idx, d in enumerate(d_list)]
#     skr_theory, skr_sim, skr_plob_theory, skr_plob_sim = [], [], [], []

#     with ProcessPoolExecutor(max_workers=max_workers) as executor:
#         futures = [executor.submit(compute_point, t) for t in tasks]
#         for f in tqdm(as_completed(futures), total=len(futures), desc="Computing SKR"):
#             skr_t, skr_s, skr_pt, skr_ps = f.result()
#             skr_theory.append(skr_t * 1e-3)
#             skr_sim.append(skr_s * 1e-3)
#             skr_plob_theory.append(skr_pt * 1e-3)
#             skr_plob_sim.append(skr_ps * 1e-3)

#     print("\nAll points computed.")

#     # Plot
#     plt.figure(figsize=(12,6))
#     plt.plot(times, skr_theory, label="SKR DW Bound (Theoretical)", color="green")
#     plt.plot(times, skr_sim, label="SKR DW Bound (Simulation)", color="orange", linestyle="--")
#     plt.plot(times, skr_plob_theory, label="SKR PLOB Bound (Theoretical)", color="blue", linestyle="-.")
#     plt.plot(times, skr_plob_sim, label="SKR PLOB Bound (Simulation)", color="red", linestyle=":")

#     plt.xlabel("Time (s)")
#     plt.ylabel("SKR (kbps)")
#     plt.grid(True)
#     plt.legend()
#     plt.savefig("skr.svg", format="svg", dpi=300, bbox_inches="tight")
#     plt.show()

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
    plt.plot(times_hours, skr_sim, label="SKR DW Bound (Simulation)", color="orange", linestyle="--", marker="o", markersize=8, markevery=8640)
    plt.plot(times_hours, skr_plob_sim, label="SKR PLOB Bound (Simulation)", color="red", linestyle=":", marker="x", markersize=8, markevery=8640)
    plt.plot(times_hours, skr_theory, label="SKR DW Bound (Theoretical)", color="green", marker="s", markersize=8, markevery=8640)
    plt.plot(times_hours, skr_plob_theory, label="SKR PLOB Bound (Theoretical)", color="blue", linestyle="-.", marker="^", markersize=8, markevery=8640)

    plt.xlabel("Time (hours)")
    plt.ylabel("Maximum SKR (kbps)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("skr_cleaned.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()

    print("✅ Plot complete (outliers replaced if detected).")