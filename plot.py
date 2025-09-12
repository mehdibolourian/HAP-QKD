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