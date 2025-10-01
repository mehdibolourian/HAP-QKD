from libraries import *

R           = 6371  # Earth's radius in km
R_TX        = 0.1
R_RX        = 0.4 # 0.3
LAMBDA      = 1550e-9
THETA       = 1.22 * LAMBDA / R_TX
H_W         = 5   # Altitude of snow/rain
H_C         = 0.5 # Altitude of fog
V_fog       = 0.5 # 0.5km visibility (Moderate fog)
V_rain_snow = 4   # 4km visibility   (Moderate rain/snow)

level     = "50"  # hPa level (~20 km altitude)
file_name = f"era5_{level}hpa_hourly.nc"

# Earth radius in km
R = 6371.0  

def latlon_to_cartesian(lat, lon):
    """Convert lat/lon (deg) to Cartesian coordinates (x,y,z)."""
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    x = R * np.cos(lat_rad) * np.cos(lon_rad)
    y = R * np.cos(lat_rad) * np.sin(lon_rad)
    z = R * np.sin(lat_rad)
    return np.array([x, y, z])

def cartesian_to_latlon(x, y, z):
    """Convert Cartesian coordinates back to lat/lon (deg)."""
    lat = np.degrees(np.arcsin(z / R))
    lon = np.degrees(np.arctan2(y, x))
    return lat, lon+360

def rotation_matrix_from_vectors(a, b):
    """Find rotation matrix that rotates vector a → b on the sphere."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)
    if np.isclose(c, 1.0):  # same vector
        return np.eye(3)
    if np.isclose(c, -1.0):  # opposite vector
        # 180° rotation around any perpendicular vector
        perp = np.array([1,0,0]) if not np.allclose(a,[1,0,0]) else np.array([0,1,0])
        v = np.cross(a, perp)
        v /= np.linalg.norm(v)
        H = np.array([[0, -v[2], v[1]],
                      [v[2], 0, -v[0]],
                      [-v[1], v[0], 0]])
        return -np.eye(3) + 2 * np.outer(v,v)
    s = np.linalg.norm(v)
    k = np.array([[0, -v[2], v[1]],
                  [v[2], 0, -v[0]],
                  [-v[1], v[0], 0]])
    Rm = np.eye(3) + k + k @ k * ((1-c) / (s**2))
    return Rm

def shift_trajectory(lats, lons, new_lat0, new_lon0):
    """
    Shift a trajectory so that its first point moves to (new_lat0, new_lon0),
    preserving the relative path on the Earth's sphere.
    """
    # Original start
    lat0, lon0 = lats[0], lons[0]

    if lat0 != new_lat0 or lon0 != new_lon0:
        # Convert all points to Cartesian
        traj_xyz = np.array([latlon_to_cartesian(lat, lon) for lat, lon in zip(lats, lons)])
        start_vec = traj_xyz[0]
        new_start_vec = latlon_to_cartesian(new_lat0, new_lon0)
        
        # Compute rotation
        Rm = rotation_matrix_from_vectors(start_vec, new_start_vec)
        
        # Rotate trajectory
        rotated_xyz = traj_xyz @ Rm.T
        
        # Convert back to lat/lon
        new_lats, new_lons = [], []
        for x, y, z in rotated_xyz:
            lat, lon = cartesian_to_latlon(x, y, z)
            new_lats.append(lat)
            new_lons.append(lon)
        
        return new_lats, new_lons
    else:
        return lats, lons


def latlon_to_tangent(lon, lat, lon0, lat0):
    """Project lat/lon (deg) into local tangent plane (x,y) in km 
       relative to reference point (lat0, lon0)."""
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    lat0_rad = np.radians(lat0)
    lon0_rad = np.radians(lon0)

    dlat = lat_rad - lat0_rad
    dlon = lon_rad - lon0_rad

    x = R * dlon * np.cos(lat0_rad)  # East
    y = R * dlat                     # North
    return x, y

def tangent_to_latlon(x, y, lon0, lat0):
    """Convert local tangent plane (x,y) back to lat/lon (deg)."""
    lat0_rad = np.radians(lat0)
    lon0_rad = np.radians(lon0)

    dlat = y / R
    dlon = x / (R * np.cos(lat0_rad))

    lat = lat0_rad + dlat
    lon = lon0_rad + dlon
    return np.degrees(lon), np.degrees(lat)






    
# Create transformer objects
# Always use `always_xy=True` so that (lon, lat) order is consistent
to_xy_transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
to_latlon_transformer = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)

def lonlat_to_xy(lon, lat):
    """Convert longitude/latitude to absolute x/y in meters (Web Mercator)."""
    x, y = to_xy_transformer.transform(lon, lat)
    return [x, y]

def xy_to_lonlat(x, y):
    """Convert absolute x/y in meters (Web Mercator) back to longitude/latitude."""
    lon, lat = to_latlon_transformer.transform(x, y)
    return [lon+360, lat]

def update_coordinates(method, hnodes, syst):
    if method == "wind":
        download_data()
    
        # Load the NetCDF file
        ds = xr.open_dataset(file_name)
    
        # Extract variables
        u = ds['u']     # u-component of wind (zonal)
        v = ds['v']     # v-component of wind (meridional)
        time = ds['valid_time']
        lat = ds['latitude']
        lon = ds['longitude']

        # Compute wind speed magnitude (m/s)
        wind_speed = np.sqrt(u**2 + v**2)
        
        u = ds.u.sel(latitude=lat, longitude=lon, method="nearest").values.flatten()
        v = ds.v.sel(latitude=lat, longitude=lon, method="nearest").values.flatten()
        time = ds.valid_time.values
    
        speed = compute_speed(u, v)
        direction = compute_direction(u, v)
        
        hnodes_xy = np.zeros((len(syst.T),2))
    elif method == "stratotegic":
        # Load CSV
        df = pd.read_csv("dataset/balloon_sim_data.csv")
        # Keep only rows where Time_s is an integer
        df_filtered = df[df["Time_s"] % 1 == 0].reset_index(drop=True)
        # Extract required columns
        result = df_filtered[["Time_s", "Longitude_deg", "Latitude_deg", "Altitude_m"]]

        df_int = result[result["Time_s"].apply(float.is_integer)]

        lons = df_int["Longitude_deg"].tolist()
        lats = df_int["Latitude_deg"].tolist()

        new_lats = {}
        new_lons = {}

        for idx_hnode, hnode in enumerate(hnodes):
            new_lats[idx_hnode], new_lons[idx_hnode] = shift_trajectory(lats, lons, hnode.la[0], hnode.lg[0])
            #print(new_lons[idx_hnode][0:int(len(syst.T))])
            #print(new_lats[idx_hnode][0:int(len(syst.T))])
    else:
        raise ValueError("Method must be either 'wind' or 'stratotegic'")

    if method == "wind":
        for idx_hnode, hnode in enumerate(hnodes):
            hnode.H[0] = 15
    elif method == "stratotegic":
        for idx_hnode, hnode in enumerate(hnodes):
            hnode.H[0] = 25
    
    for t in syst.T[1:]:
        if method == "wind":
            for idx_hnode, hnode in enumerate(hnodes):
                # Get current position in meters
                x, y = lonlat_to_xy(hnode.lg[t-1], hnode.la[t-1])
    
                # Get wind speed at that position and time
                u_val = ds.u.sel(
                    valid_time=time[t-1],
                    latitude=hnode.la[t-1],
                    longitude=hnode.lg[t-1],
                    method="nearest"
                ).values.item()
    
                v_val = ds.v.sel(
                    valid_time=time[t-1],
                    latitude=hnode.la[t-1],
                    longitude=hnode.lg[t-1],
                    method="nearest"
                ).values.item()
    
                # Update position in meters
                x_new = x + u_val * syst.THETA
                y_new = y + v_val * syst.THETA
    
                # Store new xy
                hnodes_xy[t] = [x_new, y_new]
    
                # Convert back to lon/lat
                lon_new, lat_new = xy_to_lonlat(x_new, y_new)
                hnode.lg[t] = lon_new
                hnode.la[t] = lat_new
        elif method == "stratotegic":
            # Find the row corresponding to time t
            row = result.loc[result["Time_s"] == t]
            if not row.empty:
                for idx_hnode, hnode in enumerate(hnodes):
                    lon = new_lons[idx_hnode][t] #row["Longitude_deg"].values[0]
                    lat = new_lats[idx_hnode][t] #row["Latitude_deg"].values[0]
                    alt = row["Altitude_m"].values[0]
                    
                    hnode.lg[t] = lon
                    hnode.la[t] = lat
                    hnode.H[t]  = alt / 1000

# def calculate_key_rate(method, links, fog, rain, snow, syst):
#     NUM_LINKS = len(links)
    
#     # Initialize with None for every link so indices match exactly
#     K_MAX = [None] * NUM_LINKS

#     d_list = []
#     h_list = []

#     for idx_l, l in enumerate(links):
#         # Identify which node is HAP and which is GS
#         if isinstance(l.n1, hap) and not isinstance(l.n2, hap):
#             hap_node, gs_node = l.n1, l.n2
#         elif isinstance(l.n2, hap) and not isinstance(l.n1, hap):
#             hap_node, gs_node = l.n2, l.n1
#         else:
#             # Skip links that are not HAP–GS
#             continue

#         # print(f"hap_node.la: {hap_node.la}")
#         # print(f"hap_node.lg: {hap_node.lg}")
#         # print(f"hap_node.H: {hap_node.H}")

#         # Time-varying HAP coordinates
#         la_rad_h = [math.radians(hap_node.la[t]) for t in syst.T]
#         lg_rad_h = [math.radians(hap_node.lg[t]) for t in syst.T]
#         H_h      = [hap_node.H[t] for t in syst.T]  # in km

#         # Static GS coordinates
#         la_rad_g = math.radians(gs_node.la)
#         lg_rad_g = math.radians(gs_node.lg)

#         x_g = R * math.cos(la_rad_g) * math.cos(lg_rad_g)
#         y_g = R * math.cos(la_rad_g) * math.sin(lg_rad_g)

#         # Now compute time-varying HAP coordinates
#         x_h = [R * math.cos(la_rad_h[t]) * math.cos(lg_rad_h[t]) for t in syst.T]
#         y_h = [R * math.cos(la_rad_h[t]) * math.sin(lg_rad_h[t]) for t in syst.T]

#         # Horizontal distances
#         d_los_hor = [math.sqrt((x_h[t] - x_g) ** 2 + (y_h[t] - y_g) ** 2) for t in syst.T]

#         # Elevation angles
#         alpha = [math.atan(H_h[t] / d_los_hor[t]) if d_los_hor[t] > 0 else math.pi / 2 for t in syst.T]

#         # LOS distances
#         d_los = [H_h[t] / math.sin(alpha[t]) for t in syst.T]

#         # print(f"syst.T: {syst.T}")
#         # print(f"d_los: {d_los}")
#         # print(f"alpha: {alpha}")

#         if idx_l == 0:
#             d_list = d_los
#             h_list = H_h
    
#             plot_skr("downlink", 5, d_list, h_list)

#         if method == "plob":
#             # Losses
#             L_geo = [20 * max(math.log10((R_TX + d_los[t] * 1000 * THETA) / R_RX), 0) for t in syst.T]
#             L_ma  = [0.01 * d_los[t] for t in range(len(syst.T))]
    
#             R_C = [H_C / math.sin(alpha[t]) for t in syst.T]
#             R_W = [H_W / math.sin(alpha[t]) for t in syst.T]
    
#             U = [1.6 if l.V[t] > 50 else (1.3 if 6 < l.V[t] <= 50 else 0.585 * l.V[t] ** (1 / 3)) for t in syst.T]
            
#             L_fog  = [(3.91 / l.V[t]) * ((LAMBDA / 550) ** (-U[t])) * R_C[t] for t in syst.T]
#             L_snw  = [(58   / l.V[t]) * R_W[t] for t in syst.T]
#             L_rain = [(2.8  / l.V[t]) * R_W[t] for t in syst.T]
    
#             # Total loss per time
#             L_t = [L_geo[t] + L_ma[t] + L_fog[t] * fog[t] + L_snw[t] * snow[t] + L_rain[t] * rain[t] for t in syst.T]
    
#             # Channel transmission
#             ETA = [10 ** (-L_t[t] / 10) for t in syst.T]
    
#             # Key rate over time
#             K_link = [-ts.ratesources * ts.sourceeff * math.log2(1 - ETA[t]) for t in syst.T]

#             # print(f"K_link: {K_link}")
    
#             # Save same key rate for both directions
#             K_MAX[idx_l] = K_link
#         elif method == "theoretical":
#             # start timer
#             t0 = time.perf_counter()
#             dir = ""
#             if isinstance(l.n1, gs) and isinstance(l.n2, hap): ## Uplink
#                 dir = "uplink"
#             else: ## Downlink
#                 dir = "downlink"
#             eta_theory = [ts.channel_theory(direction=dir, gs_alt=0, balloon_alt=H_h[t], distance=d_los[t], n_correction=6) for t in syst.T]
#             K_MAX[idx_l] = [ts.compute_skr(eta_theory[t]) for t in syst.T]
            
#             sys.stdout.write("\rProcessing... " + str(idx_l))
#             sys.stdout.flush()
#             #print(f"eta_theory: {eta_theory}, K_MAX[{idx_l}]: {K_MAX[idx_l]}, d_los: {d_los}")
#         elif method == "simulation":
#             eta_sim      = [ts.simulated_eff(distance=d_los[t], h_balloons=H_h[t], n=6) for t in syst.T]
#             K_MAX[idx_l] = [ts.compute_skr(eta_sim[t]) for t in syst.T]
#             #print(f"eta_sim: {eta_sim}, K_MAX[{idx_l}]: {K_MAX[idx_l]}, d_los: {d_los}")
#     return K_MAX

def _compute_key_point(task):
    """Worker for one (link, time) pair."""
    idx_l, t, l, syst, method, fog, rain, snow = task

    # Identify which node is HAP and GS
    if isinstance(l.n1, hap) and not isinstance(l.n2, hap):
        hap_node, gs_node = l.n1, l.n2
    elif isinstance(l.n2, hap) and not isinstance(l.n1, hap):
        hap_node, gs_node = l.n2, l.n1
    else:
        return idx_l, t, None  # skip non HAP–GS links

    # HAP coordinates
    la_rad_h = math.radians(hap_node.la[t])
    lg_rad_h = math.radians(hap_node.lg[t])
    H_h      = hap_node.H[t]

    # GS coordinates
    la_rad_g = math.radians(gs_node.la)
    lg_rad_g = math.radians(gs_node.lg)

    # Cartesian
    x_g = R * math.cos(la_rad_g) * math.cos(lg_rad_g)
    y_g = R * math.cos(la_rad_g) * math.sin(lg_rad_g)
    x_h = R * math.cos(la_rad_h) * math.cos(lg_rad_h)
    y_h = R * math.cos(la_rad_h) * math.sin(lg_rad_h)

    # Horizontal distance
    d_los_hor = math.sqrt((x_h - x_g) ** 2 + (y_h - y_g) ** 2)
    alpha     = math.atan(H_h / d_los_hor) if d_los_hor > 0 else math.pi / 2
    d_los     = H_h / math.sin(alpha)

    if method == "plob":
        # Geometric and medium losses
        L_geo = 20 * max(math.log10((R_TX + d_los * 1000 * THETA) / R_RX), 0)
        L_ma  = 0.01 * d_los

        R_C = H_C / math.sin(alpha)
        R_W = H_W / math.sin(alpha)

        if l.V[t] > 50:
            U = 1.6
        elif 6 < l.V[t] <= 50:
            U = 1.3
        else:
            U = 0.585 * l.V[t] ** (1 / 3)

        L_fog  = (3.91 / l.V[t]) * ((LAMBDA / 550) ** (-U)) * R_C
        L_snw  = (58   / l.V[t]) * R_W
        L_rain = (2.8  / l.V[t]) * R_W

        L_t = L_geo + L_ma + L_fog * fog[t] + L_snw * snow[t] + L_rain * rain[t]
        ETA = 10 ** (-L_t / 10)
        return idx_l, t, -ts.ratesources * ts.sourceeff * math.log2(1 - ETA)

    elif method == "theoretical":
        dir = "uplink" if isinstance(l.n1, gs) and isinstance(l.n2, hap) else "downlink"
        eta_theory = ts.channel_theory(direction=dir, gs_alt=0, balloon_alt=H_h,
                                       distance=d_los, n_correction=6)
        return idx_l, t, ts.compute_skr(eta_theory)

    elif method == "simulation":
        dir = "uplink" if isinstance(l.n1, gs) and isinstance(l.n2, hap) else "downlink"
        #eta_sim = ts.simulated_eff(distance=d_los, h_balloons=H_h, n=6)
        eta_sim = ts.channel_simulation(direction=dir, gs_alt=0, balloon_alt=H_h,
                                        distance=d_los, n_correction=6)
        # print(f"eta_sim: {eta_sim}")
        return idx_l, t, ts.compute_skr(eta_sim)

    return idx_l, t, None

def calculate_key_rate(method, links, fog, rain, snow, syst, max_workers=20):
    NUM_LINKS = len(links)
    K_MAX = [None] * NUM_LINKS

    # Build tasks
    tasks = []
    los_store = {}   # store d_los and H_h for link 0
    for idx_l, l in enumerate(links):
        if isinstance(l.n1, hap) and not isinstance(l.n2, hap):
            hap_node, gs_node = l.n1, l.n2
        elif isinstance(l.n2, hap) and not isinstance(l.n1, hap):
            hap_node, gs_node = l.n2, l.n1
        else:
            continue

        # Precompute geometry for plotting (first link only)
        if idx_l == 0:
            la_rad_h = [math.radians(hap_node.la[t]) for t in syst.T]
            lg_rad_h = [math.radians(hap_node.lg[t]) for t in syst.T]
            H_h      = [hap_node.H[t] for t in syst.T]

            la_rad_g = math.radians(gs_node.la)
            lg_rad_g = math.radians(gs_node.lg)
            x_g = R * math.cos(la_rad_g) * math.cos(lg_rad_g)
            y_g = R * math.cos(la_rad_g) * math.sin(lg_rad_g)

            x_h = [R * math.cos(la_rad_h[t]) * math.cos(lg_rad_h[t]) for t in syst.T]
            y_h = [R * math.cos(la_rad_h[t]) * math.sin(lg_rad_h[t]) for t in syst.T]
            d_los_hor = [math.sqrt((x_h[t] - x_g) ** 2 + (y_h[t] - y_g) ** 2) for t in syst.T]
            alpha = [math.atan(H_h[t] / d_los_hor[t]) if d_los_hor[t] > 0 else math.pi / 2 for t in syst.T]
            d_los = [H_h[t] / math.sin(alpha[t]) for t in syst.T]

            los_store["d_list"] = d_los
            los_store["h_list"] = H_h

        # Tasks for each time slot
        for t in syst.T:
            tasks.append((idx_l, t, l, syst, method, fog, rain, snow))

    # Parallel execution
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for idx_l, t, k_val in tqdm(executor.map(_compute_key_point, tasks, chunksize=10),
                                    total=len(tasks)):
            if k_val is None:
                continue
            if K_MAX[idx_l] is None:
                K_MAX[idx_l] = [None] * len(syst.T)
            K_MAX[idx_l][t] = k_val

    # Do the plot AFTER everything is finished
    if 0 in range(NUM_LINKS) and "d_list" in los_store:
        plot_skr("downlink", 6, los_store["d_list"], los_store["h_list"])

    return K_MAX


#######################################
####### Wind-related functions ########
#######################################

def compute_direction(u, v):
    return (np.degrees(np.arctan2(v, u)) + 360) % 360

def compute_speed(u, v):
    return np.sqrt(u**2 + v**2)

def download_data():
    lat, lon         = 49.0, 279.0  # Your location
    year, month, day = "2025", "07", "12"

    if os.path.exists(file_name):
        return

    c = cdsapi.Client()

    # Request only 4 time points: 00, 06, 12, 18
    hours = ["00:00", "06:00", "12:00", "18:00"]

    # Use a small bounding box around the location
    delta = 0.125  # One grid point in ERA5
    area = [lat + delta, lon - delta, lat - delta, lon + delta]  # N, W, S, E

    c.retrieve(
        "reanalysis-era5-pressure-levels",
        {
            "product_type": "reanalysis",
            "format": "netcdf",
            "variable": ["u_component_of_wind", "v_component_of_wind"],
            "pressure_level": [level],
            "year": year,
            "month": month,
            "day": day,
            "time": hours,
            "area": area,
        },
        file_name
    )