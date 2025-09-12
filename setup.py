from libraries import *

SYNTH_STRATO    = 0    ## 0: Synthetic, 1: Stratotegic Data

KEY_RATE_SCALE  = 1e-6
NUM_TIME_SLOTS  = 50 if SYNTH_STRATO else 4

R           = 6371  # Earth's radius in km
B           = 1e8   # Pulse rate
R_TX        = 0.1
R_RX        = 0.4
LAMBDA      = 1550e-9
THETA       = 1.22 * LAMBDA / R_TX
H_W         = 5   # Altitude of snow/rain
H_C         = 0.5 # Altitude of fog
V_fog       = 0.5 # 0.5km visibility (Moderate fog)
V_rain_snow = 4   # 4km visibility   (Moderate rain/snow)

## T     --> 12 subcarriers with 15 kHz spacing for 1 ms interval
## THETA --> Frame is 10ms --> 10, 20, 40, 80, 160  -- 3GPP TS 38.211
## G     --> 2 GSs and 2 HAPs with full connectivity
sys = system(range(NUM_TIME_SLOTS), 1, np.array([[1, 1]]))

level     = "50"  # hPa level (~20 km altitude)
file_name = f"era5_{level}hpa_hourly.nc"

def generate_keyrate_demands(
    hours=24,
    step_min=30,
    profiles=("enterprise", "residential", "critical"),
    base_mbps=(5.0, 3.0, 1.0),          # baseline secret-key rate per profile (Mb/s)
    peak_factor=(3.0, 2.5, 1.2),        # peak / base multiplier
    noise_std_frac=0.10,                # Gaussian noise as fraction of value
    floor_mbps=0.05,                    # never drop below this Mb/s
    seed=42
):
    """
    Returns:
        t_hours: array of length T (hours)
        demand_dict: {profile_name: np.array of shape (T,) in Mb/s}
        df: DataFrame with columns ['time_h', <profile...>] for plotting/debug
    """
    rng = np.random.default_rng(seed)
    T = int(hours * 60 / step_min)
    t = np.arange(T) * (step_min / 60.0)  # hours from 0..24

    demand_dict = {}
    for name, base, peak in zip(profiles, base_mbps, peak_factor):
        # Base diurnal shapes
        if name.lower() == "enterprise":
            # High 9:00–17:00, low at night
            # Sinusoid peaking ~14:00 (phase shift), scaled to [base, base*peak]
            phase = -4.0  # hours shift so peak ~14:00
            shape = 0.5 + 0.5 * np.sin(2*np.pi*(t + phase)/24)
        elif name.lower() == "residential":
            # Low daytime, peak in the evening (~20:00)
            phase = 2.0   # peak later
            shape = 0.5 + 0.5 * np.sin(2*np.pi*(t + phase)/24)
            shape = 1.0 - shape  # invert
        elif name.lower() == "critical":
            # Nearly flat (e.g., control/monitoring), mild diurnal wiggle
            shape = 0.8 + 0.2 * np.sin(2*np.pi*(t - 1.0)/24)
        else:
            # Generic diurnal
            shape = 0.5 + 0.5 * np.sin(2*np.pi*t/24)

        # Scale shape to [base, base*peak]
        shape = (shape - shape.min()) / (shape.max() - shape.min() + 1e-12)
        demand = base + (base*peak - base) * shape

        # Add noise proportional to the level
        noise = rng.normal(0.0, noise_std_frac, size=T) * demand
        demand = np.clip(demand + noise, floor_mbps, None)

        demand_dict[name] = demand

    df = pd.DataFrame({"time_h": t, **{k: v for k, v in demand_dict.items()}})
    return t, demand_dict, df

def get_neighbors(adj_matrix, node, node_type, gnodes, hnodes):
    if node_type == 0: ## GS
        return [hnodes[j] for j, connected in enumerate(adj_matrix[:,gnodes.index(node)]) if connected == 1]
    else: ## HAP
        return [gnodes[j] for j, connected in enumerate(adj_matrix[hnodes.index(node)]) if connected == 1]

def compute_direction(u, v):
    return (np.degrees(np.arctan2(v, u)) + 360) % 360

def compute_speed(u, v):
    return np.sqrt(u**2 + v**2)

def download_data():
    lat, lon         = 51.0, 4.0  # Your location
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
    return [lon, lat]

def update_coordinates(method, hnodes):
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
        
        hnodes_xy = np.zeros((len(sys.T),2))
    elif method == "stratotegic":
        # Load CSV
        df = pd.read_csv("dataset/balloon_sim_data.csv")
        # Keep only rows where Time_s is an integer
        df_filtered = df[df["Time_s"] % 1 == 0].reset_index(drop=True)
        # Extract required columns
        result = df_filtered[["Time_s", "Longitude_deg", "Latitude_deg", "Altitude_m"]]
    else:
        raise ValueError("Method must be either 'wind' or 'strategic'")
    
    for t in sys.T[1:]:
        if method == "wind":
            for hnode in hnodes:
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
                x_new = x + u_val * sys.THETA
                y_new = y + v_val * sys.THETA
    
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
                lon = row["Longitude_deg"].values[0]
                lat = row["Latitude_deg"].values[0]
                alt = row["Altitude_m"].values[0]
    
                for hnode in hnodes:
                    hnode.lg[t] = lon
                    hnode.la[t] = lat
                    hnode.H[t]  = alt

def calculate_key_rate(links, fog, rain, snow):
    NUM_LINKS = len(links)
    
    # Initialize with None for every link so indices match exactly
    K_MAX = [None] * NUM_LINKS

    for idx_l, l in enumerate(links):
        # Identify which node is HAP and which is GS
        if isinstance(l.n1, hap) and not isinstance(l.n2, hap):
            hap_node, gs_node = l.n1, l.n2
        elif isinstance(l.n2, hap) and not isinstance(l.n1, hap):
            hap_node, gs_node = l.n2, l.n1
        else:
            # Skip links that are not HAP–GS
            continue

        # Time-varying HAP coordinates
        la_rad_h = [math.radians(hap_node.la[t]) for t in sys.T]
        lg_rad_h = [math.radians(hap_node.lg[t]) for t in sys.T]
        H_h      = [hap_node.H[t] for t in sys.T]  # in km

        # Static GS coordinates
        la_rad_g = math.radians(gs_node.la)
        lg_rad_g = math.radians(gs_node.lg)

        x_g = R * math.cos(la_rad_g) * math.cos(lg_rad_g)
        y_g = R * math.cos(la_rad_g) * math.sin(lg_rad_g)

        # Now compute time-varying HAP coordinates
        x_h = [R * math.cos(la_rad_h[t]) * math.cos(lg_rad_h[t]) for t in sys.T]
        y_h = [R * math.cos(la_rad_h[t]) * math.sin(lg_rad_h[t]) for t in sys.T]

        # Horizontal distances
        d_los_hor = [math.sqrt((x_h[t] - x_g) ** 2 + (y_h[t] - y_g) ** 2) for t in sys.T]

        # Elevation angles
        alpha = [math.atan(H_h[t] / d_los_hor[t]) if d_los_hor[t] > 0 else math.pi / 2 for t in sys.T]

        # LOS distances
        d_los = [H_h[t] / math.sin(alpha[t]) for t in range(len(sys.T))]

        # Losses
        L_geo = [20 * max(math.log10((R_TX + d_los[t] * 1000 * THETA) / R_RX), 0) for t in sys.T]
        L_ma  = [0.01 * d_los[t] for t in range(len(sys.T))]

        R_C = [H_C / math.sin(alpha[t]) for t in sys.T]
        R_W = [H_W / math.sin(alpha[t]) for t in sys.T]

        U = [1.6 if l.V[t] > 50 else (1.3 if 6 < l.V[t] <= 50 else 0.585 * l.V[t] ** (1 / 3)) for t in sys.T]
        
        L_fog  = [(3.91 / l.V[t]) * ((LAMBDA / 550) ** (-U[t])) * R_C[t] for t in sys.T]
        L_snw  = [(58   / l.V[t]) * R_W[t] for t in sys.T]
        L_rain = [(2.8  / l.V[t]) * R_W[t] for t in sys.T]

        # Total loss per time
        L_t = [L_geo[t] + L_ma[t] + L_fog[t] * fog[t] + L_snw[t] * snow[t] + L_rain[t] * rain[t] for t in range(len(sys.T))]

        # Channel transmission
        ETA = [10 ** (-L_t[t] / 10) for t in range(len(sys.T))]

        # Key rate over time
        K_link = [-B * math.log2(1 - ETA[t]) for t in range(len(sys.T))]

        # Save same key rate for both directions
        K_MAX[idx_l] = K_link

    return K_MAX

def init_setup():
    # Process each node
    gnodes  = []
    hnodes  = []
    links   = []
    demands = []
    
    gnodes.append(gs(278.8, 49, 1e2, 1e2, 1e4))
    gnodes.append(gs(279.2, 49, 1e2, 1e2, 1e4))
    
    hnodes.append(hap([279]*len(sys.T), [49]*len(sys.T), [25]*len(sys.T), 1e2, 1e2, 1e4))
    
    if SYNTH_STRATO == 1:
        update_coordinates("stratotegic", hnodes)
    else:
        update_coordinates("wind", hnodes)
    
    NUM_LINKS = 8
    links.append(link(gnodes[0], hnodes[0], [100]*len(sys.T), [(0,0,0)]*len(sys.T), [1e6*KEY_RATE_SCALE]*len(sys.T)))
    links.append(link(gnodes[1], hnodes[0], [100]*len(sys.T), [(0,0,0)]*len(sys.T), [1e6*KEY_RATE_SCALE]*len(sys.T)))
    links.append(link(hnodes[0], gnodes[0], [100]*len(sys.T), [(0,0,0)]*len(sys.T), [1e6*KEY_RATE_SCALE]*len(sys.T)))
    links.append(link(hnodes[0], gnodes[1], [100]*len(sys.T), [(0,0,0)]*len(sys.T), [1e6*KEY_RATE_SCALE]*len(sys.T)))

    fog  = [0] * len(sys.T)
    rain = [0] * len(sys.T)
    snow = [0] * len(sys.T)
    K_MAX = calculate_key_rate(links, fog, rain, snow)
    
    for idx_l, l in enumerate(links):
        for t in sys.T:
            l.K_MAX[t] = K_MAX[idx_l][t]*KEY_RATE_SCALE
            
    t, demand_dict, df = generate_keyrate_demands(hours=1, step_min=1/60)

    # Pick a profile, e.g. "enterprise"
    k_req_vals = (KEY_RATE_SCALE * demand_dict["enterprise"] * 1e3).tolist()

    # Use in your demand object
    demands.append(
        demand(
            k_req_vals,
            gnodes[0],
            gnodes[1]
        )
    )

    return gnodes, hnodes, links, demands