from libraries import *

case = 1
SYNTH_STRATO    = 1    ## 0: Wind, 1: Stratotegic Data

COORDINATE_SCALE = 1
KEY_RATE_SCALE   = 1
# if case == 1:
#     NUM_TIME_SLOTS = 144 if SYNTH_STRATO else 4
# else:
#     NUM_TIME_SLOTS = 288 if SYNTH_STRATO else 4
NUM_TIME_SLOTS = 48 if SYNTH_STRATO else 4
STORAGE_SCALE    = 1

MODEL_KEY_RATE   = "theoretical" # "plob", "theoretical", "simulation"

## T     --> Set of time steps
## THETA --> 60 sec.
## G     --> 2 GSs and 2 HAPs with full connectivity
# if case == 1:
#     syst = system(range(NUM_TIME_SLOTS), 600, np.array([[1, 1]]))
# else:
#     syst = system(range(NUM_TIME_SLOTS), 300, np.array([[1, 1]]))
syst = system(range(NUM_TIME_SLOTS), 1800, np.array([[1, 1]]))

level     = "50"  # hPa level (~20 km altitude)
file_name = f"era5_{level}hpa_hourly.nc"

def compute_zenith_stats(links):
    R = 6371.0  # Earth radius (km)

    def haversine(lon1, lat1, lon2, lat2):
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
        c = 2*np.arcsin(np.sqrt(a))
        return R * c

    for idx, l in enumerate(links):

        # Only meaningful for GS-HAP links
        if not hasattr(l.n1, "H") and not hasattr(l.n2, "H"):
            continue

        zenith_list = []

        for t in syst.T:

            # determine GS and HAP
            if hasattr(l.n1, "H"):   # n1 is HAP
                hap_node = l.n1
                gs_node = l.n2
            else:
                hap_node = l.n2
                gs_node = l.n1

            lon_g = gs_node.lg
            lat_g = gs_node.la

            lon_h = hap_node.lg[t]
            lat_h = hap_node.la[t]
            alt_h = hap_node.H[t]  # km (your HAP altitude list)

            horiz_dist = haversine(lon_g, lat_g, lon_h, lat_h)

            zenith = np.degrees(np.arctan2(horiz_dist, alt_h))

            zenith_list.append(zenith)

        zenith_arr = np.array(zenith_list)

        print(f"Link {idx}: {l.n1.tag} -> {l.n2.tag}")
        print(f"   min zenith : {zenith_arr.min():.2f}°")
        print(f"   max zenith : {zenith_arr.max():.2f}°")
        print(f"   avg zenith : {zenith_arr.mean():.2f}°")
        print()

def set_link_active_by_zenith(links, link_active, ZENITH_THRESHOLD=82, p_outage=0.05, seed=42):
    R = 6371.0  # Earth radius in km
    rng = np.random.default_rng(seed)

    def haversine(lon1, lat1, lon2, lat2):
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
        c = 2*np.arcsin(np.sqrt(a))
        return R * c  # km

    for l_idx, l in enumerate(links):

        for t in syst.T:

            # Identify GS and HAP
            if hasattr(l.n1, "H"):  # n1 is HAP
                hap = l.n1
                gs  = l.n2
            elif hasattr(l.n2, "H"):  # n2 is HAP
                hap = l.n2
                gs  = l.n1
            else:
                continue

            lon_g, lat_g = gs.lg, gs.la
            lon_h, lat_h = hap.lg[t], hap.la[t]
            alt_h        = hap.H[t]

            # horizontal distance
            d = haversine(lon_g, lat_g, lon_h, lat_h)

            # zenith angle
            zenith = np.degrees(np.arctan2(d, alt_h))

            # deactivate if zenith too large
            if zenith > ZENITH_THRESHOLD:
                link_active[l_idx, t] = 0
                continue

            # 5% random outage regardless
            if rng.random() < p_outage:
                link_active[l_idx, t] = 0

    return link_active

def init_setup():
    # Process each node
    gnodes  = []
    hnodes  = []
    links   = []
    demands = []
    
    gnodes.append(gs(278.8, 49, 1e2, 1e2, 1e4))
    gnodes.append(gs(279.2, 49, 1e2, 1e2, 1e4))
    
    hnodes.append(hap([279]*len(syst.T), [49]*len(syst.T), [25]*len(syst.T), 1e2, 1e2, 1e4))
    
    if SYNTH_STRATO == 1:
        update_coordinates("stratotegic", hnodes, syst)
    else:
        update_coordinates("wind", hnodes, syst)
    
    links.append(link(gnodes[0], hnodes[0], [100]*len(syst.T), [(0,0,0)]*len(syst.T), [1e6]*len(syst.T)))
    links.append(link(gnodes[1], hnodes[0], [100]*len(syst.T), [(0,0,0)]*len(syst.T), [1e6]*len(syst.T)))
    links.append(link(hnodes[0], gnodes[0], [100]*len(syst.T), [(0,0,0)]*len(syst.T), [1e6]*len(syst.T)))
    links.append(link(hnodes[0], gnodes[1], [100]*len(syst.T), [(0,0,0)]*len(syst.T), [1e6]*len(syst.T)))

    plot_connectivity_graph(gnodes, hnodes, links)

    fog  = [0] * len(syst.T)
    rain = [0] * len(syst.T)
    snow = [0] * len(syst.T)
    K_MAX = calculate_key_rate(MODEL_KEY_RATE, links, fog, rain, snow, syst)

    # Compute efficiencies
    eta_theory = ts.theoretical_eff(distance=25, h_balloons=15, n=5)
    eta_sim    = ts.simulated_eff(distance=25, h_balloons=15, n=5)
    
    # Compute SKRs
    skr_theory = ts.compute_skr(eta_theory)
    skr_sim    = ts.compute_skr(eta_sim)
    
    # print(f"Theoretical efficiency: {eta_theory:.4f} -> SKR: {skr_theory:.2f} kbit/s")
    # print(f"Simulated  efficiency: {eta_sim:.4f} -> SKR: {skr_sim:.2f} kbit/s")
    
    for idx_l, l in enumerate(links):
        for t in syst.T:
            l.K_MAX[t] = K_MAX[idx_l][t]
            
    t, demand_dict, df = generate_keyrate_demands(hours=1, step_min=1/60)

    # Pick a profile, e.g. "enterprise"
    k_req_vals = (demand_dict["enterprise"] * 1e3).tolist()

    # Use in your demand object
    demands.append(
        demand(
            k_req_vals,
            gnodes[0],
            gnodes[1]
        )
    )

    return gnodes, hnodes, links, demands

def init_setup_real_skr_plot():
    # Process each node
    gnodes  = []
    hnodes  = []
    links   = []
    demands = []

    # Ground Stations (longitude, latitude roughly approximated in degrees)
    gnodes.append(gs(279, 49, 1, 1, 1e9, "GS")) # GS
    
    # # HAPs at 15 km altitude above Padua and Florence
    hnodes.append(hap([279]*len(syst.T), [49]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP"))  # Stratotegic coordinates
    
    # Update coordinates depending on model choice
    if SYNTH_STRATO == 1:
        update_coordinates("stratotegic", hnodes, syst)
    else:
        update_coordinates("wind", hnodes, syst)
    
    # Links: connect only GSs to HAPs
    for gs_node in gnodes:
        for hap_node in hnodes:
            links.append(link(gs_node, hap_node,
                              [100]*len(syst.T),
                              [(0,0,0)]*len(syst.T),
                              [1e6]*len(syst.T)))
            links.append(link(hap_node, gs_node,
                              [100]*len(syst.T),
                              [(0,0,0)]*len(syst.T),
                              [1e6]*len(syst.T)))

    #plot_connectivity_graph(gnodes, hnodes, links)
    plot_connectivity_graph_3d(gnodes, hnodes, links)

    return gnodes, hnodes, links, demands

def init_setup_real_planning():
    # Process each node
    gnodes  = []
    hnodes  = []
    links   = []
    demands = []

    # Area: 176.2 km x 140.9 km
    
    
    # Ground Stations (longitude, latitude roughly approximated in degrees)
    gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))         # Timmins GS
    gnodes.append(gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"))  # IroquoisFalls GS - ~70 km northeast of Timmins
    gnodes.append(gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"))     # Kapuskasing GS - ~160 km northwest of Timmins
    gnodes.append(gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"))        # Cochrane GS - ~110 km north of Timmins
    gnodes.append(gs(279.9674, 48.1512, 1, 1, 1e9, "Kirkland Lake"))   # KirklandLake GS - ~140 km southeast of Timmins
    # ---- Additional Ground Stations ----
    gnodes.append(gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay"))        # ~280 km south of Timmins
    # gnodes.append(gs(280.9906, 49.6850, 1, 1, 1e9, "Moosonee"))        # ~400 km north of Timmins
    
    # # HAPs at 15 km altitude above Padua and Florence
    hnodes.append(hap([277.06]*len(syst.T), [47.85]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_0"))  # Stratotegic coordinates
    hnodes.append(hap([277.85]*len(syst.T), [47.60]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_1"))  # Moonbeam town center
    # # HAP3 between Timmins and KirklandLake
    hnodes.append(hap([278.26]*len(syst.T), [47.22]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_2"))
    # # ---- Additional HAP ----
    hnodes.append(hap([277.61]*len(syst.T), [47.38]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_3")) # between Sudbury and North Bay
    hnodes.append(hap([278.69]*len(syst.T), [48.77]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_4")) # between Sudbury and North Bay
    hnodes.append(hap([278.04]*len(syst.T), [48.93]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_5")) # between Sudbury and North Bay
    hnodes.append(hap([278.28]*len(syst.T), [49.15]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_6")) # between Sudbury and North Bay
    hnodes.append(hap([279.48]*len(syst.T), [48.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_7")) # between Sudbury and North Bay
    hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_8")) # between Sudbury and North Bay
    hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_9")) # between Sudbury and North Bay
    hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_10")) # between Sudbury and North Bay
    hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_11")) # between Sudbury and North Bay
    hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_12")) # between Sudbury and North Bay
    hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_13")) # between Sudbury and North Bay
    hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_14")) # between Sudbury and North Bay
    hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_15")) # between Sudbury and North Bay
    
    # Update coordinates depending on model choice
    if SYNTH_STRATO == 1:
        update_coordinates("stratotegic", hnodes, syst)
    else:
        update_coordinates("wind", hnodes, syst)
    
    # Links: connect only GSs to HAPs
    for gs_node in gnodes:
        for hap_node in hnodes:
            links.append(link(gs_node, hap_node,
                              [100]*len(syst.T),
                              [(0,0,0)]*len(syst.T),
                              [1e6]*len(syst.T)))
            links.append(link(hap_node, gs_node,
                              [100]*len(syst.T),
                              [(0,0,0)]*len(syst.T),
                              [1e6]*len(syst.T)))

    for hap_node in hnodes:
        print(hap_node.H)

    # for l in links:
    #     print(f"idx: {links.index(l)}, l_n1_tag: {l.n1.tag}, l_n2_tag: {l.n2.tag}")

    plot_connectivity_graph(gnodes, hnodes, links)
    plot_connectivity_graph_3d(gnodes, hnodes, links)
    # animate_hap_trajectories(syst.T, [hnode.lg for hnode in hnodes], [hnode.la for hnode in hnodes], [f"HAP_{idx_hnode}" for idx_hnode, hnode in enumerate(hnodes)])

    fog   = [0] * len(syst.T)
    rain  = [0] * len(syst.T)
    snow  = [0] * len(syst.T)
    # N = 1
    # start = 86300
    # for n in range(N):
    #     K_MAX = calculate_key_rate(MODEL_KEY_RATE, links, fog, rain, snow, syst,
    #                                max_workers=24, chunk_size=1,
    #                                start_chunk=start + n*10, end_chunk=start + (n+1)*10,
    #                                checkpoint_file="K_MAX_checkpoint.pkl") # method: "plob", "theoretical", "simulation"

    # K_MAX, _ = calculate_key_rate_mac(MODEL_KEY_RATE, links, fog, rain, snow, syst) # method: "plob", "theoretical", "simulation"

    #print(f"K_MAX:{K_MAX}")

    # # Compute efficiencies
    # eta_theory = ts.theoretical_eff(distance=25, h_balloons=15, n=5)
    # eta_sim    = ts.simulated_eff(distance=25, h_balloons=15, n=5)
    
    # # Compute SKRs
    # skr_theory = ts.compute_skr(eta_theory)
    # skr_sim    = ts.compute_skr(eta_sim)
    
    # print(f"Theoretical efficiency: {eta_theory:.4f} -> SKR: {skr_theory:.2f} kbit/s")
    # print(f"Simulated  efficiency: {eta_sim:.4f} -> SKR: {skr_sim:.2f} kbit/s")
    
    # for idx_l, l in enumerate(links):
    #     for t in syst.T:
    #         l.K_MAX[t] = K_MAX[idx_l][t]
            
    #t, demand_dict, df = generate_keyrate_demands(hours=1, step_min=1/60)

    # Pick a profile, e.g. "enterprise"
    k_req_vals = [0.1] * len(syst.T) # 25600 bits/sec
    
    # Use in your demand object
    demands.append(
        demand(
            k_req_vals,
            gnodes[0], #gnodes[2],
            gnodes[1]  #gnodes[3]
        )
    )
    # demands.append(
    #     demand(
    #         k_req_vals,
    #         gnodes[2], #gnodes[2],
    #         gnodes[3]  #gnodes[3]
    #     )
    # )
    # demands.append(
    #     demand(
    #         k_req_vals,
    #         gnodes[1], #gnodes[2],
    #         gnodes[2]  #gnodes[3]
    #     )
    # )

    demands = generate_demands(gnodes, syst, mean_kbps=100, amp=20, noise_std=2, pattern="sinusoidal")
    # demands = generate_demands(gnodes, syst, mean_kbps=5, amp=1, noise_std=0, pattern="sinusoidal")

    # --- Plot all demands ---
    plt.figure(figsize=(8, 5))
    for d in demands:
        src_idx = gnodes.index(d.n1)
        dst_idx = gnodes.index(d.n2)
        plt.plot(
            syst.T, d.K_REQ,
            lw=1.6,
            label=f"GS{src_idx} ↔ GS{dst_idx}"
        )
    plt.xlabel("Time step (t)")
    plt.ylabel("Key Rate Demand (kb/sec)")
    plt.title("Generated GS–GS Demands over Time")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.show()
    # ------------------------

    return gnodes, hnodes, links, demands

def init_setup_real_offline(prob):
    # Process each node
    gnodes  = []
    hnodes  = []
    links   = []
    demands = []

    # Area: 176.2 km x 140.9 km
    # Ground Stations (longitude, latitude roughly approximated in degrees)
    if prob == 1: ## Only one link (To study QKP's impact directly)
        # gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))         # Timmins GS
        # gnodes.append(gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"))  # IroquoisFalls GS - ~70 km northeast of Timmins


        
        gnodes.append(gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"))     # Kapuskasing GS - ~160 km northwest of Timmins
        # ---- Additional Ground Stations ----
        gnodes.append(gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay"))        # ~280 km south of Timmins



        # gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))         # Timmins GS
        # gnodes.append(gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"))        # Cochrane GS - ~110 km north of Timmins

        
        
        #hnodes.append(hap([278.69]*len(syst.T), [48.77]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_5")) # between Sudbury and North Bay
        hnodes.append(hap([277.06]*len(syst.T), [47.85]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_0"))  # Stratotegic coordinates
    elif prob == 2:
        gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))         # Timmins GS
        gnodes.append(gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"))  # IroquoisFalls GS - ~70 km northeast of Timmins
        gnodes.append(gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"))     # Kapuskasing GS - ~160 km northwest of Timmins
        gnodes.append(gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"))        # Cochrane GS - ~110 km north of Timmins
        gnodes.append(gs(279.9674, 48.1512, 1, 1, 1e9, "Kirkland Lake"))   # KirklandLake GS - ~140 km southeast of Timmins
        # ---- Additional Ground Stations ----
        gnodes.append(gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay"))        # ~280 km south of Timmins
        # gnodes.append(gs(280.9906, 49.6850, 1, 1, 1e9, "Moosonee"))        # ~400 km north of Timmins
        
        # # # HAPs at 15 km altitude above Padua and Florence
        hnodes.append(hap([277.06]*len(syst.T), [47.85]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_0"))  # Stratotegic coordinates
        hnodes.append(hap([277.77]*len(syst.T), [47.67]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_1"))  # Moonbeam town center
        # # # HAP3 between Timmins and KirklandLake
        hnodes.append(hap([277.94]*len(syst.T), [47.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_2"))
        # # ---- Additional HAP ----
        hnodes.append(hap([278.26]*len(syst.T), [47.21]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_3")) # between Sudbury and North Bay
        hnodes.append(hap([277.61]*len(syst.T), [47.38]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_4")) # between Sudbury and North Bay
        hnodes.append(hap([278.69]*len(syst.T), [48.77]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_5")) # between Sudbury and North Bay
        # hnodes.append(hap([278.28]*len(syst.T), [49.15]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_6")) # between Sudbury and North Bay
        # hnodes.append(hap([279.48]*len(syst.T), [48.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_7")) # between Sudbury and North Bay
        # hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_8")) # between Sudbury and North Bay
        # hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_9")) # between Sudbury and North Bay
        
    # Update coordinates depending on model choice
    if SYNTH_STRATO == 1:
        update_coordinates("stratotegic", hnodes, syst)
    else:
        update_coordinates("wind", hnodes, syst)
    
    # Links: connect only GSs to HAPs
    for gs_node in gnodes:
        for hap_node in hnodes:
            links.append(link(gs_node, hap_node,
                              [100]*len(syst.T),
                              [(0,0,0)]*len(syst.T),
                              [1e6]*len(syst.T)))
            links.append(link(hap_node, gs_node,
                              [100]*len(syst.T),
                              [(0,0,0)]*len(syst.T),
                              [1e6]*len(syst.T)))

    # plot_connectivity_graph(gnodes, hnodes, links)

    fog   = [0] * len(syst.T)
    rain  = [0] * len(syst.T)
    snow  = [0] * len(syst.T)

    K_MAX, _ = calculate_key_rate_mac(MODEL_KEY_RATE, links, fog, rain, snow, syst) # method: "plob", "theoretical", "simulation"

    
    for idx_l, l in enumerate(links):
        for t in syst.T:
            l.K_MAX[t] = K_MAX[idx_l][t]
            
    #t, demand_dict, df = generate_keyrate_demands(hours=1, step_min=1/60)

    # Pick a profile, e.g. "enterprise"
    k_req_vals = [0.1] * len(syst.T) # 25600 bits/sec
    
    # Use in your demand object
    # demands.append(
    #     demand(
    #         k_req_vals,
    #         gnodes[0],
    #         gnodes[1]
    #     )
    # )

    #demands = [generate_demands(gnodes, syst, mean_kbps=100, amp=0.1, noise_std=0, pattern="sinusoidal")[0]]
    demands = generate_demands(gnodes, syst, mean_kbps=0.5, amp=0.1, noise_std=0, pattern="stadium")
    
    # --- Plot all demands ---
    plt.figure(figsize=(8, 5))
    for d in demands:
        src_idx = gnodes.index(d.n1)
        dst_idx = gnodes.index(d.n2)
        plt.plot(
            syst.T, d.K_REQ,
            lw=1.6,
            label=f"GS{src_idx} ↔ GS{dst_idx}"
        )
    plt.xlabel("Time step (t)")
    plt.ylabel("Key Rate Demand (kb/sec)")
    plt.title("Generated GS–GS Demands over Time")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(fontsize=8)
    plt.tight_layout()

    plt.savefig("demand.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()
    # ------------------------

    return gnodes, hnodes, links, demands

def init_setup_real_online(prob):
    # Process each node
    gnodes  = []
    hnodes  = []
    links   = []
    demands = []

    # Area: 176.2 km x 140.9 km
    # Ground Stations (longitude, latitude roughly approximated in degrees)
    if prob == 1: ## Only one link (To study QKP's impact directly)
        # gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))         # Timmins GS
        # gnodes.append(gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"))  # IroquoisFalls GS - ~70 km northeast of Timmins


        
        gnodes.append(gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"))     # Kapuskasing GS - ~160 km northwest of Timmins
        # ---- Additional Ground Stations ----
        gnodes.append(gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay"))        # ~280 km south of Timmins



        # gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))         # Timmins GS
        # gnodes.append(gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"))        # Cochrane GS - ~110 km north of Timmins

        
        
        #hnodes.append(hap([278.69]*len(syst.T), [48.77]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_5")) # between Sudbury and North Bay
        hnodes.append(hap([277.06]*len(syst.T), [47.85]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_0"))  # Stratotegic coordinates
    elif prob == 2:
        if case == 1:
            gnodes.extend([
                gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"),
                gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"),
                gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"),
                gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"),
                gs(279.9674, 48.1512, 1, 1, 1e9, "Kirkland Lake"),
                gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay")
            ])
        
            hnodes.extend([
                hap([277.06]*len(syst.T), [47.85]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_0"),
                hap([277.77]*len(syst.T), [47.67]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_1"),
                hap([277.94]*len(syst.T), [47.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_2"),
                hap([278.26]*len(syst.T), [47.21]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_3"),
                hap([277.61]*len(syst.T), [47.38]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_4"),
                hap([278.69]*len(syst.T), [48.77]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_5")
            ])
        elif case == 2:
            gnodes.extend([
                gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"),
                gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"),
                gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay")
            ])
        
            hnodes.extend([
                hap([277.94]*len(syst.T), [47.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_2"),
                hap([278.69]*len(syst.T), [48.77]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_5")
            ])
        else:
            gnodes.extend([
                gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"),
                gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"),
                gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"),
                gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay")
            ])
        
            hnodes.extend([
                hap([277.06]*len(syst.T), [47.85]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_0"),
                hap([277.94]*len(syst.T), [47.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_2"),
                hap([278.69]*len(syst.T), [48.77]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_5")
            ])
        
    # Update coordinates depending on model choice
    if SYNTH_STRATO == 1:
        update_coordinates("stratotegic", hnodes, syst)
    else:
        update_coordinates("wind", hnodes, syst)
    
    # Links: connect only GSs to HAPs
    for gs_node in gnodes:
        for hap_node in hnodes:
            links.append(link(gs_node, hap_node,
                              [100]*len(syst.T),
                              [(0,0,0)]*len(syst.T),
                              [1e6]*len(syst.T)))
            links.append(link(hap_node, gs_node,
                              [100]*len(syst.T),
                              [(0,0,0)]*len(syst.T),
                              [1e6]*len(syst.T)))

    # plot_connectivity_graph(gnodes, hnodes, links)

    fog   = [0] * len(syst.T)
    rain  = [0] * len(syst.T)
    snow  = [0] * len(syst.T)

    K_MAX, _ = calculate_key_rate_mac(MODEL_KEY_RATE, links, fog, rain, snow, syst) # method: "plob", "theoretical", "simulation"

    
    for idx_l, l in enumerate(links):
        for t in syst.T:
            l.K_MAX[t] = K_MAX[idx_l][t]
            
    #t, demand_dict, df = generate_keyrate_demands(hours=1, step_min=1/60)

    # Pick a profile, e.g. "enterprise"
    k_req_vals = [0.1] * len(syst.T) # 25600 bits/sec
    
    # Use in your demand object
    # demands.append(
    #     demand(
    #         k_req_vals,
    #         gnodes[0],
    #         gnodes[1]
    #     )
    # )

    #demands = [generate_demands(gnodes, syst, mean_kbps=100, amp=0.1, noise_std=0, pattern="sinusoidal")[0]]
    demands = generate_demands(gnodes, syst, mean_kbps=0.5, amp=0.5, noise_std=0.0, pattern="stadium")
    
    # --- Plot all demands ---
    plt.figure(figsize=(8, 5))
    for d in demands:
        src_idx = gnodes.index(d.n1)
        dst_idx = gnodes.index(d.n2)
        plt.plot(
            syst.T, d.K_REQ,
            lw=1.6,
            label=f"GS{src_idx} ↔ GS{dst_idx}"
        )
    plt.xlabel("Time step (t)")
    plt.ylabel("Key Rate Demand (kb/sec)")
    plt.title("Generated GS–GS Demands over Time")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(fontsize=8)
    plt.tight_layout()

    plt.savefig("demand.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()
    # ------------------------

    return gnodes, hnodes, links, demands

def generate_master_displacement(num_steps):
    """Generates a single relative movement path (delta_lon, delta_lat, delta_alt)."""
    # Assuming the first HAP's lat for degree-to-km conversion
    lat_to_km = 111.0
    lon_to_km = 111.0 * np.cos(np.radians(47.5)) # Approximation for the region
    max_radius_km = 50.0
    
    # Start at 0 displacement
    d_lons, d_lats, d_alts = [0.0], [0.0], [0.0]
    vx, vy = np.random.uniform(-0.015, 0.015), np.random.uniform(-0.015, 0.015)
    
    curr_dx, curr_dy = 0.0, 0.0
    
    for _ in range(num_steps - 1):
        vx = 0.85 * vx + 0.15 * np.random.uniform(-0.1, 0.1)
        vy = 0.85 * vy + 0.15 * np.random.uniform(-0.1, 0.1)
        
        new_dx = curr_dx + vx
        new_dy = curr_dy + vy
        
        # Check boundary relative to displacement origin
        dist_km = np.sqrt((new_dx * lon_to_km)**2 + (new_dy * lat_to_km)**2)
        
        if dist_km < max_radius_km:
            curr_dx, curr_dy = new_dx, new_dy
        else:
            vx, vy = -vx * 0.8, -vy * 0.8
            
        d_lons.append(curr_dx)
        d_lats.append(curr_dy)
        d_alts.append(1.5 * np.sin(_ / 20.0) + np.random.uniform(-0.1, 0.1))
        
    return np.array(d_lons), np.array(d_lats), np.array(d_alts)

def init_setup_training_online(plot_trajs=True):
    gnodes = []
    hnodes = []
    links = []
    
    # if case == 1:
    #     gnodes.extend([
    #         gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"),
    #         gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"),
    #         gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"),
    #         gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"),
    #         gs(279.9674, 48.1512, 1, 1, 1e9, "Kirkland Lake"),
    #         gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay")
    #     ])
    
    #     hap_seeds = [
    #         (277.06, 47.85, 15, "HAP_0"),
    #         (277.77, 47.67, 15, "HAP_1"),
    #         (277.94, 47.52, 15, "HAP_2"),
    #         (278.26, 47.21, 15, "HAP_3"),
    #         (277.61, 47.38, 15, "HAP_4"),
    #         (278.69, 48.77, 15, "HAP_5")
    #     ]
    # elif case == 2:
    #     gnodes.extend([
    #         gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"),
    #         gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"),
    #         gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay")
    #     ])
    
    #     hap_seeds = [
    #         (277.94, 47.52, 15, "HAP_2"),
    #         (278.69, 48.77, 15, "HAP_5")
    #     ]
    # else:
    #     gnodes.extend([
    #         gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"),
    #         gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"),
    #         gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"),
    #         gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay")
    #     ])
    
    #     hap_seeds = [
    #         (277.06, 47.85, 15, "HAP_0"),
    #         (277.94, 47.52, 15, "HAP_2"),
    #         (278.69, 48.77, 15, "HAP_5")
    #     ]

    gnodes.extend([])

    hap_seeds = [
        (278.69, 48.77, 15, "HAP_5")
    ]

    num_steps = len(syst.T)
    
    # GENERATE MASTER TRAJECTORY ONCE
    d_lons, d_lats, d_alts = generate_master_displacement(num_steps)

    for lon_0, lat_0, alt_0, name in hap_seeds:
        # Shift the master trajectory by the seed coordinates
        final_lons = lon_0 + d_lons
        final_lats = lat_0 + d_lats
        final_alts = alt_0 + d_alts
        
        hnodes.append(hap(final_lons.tolist(), final_lats.tolist(), final_alts.tolist(), 1, 1, 1e9, name))

    if plot_trajs:
        # ---------------- 3D Plot ----------------
        fig1 = plt.figure(figsize=(4, 3))
        ax = fig1.add_subplot(111, projection='3d')
        
        # ==========================
        # Transparent cylinder
        # ==========================
        cyl_radius=0.65,      # cylinder radius (lon/lat units)
        cyl_alpha=0.05       # transparency
        
        # Cylinder center (mean lon/lat)
        all_lons = [h.lg[0] for h in hnodes]
        all_lats = [h.la[0] for h in hnodes]
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
    
        for h in hnodes:
            ax.plot(h.lg, h.la, h.H, label=h.tag, alpha=0.7, color="red")

        ax.set_xlabel("Longitude", fontsize=14)
        ax.set_ylabel("Latitude", fontsize=14)
        ax.set_zlabel("Altitude (km)", fontsize=14)
    
        gs_lon = [g.lg for g in gnodes]
        gs_lat = [g.la for g in gnodes]

        # ==========================
        # View angle (important)
        # ==========================
        ax.view_init(elev=25, azim=100)
    
        plt.tight_layout()
        plt.savefig("gmmm_3d.pdf", format="pdf", bbox_inches="tight", pad_inches=0.01, transparent="True")
        plt.show()
    
    
        # ---------------- 2D Plot ----------------
        fig2 = plt.figure(figsize=(4, 3))
        ax2 = fig2.add_subplot(111)
    
        for h in hnodes:
            ax2.plot(h.lg, h.la, label=h.tag, linewidth=2, color="red")
            ax2.scatter(h.lg[0], h.la[0], s=40, color='red')
    
        ax2.grid(True, linestyle='--', alpha=0.5)

        plt.xlabel("Longitude", fontsize=14)
        plt.ylabel("Latitude", fontsize=14)
    
        plt.tight_layout()
        plt.savefig("gmmm_2d.pdf", format="pdf", bbox_inches="tight", pad_inches=0.01, transparent="True")
        plt.show()

    # Links: connect only GSs to HAPs
    for gs_node in gnodes:
        for hap_node in hnodes:
            links.append(link(gs_node, hap_node, [100]*num_steps, [(0,0,0)]*num_steps, [1e6]*num_steps))
            links.append(link(hap_node, gs_node, [100]*num_steps, [(0,0,0)]*num_steps, [1e6]*num_steps))

    # Calculate dynamic key rates (Crucial: distance changes every timestep)
    K_MAX, _ = calculate_key_rate_mac(MODEL_KEY_RATE, links, 0, 0, 0, syst)
    
    for idx_l, l in enumerate(links):
        for t in syst.T:
            l.K_MAX[t] = K_MAX[idx_l][t]

    demands = generate_demands(gnodes, syst, mean_kbps=0.5, amp=0.5, noise_std=0, pattern="stadium")
    
    return gnodes, hnodes, links, demands