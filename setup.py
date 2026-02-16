from libraries import *

SYNTH_STRATO    = 1    ## 0: Wind, 1: Stratotegic Data

COORDINATE_SCALE = 1
KEY_RATE_SCALE   = 1
NUM_TIME_SLOTS   = 96 if SYNTH_STRATO else 4
STORAGE_SCALE    = 1

MODEL_KEY_RATE   = "theoretical" # "plob", "theoretical", "simulation"

## T     --> Set of time steps
## THETA --> 60 sec.
## G     --> 2 GSs and 2 HAPs with full connectivity
syst = system(range(NUM_TIME_SLOTS), 900, np.array([[1, 1]]))

level     = "50"  # hPa level (~20 km altitude)
file_name = f"era5_{level}hpa_hourly.nc"

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
        #gnodes.append(gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"))  # IroquoisFalls GS - ~70 km northeast of Timmins
        gnodes.append(gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"))     # Kapuskasing GS - ~160 km northwest of Timmins
        gnodes.append(gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"))        # Cochrane GS - ~110 km north of Timmins
        #gnodes.append(gs(279.9674, 48.1512, 1, 1, 1e9, "Kirkland Lake"))   # KirklandLake GS - ~140 km southeast of Timmins
        # ---- Additional Ground Stations ----
        gnodes.append(gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay"))        # ~280 km south of Timmins
        # gnodes.append(gs(280.9906, 49.6850, 1, 1, 1e9, "Moosonee"))        # ~400 km north of Timmins
        
        # # # HAPs at 15 km altitude above Padua and Florence
        hnodes.append(hap([277.06]*len(syst.T), [47.85]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_0"))  # Stratotegic coordinates
        #hnodes.append(hap([277.77]*len(syst.T), [47.67]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_1"))  # Moonbeam town center
        # # # HAP3 between Timmins and KirklandLake
        hnodes.append(hap([277.94]*len(syst.T), [47.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_2"))
        # # ---- Additional HAP ----
        #hnodes.append(hap([278.26]*len(syst.T), [47.21]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_3")) # between Sudbury and North Bay
        #hnodes.append(hap([277.61]*len(syst.T), [47.38]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_4")) # between Sudbury and North Bay
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
    demands = generate_demands(gnodes, syst, mean_kbps=0.5, amp=0.1, noise_std=0, pattern="sinusoidal")
    
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
        gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))         # Timmins GS
        gnodes.append(gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"))  # IroquoisFalls GS - ~70 km northeast of Timmins
        gnodes.append(gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"))     # Kapuskasing GS - ~160 km northwest of Timmins
        # gnodes.append(gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"))        # Cochrane GS - ~110 km north of Timmins
        # gnodes.append(gs(279.9674, 48.1512, 1, 1, 1e9, "Kirkland Lake"))   # KirklandLake GS - ~140 km southeast of Timmins
        # # ---- Additional Ground Stations ----
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
        
        # gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))         # Timmins GS
        # gnodes.append(gs(279.3186, 48.7669, 1, 1, 1e9, "Iroquois Falls"))  # IroquoisFalls GS - ~70 km northeast of Timmins
        # gnodes.append(gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"))     # Kapuskasing GS - ~160 km northwest of Timmins
        # gnodes.append(gs(278.984,  49.0670, 1, 1, 1e9, "Cochrane"))        # Cochrane GS - ~110 km north of Timmins
        # gnodes.append(gs(279.9674, 48.1512, 1, 1, 1e9, "Kirkland Lake"))   # KirklandLake GS - ~140 km southeast of Timmins
        # # ---- Additional Ground Stations ----
        # gnodes.append(gs(276.7073, 46.3091, 1, 1, 1e9, "NorthBay"))        # ~280 km south of Timmins
        # # gnodes.append(gs(280.9906, 49.6850, 1, 1, 1e9, "Moosonee"))        # ~400 km north of Timmins
        
        # # # # HAPs at 15 km altitude above Padua and Florence
        # hnodes.append(hap([277.06]*len(syst.T), [47.85]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_0"))  # Stratotegic coordinates
        # hnodes.append(hap([277.77]*len(syst.T), [47.67]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_1"))  # Moonbeam town center
        # # # # HAP3 between Timmins and KirklandLake
        # hnodes.append(hap([277.94]*len(syst.T), [47.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_2"))
        # # # ---- Additional HAP ----
        # hnodes.append(hap([278.26]*len(syst.T), [47.21]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_3")) # between Sudbury and North Bay
        # hnodes.append(hap([277.61]*len(syst.T), [47.38]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_4")) # between Sudbury and North Bay
        # hnodes.append(hap([278.69]*len(syst.T), [48.77]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_5")) # between Sudbury and North Bay
        # # hnodes.append(hap([278.28]*len(syst.T), [49.15]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_6")) # between Sudbury and North Bay
        # # hnodes.append(hap([279.48]*len(syst.T), [48.52]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_7")) # between Sudbury and North Bay
        # # hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_8")) # between Sudbury and North Bay
        # # hnodes.append(hap([275.8]*len(syst.T), [47.6]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, "HAP_9")) # between Sudbury and North Bay
        
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
    demands = generate_demands(gnodes, syst, mean_kbps=0.5, amp=0.1, noise_std=0, pattern="sinusoidal")
    
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