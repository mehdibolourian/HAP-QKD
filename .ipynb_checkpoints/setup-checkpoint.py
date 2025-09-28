from libraries import *

SYNTH_STRATO    = 1    ## 0: Wind, 1: Stratotegic Data

COORDINATE_SCALE = 1 #1e-3
KEY_RATE_SCALE   = 1e-1
NUM_TIME_SLOTS   = 3 if SYNTH_STRATO else 4
STORAGE_SCALE    = 1

## T     --> 12 subcarriers with 15 kHz spacing for 1 ms interval
## THETA --> 1 sec.
## G     --> 2 GSs and 2 HAPs with full connectivity
syst = system(range(NUM_TIME_SLOTS), 1, np.array([[1, 1]]))

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
    K_MAX = calculate_key_rate("plob", links, fog, rain, snow, syst)

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

def init_setup_real():
    # Process each node
    gnodes  = []
    hnodes  = []
    links   = []
    demands = []
    
    # Ground Stations (longitude, latitude roughly approximated in degrees)
    # Venice, Padua, Florence, Siena
    gnodes.append(gs(278.6695, 48.4758, 1, 1, 1e9, "Timmins"))   # Timmins GS
    gnodes.append(gs(279.3186, 48.7669, 1, 1, 1e9, "IroquoisFalls"))   # IroquoisFalls GS - ~70 km northeast of Timmins
    gnodes.append(gs(277.5669, 49.4169, 1, 1, 1e9, "Kapuskasing"))   # Kapuskasing GS - ~160 km northwest of Timmins
    gnodes.append(gs(278.984, 49.0670, 1, 1, 1e9, "Cochrane"))    # Cochrane GS - ~110 km north of Timmins
    gnodes.append(gs(279.9674, 48.1512, 1, 1, 1e9, "KirklandLake"))   # KirklandLake GS - ~140 km southeast of Timmins
    
    # HAPs at 35 km altitude above Padua and Florence
    hnodes.append(hap([279]*len(syst.T), [49]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9))  # Stratotegic coordinates
    hnodes.append(hap([277.85]*len(syst.T), [49.34]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9))  # Moonbeam town center
    
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

    # for l in links:
    #     print(f"idx: {links.index(l)}, l_n1_tag: {l.n1.tag}, l_n2_tag: {l.n2.tag}")

    plot_connectivity_graph(gnodes, hnodes, links)
    animate_hap_trajectories(syst.T, [hnode.lg for hnode in hnodes], [hnode.la for hnode in hnodes], [f"HAP_{idx_hnode}" for idx_hnode, hnode in enumerate(hnodes)])

    fog   = [0] * len(syst.T)
    rain  = [0] * len(syst.T)
    snow  = [0] * len(syst.T)
    K_MAX = calculate_key_rate("plob", links, fog, rain, snow, syst) # method: "plob", "theoretical", "simulation"

    #print(f"K_MAX:{K_MAX}")

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
            
    #t, demand_dict, df = generate_keyrate_demands(hours=1, step_min=1/60)

    # Pick a profile, e.g. "enterprise"
    k_req_vals = [0.02] * len(syst.T) # 0.02 bits/sec

    # Use in your demand object
    demands.append(
        demand(
            k_req_vals,
            gnodes[0],
            gnodes[1]
        )
    )

    return gnodes, hnodes, links, demands