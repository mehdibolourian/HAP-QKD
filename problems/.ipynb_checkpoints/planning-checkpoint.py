from libraries import *

def calculate_key_rate_planning(method, link, d_los, t):
    K_MAX = 0
    if method == "plob":
        L_geo = 20 * max(math.log10((R_TX + d_los * 1000 * THETA) / R_RX), 0)
        L_ma  = 0.01 * d_los
        L_t   = L_geo + L_ma
        
        ETA = 10**(-L_t/10)
        
        K_MAX = -ts.ratesources * math.log2(1 - ETA)
    elif method == "theoretical":
        # Compute efficiencies
        eta_theory = ts.channel_theory(direction="downlink", gs_alt=0, balloon_alt=15,
                                       distance=d_los, n_correction=6)
        # Compute SKRs
        K_MAX = ts.compute_skr(eta_theory)
    elif method == "simulation":
        # Compute efficiencies
        eta_sim = ts.channel_simulation(direction="downlink", gs_alt=0, balloon_alt=15,
                                       distance=d_los, n_correction=6)
        # Compute SKRs
        K_MAX   = ts.compute_skr(eta_sim)
    
    return K_MAX

## problem: "FREE", "FIXED", "INIT"
def planning(problem, gss, haps, links):
    # Create Optimization Model
    m = gp.Model("hap-qkd")

    demands = []
    for idx_g1, g1 in enumerate(gss):
        for idx_g2, g2 in enumerate(gss):
            if idx_g1 < idx_g2:
                demands.append(demand(0, g1, g2))

    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    r, r_l, z = {}, {}, {}

    for idx_d, d in enumerate(demands):
        for t in syst.T:
            r[idx_d, t] = m.addVar(name=f"r_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

    for idx_l, l in enumerate(links):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                z[idx_l, idx_d, t]   = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}", vtype=GRB.BINARY)
                #r_l[idx_l, idx_d, t] = m.addVar(name=f"r_l_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

    # k (key rate) and d (distance)
    k, d = {}, {}
    for idx_l, l in enumerate(links):
        if problem == "FIXED" or problem == "INIT":
            d[idx_l] = m.addVar(name=f"d_{idx_l}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
        for t in syst.T:
            if problem == "FREE":
                d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                
            k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

            dpts = np.linspace(15 * COORDINATE_SCALE, 2e2 * COORDINATE_SCALE, 4)
            kpts = [calculate_key_rate_planning("plob", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

            if problem == "FREE":
                m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")
            elif problem == "FIXED" or problem == "INIT":
                m.addGenConstrPWL(d[idx_l], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

    d_min, d_max = 15.0, 200.0
    breakpoint_counts = [5, 10, 20, 40]  # numbers of breakpoints to test
    plt.figure(figsize=(8,5))
    for n_bp in breakpoint_counts:
        # breakpoints for this case
        dpts = np.linspace(d_min, d_max, n_bp)
        #kpts = [calculate_key_rate_planning(0, d) * KEY_RATE_SCALE for d in dpts]
        kpts = [calculate_key_rate_planning("plob", 0, d, 0) * KEY_RATE_SCALE for d in dpts]
    
        # plot this curve
        plt.plot(dpts, kpts, marker="o", linestyle="-", label=f"{n_bp} breakpoints")
    plt.xlabel("Distance (dpts)")
    plt.ylabel("Key Rate (kpts)")
    plt.title("Key Rate vs Distance for Different Breakpoint Counts")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
            
    # Coordinate decision variables for each HAP and time
    c1, c2 = {}, {}  # x, y in km
    for idx_h, hnode in enumerate(haps):
        if problem == "FIXED" or problem == "INIT":
            c1[idx_h] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}")
            c2[idx_h] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}")
        elif problem == "FREE":
            for t in syst.T:
                c1[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
                c2[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")

    
    # Secondary objective: maximize d
    if problem == "FIXED" or problem == "INIT":
        m.setObjective(sum(sum(r[idx_d, t]
                               for idx_d, d in enumerate(demands)
                              )
                            for t in syst.T
                          ), GRB.MAXIMIZE)

    elif problem == "FREE":
        m.setObjectiveN(-sum(sum(d[idx_l, t]
                                 for idx_l, l in enumerate(links)
                                )
                             for t in syst.T
                            ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")

    m.setParam("MIPGap", 1e-3)          # force very tight gap
    m.setParam("MIPGapAbs", 1e-3)
    m.setParam("FeasibilityTol", 1e-3)
    m.setParam("IntFeasTol", 1e-3)
    m.setParam("OptimalityTol", 1e-3)

    m.Params.Presolve = 2           # Aggressive presolve
    m.Params.Method = 2             # Barrier for root relaxation
    m.Params.Cuts = 2               # More cutting planes
    m.Params.Heuristics = 0.25      # Reasonable level of primal heuristics
    m.Params.MIPFocus = 1           # Emphasize feasible solutions first
    m.Params.NumericFocus = 1       # Reduce rounding issues due to tiny epsilons
    m.Params.Threads = 0            # Use all cores

    m.Params.NodefileStart = 0.5    # Start writing node files to disk to prevent RAM overflow
    m.Params.NoRelHeurTime = 20     # Time for heuristics before branching
    m.Params.ConcurrentMIP = 1      # Use multiple concurrent solvers for first feasible solution

    ## Constraints
    # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
    # r_h = min_{l:z_l=1}(r_1+r_2)
    m.addConstrs(
        (
            r[idx_d, t] <= k[idx_l, t] + 1e8 * KEY_RATE_SCALE * (1 - z[idx_l, idx_d, t])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="demand_link_coordination"
    )

    # Flow conservation
    m.addConstrs(
        (
            sum(z[idx_l, idx_d, t]
                for idx_l, l in enumerate(links)
                if isinstance(l.n1, gs)
                if gss.index(l.n1) == gss.index(d.n1)
               ) - sum(z[idx_l, idx_d, t]
                       for idx_l, l in enumerate(links)
                       if isinstance(l.n2, gs)
                       if gss.index(l.n2) == gss.index(d.n1)
                      ) == 1
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="flow_conservation_1"
    )
    m.addConstrs(
        (
            sum(z[idx_l, idx_d, t]
                for idx_l, l in enumerate(links)
                if isinstance(l.n2, gs)
                if gss.index(l.n2) == gss.index(d.n2)
               ) - sum(z[idx_l, idx_d, t]
                       for idx_l, l in enumerate(links)
                       if isinstance(l.n1, gs)
                       if gss.index(l.n1) == gss.index(d.n2)
                      ) == 1
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="flow_conservation_2"
    )
    m.addConstrs(
        (
            sum(z[idx_l, idx_d, t]
                for idx_l, l in enumerate(links)
                if  l.n1 == n
               ) - sum(z[idx_l, idx_d, t]
                       for idx_l, l in enumerate(links)
                       if  l.n2 == n
                      ) == 0
            for idx_d, d in enumerate(demands)
            for n in gss + haps
            if  n != d.n1 and n != d.n2
            for t        in syst.T
        ), name="flow_conservation_3"
    )
    m.addConstrs(
        (
            z[idx_l_1, idx_d, t] + z[idx_l_2, idx_d, t] <= 1
            for idx_l_1, l_1 in enumerate(links)
            for idx_l_2, l_2 in enumerate(links)
            if  idx_l_1 < idx_l_2 and (l_1.n1 == l_2.n2) and (l_1.n2 == l_2.n1)
            for idx_d, d     in enumerate(demands)
            for t            in syst.T
        ), name="loop_prevention"
    )
    
    # Maximum Tx/Rx Connection
    m.addConstrs(
        (
            sum(z[idx_l, idx_d, t]
                for idx_l, l in enumerate(links)
                if l.n1 == n
               ) <= n.N_TX
            for idx_n, n in enumerate(gss + haps)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="max_tx_connections"
    )
    
    m.addConstrs(
        (
            sum(z[idx_l, idx_d, t]
                for idx_l, l in enumerate(links)
                if l.n2 == n
               ) <= n.N_RX
            for idx_n, n in enumerate(gss + haps)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="max_rx_connections"
    )

    # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
    for idx_l, l in enumerate(links):
        # identify which endpoint is HAP and which is GS
        if isinstance(l.n1, hap) and isinstance(l.n2, gs):
            hap_idx, gs_node = haps.index(l.n1), l.n2
        elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
            hap_idx, gs_node = haps.index(l.n2), l.n1

        # precompute GS coordinates in same (x,y,z) frame and units (km)
        [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)  # ensure these are constants in km

        if problem == "FIXED" or problem == "INIT":
            # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
            dx = c1[hap_idx] - cg1*COORDINATE_SCALE
            dy = c2[hap_idx] - cg2*COORDINATE_SCALE
            m.addQConstr(d[idx_l]*d[idx_l] >= dx*dx + dy*dy + 15*15*COORDINATE_SCALE*COORDINATE_SCALE,
                         name=f"dist_cone_{idx_l}")
        elif problem == "FREE":
            for t in syst.T:
                # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
                dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
                dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
                m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + 15*15*COORDINATE_SCALE*COORDINATE_SCALE,
                             name=f"dist_cone_{idx_l}_{t}")

    ## Solve
    m.optimize()

    if m.status == GRB.OPTIMAL:
        print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        solution = {
            "r":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r.items()},
            #"r_l":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_l.items()},
            "z":   {k: v.X for k, v in z.items()},
            "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
            "d":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in d.items()},
            "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
            "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()}
        }
        

        pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
        pp.pprint(solution)

        # print(f"lonlat0:{solution["c1"][0,0]}, {solution["c2"][0,0]}, {xy_to_lonlat(solution["c1"][0,0], solution["c2"][0,0])}")

        # Actual HAP trajectories
        actual_lons = [hnode.lg for hnode in haps]
        actual_lats = [hnode.la for hnode in haps]
        actual_labels = [f"HAP_{idx_hnode}_actual" for idx_hnode, _ in enumerate(haps)]
        
        # Planned HAP trajectories
        planned_lons = []
        planned_lats = []
        for idx_hnode in range(len(haps)):
            lon_series = []
            lat_series = []
            for t in syst.T:
                x = solution["c1"].get((idx_hnode))
                y = solution["c2"].get((idx_hnode))
                if x is not None and y is not None:
                    # lon, lat = xy_to_lonlat(x, y)
                    lon, lat = tangent_to_latlon(x, y, 279, 49)
                    # lon_series.append(lon+360)   # shift if needed
                    lon_series.append(lon)   # shift if needed
                    lat_series.append(lat)
            planned_lons.append(lon_series)
            planned_lats.append(lat_series)
        planned_labels = [f"HAP_{idx_hnode}_planned" for idx_hnode in range(len(haps))]
        
        # Ground Stations (replicate coordinates across all T so they plot in animation)
        gs_lons = []
        gs_lats = []
        for gnode in gss:
            gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
            gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
        gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]
        
        # Combine everything
        # all_lons = actual_lons + planned_lons + gs_lons
        # all_lats = actual_lats + planned_lats + gs_lats
        # all_labels = actual_labels + planned_labels + gs_labels

        all_lons = planned_lons + gs_lons
        all_lats = planned_lats + gs_lats
        all_labels = planned_labels + gs_labels
        
        # Animate
        animate_hap_trajectories(
            syst.T,
            all_lons,
            all_lats,
            all_labels
        )
    else:
        print("No optimal solution found.")
        solution = None

    return solution



    
# def planning(problem, gss, haps, links, demands):
#     # Create Optimization Model
#     m = gp.Model("hap-qkd")

#     ## Decision Variables

#     # k (key rate) and d (distance)
#     k, d = {}, {}
#     for idx_l, l in enumerate(links):
#         if problem == "FIXED" or problem == "INIT":
#             d[idx_l] = m.addVar(name=f"d_{idx_l}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=5e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
#         for t in syst.T:
#             if problem == "FREE":
#                 d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=5e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                
#             k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#             dpts = np.linspace(15 * COORDINATE_SCALE, 1e3 * COORDINATE_SCALE, 5)
#             kpts = [calculate_key_rate_planning("plob", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

#             if problem == "FREE":
#                 m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")
#             elif problem == "FIXED" or problem == "INIT":
#                 m.addGenConstrPWL(d[idx_l], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

#     d_min, d_max = 15.0, 500.0
#     breakpoint_counts = [5, 10, 20, 40]  # numbers of breakpoints to test
#     plt.figure(figsize=(8,5))
#     for n_bp in breakpoint_counts:
#         # breakpoints for this case
#         dpts = np.linspace(d_min, d_max, n_bp)
#         #kpts = [calculate_key_rate_planning(0, d) * KEY_RATE_SCALE for d in dpts]
#         kpts = [calculate_key_rate_planning("plob", 0, d, 0) * KEY_RATE_SCALE for d in dpts]
    
#         # plot this curve
#         plt.plot(dpts, kpts, marker="o", linestyle="-", label=f"{n_bp} breakpoints")
#     plt.xlabel("Distance (dpts)")
#     plt.ylabel("Key Rate (kpts)")
#     plt.title("Key Rate vs Distance for Different Breakpoint Counts")
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.show()
            
#     # Coordinate decision variables for each HAP and time
#     c1, c2, c3 = {}, {}, {}  # x, y, z in km
#     for idx_h, hnode in enumerate(haps):
#         if problem == "FIXED" or problem == "INIT":
#             c1[idx_h] = m.addVar(lb=-3e2*COORDINATE_SCALE, ub=3e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}")
#             c2[idx_h] = m.addVar(lb=-3e2*COORDINATE_SCALE, ub=3e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}")
#             c3[idx_h] = m.addVar(lb=15*COORDINATE_SCALE,   ub=25*COORDINATE_SCALE,  vtype=GRB.CONTINUOUS, name=f"c3_{idx_h}")  # altitude in km
#         elif problem == "FREE":
#             for t in syst.T:
#                 c1[idx_h, t] = m.addVar(lb=-3e2*COORDINATE_SCALE, ub=3e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
#                 c2[idx_h, t] = m.addVar(lb=-3e2*COORDINATE_SCALE, ub=3e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")
#                 c3[idx_h, t] = m.addVar(lb=15*COORDINATE_SCALE,   ub=25*COORDINATE_SCALE,  vtype=GRB.CONTINUOUS, name=f"c3_{idx_h}_{t}")  # altitude in km

#     f = {}
#     for idx_g, g in enumerate(gss):
#         for t in syst.T:
#             f[t] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, name=f"f_{t}")
#     ## Objective
    
#     # Secondary objective: maximize d
#     if problem == "FIXED" or problem == "INIT":
#         # m.setObjectiveN(-sum(d[idx_l]
#         #                      for idx_l, l in enumerate(links)
#         #                     ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")
#         # m.setObjective(sum(d[idx_l] for idx_l, l in enumerate(links)), GRB.MINIMIZE)
        
#         # m.setObjective(sum(sum(k[idx_l, t]
#         #                        for idx_l, l in enumerate(links)
#         #                       )
#         #                    for t in syst.T
#         #                   ), GRB.MAXIMIZE)

#         m.setObjective(sum(f[t]
#                            for t in syst.T
#                           ), GRB.MAXIMIZE)
#     elif problem == "FREE":
#         m.setObjectiveN(-sum(sum(d[idx_l, t]
#                                  for idx_l, l in enumerate(links)
#                                 )
#                              for t in syst.T
#                             ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")

#     m.setParam("MIPGap", 1e-6)          # force very tight gap
#     m.setParam("MIPGapAbs", 1e-6)
#     m.setParam("FeasibilityTol", 1e-6)
#     m.setParam("IntFeasTol", 1e-6)
#     m.setParam("OptimalityTol", 1e-6)

#     ## Constraints
#     for t in syst.T:
#         for idx_l, l in enumerate(links):
#             m.addQConstr(f[t] <= k[idx_l, t],
#                          name=f"min_flow_{idx_l}_{t}")

#     # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
#     for idx_l, l in enumerate(links):
#         # identify which endpoint is HAP and which is GS
#         if isinstance(l.n1, hap) and isinstance(l.n2, gs):
#             hap_idx, gs_node = haps.index(l.n1), l.n2
#         elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
#             hap_idx, gs_node = haps.index(l.n2), l.n1

#         # precompute GS coordinates in same (x,y,z) frame and units (km)
#         #[cg1, cg2] = lonlat_to_xy(gs_node.lg, gs_node.la)  # ensure these are constants in km
#         [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)  # ensure these are constants in km

#         #print(cg1, cg2)
#         # print(xy_to_lonlat(cg1, cg2))
#         #print(tangent_to_latlon(cg1, cg2, 279, 49))

#         if problem == "FIXED" or problem == "INIT":
#             # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
#             dx = c1[hap_idx] - cg1*COORDINATE_SCALE
#             dy = c2[hap_idx] - cg2*COORDINATE_SCALE
#             dz = c3[hap_idx]
#             m.addQConstr(d[idx_l]*d[idx_l] >= dx*dx + dy*dy + dz*dz*COORDINATE_SCALE*COORDINATE_SCALE,
#                          name=f"dist_cone_{idx_l}")
#         elif problem == "FREE":
#             for t in syst.T:
#                 # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
#                 dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
#                 dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
#                 dz = c3[hap_idx, t]
#                 m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + dz*dz*COORDINATE_SCALE*COORDINATE_SCALE,
#                              name=f"dist_cone_{idx_l}_{t}")

#     ## Solve
#     m.optimize()

#     if m.status == GRB.OPTIMAL:
#         print("\n=========== OPTIMAL SOLUTION FOUND ===========")

#         solution = {
#             "k":  {k: round(v.X, 3) for k, v in k.items()},
#             "d":  {k: round(v.X, 3) for k, v in d.items()},
#             "c1": {k: round(v.X, 3) for k, v in c1.items()},
#             "c2": {k: round(v.X, 3) for k, v in c2.items()},
#             "c3": {k: round(v.X, 3) for k, v in c3.items()}
#         }
#         print(solution["c1"], solution["c2"], solution["c3"])

#         # Actual HAP trajectories
#         actual_lons = [hnode.lg for hnode in haps]
#         actual_lats = [hnode.la for hnode in haps]
#         actual_labels = [f"HAP_{idx_hnode}_actual" for idx_hnode, _ in enumerate(haps)]
        
#         # Planned HAP trajectories
#         planned_lons = []
#         planned_lats = []
#         for idx_hnode in range(len(haps)):
#             lon_series = []
#             lat_series = []
#             for t in syst.T:
#                 x = solution["c1"].get(idx_hnode)
#                 y = solution["c2"].get(idx_hnode)
#                 if x is not None and y is not None:
#                     # lon, lat = xy_to_lonlat(x, y)
#                     lon, lat = tangent_to_latlon(x, y, 279, 49)
#                     # lon_series.append(lon+360)   # shift if needed
#                     lon_series.append(lon)   # shift if needed
#                     lat_series.append(lat)
#             planned_lons.append(lon_series)
#             planned_lats.append(lat_series)
#         planned_labels = [f"HAP_{idx_hnode}_planned" for idx_hnode in range(len(haps))]
        
#         # Ground Stations (replicate coordinates across all T so they plot in animation)
#         gs_lons = []
#         gs_lats = []
#         for gnode in gss:
#             gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
#             gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
#         gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]
        
#         # Combine everything
#         all_lons = actual_lons + planned_lons + gs_lons
#         all_lats = actual_lats + planned_lats + gs_lats
#         all_labels = actual_labels + planned_labels + gs_labels

#         # Animate
#         animate_hap_trajectories(
#             syst.T,
#             all_lons,
#             all_lats,
#             all_labels
#         )
#     else:
#         print("No optimal solution found.")
#         solution = None

#     return solution




    
# def planning(problem, gss, haps, links, demands):
#     # Create Optimization Model
#     m = gp.Model("hap-qkd")

#     ## Decision Variables
#     # # Dictionaries of decision variables instead of MVar arrays
#     # r_1, r_2, r_h, a, z = {}, {}, {}, {}, {}

#     # for idx_l, l in enumerate(links):
#     #     for idx_d, d in enumerate(demands):
#     #         for t in syst.T:
#     #             r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
#     #             r_2[idx_l, idx_d, t] = m.addVar(name=f"r_2_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
#     #             z[idx_l, idx_d, t]   = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}",   vtype=GRB.BINARY)

                
#     # for idx_d, d in enumerate(demands):
#     #     for t in syst.T:
#     #         r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
            
#     # for idx_l, l in enumerate(links):
#     #     for t in syst.T:
#     #         a[idx_l, t] = m.addVar(name=f"a_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#     # k (key rate) and d (distance)
#     k, d = {}, {}
#     for idx_l, l in enumerate(links):
#         if problem == "FIXED" or problem == "INIT":
#             d[idx_l] = m.addVar(name=f"d_{idx_l}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=5e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
#         for t in syst.T:
#             if problem == "FREE":
#                 d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=5e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                
#             k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#             dpts = np.linspace(15 * COORDINATE_SCALE, 1e3 * COORDINATE_SCALE, 5)
#             kpts = [calculate_key_rate_planning("plob", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

#             if problem == "FREE":
#                 m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")
#             elif problem == "FIXED" or problem == "INIT":
#                 m.addGenConstrPWL(d[idx_l], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

#     d_min, d_max = 15.0, 500.0
#     breakpoint_counts = [5, 10, 20, 40]  # numbers of breakpoints to test
#     plt.figure(figsize=(8,5))
#     for n_bp in breakpoint_counts:
#         # breakpoints for this case
#         dpts = np.linspace(d_min, d_max, n_bp)
#         #kpts = [calculate_key_rate_planning(0, d) * KEY_RATE_SCALE for d in dpts]
#         kpts = [calculate_key_rate_planning("plob", 0, d, 0) * KEY_RATE_SCALE for d in dpts]
    
#         # plot this curve
#         plt.plot(dpts, kpts, marker="o", linestyle="-", label=f"{n_bp} breakpoints")
#     plt.xlabel("Distance (dpts)")
#     plt.ylabel("Key Rate (kpts)")
#     plt.title("Key Rate vs Distance for Different Breakpoint Counts")
#     plt.grid(True)
#     plt.legend()
#     plt.tight_layout()
#     plt.show()
            
#     # Coordinate decision variables for each HAP and time
#     c1, c2, c3 = {}, {}, {}  # x, y, z in km
#     for idx_h, hnode in enumerate(haps):
#         if problem == "FIXED" or problem == "INIT":
#             c1[idx_h] = m.addVar(lb=-3e2*COORDINATE_SCALE, ub=3e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}")
#             c2[idx_h] = m.addVar(lb=-3e2*COORDINATE_SCALE, ub=3e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}")
#             c3[idx_h] = m.addVar(lb=15*COORDINATE_SCALE,   ub=25*COORDINATE_SCALE,  vtype=GRB.CONTINUOUS, name=f"c3_{idx_h}")  # altitude in km
#         elif problem == "FREE":
#             for t in syst.T:
#                 c1[idx_h, t] = m.addVar(lb=-3e2*COORDINATE_SCALE, ub=3e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
#                 c2[idx_h, t] = m.addVar(lb=-3e2*COORDINATE_SCALE, ub=3e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")
#                 c3[idx_h, t] = m.addVar(lb=15*COORDINATE_SCALE,   ub=25*COORDINATE_SCALE,  vtype=GRB.CONTINUOUS, name=f"c3_{idx_h}_{t}")  # altitude in km

#     ## Objective
#     # m.ModelSense = GRB.MAXIMIZE

#     # # Primary objective: maximize r_h
#     # m.setObjectiveN(sum(sum(r_h[idx_d, t]
#     #                         for idx_d, d in enumerate(demands)
#     #                        )
#     #                    for t in syst.T
#     #                   ) * syst.THETA, index=0, priority=2, weight=1.0, abstol=1e-6, reltol=1e-6, name="Primary")

#     # # Secondary objective: maximize a
#     # m.setObjectiveN(sum(sum(-a[idx_l, t]
#     #                        for idx_l, l in enumerate(links)
#     #                       )
#     #                    for t in syst.T
#     #                   ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")

#     # # Secondary objective: maximize k
#     # m.setObjectiveN(-sum(sum(d[idx_l]
#     #                          for idx_l, l in enumerate(links)
#     #                         )
#     #                      for t in syst.T
#     #                     ) * syst.THETA, index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")

    
#     # Secondary objective: maximize d
#     if problem == "FIXED" or problem == "INIT":
#         # m.setObjectiveN(-sum(d[idx_l]
#         #                      for idx_l, l in enumerate(links)
#         #                     ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")
#         # m.setObjective(sum(d[idx_l] for idx_l, l in enumerate(links)), GRB.MINIMIZE)
#         m.setObjective(sum(sum(k[idx_l, t]
#                                for idx_l, l in enumerate(links)
#                               )
#                            for t in syst.T
#                           ), GRB.MAXIMIZE)
#     elif problem == "FREE":
#         m.setObjectiveN(-sum(sum(d[idx_l, t]
#                                  for idx_l, l in enumerate(links)
#                                 )
#                              for t in syst.T
#                             ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")

#     m.setParam("MIPGap", 1e-6)          # force very tight gap
#     m.setParam("MIPGapAbs", 1e-6)
#     m.setParam("FeasibilityTol", 1e-6)
#     m.setParam("IntFeasTol", 1e-6)
#     m.setParam("OptimalityTol", 1e-6)

#     ## Constraints
#     # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
#     # r_h = min_{l:z_l=1}(r_1+r_2)
#     # m.addConstrs(
#     #     (
#     #         r_h[idx_d, t] <= r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] + d.K_REQ[t] * KEY_RATE_SCALE * (1 - z[idx_l, idx_d, t])
#     #         for idx_l, l in enumerate(links)
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="demand_link_coordination"
#     # )
    
#     # EPSILON = 1e-3
#     # # Key rate and routing coordination (1)
#     # m.addConstrs(
#     #     (
#     #         r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] >= EPSILON * z[idx_l, idx_d, t]
#     #         for idx_l, l in enumerate(links)
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="demand_link_coordination_1"
#     # )
    
#     # # Key rate and routing coordination (2)
#     # m.addConstrs(
#     #     (
#     #         r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] <= d.K_REQ[t] * KEY_RATE_SCALE * z[idx_l, idx_d, t]
#     #         for idx_l, l in enumerate(links)
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="key_rate_routing_coordination_2"
#     # )
    
#     # # Max key rate
#     # m.addConstrs(
#     #     (
#     #         sum(r_1[idx_l, idx_d, t]
#     #             for idx_d, d in enumerate(demands)
#     #            ) <= k[idx_l, t]
#     #         for idx_l, l in enumerate(links)
#     #         for t        in syst.T
#     #     ), name="max_key_rate"
#     # )
    
#     # # Flow conservation
#     # m.addConstrs(
#     #     (
#     #         sum(z[idx_l, idx_d, t]
#     #             for idx_l, l in enumerate(links)
#     #             if isinstance(l.n1, gs)
#     #             if gss.index(l.n1) == gss.index(d.n1)
#     #            ) - sum(z[idx_l, idx_d, t]
#     #                    for idx_l, l in enumerate(links)
#     #                    if isinstance(l.n2, gs)
#     #                    if gss.index(l.n2) == gss.index(d.n1)
#     #                   ) == 1
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="flow_conservation_1"
#     # )
#     # m.addConstrs(
#     #     (
#     #         sum(z[idx_l, idx_d, t]
#     #             for idx_l, l in enumerate(links)
#     #             if isinstance(l.n2, gs)
#     #             if gss.index(l.n2) == gss.index(d.n2)
#     #            ) - sum(z[idx_l, idx_d, t]
#     #                    for idx_l, l in enumerate(links)
#     #                    if isinstance(l.n1, gs)
#     #                    if gss.index(l.n1) == gss.index(d.n2)
#     #                   ) == 1
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="flow_conservation_2"
#     # )
#     # m.addConstrs(
#     #     (
#     #         sum(z[idx_l, idx_d, t]
#     #             for idx_l, l in enumerate(links)
#     #             if  l.n1 == n
#     #            ) - sum(z[idx_l, idx_d, t]
#     #                    for idx_l, l in enumerate(links)
#     #                    if  l.n2 == n
#     #                   ) == 0
#     #         for idx_d, d in enumerate(demands)
#     #         for n in gss + haps
#     #         if  n != d.n1 and n != d.n2
#     #         for t        in syst.T
#     #     ), name="flow_conservation_3"
#     # )
#     # m.addConstrs(
#     #     (
#     #         z[idx_l_1, idx_d, t] + z[idx_l_2, idx_d, t] <= 1
#     #         for idx_l_1, l_1 in enumerate(links)
#     #         for idx_l_2, l_2 in enumerate(links)
#     #         if  idx_l_1 < idx_l_2 and (l_1.n1 == l_2.n2) and (l_1.n2 == l_2.n1)
#     #         for idx_d, d     in enumerate(demands)
#     #         for t            in syst.T
#     #     ), name="loop_prevention"
#     # )
    
#     # # Maximum Tx/Rx Connection
#     # m.addConstrs(
#     #     (
#     #         sum(z[idx_l, idx_d, t]
#     #             for idx_l, l in enumerate(links)
#     #             if l.n1 == n
#     #            ) <= n.N_TX
#     #         for idx_n, n in enumerate(gss + haps)
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="max_tx_connections"
#     # )
    
#     # m.addConstrs(
#     #     (
#     #         sum(z[idx_l, idx_d, t]
#     #             for idx_l, l in enumerate(links)
#     #             if l.n2 == n
#     #            ) <= n.N_RX
#     #         for idx_n, n in enumerate(gss + haps)
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="max_rx_connections"
#     # )
    
#     # # QKP on HAPs and GSs
#     # m.addConstrs(
#     #     (
#     #         sum(a[idx_l, tp]
#     #             for tp in range(t)
#     #            ) >= syst.THETA * sum(r_2[idx_l, idx_d, t]
#     #                                  for idx_d, d in enumerate(demands)
#     #                                 ) * STORAGE_SCALE
#     #         for idx_l, l in enumerate(links)
#     #         for t        in syst.T
#     #     ), name="qkp_min_capacity"
#     # )
    
#     # m.addConstrs(
#     #     (
#     #         a[idx_l, t] == syst.THETA * (k[idx_l, t] - sum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
#     #                                                        for idx_d, d in enumerate(demands)
#     #                                                       )
#     #                                     ) * STORAGE_SCALE
#     #         for idx_l, l in enumerate(links)
#     #         for t        in syst.T
#     #     ), name="qkp_sequence"
#     # )
    
#     # m.addConstrs(
#     #     (
#     #         sum(a[idx_l, tp]
#     #             for tp in range(t)
#     #            ) <= min(l.n1.A_MAX, l.n2.A_MAX) * STORAGE_SCALE
#     #         for idx_l, l in enumerate(links)
#     #         for t        in syst.T
#     #     ), name="qkp_max_capacity"
#     # )

#     # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
#     for idx_l, l in enumerate(links):
#         # identify which endpoint is HAP and which is GS
#         if isinstance(l.n1, hap) and isinstance(l.n2, gs):
#             hap_idx, gs_node = haps.index(l.n1), l.n2
#         elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
#             hap_idx, gs_node = haps.index(l.n2), l.n1

#         # precompute GS coordinates in same (x,y,z) frame and units (km)
#         #[cg1, cg2] = lonlat_to_xy(gs_node.lg, gs_node.la)  # ensure these are constants in km
#         [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)  # ensure these are constants in km

#         #print(cg1, cg2)
#         # print(xy_to_lonlat(cg1, cg2))
#         #print(tangent_to_latlon(cg1, cg2, 279, 49))

#         if problem == "FIXED" or problem == "INIT":
#             # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
#             dx = c1[hap_idx] - cg1*COORDINATE_SCALE
#             dy = c2[hap_idx] - cg2*COORDINATE_SCALE
#             dz = c3[hap_idx]
#             m.addQConstr(d[idx_l]*d[idx_l] >= dx*dx + dy*dy + dz*dz*COORDINATE_SCALE*COORDINATE_SCALE,
#                          name=f"dist_cone_{idx_l}")
#         elif problem == "FREE":
#             for t in syst.T:
#                 # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
#                 dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
#                 dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
#                 dz = c3[hap_idx, t]
#                 m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + dz*dz*COORDINATE_SCALE*COORDINATE_SCALE,
#                              name=f"dist_cone_{idx_l}_{t}")

#     # for idx_l, l in enumerate(links):
#     #     for t in syst.T:
#     #         m.addConstr(d[idx_l, t] == d[idx_l, 0], name=f"fix_dist_{idx_l}_{t}")

#     ## Solve
#     m.optimize()

#     if m.status == GRB.OPTIMAL:
#         print("\n=========== OPTIMAL SOLUTION FOUND ===========")

#         # solution = {
#         #     "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
#         #     "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
#         #     "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
#         #     "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items()},
#         #     "z":   {k: v.X for k, v in z.items()},
#         #     "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
#         #     "d":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in d.items()},
#         #     "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
#         #     "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()},
#         #     "c3":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c3.items()}
#         # }
        
#         # solution = {
#         #     "r_1": {k: round(v.X, 3) for k, v in r_1.items()},
#         #     "r_2": {k: round(v.X, 3) for k, v in r_2.items()},
#         #     "r_h": {k: round(v.X, 3) for k, v in r_h.items()},
#         #     "a":   {k: round(v.X, 3) for k, v in a.items()},
#         #     "z":   {k: v.X for k, v in z.items()},
#         #     "k":   {k: round(v.X, 3) for k, v in k.items()}
#         # }

#         solution = {
#             "k":  {k: round(v.X, 3) for k, v in k.items()},
#             "d":  {k: round(v.X, 3) for k, v in d.items()},
#             "c1": {k: round(v.X, 3) for k, v in c1.items()},
#             "c2": {k: round(v.X, 3) for k, v in c2.items()},
#             "c3": {k: round(v.X, 3) for k, v in c3.items()}
#         }

#         for hap in haps:
#             hap.
        
#         # if problem == "FREE":
#         #     extra = {
#         #         "d":  {k: round(v.X, 3) for k, v in d.items()},
#         #         "c1": {k: round(v.X, 3) for k, v in c1.items()},
#         #         "c2": {k: round(v.X, 3) for k, v in c2.items()},
#         #         "c3": {k: round(v.X, 3) for k, v in c3.items()}
#         #     }
#         #     solution.update(extra)

#         plot_connectivity_graph(gnodes, hnodes, links)

#         pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
#         pp.pprint(solution)

#         # print(f"lonlat0:{solution["c1"][0,0]}, {solution["c2"][0,0]}, {xy_to_lonlat(solution["c1"][0,0], solution["c2"][0,0])}")

#         # # Actual HAP trajectories
#         # actual_lons = [hnode.lg for hnode in haps]
#         # actual_lats = [hnode.la for hnode in haps]
#         # actual_labels = [f"HAP_{idx_hnode}_actual" for idx_hnode, _ in enumerate(haps)]
        
#         # # Planned HAP trajectories
#         # planned_lons = []
#         # planned_lats = []
#         # for idx_hnode in range(len(haps)):
#         #     lon_series = []
#         #     lat_series = []
#         #     for t in syst.T:
#         #         x = solution["c1"].get((idx_hnode, t))
#         #         y = solution["c2"].get((idx_hnode, t))
#         #         if x is not None and y is not None:
#         #             # lon, lat = xy_to_lonlat(x, y)
#         #             lon, lat = tangent_to_latlon(x, y, 279, 49)
#         #             # lon_series.append(lon+360)   # shift if needed
#         #             lon_series.append(lon)   # shift if needed
#         #             lat_series.append(lat)
#         #     planned_lons.append(lon_series)
#         #     planned_lats.append(lat_series)
#         # planned_labels = [f"HAP_{idx_hnode}_planned" for idx_hnode in range(len(haps))]
        
#         # # Ground Stations (replicate coordinates across all T so they plot in animation)
#         # gs_lons = []
#         # gs_lats = []
#         # for gnode in gss:
#         #     gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
#         #     gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
#         # gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]
        
#         # # Combine everything
#         # all_lons = actual_lons + planned_lons + gs_lons
#         # all_lats = actual_lats + planned_lats + gs_lats
#         # all_labels = actual_labels + planned_labels + gs_labels
        
#         # # Animate
#         # animate_hap_trajectories(
#         #     syst.T,
#         #     all_lons,
#         #     all_lats,
#         #     all_labels
#         # )
#     else:
#         print("No optimal solution found.")
#         solution = None

#     return solution