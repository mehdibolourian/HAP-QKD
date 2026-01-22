from libraries import *

lambda_1 = 100

def placement_greedy(gss, haps, links):
    lon0, lat0 = 279, 49

    solution = {"c1": {}, "c2": {}}

    demands = []
    for idx_g1, g1 in enumerate(gss):
        for idx_g2, g2 in enumerate(gss):
            if idx_g2 > idx_g1:
                demands.append(
                    demand(
                        100,
                        gss[idx_g1],
                        gss[idx_g2]
                    )
                )

    # Reference HAP trajectory
    c1_list_ref, c2_list_ref = [], []
    for lon, lat in zip(haps[0].lg, haps[0].la):
        c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
        c1_list_ref.append(c1)
        c2_list_ref.append(c2)

    # Mean altitude
    H = haps[0].H
    H_mean = sum(H) / len(H) if H else 0.0

    # Compute distances
    dist = {}
    for idx_d, d in enumerate(demands):
        src, dst = d.n1, d.n2

        src_c1, src_c2 = latlon_to_tangent(src.lg, src.la, lon0, lat0)

        dst_c1, dst_c2 = latlon_to_tangent(dst.lg, dst.la, lon0, lat0)

        dist[d] = math.sqrt(
            (src_c1 - dst_c1)**2 +
            (src_c2 - dst_c2)**2
        )

    # Sort links
    ordered_demands = sorted(demands, key=lambda l: dist[l], reverse=True)

    # Place HAPs
    c1_ref_mean = sum(c1_list_ref) / len(c1_list_ref)
    c2_ref_mean = sum(c2_list_ref) / len(c2_list_ref)
    for hap_idx, d in enumerate(ordered_demands[:len(haps)]):
        src, dst = d.n1, d.n2

        src_c1, src_c2 = latlon_to_tangent(src.lg, src.la, lon0, lat0)
        dst_c1, dst_c2 = latlon_to_tangent(dst.lg, dst.la, lon0, lat0)
    
        # Desired midpoint
        mid_c1 = (src_c1 + dst_c1) / 2
        mid_c2 = (src_c2 + dst_c2) / 2
    
        # Shift so that MEAN of trajectory equals midpoint
        shift_c1 = mid_c1 - c1_ref_mean
        shift_c2 = mid_c2 - c2_ref_mean
    
        for t in syst.T:
            solution["c1"][(hap_idx, t)] = c1_list_ref[t] + shift_c1
            solution["c2"][(hap_idx, t)] = c2_list_ref[t] + shift_c2
            haps[hap_idx].lg[t], haps[hap_idx].la[t] = tangent_to_latlon(solution["c1"][(hap_idx, t)], solution["c2"][(hap_idx, t)], lon0, lat0)

    dist_bottle = 0
    dist_sum    = 0
    for g in gss:
        dist_min = math.inf
        for h in haps:
            gs_c1,  gs_c2  = latlon_to_tangent(g.lg, g.la, lon0, lat0)
            hap_c1, hap_c2 = latlon_to_tangent(sum(h.lg)/len(h.lg), sum(h.la)/len(h.la), lon0, lat0)
            
            dist_min = min(dist_min, math.sqrt((gs_c1 - hap_c1)**2 + (gs_c2 - hap_c2)**2))

            print(f"dist_min: {dist_min}")

        dist_sum    = dist_sum + dist_min
        dist_bottle = max(dist_bottle, dist_min)

    print(f"dist_bottle: {dist_bottle}")
    print(f"dist_avg: {dist_sum/len(gss)}")

    # Planned trajectories
    planned_lons, planned_lats = [], []
    for idx in range(len(haps)):
        lon_series, lat_series = [], []
        for t in syst.T:
            x = solution["c1"].get((idx, t))
            y = solution["c2"].get((idx, t))
            if x is not None and y is not None:
                lon, lat = tangent_to_latlon(x, y, lon0, lat0)
                lon_series.append(lon)
                lat_series.append(lat)
        planned_lons.append(lon_series)
        planned_lats.append(lat_series)

    planned_labels = [f"HAP_{i}*" for i in range(len(haps))]

    plot_connectivity_graph_planning(
                        gss, haps, links,
                        planned_lons=planned_lons,
                        planned_lats=planned_lats,
                        planned_labels=planned_labels,
                        alg="grd"
                       )

    plot_connectivity_graph_planning_3d(
                        gss, haps, links,
                        planned_lons=planned_lons,
                        planned_lats=planned_lats,
                        planned_alts=haps[0].H,
                        planned_labels=planned_labels,
                        alg="grd"
                       )

    # plot_connectivity_graph_planning_3d(
    #     gss, haps, links,
    #     planned_lons=planned_lons,
    #     planned_lats=planned_lats,
    #     planned_labels=planned_labels
    # )

    return solution, calculate_key_rate_planning("theoretical", 0, dist_bottle, 0), calculate_key_rate_planning("theoretical", 0, dist_sum/len(gss), 0)

def placement_kmeans(gss, haps, links):
    lon0, lat0 = 279, 49

    solution = {"c1": {}, "c2": {}}

    #print(f"len(haps): {len(haps)}")

    demands     = []
    gss_demands = []
    for idx_g1, g1 in enumerate(gss):
        for idx_g2, g2 in enumerate(gss):
            if idx_g2 > idx_g1:
                demands.append(
                    demand(
                        100,
                        gss[idx_g1],
                        gss[idx_g2]
                    )
                )
                if gss[idx_g1] not in gss_demands:
                    gss_demands.append(gss[idx_g1])
                if gss[idx_g2] not in gss_demands:
                    gss_demands.append(gss[idx_g2])

    # Reference HAP trajectory
    c1_list_ref, c2_list_ref = [], []
    for lon, lat in zip(haps[0].lg, haps[0].la):
        c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
        c1_list_ref.append(c1)
        c2_list_ref.append(c2)

    c1_ref_mean = sum(c1_list_ref) / len(c1_list_ref)
    c2_ref_mean = sum(c2_list_ref) / len(c2_list_ref)

    # Mean altitude
    H = haps[0].H
    H_mean = sum(H) / len(H) if H else 0.0

    # GSs coordinates
    gs_coordinates = []
    for g in gss_demands:
        c1, c2 = latlon_to_tangent(g.lg, g.la, lon0, lat0)
        gs_coordinates.append([c1, c2])

    if len(haps) == 1:
        k = 2
    elif len(haps) <= 3:
        k = 3
    elif len(haps) <= 6:
        k = 4
    elif len(haps) <= 10:
        k = 5
    else:
        k = 6
        
    kmeans = KMeans(n_clusters=k, random_state=0, n_init=10)
    
    labels = kmeans.fit_predict(gs_coordinates)
    centers = kmeans.cluster_centers_
    
    #print(labels)
    #print(centers)

    # --- ordered center pairs (farthest first) ---
    center_pairs = []
    
    for (i, c1), (j, c2) in combinations(enumerate(centers), 2):
        dist = np.linalg.norm(c1 - c2)
        center_pairs.append((dist, i, j))
    
    # sort by distance descending
    center_pairs.sort(reverse=True, key=lambda x: x[0])

    num_cluster_conns = k * (k-1) / 2
    hap_idx = 0
    for dist, idx_gc1, idx_gc2 in center_pairs:
        mid_c1 = (centers[idx_gc1][0] + centers[idx_gc2][0])/2
        mid_c2 = (centers[idx_gc1][1] + centers[idx_gc2][1])/2
        
        # Shift so that MEAN of trajectory equals midpoint
        shift_c1 = mid_c1 - c1_ref_mean
        shift_c2 = mid_c2 - c2_ref_mean
    
        for t in syst.T:
            solution["c1"][(hap_idx, t)] = c1_list_ref[t] + shift_c1
            solution["c2"][(hap_idx, t)] = c2_list_ref[t] + shift_c2
            haps[hap_idx].lg[t], haps[hap_idx].la[t] = tangent_to_latlon(solution["c1"][(hap_idx, t)], solution["c2"][(hap_idx, t)], lon0, lat0)

        hap_idx += 1
        if hap_idx == len(haps):
            break

    dist_bottle = 0
    dist_sum    = 0
    for g in gss:
        dist_min = math.inf
        for h in haps:
            gs_c1,  gs_c2  = latlon_to_tangent(g.lg, g.la, lon0, lat0)
            hap_c1, hap_c2 = latlon_to_tangent(sum(h.lg)/len(h.lg), sum(h.la)/len(h.la), lon0, lat0)
            
            dist_min = min(dist_min, math.sqrt((gs_c1 - hap_c1)**2 + (gs_c2 - hap_c2)**2))

            print(f"dist_min: {dist_min}")

        dist_sum    = dist_sum + dist_min
        dist_bottle = max(dist_bottle, dist_min)

    print(f"dist_bottle: {dist_bottle}")
    print(f"dist_avg: {dist_sum/len(gss)}")

    # Planned trajectories
    planned_lons, planned_lats = [], []
    for idx in range(len(haps)):
        lon_series, lat_series = [], []
        for t in syst.T:
            x = solution["c1"].get((idx, t))
            y = solution["c2"].get((idx, t))
            if x is not None and y is not None:
                lon, lat = tangent_to_latlon(x, y, lon0, lat0)
                lon_series.append(lon)
                lat_series.append(lat)
        planned_lons.append(lon_series)
        planned_lats.append(lat_series)

    planned_labels = [f"HAP_{i}*" for i in range(len(haps))]

    plot_connectivity_graph_planning(
                        gss, haps, links,
                        planned_lons=planned_lons,
                        planned_lats=planned_lats,
                        planned_labels=planned_labels,
                        alg="kmn"
                       )
    plot_connectivity_graph_planning_3d(
                        gss, haps, links,
                        planned_lons=planned_lons,
                        planned_lats=planned_lats,
                        planned_alts=haps[0].H,
                        planned_labels=planned_labels,
                        alg="kmn"
                       )

    # plot_connectivity_graph_planning_3d(
    #     gss, haps, links,
    #     planned_lons=planned_lons,
    #     planned_lats=planned_lats,
    #     planned_labels=planned_labels
    # )

    return solution, calculate_key_rate_planning("theoretical", 0, dist_bottle, 0), calculate_key_rate_planning("theoretical", 0, dist_sum/len(gss), 0)






    
## Find the optimal placement of the HAPs to reach the maximum key generation for all the end-to-end paths between GS pairs. 
def placement(gss, haps, links):
    c1_list_ref, c2_list_ref, c3_list_ref = [], [], []
    lon0, lat0 = 279, 49

    demands = []
    for idx_g1, g1 in enumerate(gss):
        for idx_g2, g2 in enumerate(gss):
            if idx_g2 > idx_g1:
                demands.append(
                    demand(
                        100,
                        gss[idx_g1],
                        gss[idx_g2]
                    )
                )

    # Create Optimization Model
    m = gp.Model("hap-qkd")

    for lon, lat in zip(haps[0].lg, haps[0].la):
        c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
        c1_list_ref.append(c1)
        c2_list_ref.append(c2)
    c3_list_ref = haps[0].H

    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    z, o = {}, {}

    for idx_l, l in enumerate(links):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                z[idx_l, idx_d, t] = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}", vtype=GRB.BINARY)

    nodes = gss + haps

    dpts = np.linspace(15 * COORDINATE_SCALE, 3e3 * COORDINATE_SCALE, 4)
    kpts = [calculate_key_rate_planning("theoretical", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

    # k (key rate) and d (distance)
    k, di, ke = {}, {}, {}
    for idx_l, l in enumerate(links):
        for t in syst.T:
            di[idx_l, t] = m.addVar(name=f"di_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=3e3 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                
            k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

            #print(f"kpts: {kpts}")
            m.addGenConstrPWL(di[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

    plt.plot(dpts, kpts)
    plt.show()

    for idx_d, d in enumerate(demands):
        for t in syst.T:
            ke[idx_d, t] = m.addVar(name=f"ke_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

    # for t in syst.T:
    #     ke[t] = m.addVar(name=f"ke_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

    # Coordinate decision variables for each HAP and time
    c1, c2 = {}, {}  # x, y in km
    for idx_h, hnode in enumerate(haps):
        for t in syst.T:
            c1[idx_h, t] = m.addVar(lb=-1.5e3*COORDINATE_SCALE, ub=1.5e3*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
            c2[idx_h, t] = m.addVar(lb=-1.5e3*COORDINATE_SCALE, ub=1.5e3*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")

    # Secondary objective: maximize d
    # m.setObjective(gp.quicksum(ke[t]
    #                            for t in syst.T
    #                           )
    #                , GRB.MAXIMIZE)

    m.setObjective(gp.quicksum(gp.quicksum(ke[idx_d, t]
                                           for t in syst.T
                                          )
                               for idx_d, d in enumerate(demands)
                              )
                   , GRB.MAXIMIZE)

    # m.setObjective(gp.quicksum(gp.quicksum(ke[idx_d, t]
    #                                        for t in syst.T
    #                                       )
    #                            for idx_d, d in enumerate(demands)
    #                           ) - lambda_1 * gp.quicksum(gp.quicksum((gp.quicksum(gp.quicksum(z[idx_l, idx_d, t]
    #                                                                                           for idx_l, l in enumerate(links)
    #                                                                                           if l.n1 == h or l.n2 == h
    #                                                                                          )
    #                                                                               for idx_d, d in enumerate(demands)
    #                                                                              ))**2
    #                                                                  for idx_h, h in enumerate(haps)
    #                                                                 )
    #                                                      for t in syst.T
    #                                                     )
    #                , GRB.MAXIMIZE)

     # - gp.quicksum(gp.quicksum(gp.quicksum(z[idx_l, idx_d, t]
     #                                                                  for idx_l, l in enumerate(links)
     #                                                                 )
     #                                                      for idx_d, d in enumerate(demands)
     #                                                     )
     #                                          for t in syst.T
     #                                         )
    

    m.setParam("MIPGap", 1e-4)
    m.setParam("MIPGapAbs", 1e-4)
    m.setParam("FeasibilityTol", 1e-4)
    m.setParam("IntFeasTol", 1e-4)
    m.setParam("OptimalityTol", 1e-4)
    m.setParam("TimeLimit", 60)  # seconds
    m.Params.Presolve = 2
    m.Params.Method = 2
    m.Params.Cuts = 2
    m.Params.Heuristics = 1
    m.Params.MIPFocus = 1
    m.Params.NumericFocus = 1
    m.Params.Threads = 0
    m.Params.NodefileStart = 0.5
    m.Params.NoRelHeurTime = 60
    m.Params.ConcurrentMIP = 1

    # One-hop paths
    m.addConstrs(
        (gp.quicksum(z[idx_l, idx_d, t]
                     for idx_l, l in enumerate(links)
                    )
         == 2
         for idx_d, d in enumerate(demands)
         for t in syst.T),
        name="one_hop"
    )

    # Flow conservation
    m.addConstrs(
        (gp.quicksum(z[idx_l, idx_d, t]
                     for idx_l, l in enumerate(links)
                     if l.n1 == d.n1)
         - gp.quicksum(z[idx_l, idx_d, t]
                       for idx_l, l in enumerate(links)
                       if l.n2 == d.n1)
         == 1
         for idx_d, d in enumerate(demands)
         for t in syst.T),
        name="flow_conservation_1"
    )

    m.addConstrs(
        (gp.quicksum(z[idx_l, idx_d, t]
                     for idx_l, l in enumerate(links)
                     if l.n2 == d.n2)
         - gp.quicksum(z[idx_l, idx_d, t]
                       for idx_l, l in enumerate(links)
                       if l.n1 == d.n2)
         == 1
         for idx_d, d in enumerate(demands)
         for t in syst.T),
        name="flow_conservation_2"
    )

    m.addConstrs(
        (gp.quicksum(z[idx_l, idx_d, t]
                     for idx_l, l in enumerate(links)
                     if l.n1 == n
                    )
         - gp.quicksum(z[idx_l, idx_d, t]
                       for idx_l, l in enumerate(links)
                       if l.n2 == n
                      )
         == 0
         for idx_d, d in enumerate(demands)
         for n in gss + haps
         if n != d.n1 and n != d.n2
         for t in syst.T),
        name="flow_conservation_3"
    )

    ########### Exclusive constraints ###########

    m.addConstrs(
        (
            ke[idx_d, t] <= k[idx_l, t] + kpts[0] * (1 - z[idx_l, idx_d, t])
            for idx_d, d in enumerate(demands)
            for idx_l, l in enumerate(links)
            for t in syst.T
        ), name="keyrate_active_link_3"
    )
    # m.addConstrs(
    #     (
    #         ke[t] <= k[idx_l, t] + kpts[0] * (1 - z[idx_l, idx_d, t])
    #         for idx_d, d in enumerate(demands)
    #         for idx_l, l in enumerate(links)
    #         for t in syst.T
    #     ), name="keyrate_active_link_3"
    # )

    m.addConstrs(
        (c1[idx_h, t] == c1_list_ref[t] + (c1[idx_h, 0] - c1_list_ref[0])
         for idx_h, h in enumerate(haps)
         for t in syst.T if t >= 1),
        name="shift_trajectory_1"
    )

    m.addConstrs(
        (c2[idx_h, t] == c2_list_ref[t] + (c2[idx_h, 0] - c2_list_ref[0])
         for idx_h, h in enumerate(haps)
         for t in syst.T if t >= 1),
        name="shift_trajectory_2"
    )

    # m.addConstrs(
    #     (c1[idx_h, t] > c1[idx_h, 0]
    #      for idx_h, h in enumerate(haps)
    #      for t in syst.T if t >= 1),
    #     name="placement_1"
    # )

    # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
    for idx_l, l in enumerate(links):
        # identify which endpoint is HAP and which is GS
        if isinstance(l.n1, hap) and isinstance(l.n2, gs):
            hap_idx, gs_node = haps.index(l.n1), l.n2
        elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
            hap_idx, gs_node = haps.index(l.n2), l.n1

        [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)

        for t in syst.T:
            dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
            dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
            m.addQConstr(di[idx_l, t]*di[idx_l, t] >= dx*dx + dy*dy + haps[hap_idx].H[t]*haps[hap_idx].H[t]*COORDINATE_SCALE*COORDINATE_SCALE,
                         name=f"dist_cone_{idx_l}_{t}")
                         
    ## Solve
    m.optimize()

    if m.status in (GRB.OPTIMAL, GRB.TIME_LIMIT) and m.SolCount > 0:
        print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        solution = {
            "z":   {k: v.X for k, v in z.items()},
            #"o":   {k: round(v.X, 3) for k, v in o.items()},
            "ke":  {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in ke.items()},
            "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
            "di":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in di.items()},
            "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
            "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()}
        }

        key_rate_m = min(solution["k"][idx_l, t]
                         for idx_l, l in enumerate(links)
                         for t in syst.T
                        )
        #key_rate_s = sum(solution["k"])
        
        print(f"Minimum link key rate: {key_rate_m}") #, Sum: {key_rate_s}")

        # for idx_l, l in enumerate(links):
        #     for idx_d, d in enumerate(demands):
        #         for t in syst.T:
        #             Z = solution["z"][idx_l, idx_d, t]
        #             if Z >= 0.9:
        #                 print(f"SELECTED: link.n1: {l.n1.tag}, link.n2: {l.n2.tag}, demand: {idx_d}, t: {t}")
                    
        # plot_z_timeline(solution["z"], links, demands, figsize=(10,6))
        
        # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
        # pp.pprint(solution)

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
                x = solution["c1"].get((idx_hnode,t))
                y = solution["c2"].get((idx_hnode,t))
                haps[idx_hnode].lg[t], haps[idx_hnode].la[t] = tangent_to_latlon(x, y, lon0, lat0)
                if x is not None and y is not None:
                    # lon, lat = xy_to_lonlat(x, y)
                    lon, lat = tangent_to_latlon(x, y, 279, 49)
                    lon_series.append(lon)   # shift if needed
                    lat_series.append(lat)
            planned_lons.append(lon_series)
            planned_lats.append(lat_series)
        planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]
        
        # Ground Stations (replicate coordinates across all T so they plot in animation)
        gs_lons = []
        gs_lats = []
        for gnode in gss:
            gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
            gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
        gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]

        all_lons = planned_lons + gs_lons
        all_lats = planned_lats + gs_lats
        all_labels = planned_labels + gs_labels

        # plot_connectivity_graph_planning(gss, haps, links, 
        #                 planned_lons=planned_lons,
        #                 planned_lats=planned_lats,
        #                 planned_labels=planned_labels,
        #                 alg="opt"
        #                                 )

        # plot_connectivity_graph_planning_3d(gss, haps, links, 
        #                 planned_lons=planned_lons,
        #                 planned_lats=planned_lats,
        #                 planned_labels=planned_labels)

        # return solution
    else:
        solution = None

    dist_bottle = 0
    dist_sum    = 0
    for g in gss:
        dist_min = math.inf
        for h in haps:
            gs_c1,  gs_c2  = latlon_to_tangent(g.lg, g.la, lon0, lat0)
            hap_c1, hap_c2 = latlon_to_tangent(sum(h.lg)/len(h.lg), sum(h.la)/len(h.la), lon0, lat0)
            
            dist_min = min(dist_min, math.sqrt((gs_c1 - hap_c1)**2 + (gs_c2 - hap_c2)**2))

            print(f"dist_min: {dist_min}")

        dist_sum    = dist_sum + dist_min
        dist_bottle = max(dist_bottle, dist_min)

    print(f"dist_bottle: {dist_bottle}")
    print(f"dist_avg: {dist_sum/len(gss)}")

    # Planned HAP trajectories
    planned_lons = []
    planned_lats = []
    for idx_hnode in range(len(haps)):
        lon_series = []
        lat_series = []
        for t in syst.T:
            x = solution["c1"].get((idx_hnode,t))
            y = solution["c2"].get((idx_hnode,t))
            haps[idx_hnode].lg[t], haps[idx_hnode].la[t] = tangent_to_latlon(x, y, lon0, lat0)
            if x is not None and y is not None:
                # lon, lat = xy_to_lonlat(x, y)
                lon, lat = tangent_to_latlon(x, y, 279, 49)
                lon_series.append(lon)   # shift if needed
                lat_series.append(lat)
        planned_lons.append(lon_series)
        planned_lats.append(lat_series)
    planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]

    plot_connectivity_graph_planning(gss, haps, links, 
                    planned_lons=planned_lons,
                    planned_lats=planned_lats,
                    planned_labels=planned_labels,
                    alg="opt"
                   )
    plot_connectivity_graph_planning_3d(gss, haps, links, 
                    planned_lons=planned_lons,
                    planned_lats=planned_lats,
                    planned_alts=haps[0].H,
                    planned_labels=planned_labels,
                    alg="opt"
                   )
    
    return solution, calculate_key_rate_planning("theoretical", 0, dist_bottle, 0), calculate_key_rate_planning("theoretical", 0, dist_sum/len(gss), 0)

def plot_connectivity_graph_planning(gnodes, hnodes, links, 
                                    planned_lons=None, planned_lats=None, planned_labels=None, alg=""):
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
    plt.figure(figsize=(5, 4))
    
    # --- Plot GS nodes ---
    for gs_node in gnodes:
        plt.scatter(gs_node.lg, gs_node.la, color='skyblue', s=80, zorder=5, marker='^')
        # Optional: label the GS
        if hasattr(gs_node, 'tag'):
            plt.text(gs_node.lg + 0.04, gs_node.la + 0.04, gs_node.tag, fontsize=9)
    
    # # --- Plot HAP nodes (initial position) ---
    # for hap_node in hnodes:
    #     plt.scatter(hap_node.lg[0], hap_node.la[0], color='orange', s=5, zorder=5)
    #     if hasattr(hap_node, 'tag'):
    #         plt.text(hap_node.lg[0] - 0.4, hap_node.la[0] - 0.2, hap_node.tag, fontsize=9)
    
    # # --- Plot edges without duplicates ---
    # plotted_edges = set()
    # for l in links:
    #     # Use frozenset to make the edge unordered (A-B same as B-A)
    #     edge_key = frozenset([l.n1, l.n2])
    #     if edge_key in plotted_edges:
    #         continue  # already plotted
    #     plotted_edges.add(edge_key)
    
    #     # Determine coordinates for nodes
    #     x = [l.n1.lg[0] if isinstance(l.n1.lg, list) else l.n1.lg,
    #          l.n2.lg[0] if isinstance(l.n2.lg, list) else l.n2.lg]
    #     y = [l.n1.la[0] if isinstance(l.n1.la, list) else l.n1.la,
    #          l.n2.la[0] if isinstance(l.n2.la, list) else l.n2.la]
        
    #     # Decide line style
    #     plt.plot(x, y, color='grey', linestyle='--', alpha=0.6, linewidth=0.5)
    
    # # --- Plot HAP trajectories ---
    # for hap_node in hnodes:
    #     plt.plot(hap_node.lg, hap_node.la, color='orange', linewidth=2, alpha=0.8)
    
    
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
        Line2D([], [], marker='o', color='red', linestyle='None', markersize=6, label='HAP')
    ]
    plt.legend(handles=custom_handles, loc='upper right', frameon=True, fontsize=11)
    
    plt.grid(True, alpha=0.3)

    # --- Plot optimal trajectories (if provided) ---
    if planned_lons and planned_lats:
        for idx, (plon, plat) in enumerate(zip(planned_lons, planned_lats)):
            label = (planned_labels[idx] if planned_labels and idx < len(planned_labels)
                     else f"Optimal_{idx}")
            
            if len(plon) == 0 or len(plat) == 0:
                continue

            # Plot connecting red line (trajectory)
            plt.plot(plon, plat, color='red', linewidth=1.1, alpha=1, label=None)

            # Mark initial point
            plt.scatter(plon[0], plat[0], color='red', s=5, marker='o', zorder=6)

            # Add text label near the last planned point with coordinates
            lon_last, lat_last = plon[-1], plat[-1]
            lon_first, lat_first = plon[0], plat[0]
            coord_text = f"{idx}:({lon_first:.2f}, {lat_first:.2f})"
            print(coord_text)
            plt.text(lon_last - 0.4, lat_last - 0.2, coord_text,
                     fontsize=9, fontstyle='italic', color='red')

    # ==============================
    # Zoomed-in inset for one HAP
    # ==============================
    hap_zoom = hnodes[0]   # choose which HAP to zoom

    # Create inset axis
    ax = plt.gca()
    axins = inset_axes(
        ax,
        width="25%",   # relative size
        height="25%",
        loc="lower right",
        borderpad=1.2
    )

    # Plot trajectory inside inset
    axins.plot(hap_zoom.lg, hap_zoom.la,
               color='red', linewidth=1.1)

    # Optional: mark start point
    axins.scatter(hap_zoom.lg[0], hap_zoom.la[0],
                  color='red', s=2, zorder=5)

    # Set zoom window (tight bounds)
    margin = 0.03
    axins.set_xlim(min(hap_zoom.lg) - margin, max(hap_zoom.lg) + margin)
    axins.set_ylim(min(hap_zoom.la) - margin, max(hap_zoom.la) + margin)

    # Clean inset appearance
    axins.set_xticks([])
    axins.set_yticks([])
    axins.grid(True, alpha=0.3)

    # Draw rectangle on main plot to show zoomed region
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
    plt.savefig(f"hap_qkd_network_{len(hnodes)}_{alg}.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()






















    
## For a given set of demands what is the minimum number of HAPs and where to place them to satisfy all the demands.
def demand_feasibility(gss, demands):
    MAX_HAPS = 10
    c1_list_ref, c2_list_ref, c3_list_ref = [], [], []
    lon0, lat0 = 279, 49

    for num_h in range(1, MAX_HAPS):
        haps = []
        links = []

        # Create Optimization Model
        m = gp.Model("hap-qkd")

        for n in range(num_h):
            haps.append(hap([279]*len(syst.T), [49]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, f"HAP_{n}"))

        # Update coordinates depending on model choice
        if SYNTH_STRATO == 1:
            update_coordinates("stratotegic", haps, syst)
        else:
            update_coordinates("wind", haps, syst)

        ## Only once is enough
        if num_h == 1:
            for lon, lat in zip(haps[0].lg, haps[0].la):
                c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
                c1_list_ref.append(c1)
                c2_list_ref.append(c2)
            c3_list_ref = haps[0].H

        # Links: connect only GSs to HAPs
        for gs_node in gss:
            for hap_node in haps:
                links.append(link(gs_node, hap_node, [100]*len(syst.T), [(0,0,0)]*len(syst.T), [1]*len(syst.T)))
                links.append(link(hap_node, gs_node, [100]*len(syst.T), [(0,0,0)]*len(syst.T), [1]*len(syst.T)))

        ## Decision Variables
        # Dictionaries of decision variables instead of MVar arrays
        r_1, r_2, r_h, a, z = {}, {}, {}, {}, {}

        for idx_l, l in enumerate(links):
            for idx_d, d in enumerate(demands):
                for t in syst.T:
                    r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
                    r_2[idx_l, idx_d, t] = m.addVar(name=f"r_2_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
                    z[idx_l, idx_d, t] = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}", vtype=GRB.BINARY)

        for idx_d, d in enumerate(demands):
            for t in syst.T:
                r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)

        for idx_l, l in enumerate(links):
            for t in syst.T:
                a[idx_l, t] = m.addVar(name=f"a_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY)

        nodes = gss + haps

        # k (key rate) and d (distance)
        k, d, kz = {}, {}, {}
        for idx_l, l in enumerate(links):
            for t in syst.T:
                
                d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                    
                k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
                kz[idx_l, t] = m.addVar(name=f"kz_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

                dpts = np.linspace(15 * COORDINATE_SCALE, 2e2 * COORDINATE_SCALE, 4)
                kpts = [calculate_key_rate_planning("theoretical", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

                #print(f"kpts: {kpts}")
                m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

        # Coordinate decision variables for each HAP and time
        c1, c2 = {}, {}  # x, y in km
        for idx_h, hnode in enumerate(haps):
            for t in syst.T:
                c1[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
                c2[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")

        # Secondary objective: maximize d
        m.setObjective(1, GRB.MAXIMIZE)

        m.setParam("MIPGap", 1e-4)
        m.setParam("MIPGapAbs", 1e-4)
        m.setParam("FeasibilityTol", 1e-4)
        m.setParam("IntFeasTol", 1e-4)
        m.setParam("OptimalityTol", 1e-4)
        m.Params.Presolve = 2
        m.Params.Method = 2
        m.Params.Cuts = 2
        m.Params.Heuristics = 0.25
        m.Params.MIPFocus = 1
        m.Params.NumericFocus = 1
        m.Params.Threads = 0
        m.Params.NodefileStart = 0.5
        m.Params.NoRelHeurTime = 20
        m.Params.ConcurrentMIP = 1

        ## Constraints     
        # Maximum Link Capacity --> Enforces only r_1
        m.addConstrs(
            (
                gp.quicksum(r_1[idx_l, idx_d, t]
                            for idx_d, d in enumerate(demands)
                           ) <= k[idx_l, t]
             for idx_l, l in enumerate(links)
             for t in syst.T
            ), name="max_link_capacity"
        )

        # Maximum Tx/Rx Connection
        m.addConstrs(
            (gp.quicksum(z[idx_l, idx_d, t]
                         for idx_l, l in enumerate(links)
                         for idx_d, d in enumerate(demands)
                         if l.n1 == n
                        ) <= n.N_TX
             for idx_n, n in enumerate(nodes)
             for t in syst.T),
            name="max_tx_connections"
        )

        m.addConstrs(
            (gp.quicksum(z[idx_l, idx_d, t]
                         for idx_l, l in enumerate(links)
                         for idx_d, d in enumerate(demands)
                         if l.n2 == n
                        ) <= n.N_RX
             for idx_n, n in enumerate(nodes)
             for t in syst.T),
            name="max_rx_connections"
        )

        # Flow conservation
        m.addConstrs(
            (gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                         for idx_l, l in enumerate(links)
                         if l.n1 == d.n1)
             - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                           for idx_l, l in enumerate(links)
                           if l.n2 == d.n1)
             == r_h[idx_d, t]
             for idx_d, d in enumerate(demands)
             for t in syst.T),
            name="flow_conservation_1"
        )

        m.addConstrs(
            (gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                         for idx_l, l in enumerate(links)
                         if l.n2 == d.n2)
             - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                           for idx_l, l in enumerate(links)
                           if l.n1 == d.n2)
             == r_h[idx_d, t]
             for idx_d, d in enumerate(demands)
             for t in syst.T),
            name="flow_conservation_2"
        )

        m.addConstrs(
            (gp.quicksum((r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t])
                         for idx_l, l in enumerate(links)
                         if l.n1 == n
                        )
             - gp.quicksum((r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t])
                           for idx_l, l in enumerate(links)
                           if l.n2 == n
                          )
             == 0
             for idx_d, d in enumerate(demands)
             for n in gss + haps
             if n != d.n1 and n != d.n2
             for t in syst.T),
            name="flow_conservation_3"
        )
        
        # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
        m.addConstrs(
            (r_h[idx_d, t] >= r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
             for idx_l, l in enumerate(links)
             for idx_d, d in enumerate(demands)
             for t in syst.T),
            name="demand_link_coordination_1"
        )

        # # Key rate and routing coordination (1)
        # m.addConstrs(
        #     (r_1[idx_l, idx_d, t] >= 1e-2 * z[idx_l, idx_d, t]
        #      for idx_l, l in enumerate(links)
        #      for idx_d, d in enumerate(demands)
        #      for t in syst.T),
        #     name="key_rate_routing_coordination_1"
        # )

        # Key rate and routing coordination (2)
        m.addConstrs(
            (r_1[idx_l, idx_d, t] <= d.K_REQ[t] * z[idx_l, idx_d, t]
             for idx_l, l in enumerate(links)
             for idx_d, d in enumerate(demands)
             for t in syst.T),
            name="key_rate_routing_coordination_2"
        )

        # QKP on HAPs and GSs
        m.addConstrs(
            (r_2[idx_l, idx_d, 0] == 0
             for idx_l, l in enumerate(links)
             for idx_d, d in enumerate(demands)),
            name="initial_empty_QKP"
        )
        
        m.addConstrs(
            (gp.quicksum(a[idx_l, tp] for tp in range(t))
             >= syst.THETA * gp.quicksum(r_2[idx_l, idx_d, t] for idx_d, d in enumerate(demands)) * STORAGE_SCALE
             for idx_l, l in enumerate(links)
             for t in syst.T[1:]),
            name="qkp_min_capacity"
        )

        m.addConstrs(
            (a[idx_l, t] == syst.THETA * (kz[idx_l, t]
                                          - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                                                        for idx_d, d in enumerate(demands))) * STORAGE_SCALE
             for idx_l, l in enumerate(links)
             for t in syst.T),
            name="qkp_sequence"
        )

        m.addConstrs(
            (
                gp.quicksum(a[idx_l, tp]
                            for tp in range(t+1)
                           ) >= 0
                for idx_l, l in enumerate(links)
                for t        in syst.T
            ), name="positive_storage"
        )

        ########### Exclusive constraints ###########

        # Demand satisfaction guarantee
        m.addConstrs(
            (gp.quicksum(r_h[idx_d, t] for t in syst.T) == sum(d.K_REQ[t] for t in syst.T) * KEY_RATE_SCALE
             for idx_d, d in enumerate(demands)),
            name="demand_satisfaction_guarantee"
        )

        # Active link check
        m.addConstrs(
            (
                kz[idx_l, t] <= k[idx_l, t]
                for idx_l, l in enumerate(links)
                for t in syst.T
            ), name="keyrate_active_link_1"
        )

        m.addConstrs(
            (
                kz[idx_l, t] <= 1e7 * gp.quicksum(z[idx_l, idx_d, t]
                                                   for idx_d, d in enumerate(demands)
                                                  )
                for idx_l, l in enumerate(links)
                for t in syst.T
            ), name="keyrate_active_link_2"
        )

        m.addConstrs(
            (
                kz[idx_l, t] >= k[idx_l, t] - 1e7 * (1 - gp.quicksum(z[idx_l, idx_d, t]
                                                                      for idx_d, d in enumerate(demands)
                                                                     )
                                                     )
                for idx_l, l in enumerate(links)
                for t in syst.T
            ), name="keyrate_active_link_3"
        )

        m.addConstrs(
            (c1[idx_h, t] == c1_list_ref[t] + (c1[idx_h, 0] - c1_list_ref[0])
             for idx_h, h in enumerate(haps)
             for t in syst.T if t >= 1),
            name="shift_trajectory_1"
        )

        m.addConstrs(
            (c2[idx_h, t] == c2_list_ref[t] + (c2[idx_h, 0] - c2_list_ref[0])
             for idx_h, h in enumerate(haps)
             for t in syst.T if t >= 1),
            name="shift_trajectory_2"
        )

        # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
        for idx_l, l in enumerate(links):
            # identify which endpoint is HAP and which is GS
            if isinstance(l.n1, hap) and isinstance(l.n2, gs):
                hap_idx, gs_node = haps.index(l.n1), l.n2
            elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
                hap_idx, gs_node = haps.index(l.n2), l.n1

            [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)

            for t in syst.T:
                dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
                dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
                m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + haps[hap_idx].H[t]*haps[hap_idx].H[t]*COORDINATE_SCALE*COORDINATE_SCALE,
                             name=f"dist_cone_{idx_l}_{t}")
                             
        ## Solve
        m.optimize()

        if m.status == GRB.OPTIMAL:
            print("\n=========== OPTIMAL SOLUTION FOUND ===========")
    
            solution = {
                "r_h":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
                "r_1":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
                "r_2":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
                "a":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in a.items()},
                "z":   {k: v.X for k, v in z.items()},
                "kz":  {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in kz.items()},
                "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
                "d":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in d.items()},
                "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
                "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()}
            }

            key_rate_m = min(solution["k"][idx_l, t]
                             for idx_l, l in enumerate(links)
                             for t in syst.T
                            )
            #key_rate_s = sum(solution["k"])
            
            print(f"Minimum link key rate: {key_rate_m}") #, Sum: {key_rate_s}")

            for idx_l, l in enumerate(links):
                for idx_d, d in enumerate(demands):
                    Z = sum(solution["z"][idx_l, idx_d, t]
                            for t in syst.T
                           )
                    if Z > 0:
                        print(f"SELECTED: link.n1: {l.n1.tag}, link.n2: {l.n2.tag}")
                        
            plot_z_timeline(solution["z"], links, demands, figsize=(10,6))
            
            # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
            # pp.pprint(solution)
    
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
                    x = solution["c1"].get((idx_hnode,t))
                    y = solution["c2"].get((idx_hnode,t))
                    if x is not None and y is not None:
                        # lon, lat = xy_to_lonlat(x, y)
                        lon, lat = tangent_to_latlon(x, y, 279, 49)
                        lon_series.append(lon)   # shift if needed
                        lat_series.append(lat)
                planned_lons.append(lon_series)
                planned_lats.append(lat_series)
            planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]
            
            # Ground Stations (replicate coordinates across all T so they plot in animation)
            gs_lons = []
            gs_lats = []
            for gnode in gss:
                gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
                gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
            gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]
    
            all_lons = planned_lons + gs_lons
            all_lats = planned_lats + gs_lats
            all_labels = planned_labels + gs_labels
    
            plot_connectivity_graph_planning(gss, haps, links, 
                            planned_lons=planned_lons,
                            planned_lats=planned_lats,
                            planned_labels=planned_labels)

            plot_connectivity_graph_planning_3d(gss, haps, links, 
                            planned_lons=planned_lons,
                            planned_lats=planned_lats,
                            planned_labels=planned_labels)

            print(f"The minimum number of required HAPs: {num_h}")

            # return solution
        else:
            print(f"No optimal solution found for {num_h} HAPs.")
            solution = None

    return solution

# ## Find the optimal placement of the HAPs to reach the maximum key generation for all the end-to-end paths between GS pairs. 
# def placement(gss, haps, links):
#     c1_list_ref, c2_list_ref, c3_list_ref = [], [], []
#     lon0, lat0 = 279, 49

#     demands = []
#     for idx_g1, g1 in enumerate(gss):
#         for idx_g2, g2 in enumerate(gss):
#             if idx_g2 != idx_g1:
#                 demands.append(
#                     demand(
#                         100,
#                         gss[idx_g1],
#                         gss[idx_g2]
#                     )
#                 )

#     # Create Optimization Model
#     m = gp.Model("hap-qkd")

#     for lon, lat in zip(haps[0].lg, haps[0].la):
#         c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
#         c1_list_ref.append(c1)
#         c2_list_ref.append(c2)
#     c3_list_ref = haps[0].H

#     ## Decision Variables
#     # Dictionaries of decision variables instead of MVar arrays
#     z, o = {}, {}

#     for idx_l, l in enumerate(links):
#         for idx_d, d in enumerate(demands):
#             for t in syst.T:
#                 z[idx_l, idx_d, t] = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}", vtype=GRB.BINARY)

#     nodes = gss + haps

#     # ## Node order variable --> To prevent subtours
#     # for idx_n, n in enumerate(nodes):
#     #     for idx_d, d in enumerate(demands):
#     #         for t in syst.T:
#     #             o[idx_n, idx_d, t] = m.addVar(name=f"o_{idx_n}_{idx_d}_{t}", vtype=GRB.CONTINUOUS)

#     # k (key rate) and d (distance)
#     k, d, kz = {}, {}, {}
#     for idx_l, l in enumerate(links):
#         for t in syst.T:
#             d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                
#             k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
#             kz[idx_l, t] = m.addVar(name=f"kz_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#             dpts = np.linspace(15 * COORDINATE_SCALE, 2e2 * COORDINATE_SCALE, 2)
#             kpts = [calculate_key_rate_planning("theoretical", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

#             #print(f"kpts: {kpts}")
#             m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

#     # Coordinate decision variables for each HAP and time
#     c1, c2 = {}, {}  # x, y in km
#     for idx_h, hnode in enumerate(haps):
#         for t in syst.T:
#             c1[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
#             c2[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")

#     # Secondary objective: maximize d
#     m.setObjective(gp.quicksum(gp.quicksum(kz[idx_l, t]
#                                            for idx_l, l in enumerate(links)
#                                           )
#                                for t in syst.T
#                               ), GRB.MAXIMIZE)

#     m.setParam("MIPGap", 1e-2)
#     m.setParam("MIPGapAbs", 1e-2)
#     m.setParam("FeasibilityTol", 1e-2)
#     m.setParam("IntFeasTol", 1e-2)
#     m.setParam("OptimalityTol", 1e-2)
#     m.Params.Presolve = 2
#     m.Params.Method = 2
#     m.Params.Cuts = 2
#     m.Params.Heuristics = 0.25
#     m.Params.MIPFocus = 1
#     m.Params.NumericFocus = 1
#     m.Params.Threads = 0
#     m.Params.NodefileStart = 0.5
#     m.Params.NoRelHeurTime = 20
#     m.Params.ConcurrentMIP = 1

#     ## Constraints
#     # # Maximum Tx/Rx Connection
#     # m.addConstrs(
#     #     (gp.quicksum(z[idx_l, idx_d, t]
#     #                  for idx_l, l in enumerate(links)
#     #                  for idx_d, d in enumerate(demands)
#     #                  if l.n1 == n
#     #                 ) <= n.N_TX
#     #      for idx_n, n in enumerate(nodes)
#     #      for t in syst.T),
#     #     name="max_tx_connections"
#     # )

#     # m.addConstrs(
#     #     (gp.quicksum(z[idx_l, idx_d, t]
#     #                  for idx_l, l in enumerate(links)
#     #                  for idx_d, d in enumerate(demands)
#     #                  if l.n2 == n
#     #                 ) <= n.N_RX
#     #      for idx_n, n in enumerate(nodes)
#     #      for t in syst.T),
#     #     name="max_rx_connections"
#     # )

#     # Flow conservation
#     m.addConstrs(
#         (gp.quicksum(z[idx_l, idx_d, t]
#                      for idx_l, l in enumerate(links)
#                      if l.n1 == d.n1)
#          - gp.quicksum(z[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if l.n2 == d.n1)
#          == 1
#          for idx_d, d in enumerate(demands)
#          for t in syst.T),
#         name="flow_conservation_1"
#     )

#     m.addConstrs(
#         (gp.quicksum(z[idx_l, idx_d, t]
#                      for idx_l, l in enumerate(links)
#                      if l.n2 == d.n2)
#          - gp.quicksum(z[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if l.n1 == d.n2)
#          == 1
#          for idx_d, d in enumerate(demands)
#          for t in syst.T),
#         name="flow_conservation_2"
#     )

#     m.addConstrs(
#         (gp.quicksum(z[idx_l, idx_d, t]
#                      for idx_l, l in enumerate(links)
#                      if l.n1 == n
#                     )
#          - gp.quicksum(z[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if l.n2 == n
#                       )
#          == 0
#          for idx_d, d in enumerate(demands)
#          for n in gss + haps
#          if n != d.n1 and n != d.n2
#          for t in syst.T),
#         name="flow_conservation_3"
#     )

#     # # MTZ subtour elimination --> Eliminates pointless single/multi-hop loops in the flows --> Uses an ordering values for all nodes
#     # # --> The order values should only increase on the path --> A decrease in order value == a loop (X)
#     # M = len(nodes)
#     # m.addConstrs(
#     #     (
#     #         o[nodes.index(l.n2), idx_d, t] >= o[nodes.index(l.n1), idx_d, t] + 1 - M * (1 - z[idx_l, idx_d, t])
#     #         for idx_l, l in enumerate(links)
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="ordering_constraint_1"
#     # )
#     # m.addConstrs(
#     #     (
#     #         o[nodes.index(d.n1), idx_d, t] == 1
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="ordering_constraint_2"
#     # )
#     # m.addConstrs(
#     #     (
#     #         o[nodes.index(d.n2), idx_d, t] == M
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="ordering_constraint_2"
#     # )

#     ########### Exclusive constraints ###########

#     # Active link check
#     m.addConstrs(
#         (
#             kz[idx_l, t] <= k[idx_l, t]
#             for idx_l, l in enumerate(links)
#             for t in syst.T
#         ), name="keyrate_active_link_1"
#     )

#     m.addConstrs(
#         (
#             kz[idx_l, t] <= 1e7 * gp.quicksum(z[idx_l, idx_d, t]
#                                                for idx_d, d in enumerate(demands)
#                                               )
#             for idx_l, l in enumerate(links)
#             for t in syst.T
#         ), name="keyrate_active_link_2"
#     )

#     m.addConstrs(
#         (
#             kz[idx_l, t] >= k[idx_l, t] - 1e7 * (1 - gp.quicksum(z[idx_l, idx_d, t]
#                                                                   for idx_d, d in enumerate(demands)
#                                                                  )
#                                                  )
#             for idx_l, l in enumerate(links)
#             for t in syst.T
#         ), name="keyrate_active_link_3"
#     )

#     m.addConstrs(
#         (c1[idx_h, t] == c1_list_ref[t] + (c1[idx_h, 0] - c1_list_ref[0])
#          for idx_h, h in enumerate(haps)
#          for t in syst.T if t >= 1),
#         name="shift_trajectory_1"
#     )

#     m.addConstrs(
#         (c2[idx_h, t] == c2_list_ref[t] + (c2[idx_h, 0] - c2_list_ref[0])
#          for idx_h, h in enumerate(haps)
#          for t in syst.T if t >= 1),
#         name="shift_trajectory_2"
#     )

#     # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
#     for idx_l, l in enumerate(links):
#         # identify which endpoint is HAP and which is GS
#         if isinstance(l.n1, hap) and isinstance(l.n2, gs):
#             hap_idx, gs_node = haps.index(l.n1), l.n2
#         elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
#             hap_idx, gs_node = haps.index(l.n2), l.n1

#         [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)

#         for t in syst.T:
#             dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
#             dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
#             m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + haps[hap_idx].H[t]*haps[hap_idx].H[t]*COORDINATE_SCALE*COORDINATE_SCALE,
#                          name=f"dist_cone_{idx_l}_{t}")
                         
#     ## Solve
#     m.optimize()

#     if m.status == GRB.OPTIMAL:
#         print("\n=========== OPTIMAL SOLUTION FOUND ===========")

#         solution = {
#             "z":   {k: v.X for k, v in z.items()},
#             #"o":   {k: round(v.X, 3) for k, v in o.items()},
#             "kz":  {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in kz.items()},
#             "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
#             "d":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in d.items()},
#             "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
#             "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()}
#         }

#         key_rate_m = min(solution["k"][idx_l, t]
#                          for idx_l, l in enumerate(links)
#                          for t in syst.T
#                         )
#         #key_rate_s = sum(solution["k"])
        
#         print(f"Minimum link key rate: {key_rate_m}") #, Sum: {key_rate_s}")

#         for idx_l, l in enumerate(links):
#             for idx_d, d in enumerate(demands):
#                 Z = sum(solution["z"][idx_l, idx_d, t]
#                         for t in syst.T
#                        )
#                 if Z > 0:
#                     print(f"SELECTED: link.n1: {l.n1.tag}, link.n2: {l.n2.tag}")
                    
#         plot_z_timeline(solution["z"], links, demands, figsize=(10,6))
        
#         # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
#         # pp.pprint(solution)

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
#                 x = solution["c1"].get((idx_hnode,t))
#                 y = solution["c2"].get((idx_hnode,t))
#                 if x is not None and y is not None:
#                     # lon, lat = xy_to_lonlat(x, y)
#                     lon, lat = tangent_to_latlon(x, y, 279, 49)
#                     lon_series.append(lon)   # shift if needed
#                     lat_series.append(lat)
#             planned_lons.append(lon_series)
#             planned_lats.append(lat_series)
#         planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]
        
#         # Ground Stations (replicate coordinates across all T so they plot in animation)
#         gs_lons = []
#         gs_lats = []
#         for gnode in gss:
#             gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
#             gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
#         gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]

#         all_lons = planned_lons + gs_lons
#         all_lats = planned_lats + gs_lats
#         all_labels = planned_labels + gs_labels

#         plot_connectivity_graph_planning(gss, haps, links, 
#                         planned_lons=planned_lons,
#                         planned_lats=planned_lats,
#                         planned_labels=planned_labels)

#         plot_connectivity_graph_planning_3d(gss, haps, links, 
#                         planned_lons=planned_lons,
#                         planned_lats=planned_lats,
#                         planned_labels=planned_labels)

#         # return solution
#     else:
#         solution = None

#     return solution

# ## Find the optimal placement of the HAPs to reach the maximum key generation for all the end-to-end paths between GS pairs. 
# def placement(gss, haps, links):
#     c1_list_ref, c2_list_ref, c3_list_ref = [], [], []
#     lon0, lat0 = 279, 49

#     demands = []
#     for idx_g1, g1 in enumerate(gss):
#         for idx_g2, g2 in enumerate(gss):
#             if idx_g2 != idx_g1:
#                 demands.append(
#                     demand(
#                         100,
#                         gss[idx_g1],
#                         gss[idx_g2]
#                     )
#                 )

#     # Create Optimization Model
#     m = gp.Model("hap-qkd")

#     for lon, lat in zip(haps[0].lg, haps[0].la):
#         c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
#         c1_list_ref.append(c1)
#         c2_list_ref.append(c2)
#     c3_list_ref = haps[0].H

#     ## Decision Variables
#     # Dictionaries of decision variables instead of MVar arrays
#     r_1, r_h, a, z, o = {}, {}, {}, {}, {}

#     for idx_l, l in enumerate(links):
#         for idx_d, d in enumerate(demands):
#             for t in syst.T:
#                 r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
#                 z[idx_l, idx_d, t] = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}", vtype=GRB.BINARY)

#     for idx_d, d in enumerate(demands):
#         for t in syst.T:
#             r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#     nodes = gss + haps

#     # ## Node order variable --> To prevent subtours
#     # for idx_n, n in enumerate(nodes):
#     #     for idx_d, d in enumerate(demands):
#     #         for t in syst.T:
#     #             o[idx_n, idx_d, t] = m.addVar(name=f"o_{idx_n}_{idx_d}_{t}", vtype=GRB.CONTINUOUS)

#     # k (key rate) and d (distance)
#     k, d, kz = {}, {}, {}
#     for idx_l, l in enumerate(links):
#         for t in syst.T:
#             d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                
#             k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
#             kz[idx_l, t] = m.addVar(name=f"kz_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#             dpts = np.linspace(15 * COORDINATE_SCALE, 2e2 * COORDINATE_SCALE, 4)
#             kpts = [calculate_key_rate_planning("theoretical", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

#             #print(f"kpts: {kpts}")
#             m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

#     # Coordinate decision variables for each HAP and time
#     c1, c2 = {}, {}  # x, y in km
#     for idx_h, hnode in enumerate(haps):
#         for t in syst.T:
#             c1[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
#             c2[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")

#     # Secondary objective: maximize d
#     m.setObjective(gp.quicksum(gp.quicksum(r_h[idx_d, t]
#                                            for idx_d, d in enumerate(demands)
#                                           )
#                                for t in syst.T
#                               ), GRB.MAXIMIZE)

#     m.setParam("MIPGap", 1e-2)
#     m.setParam("MIPGapAbs", 1e-2)
#     m.setParam("FeasibilityTol", 1e-2)
#     m.setParam("IntFeasTol", 1e-2)
#     m.setParam("OptimalityTol", 1e-2)
#     m.Params.Presolve = 2
#     m.Params.Method = 2
#     m.Params.Cuts = 2
#     m.Params.Heuristics = 0.25
#     m.Params.MIPFocus = 1
#     m.Params.NumericFocus = 1
#     m.Params.Threads = 0
#     m.Params.NodefileStart = 0.5
#     m.Params.NoRelHeurTime = 20
#     m.Params.ConcurrentMIP = 1

#     ## Constraints
#     # Maximum Tx/Rx Connection
#     m.addConstrs(
#         (gp.quicksum(z[idx_l, idx_d, t]
#                      for idx_l, l in enumerate(links)
#                      for idx_d, d in enumerate(demands)
#                      if l.n1 == n
#                     ) <= n.N_TX
#          for idx_n, n in enumerate(nodes)
#          for t in syst.T),
#         name="max_tx_connections"
#     )

#     m.addConstrs(
#         (gp.quicksum(z[idx_l, idx_d, t]
#                      for idx_l, l in enumerate(links)
#                      for idx_d, d in enumerate(demands)
#                      if l.n2 == n
#                     ) <= n.N_RX
#          for idx_n, n in enumerate(nodes)
#          for t in syst.T),
#         name="max_rx_connections"
#     )

#     # Flow conservation
#     m.addConstrs(
#         (gp.quicksum(r_1[idx_l, idx_d, t]
#                      for idx_l, l in enumerate(links)
#                      if l.n1 == d.n1)
#          - gp.quicksum(r_1[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if l.n2 == d.n1)
#          == r_h[idx_d, t]
#          for idx_d, d in enumerate(demands)
#          for t in syst.T),
#         name="flow_conservation_1"
#     )

#     m.addConstrs(
#         (gp.quicksum(r_1[idx_l, idx_d, t]
#                      for idx_l, l in enumerate(links)
#                      if l.n2 == d.n2)
#          - gp.quicksum(r_1[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if l.n1 == d.n2)
#          == r_h[idx_d, t]
#          for idx_d, d in enumerate(demands)
#          for t in syst.T),
#         name="flow_conservation_2"
#     )

#     m.addConstrs(
#         (gp.quicksum(r_1[idx_l, idx_d, t]
#                      for idx_l, l in enumerate(links)
#                      if l.n1 == n
#                     )
#          - gp.quicksum(r_1[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if l.n2 == n
#                       )
#          == 0
#          for idx_d, d in enumerate(demands)
#          for n in gss + haps
#          if n != d.n1 and n != d.n2
#          for t in syst.T),
#         name="flow_conservation_3"
#     )

#     # # MTZ subtour elimination --> Eliminates pointless single/multi-hop loops in the flows --> Uses an ordering values for all nodes
#     # # --> The order values should only increase on the path --> A decrease in order value == a loop (X)
#     # M = len(nodes)
#     # m.addConstrs(
#     #     (
#     #         o[nodes.index(l.n2), idx_d, t] >= o[nodes.index(l.n1), idx_d, t] + 1 - M * (1 - z[idx_l, idx_d, t])
#     #         for idx_l, l in enumerate(links)
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="ordering_constraint_1"
#     # )
#     # m.addConstrs(
#     #     (
#     #         o[nodes.index(d.n1), idx_d, t] == 1
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="ordering_constraint_2"
#     # )
#     # m.addConstrs(
#     #     (
#     #         o[nodes.index(d.n2), idx_d, t] == M
#     #         for idx_d, d in enumerate(demands)
#     #         for t        in syst.T
#     #     ), name="ordering_constraint_2"
#     # )
    
#     # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
#     m.addConstrs(
#         (r_h[idx_d, t] >= r_1[idx_l, idx_d, t]
#          for idx_l, l in enumerate(links)
#          for idx_d, d in enumerate(demands)
#          for t in syst.T),
#         name="demand_link_coordination_1"
#     )

#     # Key rate and routing coordination (1)
#     m.addConstrs(
#         (r_1[idx_l, idx_d, t] >= 1 * z[idx_l, idx_d, t]
#          for idx_l, l in enumerate(links)
#          for idx_d, d in enumerate(demands)
#          for t in syst.T),
#         name="key_rate_routing_coordination_1"
#     )

#     # Key rate and routing coordination (2)
#     # We do not limit r_1 to d.K_REQ to enforce the optimal placement as well as demand satisfaction.
#     # Note that l.K_MAX[t] is not a fixed value in this problem.
#     m.addConstrs(
#         (r_1[idx_l, idx_d, t] <= kz[idx_l, t]
#          for idx_l, l in enumerate(links)
#          for idx_d, d in enumerate(demands)
#          for t in syst.T),
#         name="key_rate_routing_coordination_2"
#     )

#     ########### Exclusive constraints ###########

#     # Active link check
#     m.addConstrs(
#         (
#             kz[idx_l, t] <= k[idx_l, t]
#             for idx_l, l in enumerate(links)
#             for t in syst.T
#         ), name="keyrate_active_link_1"
#     )

#     m.addConstrs(
#         (
#             kz[idx_l, t] <= 1e7 * gp.quicksum(z[idx_l, idx_d, t]
#                                                for idx_d, d in enumerate(demands)
#                                               )
#             for idx_l, l in enumerate(links)
#             for t in syst.T
#         ), name="keyrate_active_link_2"
#     )

#     m.addConstrs(
#         (
#             kz[idx_l, t] >= k[idx_l, t] - 1e7 * (1 - gp.quicksum(z[idx_l, idx_d, t]
#                                                                   for idx_d, d in enumerate(demands)
#                                                                  )
#                                                  )
#             for idx_l, l in enumerate(links)
#             for t in syst.T
#         ), name="keyrate_active_link_3"
#     )

#     m.addConstrs(
#         (c1[idx_h, t] == c1_list_ref[t] + (c1[idx_h, 0] - c1_list_ref[0])
#          for idx_h, h in enumerate(haps)
#          for t in syst.T if t >= 1),
#         name="shift_trajectory_1"
#     )

#     m.addConstrs(
#         (c2[idx_h, t] == c2_list_ref[t] + (c2[idx_h, 0] - c2_list_ref[0])
#          for idx_h, h in enumerate(haps)
#          for t in syst.T if t >= 1),
#         name="shift_trajectory_2"
#     )

#     # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
#     for idx_l, l in enumerate(links):
#         # identify which endpoint is HAP and which is GS
#         if isinstance(l.n1, hap) and isinstance(l.n2, gs):
#             hap_idx, gs_node = haps.index(l.n1), l.n2
#         elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
#             hap_idx, gs_node = haps.index(l.n2), l.n1

#         [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)

#         for t in syst.T:
#             dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
#             dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
#             m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + haps[hap_idx].H[t]*haps[hap_idx].H[t]*COORDINATE_SCALE*COORDINATE_SCALE,
#                          name=f"dist_cone_{idx_l}_{t}")
                         
#     ## Solve
#     m.optimize()

#     if m.status == GRB.OPTIMAL:
#         print("\n=========== OPTIMAL SOLUTION FOUND ===========")

#         solution = {
#             "r_h":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
#             "r_1":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
#             "z":   {k: v.X for k, v in z.items()},
#             #"o":   {k: round(v.X, 3) for k, v in o.items()},
#             "kz":  {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in kz.items()},
#             "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
#             "d":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in d.items()},
#             "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
#             "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()}
#         }

#         key_rate_m = min(solution["k"][idx_l, t]
#                          for idx_l, l in enumerate(links)
#                          for t in syst.T
#                         )
#         #key_rate_s = sum(solution["k"])
        
#         print(f"Minimum link key rate: {key_rate_m}") #, Sum: {key_rate_s}")

#         for idx_l, l in enumerate(links):
#             for idx_d, d in enumerate(demands):
#                 Z = sum(solution["z"][idx_l, idx_d, t]
#                         for t in syst.T
#                        )
#                 if Z > 0:
#                     print(f"SELECTED: link.n1: {l.n1.tag}, link.n2: {l.n2.tag}")
                    
#         plot_z_timeline(solution["z"], links, demands, figsize=(10,6))
        
#         # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
#         # pp.pprint(solution)

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
#                 x = solution["c1"].get((idx_hnode,t))
#                 y = solution["c2"].get((idx_hnode,t))
#                 if x is not None and y is not None:
#                     # lon, lat = xy_to_lonlat(x, y)
#                     lon, lat = tangent_to_latlon(x, y, 279, 49)
#                     lon_series.append(lon)   # shift if needed
#                     lat_series.append(lat)
#             planned_lons.append(lon_series)
#             planned_lats.append(lat_series)
#         planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]
        
#         # Ground Stations (replicate coordinates across all T so they plot in animation)
#         gs_lons = []
#         gs_lats = []
#         for gnode in gss:
#             gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
#             gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
#         gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]

#         all_lons = planned_lons + gs_lons
#         all_lats = planned_lats + gs_lats
#         all_labels = planned_labels + gs_labels

#         plot_connectivity_graph_planning(gss, haps, links, 
#                         planned_lons=planned_lons,
#                         planned_lats=planned_lats,
#                         planned_labels=planned_labels)

#         plot_connectivity_graph_planning_3d(gss, haps, links, 
#                         planned_lons=planned_lons,
#                         planned_lats=planned_lats,
#                         planned_labels=planned_labels)

#         print(f"The minimum number of required HAPs: {num_h}")

#         # return solution
#     else:
#         print(f"No optimal solution found for {num_h} HAPs.")
#         solution = None

#     return solution



import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (needed for 3D projection)


def plot_connectivity_graph_planning_3d(gnodes, hnodes, links,
                                        planned_lons=None, planned_lats=None, planned_alts=None,
                                        planned_labels=None, alg=""):
    """
    Plot 3D connectivity graph of GS and HAP nodes, with optional planned 3D HAP trajectories.

    Parameters:
    -----------
    gnodes : list
        Ground stations with attributes 'la' (latitude), 'lg' (longitude), and optional 'tag'.
        Altitude is assumed to be 0.
    hnodes : list
        HAP nodes with trajectory lists ('la', 'lg', 'H') and optional 'tag'.
    links : list
        Link objects with attributes 'n1', 'n2'.
    planned_lons, planned_lats, planned_alts : list[list[float]], optional
        Planned longitude, latitude, and altitude trajectories for each HAP (same order as hnodes).
    planned_labels : list[str], optional
        Labels for planned trajectories.
    """

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111, projection='3d')

    # --- Plot GS nodes ---
    for gs_node in gnodes:
        ax.scatter(gs_node.lg, gs_node.la, 0, color='skyblue', s=80, zorder=5, marker='^')
        if hasattr(gs_node, 'tag'):
            if gs_node.tag == "Timmins":
                ax.text(gs_node.lg - 0.8, gs_node.la + 0.04, 0.05, gs_node.tag,
                    fontsize=9)
            else:
                ax.text(gs_node.lg + 0.04, gs_node.la + 0.04, 0.05, gs_node.tag,
                    fontsize=9)

    # --- Plot optimal trajectories (if provided) ---
    if planned_lons and planned_lats:
        for idx, (plon, plat) in enumerate(zip(planned_lons, planned_lats)):

            # print(f"idx: {idx}, hnodes[0].H[idx]: {hnodes[0].H[idx]}")
            # print(f"plon:{plon}")
            # Plot trajectory
            ax.plot(plon, plat, hnodes[0].H, color='red', linewidth=1.1, alpha=1)

            # Mark initial point
            ax.scatter(plon[0], plat[0], hnodes[0].H[0], color='red', s=5, marker='o', zorder=6)

            # Label near last planned point
            lon_last, lat_last, alt_last = plon[-1], plat[-1], hnodes[0].H[-1]
            coord_text = f"{idx}"
            ax.text(lon_last - 0.4, lat_last - 0.2, alt_last + 0.2,
                    coord_text, fontsize=9, fontstyle='italic', color='red')

    # --- Axis setup ---
    all_lons = [gs.lg for gs in gnodes] + [hap.lg[0] for hap in hnodes]
    all_lats = [gs.la for gs in gnodes] + [hap.la[0] for hap in hnodes]
    all_alts = [0] + [hap.H[0] for hap in hnodes]

    ax.set_xlabel("Longitude", fontsize=13, labelpad=10)
    ax.set_ylabel("Latitude", fontsize=13, labelpad=10)
    ax.set_zlabel("Altitude (km)", fontsize=13, labelpad=10)

    ax.set_xlim(min(all_lons) - 0.2, max(all_lons) + 1.8)
    ax.set_ylim(min(all_lats) - 0.2, max(all_lats) + 0.2)
    ax.set_zlim(min(all_alts) - 0.5, max(all_alts) + 1)

    # --- Legend ---
    custom_handles = [
        Line2D([], [], marker='o', color='skyblue', linestyle='None', markersize=6, label='GS'),
        Line2D([], [], marker='o', color='red', linestyle='None', markersize=6, label='HAP')
    ]
    ax.legend(handles=custom_handles, loc='upper left', frameon=True)

    # --- Final layout ---
    ax.grid(True, alpha=0.3)
    ax.view_init(elev=25, azim=-60)  # good default viewing angle
    plt.tight_layout()

    plt.savefig(f"hap_qkd_network_3d_{len(hnodes)}_{alg}.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.show()

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
    

def _compute_kpt(args):
    """Helper for parallel computation."""
    d, t = args
    k_val = calculate_key_rate_planning("plob", 0, d / COORDINATE_SCALE, t)
    return k_val * KEY_RATE_SCALE

def compute_kpts_parallel(dpts, t, max_workers=12):
    """
    Compute key rates for all distances in dpts in parallel using calculate_key_rate_planning().
    """
    print(f"⏱️ Computing key rates for {len(dpts)} distance points at t={t}...")

    tasks = [(d, t) for d in dpts]

    # Run parallel computation
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(
            tqdm(
                executor.map(_compute_kpt, tasks, chunksize=10),
                total=len(tasks),
                desc="Computing kpts"
            )
        )

    print(f"✅ Finished computing {len(dpts)} key rate points.")
    return results

def plot_z_timeline(z, links, demands, figsize=(10,6)):
    """
    Visualize z[idx_l, idx_d, t] binary results as a timeline.
    
    Parameters
    ----------
    z : 3D array-like (n_links x n_demands x n_time)
        Binary or 0/1 values.
    links : list
        List of link objects with attributes n1.tag and n2.tag
    demands : list
        List of demand objects with attributes n1.tag and n2.tag
    """
    n_links = len(links)
    n_demands = len(demands)
    n_time = len(syst.T)

    # Assign each demand a color
    cmap = plt.cm.get_cmap("tab10", n_demands)
    demand_colors = [cmap(i) for i in range(n_demands)]

    fig, ax = plt.subplots(figsize=figsize)

    for idx_l, link in enumerate(links):
        for idx_d, d in enumerate(demands):
            active_times = [t for t in syst.T if z[idx_l, idx_d, t] > 0.5]
            if active_times:
                # Group consecutive time steps into intervals
                segments = []
                start = active_times[0]
                for i in range(1, len(active_times)):
                    if active_times[i] != active_times[i-1] + 1:
                        segments.append((start, active_times[i-1]))
                        start = active_times[i]
                segments.append((start, active_times[-1]))

                # Draw horizontal bars for each active interval
                for (t1, t2) in segments:
                    ax.barh(
                        y=idx_l,
                        width=t2 - t1 + 1,
                        left=t1,
                        height=0.6,
                        color=demand_colors[idx_d],
                        alpha=0.7,
                        edgecolor='k'
                    )

    # --- Y-axis: link labels ---
    y_labels = [f"{l.n1.tag}→{l.n2.tag}" if hasattr(l.n1, "tag") and hasattr(l.n2, "tag")
                else f"L{idx}" for idx, l in enumerate(links)]
    ax.set_yticks(np.arange(n_links))
    ax.set_yticklabels(y_labels, fontsize=9)

    # --- X-axis: time steps ---
    ax.set_xticks(syst.T)
    ax.set_xlabel("Time Step", fontsize=11)
    ax.set_ylabel("Links", fontsize=11)

    # --- Legend: demands ---
    patches = [mpatches.Patch(color=demand_colors[i], label=f"Demand {i}: {demands[i].n1.tag}→{demands[i].n2.tag}")
               for i in range(n_demands)]
    ax.legend(handles=patches, loc="upper right", fontsize=9)

    ax.grid(True, axis='x', alpha=0.3)
    plt.title("Link Utilization Timeline by Demand (z variables)", fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.show()



































    



# ## problem: "FIXED", "INIT"
# ## For a given set of demands what is the minimum number of HAPs and where to place them to satisfy all the demands.
# def demand_feasibility(problem, gss, demands):
#     MAX_HAPS = 10

#     c1_list_ref, c2_list_ref, c3_list_ref = [], [], []
#     lon0, lat0 = 279, 49

#     for num_h in range(1, MAX_HAPS):
#         haps  = []
#         links = []

#         # Create Optimization Model
#         m = gp.Model("hap-qkd")
        
#         for n in range(num_h):
#             haps.append(hap([279]*len(syst.T), [49]*len(syst.T), [15]*len(syst.T), 1, 1, 1e9, f"HAP_{n}"))

#         # Update coordinates depending on model choice
#         if SYNTH_STRATO == 1:
#             update_coordinates("stratotegic", haps, syst)
#         else:
#             update_coordinates("wind", haps, syst)

#         ## Only once is enough
#         if num_h == 1:
#             for lon, lat in zip(haps[0].lg, haps[0].la):
#                 c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
#                 c1_list_ref.append(c1)
#                 c2_list_ref.append(c2)
#             c3_list_ref = haps[0].H

#         print(f"haps[0].lg: {haps[0].lg}")
#         print(f"haps[0].la: {haps[0].la}")
#         print(f"c1_list_ref: {c1_list_ref}")
#         print(f"c2_list_ref: {c2_list_ref}")
#         print(f"c3_list_ref: {c3_list_ref}")
    
#         # Links: connect only GSs to HAPs
#         for gs_node in gss:
#             for hap_node in haps:
#                 links.append(link(gs_node, hap_node,
#                                   [100]*len(syst.T),
#                                   [(0,0,0)]*len(syst.T),
#                                   [1]*len(syst.T)))
#                 links.append(link(hap_node, gs_node,
#                                   [100]*len(syst.T),
#                                   [(0,0,0)]*len(syst.T),
#                                   [1]*len(syst.T)))

#         ## Decision Variables
#         # Dictionaries of decision variables instead of MVar arrays
#         r_1, r_2, r_h, a, z, u = {}, {}, {}, {}, {}, {}
    
#         for idx_l, l in enumerate(links):
#             for idx_d, d in enumerate(demands):
#                 for t in syst.T:
#                     r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
#                     r_2[idx_l, idx_d, t] = m.addVar(name=f"r_2_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
#                     z[idx_l, idx_d, t]   = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}",   vtype=GRB.BINARY)
                    
#         for idx_d, d in enumerate(demands):
#             for t in syst.T:
#                 r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
#                 # u[idx_d, t]   = m.addVar(name=f"u_{idx_d}_{t}",   vtype=GRB.CONTINUOUS, lb=0.0)
                
#         for idx_l, l in enumerate(links):
#             for t in syst.T:
#                 a[idx_l, t] = m.addVar(name=f"a_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
    
#         # k (key rate) and d (distance)
#         k, d, kz = {}, {}, {}
#         for idx_l, l in enumerate(links):
#             for t in syst.T:
#                 d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                    
#                 k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
#                 kz[idx_l, t] = m.addVar(name=f"kz_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
    
#                 dpts = np.linspace(15 * COORDINATE_SCALE, 2e2 * COORDINATE_SCALE, 4)
#                 kpts = [calculate_key_rate_planning("plob", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]
#                 # kpts = compute_kpts_parallel(
#                 #            dpts=dpts,
#                 #            t=t,
#                 #            max_workers=24
#                 #        )
    
#                 m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")
    
#         # d_min, d_max = 15.0, 200.0
#         # breakpoint_counts = [5, 10, 20, 40]  # numbers of breakpoints to test
#         # plt.figure(figsize=(8,5))
#         # for n_bp in breakpoint_counts:
#         #     # breakpoints for this case
#         #     dpts = np.linspace(d_min, d_max, n_bp)
#         #     #kpts = [calculate_key_rate_planning(0, d) * KEY_RATE_SCALE for d in dpts]
#         #     kpts = [calculate_key_rate_planning("theoretical", 0, d, 0) * KEY_RATE_SCALE for d in dpts]
        
#         #     # plot this curve
#         #     plt.plot(dpts, kpts, marker="o", linestyle="-", label=f"{n_bp} breakpoints")
#         # plt.xlabel("Distance (dpts)")
#         # plt.ylabel("Key Rate (kpts)")
#         # plt.title("Key Rate vs Distance for Different Breakpoint Counts")
#         # plt.grid(True)
#         # plt.legend()
#         # plt.tight_layout()
#         # plt.show()
                
#         # Coordinate decision variables for each HAP and time
#         c1, c2 = {}, {}  # x, y in km
#         for idx_h, hnode in enumerate(haps):
#             for t in syst.T:
#                 c1[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
#                 c2[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")
    
        
#         # Secondary objective: maximize d
#         if problem == "FIXED" or problem == "INIT":
#             m.setObjective(gp.quicksum(gp.quicksum(r_h[idx_d, t]
#                                                    for idx_d, d in enumerate(demands)
#                                                   )
#                                        for t in syst.T
#                                       ), GRB.MAXIMIZE)
    
#         elif problem == "FREE":
#             m.setObjectiveN(-gp.quicksum(gp.quicksum(d[idx_l, t]
#                                                      for idx_l, l in enumerate(links)
#                                                     )
#                                          for t in syst.T
#                                         ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")
    
#         m.setParam("MIPGap", 1e-3)          # force very tight gap
#         m.setParam("MIPGapAbs", 1e-3)
#         m.setParam("FeasibilityTol", 1e-3)
#         m.setParam("IntFeasTol", 1e-3)
#         m.setParam("OptimalityTol", 1e-3)
    
#         m.Params.Presolve = 2           # Aggressive presolve
#         m.Params.Method = 2             # Barrier for root relaxation
#         m.Params.Cuts = 2               # More cutting planes
#         m.Params.Heuristics = 0.25      # Reasonable level of primal heuristics
#         m.Params.MIPFocus = 1           # Emphasize feasible solutions first
#         m.Params.NumericFocus = 1       # Reduce rounding issues due to tiny epsilons
#         m.Params.Threads = 0            # Use all cores
    
#         m.Params.NodefileStart = 0.5    # Start writing node files to disk to prevent RAM overflow
#         m.Params.NoRelHeurTime = 20     # Time for heuristics before branching
#         m.Params.ConcurrentMIP = 1      # Use multiple concurrent solvers for first feasible solution
    
#         ## Constraints
#         # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
#         # r_h = min_{l:z_l=1}(r_1+r_2)
#         m.addConstrs(
#             (
#                 r_h[idx_d, t] <= r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] + 1e10 * KEY_RATE_SCALE * (1 - z[idx_l, idx_d, t])
#                 for idx_l, l in enumerate(links)
#                 for idx_d, d in enumerate(demands)
#                 for t        in syst.T
#             ), name="demand_link_coordination"
#         )
        
#         m.addConstrs(
#             (
#                 gp.quicksum(r_h[idx_d, t]
#                             for t in syst.T
#                            ) >= sum(d.K_REQ[t]
#                                     for t in syst.T
#                                    ) * KEY_RATE_SCALE
#                 for idx_d, d in enumerate(demands)
#             ), name="demand_satisfaction_guarantee"
#         )

        
#         m.addConstrs(
#             (
#                 gp.quicksum(r_1[idx_l, idx_d, t]
#                             for idx_d, d in enumerate(demands)
#                            ) <= k[idx_l, t]
#                 for idx_l, l in enumerate(links)
#                 for t        in syst.T
#             ), name="max_key_rate"
#         )

#         m.addConstrs(
#             (
#                 k[idx_l, t] <= 1e10 * gp.quicksum(z[idx_l, idx_d, t]
#                                                   for idx_d, d in enumerate(demands)
#                                                  )
#                 for idx_l, l in enumerate(links)
#                 for t        in syst.T
#             ), name="max_key_rate"
#         )

#         EPSILON = 1e-5
#         # Key rate and routing coordination (1)
#         m.addConstrs(
#             (
#                 r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] >= EPSILON * z[idx_l, idx_d, t]
#                 for idx_l, l in enumerate(links)
#                 for idx_d, d in enumerate(demands)
#                 for t        in syst.T
#             ), name="demand_link_coordination_1"
#         )
#         # Key rate and routing coordination (2)
#         m.addConstrs(
#             (
#                 r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] <= 1e10 * KEY_RATE_SCALE * z[idx_l, idx_d, t]
#                 for idx_l, l in enumerate(links)
#                 for idx_d, d in enumerate(demands)
#                 for t        in syst.T
#             ), name="key_rate_routing_coordination_2"
#         )
    
#         # Flow conservation
#         m.addConstrs(
#             (
#                 gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
#                             for idx_l, l in enumerate(links)
#                             if isinstance(l.n1, gs)
#                             if gss.index(l.n1) == gss.index(d.n1)
#                            ) - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
#                                            for idx_l, l in enumerate(links)
#                                            if isinstance(l.n2, gs)
#                                            if gss.index(l.n2) == gss.index(d.n1)
#                                           ) == r_h[idx_d, t]
#                 for idx_d, d in enumerate(demands)
#                 for t        in syst.T
#             ), name="flow_conservation_1"
#         )
#         m.addConstrs(
#             (
#                 gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
#                             for idx_l, l in enumerate(links)
#                             if isinstance(l.n2, gs)
#                             if gss.index(l.n2) == gss.index(d.n2)
#                            ) - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
#                                            for idx_l, l in enumerate(links)
#                                            if isinstance(l.n1, gs)
#                                            if gss.index(l.n1) == gss.index(d.n2)
#                                           ) == r_h[idx_d, t]
#                 for idx_d, d in enumerate(demands)
#                 for t        in syst.T
#             ), name="flow_conservation_2"
#         )
#         m.addConstrs(
#             (
#                 gp.quicksum((r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]) * (l.n1 == n)
#                             for idx_l, l in enumerate(links)
#                            ) - gp.quicksum((r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]) * (l.n2 == n)
#                                            for idx_l, l in enumerate(links) 
#                                           ) == 0
#                 for idx_d, d in enumerate(demands)
#                 for n in gss + haps
#                 if  n != d.n1 and n != d.n2
#                 for t in syst.T
#             ), name="flow_conservation_3"
#         )
#         m.addConstrs(
#             (
#                 z[idx_l_1, idx_d, t] + z[idx_l_2, idx_d, t] <= 1
#                 for idx_l_1, l_1 in enumerate(links)
#                 for idx_l_2, l_2 in enumerate(links)
#                 if  idx_l_1 < idx_l_2 and (l_1.n1 == l_2.n2) and (l_1.n2 == l_2.n1)
#                 for idx_d, d     in enumerate(demands)
#                 for t            in syst.T
#             ), name="loop_prevention"
#         )
        
#         # Maximum Tx/Rx Connection
#         m.addConstrs(
#             (
#                 gp.quicksum(z[idx_l, idx_d, t] * (l.n1 == n)
#                             for idx_l, l in enumerate(links)
#                             for idx_d, d in enumerate(demands)
#                            ) <= n.N_TX
#                 for idx_n, n in enumerate(gss + haps)
#                 for t        in syst.T
#             ), name="max_tx_connections"
#         )
        
#         m.addConstrs(
#             (
#                 gp.quicksum(z[idx_l, idx_d, t] * (l.n2 == n)
#                             for idx_l, l in enumerate(links)
#                             for idx_d, d in enumerate(demands)
#                            ) <= n.N_RX
#                 for idx_n, n in enumerate(gss + haps)
#                 for t        in syst.T
#             ), name="max_rx_connections"
#         )

#         # QKP on HAPs and GSs
#         m.addConstrs(
#             (
#                 gp.quicksum(a[idx_l, tp]
#                             for tp in range(t-1)
#                            ) >= syst.THETA * gp.quicksum(r_2[idx_l, idx_d, t]
#                                                          for idx_d, d in enumerate(demands)
#                                                         ) * STORAGE_SCALE
#                 for idx_l, l in enumerate(links)
#                 for t        in syst.T[1:]
#             ), name="qkp_min_capacity"
#         )

#         m.addConstrs(
#             (
#                 r_2[idx_l, idx_d, 0] == 0
#                 for idx_l, l in enumerate(links)
#                 for idx_d, d in enumerate(demands)
#             ), name="initial_empty_QKP"
#         )

#         # m.addConstrs(
#         #     (
#         #         kz[idx_l, t] <= k[idx_l, t]
#         #         for idx_l, l in enumerate(links)
#         #         for idx_d, d in enumerate(demands)
#         #     ), name="keyrate_active_link_1"
#         # )

#         # m.addConstrs(
#         #     (
#         #         kz[idx_l, t] <= 1e10 * gp.quicksum(z[idx_l, idx_d, t]
#         #                                            for idx_d, d in enumerate(demands)
#         #                                           )
#         #         for idx_l, l in enumerate(links)
#         #         for idx_d, d in enumerate(demands)
#         #     ), name="keyrate_active_link_2"
#         # )

#         # m.addConstrs(
#         #     (
#         #         kz[idx_l, t] >= k[idx_l, t] - 1e10 * (1 - gp.quicksum(z[idx_l, idx_d, t]
#         #                                                               for idx_d, d in enumerate(demands)
#         #                                                              )
#         #                                              )
#         #         for idx_l, l in enumerate(links)
#         #         for idx_d, d in enumerate(demands)
#         #     ), name="keyrate_active_link_3"
#         # )
        
#         m.addConstrs(
#             (
#                 a[idx_l, t] == syst.THETA * (k[idx_l, t] * KEY_RATE_SCALE - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
#                                                                                         for idx_d, d in enumerate(demands)
#                                                                                        )
#                                             ) * STORAGE_SCALE
#                 for idx_l, l in enumerate(links)
#                 for t        in syst.T
#             ), name="qkp_sequence"
#         )

#         m.addConstrs(
#             (
#                 c1[idx_h, t] == c1_list_ref[t] + (c1[idx_h, 0] - c1_list_ref[0])
#                 for idx_h, h in enumerate(haps)
#                 for t in syst.T
#                 if  t >= 1
#             ), name="shift_trajectory_1"
#         )

#         m.addConstrs(
#             (
#                 c2[idx_h, t] == c2_list_ref[t] + (c2[idx_h, 0] - c2_list_ref[0])
#                 for idx_h, h in enumerate(haps)
#                 for t in syst.T
#                 if  t >= 1
#             ), name="shift_trajectory_2"
#         )
    
#         # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
#         for idx_l, l in enumerate(links):
#             # identify which endpoint is HAP and which is GS
#             if isinstance(l.n1, hap) and isinstance(l.n2, gs):
#                 hap_idx, gs_node = haps.index(l.n1), l.n2
#             elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
#                 hap_idx, gs_node = haps.index(l.n2), l.n1
    
#             # precompute GS coordinates in same (x,y,z) frame and units (km)
#             [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)  # ensure these are constants in km
    
#             # if problem == "FIXED" or problem == "INIT":
#             #     # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
#             #     dx = c1[hap_idx] - cg1*COORDINATE_SCALE
#             #     dy = c2[hap_idx] - cg2*COORDINATE_SCALE
#             #     m.addQConstr(d[idx_l]*d[idx_l] >= dx*dx + dy*dy + 15*15*COORDINATE_SCALE*COORDINATE_SCALE,
#             #                  name=f"dist_cone_{idx_l}")
#             # elif problem == "FREE":
#             for t in syst.T:
#                 # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
#                 dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
#                 dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
#                 m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + c3_list_ref[t]*c3_list_ref[t]*COORDINATE_SCALE*COORDINATE_SCALE,
#                              name=f"dist_cone_{idx_l}_{t}")
    
#         ## Solve
#         m.optimize()
    
#         if m.status == GRB.OPTIMAL:
#             print("\n=========== OPTIMAL SOLUTION FOUND ===========")
    
#             solution = {
#                 "r_h":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
#                 "r_1":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
#                 "r_2":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
#                 "a":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in a.items()},
#                 "z":   {k: v.X for k, v in z.items()},
#                 "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
#                 "d":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in d.items()},
#                 "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
#                 "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()}
#             }

#             key_rate_m = min(solution["k"][idx_l, t]
#                              for idx_l, l in enumerate(links)
#                              for t in syst.T
#                             )
#             #key_rate_s = sum(solution["k"])
            
#             print(f"Minimum link key rate: {key_rate_m}") #, Sum: {key_rate_s}")

#             for idx_l, l in enumerate(links):
#                 for idx_d, d in enumerate(demands):
#                     Z = sum(solution["z"][idx_l, idx_d, t]
#                             for t in syst.T
#                            )
#                     if Z > 0:
#                         print(f"SELECTED: link.n1: {l.n1.tag}, link.n2: {l.n2.tag}")
                        
#             plot_z_timeline(solution["z"], links, demands, figsize=(10,6))
            
#             pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
#             pp.pprint(solution)
    
#             # Actual HAP trajectories
#             actual_lons = [hnode.lg for hnode in haps]
#             actual_lats = [hnode.la for hnode in haps]
#             actual_labels = [f"HAP_{idx_hnode}_actual" for idx_hnode, _ in enumerate(haps)]
            
#             # Planned HAP trajectories
#             planned_lons = []
#             planned_lats = []
#             for idx_hnode in range(len(haps)):
#                 lon_series = []
#                 lat_series = []
#                 for t in syst.T:
#                     x = solution["c1"].get((idx_hnode,t))
#                     y = solution["c2"].get((idx_hnode,t))
#                     if x is not None and y is not None:
#                         # lon, lat = xy_to_lonlat(x, y)
#                         lon, lat = tangent_to_latlon(x, y, 279, 49)
#                         lon_series.append(lon)   # shift if needed
#                         lat_series.append(lat)
#                 planned_lons.append(lon_series)
#                 planned_lats.append(lat_series)
#             planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]
            
#             # Ground Stations (replicate coordinates across all T so they plot in animation)
#             gs_lons = []
#             gs_lats = []
#             for gnode in gss:
#                 gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
#                 gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
#             gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]
    
#             all_lons = planned_lons + gs_lons
#             all_lats = planned_lats + gs_lats
#             all_labels = planned_labels + gs_labels
    
#             plot_connectivity_graph_planning(gss, haps, links, 
#                             planned_lons=planned_lons,
#                             planned_lats=planned_lats,
#                             planned_labels=planned_labels)

#             print(f"The minimum number of required HAPs: {num_h}")

#             # return solution
#         else:
#             print(f"No optimal solution found for {num_h} HAPs.")
#             solution = None

#     return solution
    





























# ## problem: "FREE", "FIXED", "INIT"
# def placement(problem, gss, haps, links):
#     # Create Optimization Model
#     m = gp.Model("hap-qkd")

#     demands = []
#     for idx_g1, g1 in enumerate(gss):
#         for idx_g2, g2 in enumerate(gss):
#             if idx_g1 < idx_g2:
#                 demands.append(demand(0, g1, g2))

#     ## Decision Variables
#     # Dictionaries of decision variables instead of MVar arrays
#     r, r_l, z = {}, {}, {}

#     for idx_d, d in enumerate(demands):
#         for t in syst.T:
#             r[idx_d, t] = m.addVar(name=f"r_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#     for idx_l, l in enumerate(links):
#         for idx_d, d in enumerate(demands):
#             for t in syst.T:
#                 z[idx_l, idx_d, t]   = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}", vtype=GRB.BINARY)
#                 #r_l[idx_l, idx_d, t] = m.addVar(name=f"r_l_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#     # k (key rate) and d (distance)
#     k, d = {}, {}
#     for idx_l, l in enumerate(links):
#         if problem == "FIXED" or problem == "INIT":
#             d[idx_l] = m.addVar(name=f"d_{idx_l}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
#         for t in syst.T:
#             if problem == "FREE":
#                 d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e2 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                
#             k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

#             dpts = np.linspace(15 * COORDINATE_SCALE, 2e2 * COORDINATE_SCALE, 4)
#             kpts = [calculate_key_rate_planning("plob", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

#             if problem == "FREE":
#                 m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")
#             elif problem == "FIXED" or problem == "INIT":
#                 m.addGenConstrPWL(d[idx_l], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

#     d_min, d_max = 15.0, 200.0
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
#     c1, c2 = {}, {}  # x, y in km
#     for idx_h, hnode in enumerate(haps):
#         if problem == "FIXED" or problem == "INIT":
#             c1[idx_h] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}")
#             c2[idx_h] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}")
#         elif problem == "FREE":
#             for t in syst.T:
#                 c1[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
#                 c2[idx_h, t] = m.addVar(lb=-1e2*COORDINATE_SCALE, ub=1e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")

    
#     # Secondary objective: maximize d
#     if problem == "FIXED" or problem == "INIT":
#         m.setObjective(sum(sum(r[idx_d, t]
#                                for idx_d, d in enumerate(demands)
#                               )
#                             for t in syst.T
#                           ), GRB.MAXIMIZE)

#     elif problem == "FREE":
#         m.setObjectiveN(-sum(sum(d[idx_l, t]
#                                  for idx_l, l in enumerate(links)
#                                 )
#                              for t in syst.T
#                             ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")

#     m.setParam("MIPGap", 1e-3)          # force very tight gap
#     m.setParam("MIPGapAbs", 1e-3)
#     m.setParam("FeasibilityTol", 1e-3)
#     m.setParam("IntFeasTol", 1e-3)
#     m.setParam("OptimalityTol", 1e-3)

#     m.Params.Presolve = 2           # Aggressive presolve
#     m.Params.Method = 2             # Barrier for root relaxation
#     m.Params.Cuts = 2               # More cutting planes
#     m.Params.Heuristics = 0.25      # Reasonable level of primal heuristics
#     m.Params.MIPFocus = 1           # Emphasize feasible solutions first
#     m.Params.NumericFocus = 1       # Reduce rounding issues due to tiny epsilons
#     m.Params.Threads = 0            # Use all cores

#     m.Params.NodefileStart = 0.5    # Start writing node files to disk to prevent RAM overflow
#     m.Params.NoRelHeurTime = 20     # Time for heuristics before branching
#     m.Params.ConcurrentMIP = 1      # Use multiple concurrent solvers for first feasible solution

#     ## Constraints
#     m.addConstrs(
#         (
#             r[idx_d, t] <= k[idx_l, t] + 1e8 * KEY_RATE_SCALE * (1 - z[idx_l, idx_d, t])
#             for idx_l, l in enumerate(links)
#             for idx_d, d in enumerate(demands)
#             for t        in syst.T
#         ), name="end_to_end_key_rate"
#     )

#     # Flow conservation
#     m.addConstrs(
#         (
#             sum(z[idx_l, idx_d, t]
#                 for idx_l, l in enumerate(links)
#                 if isinstance(l.n1, gs)
#                 if gss.index(l.n1) == gss.index(d.n1)
#                ) - sum(z[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if isinstance(l.n2, gs)
#                        if gss.index(l.n2) == gss.index(d.n1)
#                       ) == 1
#             for idx_d, d in enumerate(demands)
#             for t        in syst.T
#         ), name="flow_conservation_1"
#     )
#     m.addConstrs(
#         (
#             sum(z[idx_l, idx_d, t]
#                 for idx_l, l in enumerate(links)
#                 if isinstance(l.n2, gs)
#                 if gss.index(l.n2) == gss.index(d.n2)
#                ) - sum(z[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if isinstance(l.n1, gs)
#                        if gss.index(l.n1) == gss.index(d.n2)
#                       ) == 1
#             for idx_d, d in enumerate(demands)
#             for t        in syst.T
#         ), name="flow_conservation_2"
#     )
#     m.addConstrs(
#         (
#             sum(z[idx_l, idx_d, t]
#                 for idx_l, l in enumerate(links)
#                 if  l.n1 == n
#                ) - sum(z[idx_l, idx_d, t]
#                        for idx_l, l in enumerate(links)
#                        if  l.n2 == n
#                       ) == 0
#             for idx_d, d in enumerate(demands)
#             for n in gss + haps
#             if  n != d.n1 and n != d.n2
#             for t        in syst.T
#         ), name="flow_conservation_3"
#     )
#     m.addConstrs(
#         (
#             z[idx_l_1, idx_d, t] + z[idx_l_2, idx_d, t] <= 1
#             for idx_l_1, l_1 in enumerate(links)
#             for idx_l_2, l_2 in enumerate(links)
#             if  idx_l_1 < idx_l_2 and (l_1.n1 == l_2.n2) and (l_1.n2 == l_2.n1)
#             for idx_d, d     in enumerate(demands)
#             for t            in syst.T
#         ), name="loop_prevention"
#     )
    
#     # Maximum Tx/Rx Connection
#     m.addConstrs(
#         (
#             sum(z[idx_l, idx_d, t]
#                 for idx_l, l in enumerate(links)
#                 if l.n1 == n
#                ) <= n.N_TX
#             for idx_n, n in enumerate(gss + haps)
#             for idx_d, d in enumerate(demands)
#             for t        in syst.T
#         ), name="max_tx_connections"
#     )
    
#     m.addConstrs(
#         (
#             sum(z[idx_l, idx_d, t]
#                 for idx_l, l in enumerate(links)
#                 if l.n2 == n
#                ) <= n.N_RX
#             for idx_n, n in enumerate(gss + haps)
#             for idx_d, d in enumerate(demands)
#             for t        in syst.T
#         ), name="max_rx_connections"
#     )

#     # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
#     for idx_l, l in enumerate(links):
#         # identify which endpoint is HAP and which is GS
#         if isinstance(l.n1, hap) and isinstance(l.n2, gs):
#             hap_idx, gs_node = haps.index(l.n1), l.n2
#         elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
#             hap_idx, gs_node = haps.index(l.n2), l.n1

#         # precompute GS coordinates in same (x,y,z) frame and units (km)
#         [cg1, cg2] = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)  # ensure these are constants in km

#         if problem == "FIXED" or problem == "INIT":
#             # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
#             dx = c1[hap_idx] - cg1*COORDINATE_SCALE
#             dy = c2[hap_idx] - cg2*COORDINATE_SCALE
#             m.addQConstr(d[idx_l]*d[idx_l] >= dx*dx + dy*dy + 15*15*COORDINATE_SCALE*COORDINATE_SCALE,
#                          name=f"dist_cone_{idx_l}")
#         elif problem == "FREE":
#             for t in syst.T:
#                 # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
#                 dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
#                 dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
#                 m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + 15*15*COORDINATE_SCALE*COORDINATE_SCALE,
#                              name=f"dist_cone_{idx_l}_{t}")

#     ## Solve
#     m.optimize()

#     if m.status == GRB.OPTIMAL:
#         print("\n=========== OPTIMAL SOLUTION FOUND ===========")

#         solution = {
#             "r":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r.items()},
#             "z":   {k: v.X for k, v in z.items()},
#             "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
#             "d":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in d.items()},
#             "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
#             "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()}
#         }
        

#         # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
#         # pp.pprint(solution)

#         # print(f"lonlat0:{solution["c1"][0,0]}, {solution["c2"][0,0]}, {xy_to_lonlat(solution["c1"][0,0], solution["c2"][0,0])}")

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
#                 x = solution["c1"].get((idx_hnode))
#                 y = solution["c2"].get((idx_hnode))
#                 if x is not None and y is not None:
#                     # lon, lat = xy_to_lonlat(x, y)
#                     lon, lat = tangent_to_latlon(x, y, 279, 49)
#                     # lon_series.append(lon+360)   # shift if needed
#                     lon_series.append(lon)   # shift if needed
#                     lat_series.append(lat)
#             planned_lons.append(lon_series)
#             planned_lats.append(lat_series)
#         planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]
        
#         # Ground Stations (replicate coordinates across all T so they plot in animation)
#         gs_lons = []
#         gs_lats = []
#         for gnode in gss:
#             gs_lons.append([gnode.lg] * len(syst.T))   # repeat longitude for all time steps
#             gs_lats.append([gnode.la] * len(syst.T))   # repeat latitude for all time steps
#         gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]
        
#         # Combine everything
#         # all_lons = actual_lons + planned_lons + gs_lons
#         # all_lats = actual_lats + planned_lats + gs_lats
#         # all_labels = actual_labels + planned_labels + gs_labels

#         all_lons = planned_lons + gs_lons
#         all_lats = planned_lats + gs_lats
#         all_labels = planned_labels + gs_labels

#         plot_connectivity_graph_planning(gss, haps, links, 
#                         planned_lons=planned_lons,
#                         planned_lats=planned_lats,
#                         planned_labels=planned_labels)
        
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
#     #         gp.quicksum(a[idx_l, tp]
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