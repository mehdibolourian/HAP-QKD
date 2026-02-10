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

        dist_bottle = 0.0
        dist_sum    = 0.0
        num_pairs   = 0
        
        mean_alt = sum(haps[0].H) / len(haps[0].H)
        
        # Precompute HAP mean positions
        hap_positions = []
        for h in haps:
            hap_c1, hap_c2 = latlon_to_tangent(
                sum(h.lg) / len(h.lg),
                sum(h.la) / len(h.la),
                lon0, lat0
            )
            hap_positions.append((hap_c1, hap_c2))

        #print(f"hap_positions: {hap_positions}")
        
        # GS pairs
        for i, gs1 in enumerate(gss):
            gs1_c1, gs1_c2 = latlon_to_tangent(gs1.lg, gs1.la, lon0, lat0)
        
            for j, gs2 in enumerate(gss):
                if j <= i:
                    continue
        
                gs2_c1, gs2_c2 = latlon_to_tangent(gs2.lg, gs2.la, lon0, lat0)
        
                # Best HAP for this pair
                pair_best = math.inf
        
                for (hap_c1, hap_c2) in hap_positions:
        
                    d1 = math.sqrt((gs1_c1 - hap_c1)**2 +
                                    (gs1_c2 - hap_c2)**2 +
                                    mean_alt**2)
        
                    d2 = math.sqrt((gs2_c1 - hap_c1)**2 +
                                    (gs2_c2 - hap_c2)**2 +
                                    mean_alt**2)

                    #print(f"(d1, d2): ({d1}, {d2})")
                    pair_best = min(pair_best, max(d1, d2))
        
                dist_sum    += pair_best
                dist_bottle = max(dist_bottle, pair_best)
                num_pairs  += 1

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

    print(f"dist_bottle: {dist_bottle}")
    print(f"dist_avg: {dist_sum/num_pairs}")
    print(f"compare: {calculate_key_rate_planning('theoretical', 0, 177.378, sum(haps[0].H) / len(haps[0].H))}")

    return (
        solution,
        calculate_key_rate_planning("theoretical", 0, dist_bottle,         sum(haps[0].H) / len(haps[0].H)),
        calculate_key_rate_planning("theoretical", 0, dist_sum / num_pairs, sum(haps[0].H) / len(haps[0].H))
    )

def placement_kmeans(gss, haps, links):
    lon0, lat0 = 279, 49

    solution = {"c1": {}, "c2": {}}

    # -------------------------------
    # Build GS–GS demands
    # -------------------------------
    demands     = []
    gss_demands = []
    for i, g1 in enumerate(gss):
        for j, g2 in enumerate(gss):
            if j > i:
                demands.append(demand(100, g1, g2))
                if g1 not in gss_demands:
                    gss_demands.append(g1)
                if g2 not in gss_demands:
                    gss_demands.append(g2)

    # -------------------------------
    # Reference HAP trajectory
    # -------------------------------
    c1_list_ref, c2_list_ref = [], []
    for lon, lat in zip(haps[0].lg, haps[0].la):
        c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
        c1_list_ref.append(c1)
        c2_list_ref.append(c2)

    c1_ref_mean = sum(c1_list_ref) / len(c1_list_ref)
    c2_ref_mean = sum(c2_list_ref) / len(c2_list_ref)

    # -------------------------------
    # GS coordinates
    # -------------------------------
    gs_coordinates = []
    for g in gss_demands:
        c1, c2 = latlon_to_tangent(g.lg, g.la, lon0, lat0)
        gs_coordinates.append([c1, c2])

    num_gs = len(gs_coordinates)

    # -------------------------------
    # Choose k safely
    # -------------------------------
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

    k = min(k, num_gs)          # CRITICAL FIX

    # -------------------------------
    # KMeans (safe for small GS count)
    # -------------------------------
    if num_gs == 1:
        centers = np.array(gs_coordinates)
    else:
        kmeans  = KMeans(n_clusters=k, random_state=0, n_init=10)
        labels  = kmeans.fit_predict(gs_coordinates)
        centers = kmeans.cluster_centers_

    # -------------------------------
    # Ordered center pairs (farthest first)
    # -------------------------------
    center_pairs = []
    for (i, c1), (j, c2) in combinations(enumerate(centers), 2):
        dist = np.linalg.norm(c1 - c2)
        center_pairs.append((dist, i, j))

    if not center_pairs:
        center_pairs = [(0.0, 0, 0)]

    center_pairs.sort(reverse=True, key=lambda x: x[0])

    # -------------------------------
    # Assign HAP trajectories (with reuse)
    # -------------------------------
    hap_idx   = 0
    num_pairs = len(center_pairs)

    while hap_idx < len(haps):
        dist, idx_gc1, idx_gc2 = center_pairs[hap_idx % num_pairs]

        mid_c1 = (centers[idx_gc1][0] + centers[idx_gc2][0]) / 2
        mid_c2 = (centers[idx_gc1][1] + centers[idx_gc2][1]) / 2

        shift_c1 = mid_c1 - c1_ref_mean
        shift_c2 = mid_c2 - c2_ref_mean

        for t in syst.T:
            solution["c1"][(hap_idx, t)] = c1_list_ref[t] + shift_c1
            solution["c2"][(hap_idx, t)] = c2_list_ref[t] + shift_c2

            haps[hap_idx].lg[t], haps[hap_idx].la[t] = tangent_to_latlon(
                solution["c1"][(hap_idx, t)],
                solution["c2"][(hap_idx, t)],
                lon0, lat0
            )

        hap_idx += 1

    # -------------------------------
    # Distance metrics
    # -------------------------------
    dist_bottle = 0.0
    dist_sum    = 0.0
    num_pairs   = 0
    
    mean_alt = sum(haps[0].H) / len(haps[0].H)
    
    # Precompute HAP mean positions
    hap_positions = []
    for h in haps:
        hap_c1, hap_c2 = latlon_to_tangent(
            sum(h.lg) / len(h.lg),
            sum(h.la) / len(h.la),
            lon0, lat0
        )
        hap_positions.append((hap_c1, hap_c2))
    
    # GS pairs
    for i, gs1 in enumerate(gss):
        gs1_c1, gs1_c2 = latlon_to_tangent(gs1.lg, gs1.la, lon0, lat0)
    
        for j, gs2 in enumerate(gss):
            if j <= i:
                continue
    
            gs2_c1, gs2_c2 = latlon_to_tangent(gs2.lg, gs2.la, lon0, lat0)
    
            # Best HAP for this pair
            pair_best = math.inf
    
            for (hap_c1, hap_c2) in hap_positions:
    
                d1 = math.sqrt((gs1_c1 - hap_c1)**2 +
                                (gs1_c2 - hap_c2)**2 +
                                mean_alt**2)
    
                d2 = math.sqrt((gs2_c1 - hap_c1)**2 +
                                (gs2_c2 - hap_c2)**2 +
                                mean_alt**2)
    
                pair_best = min(pair_best, max(d1, d2))
    
            dist_sum    += pair_best
            dist_bottle = max(dist_bottle, pair_best)
            num_pairs  += 1

    # -------------------------------
    # Planned trajectories (for plots)
    # -------------------------------
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

    return (
        solution,
        calculate_key_rate_planning("theoretical", 0, dist_bottle,         sum(haps[0].H) / len(haps[0].H)),
        calculate_key_rate_planning("theoretical", 0, dist_sum / num_pairs, sum(haps[0].H) / len(haps[0].H))
    )

def placement_path_lower(gss, haps, links, prob):
    """
    Optimize HAP placement for maximum QKD key rates.

    gss: list of GS nodes
    haps: list of HAP nodes
    links: list of link objects
    prob: "btl" for bottleneck, "avg" for average case
    """
    c1_list_ref, c2_list_ref, c3_list_ref = [], [], []
    lon0, lat0 = 279, 49  # reference for tangent conversion

    # Create demands (GS-GS pairs)
    demands = []
    for idx_g1, g1 in enumerate(gss):
        for idx_g2, g2 in enumerate(gss):
            if idx_g2 > idx_g1:
                demands.append(demand(100, g1, g2))

    # Generate paths (GS -> HAP -> GS)
    paths = []
    for gs_src in gss:
        for gs_dst in gss:
            if gs_src == gs_dst:
                continue
            for l1, l2 in product(links, links):
                if (
                    isinstance(l1.n1, gs) and isinstance(l1.n2, hap)
                    and isinstance(l2.n1, hap) and isinstance(l2.n2, gs)
                    and l1.n1 == gs_src and l2.n2 == gs_dst
                    and l1.n2 == l2.n1
                ):
                    paths.append(path(l1, l2))

    # -----------------------------------------
    # Optimization Model
    # -----------------------------------------
    m = gp.Model("hap-qkd")

    # Reference coordinates
    for lon, lat in zip(haps[0].lg, haps[0].la):
        c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
        c1_list_ref.append(c1)
        c2_list_ref.append(c2)
    c3_list_ref = haps[0].H

    # Decision variables
    z = {}
    for idx_p, p in enumerate(paths):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                z[idx_p, idx_d, t] = m.addVar(
                    vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name=f"z_{idx_p}_{idx_d}_{t}"
                )

    di, dm = {}, {}
    for idx_l, l in enumerate(links):
        for t in syst.T:
            di[idx_l, t] = m.addVar(
                vtype=GRB.CONTINUOUS,
                lb=15 * COORDINATE_SCALE,
                ub=2e3 * COORDINATE_SCALE,
                name=f"di_{idx_l}_{t}"
            )

    if prob == "btl":
        for t in syst.T:
            dm[t] = m.addVar(
                vtype=GRB.CONTINUOUS, lb=15.0, ub=2e3 * COORDINATE_SCALE, name=f"dm_{t}"
            )
    elif prob == "avg":
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                dm[idx_d, t] = m.addVar(
                    vtype=GRB.CONTINUOUS, lb=15.0, ub=2e3 * COORDINATE_SCALE, name=f"dm_{idx_d}_{t}"
                )

    # HAP coordinates
    c1, c2 = {}, {}
    for idx_h, hnode in enumerate(haps):
        for t in syst.T:
            c1[idx_h, t] = m.addVar(lb=-1e3*COORDINATE_SCALE, ub=1e3*COORDINATE_SCALE,
                                    vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
            c2[idx_h, t] = m.addVar(lb=-1e3*COORDINATE_SCALE, ub=1e3*COORDINATE_SCALE,
                                    vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")

    m.ModelSense = GRB.MINIMIZE

    # Objective
    if prob == "btl":
        m.setObjectiveN(gp.quicksum(dm[t] for t in syst.T)/len(syst.T),
                        index=0, priority=3, weight=1.0)
    elif prob == "avg":
        m.setObjectiveN(
            gp.quicksum(gp.quicksum(dm[idx_d, t] for idx_d, d in enumerate(demands))
                        for t in syst.T) / len(demands) / len(syst.T),
            index=0, priority=3, weight=1.0
        )

    # Singular path selection
    m.addConstrs(
        (gp.quicksum(z[idx_p, idx_d, t] for idx_p, p in enumerate(paths)
                     if p.l1.n1 == d.n1 and p.l2.n2 == d.n2) == 1
         for idx_d, d in enumerate(demands)
         for t in syst.T),
        name="singular_path_selection"
    )

    # Max distance constraints
    if prob == "btl":
        m.addConstrs(
            (dm[t] >= di[links.index(p.l1), t] - 1e3*(1 - z[idx_p, idx_d, t])
             for idx_d, d in enumerate(demands)
             for idx_p, p in enumerate(paths)
             for t in syst.T),
            name="max_distance_1"
        )
        m.addConstrs(
            (dm[t] >= di[links.index(p.l2), t] - 1e3*(1 - z[idx_p, idx_d, t])
             for idx_d, d in enumerate(demands)
             for idx_p, p in enumerate(paths)
             for t in syst.T),
            name="max_distance_2"
        )
    elif prob == "avg":
        m.addConstrs(
            (dm[idx_d, t] >= di[links.index(p.l1), t] - 1e3*(1 - z[idx_p, idx_d, t])
             for idx_d, d in enumerate(demands)
             for idx_p, p in enumerate(paths)
             for t in syst.T),
            name="max_distance_1"
        )
        m.addConstrs(
            (dm[idx_d, t] >= di[links.index(p.l2), t] - 1e3*(1 - z[idx_p, idx_d, t])
             for idx_d, d in enumerate(demands)
             for idx_p, p in enumerate(paths)
             for t in syst.T),
            name="max_distance_2"
        )

    # Trajectory shift constraints
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

    # Distance cone constraints
    for idx_l, l in enumerate(links):
        if isinstance(l.n1, hap) and isinstance(l.n2, gs):
            hap_idx, gs_node = haps.index(l.n1), l.n2
        elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
            hap_idx, gs_node = haps.index(l.n2), l.n1
        else:
            continue

        cg1, cg2 = latlon_to_tangent(gs_node.lg, gs_node.la, lon0, lat0)
        for t in syst.T:
            dx = c1[hap_idx, t] - cg1*COORDINATE_SCALE
            dy = c2[hap_idx, t] - cg2*COORDINATE_SCALE
            m.addQConstr(
                di[idx_l, t]*di[idx_l, t] >= dx*dx + dy*dy + haps[hap_idx].H[t]**2*COORDINATE_SCALE**2,
                name=f"dist_cone_{idx_l}_{t}"
            )

    # -----------------------------------------
    # Solve
    # -----------------------------------------
    m.optimize()

    if m.status not in (GRB.OPTIMAL, GRB.TIME_LIMIT) or m.SolCount == 0:
        return None, None, None

    # Extract solution
    solution = {
        "z": {k: v.X for k, v in z.items()},
        "di": {k: v.X for k, v in di.items()},
        "dm": {k: v.X for k, v in dm.items()},
        "c1": {k: v.X for k, v in c1.items()},
        "c2": {k: v.X for k, v in c2.items()},
    }

    # -----------------------------------------
    # Deterministic rounding of z
    # -----------------------------------------
    z_rounded = {}
    for idx_d, d in enumerate(demands):
        for t in syst.T:
            max_val = -1
            selected_path = None
            for idx_p, p in enumerate(paths):
                val = solution["z"].get((idx_p, idx_d, t), 0)
                if val > max_val:
                    max_val = val
                    selected_path = idx_p
            for idx_p, p in enumerate(paths):
                z_rounded[idx_p, idx_d, t] = 1.0 if idx_p == selected_path else 0.0

    #print("z_rounded (only 1s):", [key for key, val in z_rounded.items() if val == 1.0])

    # -----------------------------------------
    # Recompute bottleneck and average distances
    # -----------------------------------------
    dist_bottle = 0.0
    dist_sum    = 0.0
    num_pairs   = 0
    mean_alt    = sum(haps[0].H)/len(haps[0].H)

    # Precompute HAP mean positions
    hap_positions = []
    for idx_h, h in enumerate(haps):
        assigned_c1 = []
        assigned_c2 = []
        for t in syst.T:
            for idx_p, p in enumerate(paths):
                for idx_d, d in enumerate(demands):
                    if z_rounded[idx_p, idx_d, t] > 0 and (p.l1.n2 == h or p.l2.n1 == h):
                        c1_val = solution["c1"].get((idx_h, t))
                        c2_val = solution["c2"].get((idx_h, t))
                        if c1_val is not None: assigned_c1.append(c1_val)
                        if c2_val is not None: assigned_c2.append(c2_val)
        hap_positions.append((
            sum(assigned_c1)/len(assigned_c1) if assigned_c1 else sum(h.lg)/len(h.lg),
            sum(assigned_c2)/len(assigned_c2) if assigned_c2 else sum(h.la)/len(h.la)
        ))

    # GS pairs distances
    for i, gs1 in enumerate(gss):
        gs1_c1, gs1_c2 = latlon_to_tangent(gs1.lg, gs1.la, lon0, lat0)
        for j, gs2 in enumerate(gss):
            if j <= i:
                continue
            gs2_c1, gs2_c2 = latlon_to_tangent(gs2.lg, gs2.la, lon0, lat0)
            pair_best = math.inf
            #print(f"hap_positions: {hap_positions}")
            for hap_c1, hap_c2 in hap_positions:
                d1 = math.sqrt((gs1_c1 - hap_c1)**2 + (gs1_c2 - hap_c2)**2 + mean_alt**2)
                d2 = math.sqrt((gs2_c1 - hap_c1)**2 + (gs2_c2 - hap_c2)**2 + mean_alt**2)
                #print(f"(d1, d2): ({d1}, {d2})")
                pair_best = min(pair_best, max(d1, d2))
            dist_sum += pair_best
            dist_bottle = max(dist_bottle, pair_best)
            num_pairs += 1

    dist_avg = dist_sum / num_pairs if num_pairs > 0 else float("inf")

    print(f"dist_bottle: {dist_bottle}, dist_avg: {dist_avg}")
    print(f"solution['c1']: {solution['c1'].get((0, 0))}")
    print(f"solution['c2']: {solution['c2'].get((0, 0))}")
    print(f"{tangent_to_latlon(solution['c1'].get((0, 0)), solution['c2'].get((0, 0)), lon0, lat0)}")

    # -----------------------------------------
    # Return rounded solution + key rates
    # -----------------------------------------
    return (
        z_rounded,
        calculate_key_rate_planning("theoretical", 0, dist_bottle, mean_alt),
        calculate_key_rate_planning("theoretical", 0, dist_avg, mean_alt)
    )


## Find the optimal placement of the HAPs to reach the maximum key generation for all the end-to-end paths between GS pairs.
def placement_path(gss, haps, links, prob):
    c1_list_ref, c2_list_ref, c3_list_ref = [], [], []
    lon0, lat0 = 279, 49

    alpha = 0.1

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

    paths = []
    for gs_src in gss:
        for gs_dst in gss:
            if gs_src == gs_dst:
                continue

            for l1, l2 in product(links, links):
                if (
                    isinstance(l1.n1, gs)
                    and isinstance(l1.n2, hap)
                    and isinstance(l2.n1, hap)
                    and isinstance(l2.n2, gs)
                    and l1.n1 == gs_src
                    and l2.n2 == gs_dst
                    and l1.n2 == l2.n1
                ):
                    paths.append(path(l1, l2))
                    #print(f"l1: {links.index(l1)}, l2: {links.index(l2)}")

    # Create Optimization Model
    m = gp.Model("hap-qkd")

    for lon, lat in zip(haps[0].lg, haps[0].la):
        c1, c2 = latlon_to_tangent(lon, lat, lon0, lat0)
        c1_list_ref.append(c1)
        c2_list_ref.append(c2)
    c3_list_ref = haps[0].H

    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    z = {}
    for idx_p, p in enumerate(paths):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                z[idx_p, idx_d, t] = m.addVar(
                    name=f"z_{idx_p}_{idx_d}_{t}",
                    vtype=GRB.BINARY,
                    lb=0.0,
                    ub=1.0
                )

    nodes = gss + haps

    # di (distance)
    di, dm = {}, {}
    for idx_l, l in enumerate(links):
        for t in syst.T:
            di[idx_l, t] = m.addVar(
                name=f"di_{idx_l}_{t}",
                vtype=GRB.CONTINUOUS,
                lb=15 * COORDINATE_SCALE,
                ub=2e3 * COORDINATE_SCALE
            )  # LoS distance (Min height in strat.)

    if prob == "btl":
        for t in syst.T:
            dm[t] = m.addVar(
                name=f"dm_{t}",
                vtype=GRB.CONTINUOUS,
                lb=15.0,
                ub=2e3 * COORDINATE_SCALE
            )
    elif prob == "avg":
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                dm[idx_d, t] = m.addVar(
                    name=f"dm_{idx_d}_{t}",
                    vtype=GRB.CONTINUOUS,
                    lb=15.0,
                    ub=2e3 * COORDINATE_SCALE
                )

    # Coordinate decision variables for each HAP and time
    c1, c2 = {}, {}  # x, y in km
    for idx_h, hnode in enumerate(haps):
        for t in syst.T:
            c1[idx_h, t] = m.addVar(
                lb=-1e3 * COORDINATE_SCALE,
                ub=1e3 * COORDINATE_SCALE,
                vtype=GRB.CONTINUOUS,
                name=f"c1_{idx_h}_{t}"
            )
            c2[idx_h, t] = m.addVar(
                lb=-1e3 * COORDINATE_SCALE,
                ub=1e3 * COORDINATE_SCALE,
                vtype=GRB.CONTINUOUS,
                name=f"c2_{idx_h}_{t}"
            )

    # Objective
    # if prob == "btl":
    #     m.setObjective(
    #         gp.quicksum(dm[t] for t in syst.T), GRB.MINIMIZE
    #     )
    # elif prob == "avg":
    #     m.setObjective(
    #         gp.quicksum(gp.quicksum(dm[idx_d, t] for idx_d, d in enumerate(demands)) for t in syst.T), GRB.MINIMIZE
    #     )

    # m.setParam("MIPGap", 1e-4)
    # m.setParam("MIPGapAbs", 1e-4)
    # m.setParam("FeasibilityTol", 1e-4)
    # m.setParam("IntFeasTol", 1e-4)
    # m.setParam("OptimalityTol", 1e-4)
    # m.setParam("TimeLimit", 60)  # seconds
    # m.Params.Presolve = 2
    # m.Params.Method = 2
    # m.Params.Cuts = 2
    # m.Params.Heuristics = 0.5
    # m.Params.MIPFocus = 1
    # m.Params.NumericFocus = 1
    # m.Params.Threads = 0
    # m.Params.NodefileStart = 0.5
    # m.Params.NoRelHeurTime = 60
    # m.Params.ConcurrentMIP = 1

    # m.setObjective(
    #     gp.quicksum(dm[t] for t in syst.T) + alpha * gp.quicksum(gp.quicksum(di[idx_l, t]
    #                                                                          for idx_l, l in enumerate(links)
    #                                                                         )
    #                                                              for t in syst.T),
    #     GRB.MINIMIZE
    # )

    m.ModelSense = GRB.MINIMIZE

    # Primary objective: maximize demand satisfaction
    if prob == "btl":
        m.setObjectiveN(gp.quicksum(dm[t] for t in syst.T) / len(syst.T), index=0, priority=3, weight=1.0, abstol=1e-5, reltol=1e-5, name="Primary")
    elif prob == "avg":
        m.setObjectiveN(gp.quicksum(gp.quicksum(dm[idx_d, t] for idx_d, d in enumerate(demands)) for t in syst.T) / len(demands) / len(syst.T), index=0, priority=3, weight=1.0, abstol=1e-5, reltol=1e-5, name="Primary")
    
    # Secondary objective: maximize keys served
    m.setObjectiveN(gp.quicksum(gp.quicksum(di[idx_l, t]
                                            for idx_l, l in enumerate(links)
                                           )
                                for t in syst.T
                               ), index=1, priority=2, weight=1, abstol=1e-5, reltol=1e-5, name="Secondary")

    ########### Exclusive constraints ###########

    m.addConstrs(
        (
            gp.quicksum(
                z[idx_p, idx_d, t]
                for idx_p, p in enumerate(paths)
                if p.l1.n1 == d.n1 and p.l2.n2 == d.n2
            ) == 1
            for idx_d, d in enumerate(demands)
            for t in syst.T
        ),
        name="singular_path_selection"
    )

    if prob == "btl":
        m.addConstrs(
            (
                dm[t] >= di[links.index(p.l1), t]
                - 1e3 * (1 - z[idx_p, idx_d, t])
                for idx_d, d in enumerate(demands)
                for idx_p, p in enumerate(paths)
                for t in syst.T
            ),
            name="max_distance_1"
        )
    elif prob == "avg":
        m.addConstrs(
            (
                dm[idx_d, t] >= di[links.index(p.l1), t]
                - 1e3 * (1 - z[idx_p, idx_d, t])
                for idx_d, d in enumerate(demands)
                for idx_p, p in enumerate(paths)
                for t in syst.T
            ),
            name="max_distance_1"
        )

    if prob == "btl":
        m.addConstrs(
            (
                dm[t] >= di[links.index(p.l2), t]
                - 1e3 * (1 - z[idx_p, idx_d, t])
                for idx_d, d in enumerate(demands)
                for idx_p, p in enumerate(paths)
                for t in syst.T
            ),
            name="max_distance_2"
        )
    elif prob == "avg":
        m.addConstrs(
            (
                dm[idx_d, t] >= di[links.index(p.l2), t]
                - 1e3 * (1 - z[idx_p, idx_d, t])
                for idx_d, d in enumerate(demands)
                for idx_p, p in enumerate(paths)
                for t in syst.T
            ),
            name="max_distance_2"
        )

    m.addConstrs(
        (
            c1[idx_h, t]
            == c1_list_ref[t]
            + (c1[idx_h, 0] - c1_list_ref[0])
            for idx_h, h in enumerate(haps)
            for t in syst.T
            if t >= 1
        ),
        name="shift_trajectory_1"
    )

    m.addConstrs(
        (
            c2[idx_h, t]
            == c2_list_ref[t]
            + (c2[idx_h, 0] - c2_list_ref[0])
            for idx_h, h in enumerate(haps)
            for t in syst.T
            if t >= 1
        ),
        name="shift_trajectory_2"
    )

    # Distance cone constraints
    for idx_l, l in enumerate(links):
        if isinstance(l.n1, hap) and isinstance(l.n2, gs):
            hap_idx, gs_node = haps.index(l.n1), l.n2
        elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
            hap_idx, gs_node = haps.index(l.n2), l.n1
        else:
            continue

        cg1, cg2 = latlon_to_tangent(gs_node.lg, gs_node.la, 279, 49)

        for t in syst.T:
            dx = c1[hap_idx, t] - cg1 * COORDINATE_SCALE
            dy = c2[hap_idx, t] - cg2 * COORDINATE_SCALE
            m.addQConstr(
                di[idx_l, t] * di[idx_l, t]
                >= dx * dx
                + dy * dy
                + haps[hap_idx].H[t]
                * haps[hap_idx].H[t]
                * COORDINATE_SCALE
                * COORDINATE_SCALE,
                name=f"dist_cone_{idx_l}_{t}"
            )

    ## Solve
    m.optimize()

    if m.status in (GRB.OPTIMAL, GRB.TIME_LIMIT) and m.SolCount > 0:
        print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        print(f"Optimal value: {m.ObjVal}")

        solution = {
            "z": {k: v.X for k, v in z.items()},
            "di": {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in di.items()},
            "dm": {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in dm.items()},
            "c1": {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
            "c2": {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()},
        }

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
                    lon_series.append(lon) # shift if needed
                    lat_series.append(lat)
                    planned_lons.append(lon_series)
                    planned_lats.append(lat_series)
                    planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]
                    
                    # Ground Stations (replicate coordinates across all T so they plot in animation)
                    gs_lons = []
                    gs_lats = []
                    
                    for gnode in gss:
                        gs_lons.append([gnode.lg] * len(syst.T)) # repeat longitude for all time steps
                        gs_lats.append([gnode.la] * len(syst.T)) # repeat latitude for all time steps
                        gs_labels = [f"GS_{idx_gs}" for idx_gs, _ in enumerate(gss)]
                        all_lons = planned_lons + gs_lons
                        all_lats = planned_lats + gs_lats
                        all_labels = planned_labels + gs_labels
                        
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
                                    lon_series.append(lon) # shift if needed
                                    lat_series.append(lat)
                                    planned_lons.append(lon_series)
                                    planned_lats.append(lat_series)
                                    planned_labels = [f"HAP_{idx_hnode}*" for idx_hnode in range(len(haps))]

        plot_connectivity_graph_planning(gss, haps, links, planned_lons=planned_lons, planned_lats=planned_lats, planned_labels=planned_labels, alg="opt")
        plot_connectivity_graph_planning_3d(gss, haps, links, planned_lons=planned_lons, planned_lats=planned_lats, planned_alts=haps[0].H, planned_labels=planned_labels, alg="opt")

    else:
        solution = None

    for idx_l, l in enumerate(links):
        for t in syst.T:
            if isinstance(l.n1, gs):
                #print(f"di: {idx_l}, {solution['di'].get((idx_l,t))}")
                gc1, gc2 = latlon_to_tangent(l.n1.lg, l.n1.la, lon0, lat0)
                hc1, hc2 = solution["c1"].get((haps.index(l.n2),t)), solution["c2"].get((haps.index(l.n2),t))

                distance = math.sqrt((gc1-hc1)**2 + (gc2-hc2)**2 + (l.n2.H[t])**2)
                #print(f"calc di: {idx_l}, {distance}")

    dist_bottle = 0
    dist_sum = 0
    for g in gss:
        dist_min = math.inf
        for h in haps:
            gs_c1, gs_c2 = latlon_to_tangent(g.lg, g.la, lon0, lat0)
            hap_c1, hap_c2 = latlon_to_tangent(
                sum(h.lg) / len(h.lg),
                sum(h.la) / len(h.la),
                lon0,
                lat0,
            )
            dist_min = min(
                dist_min,
                math.sqrt((gs_c1 - hap_c1) ** 2 + (gs_c2 - hap_c2) ** 2 + (sum(h.H) / len(h.H))**2),
            )

        dist_sum += dist_min
        dist_bottle = max(dist_bottle, dist_min)

    if prob == "btl":
        dist_avg_obj    = sum(sum(min(max(solution["di"].get((links.index(p.l1), t)), solution["di"].get((links.index(p.l2), t)))
                                      for idx_p, p in enumerate(paths)
                                      if p.l1.n1 == d.n1 and p.l2.n2 == d.n2
                                     )
                                  for idx_d, d in enumerate(demands)
                                 )
                              for t in syst.T
                             ) / len(demands) / len(syst.T)
        dist_bottle_obj = m.ObjVal
    elif prob == "avg":
        dist_avg_obj    = m.ObjVal
        dist_bottle_obj = sum(max(solution["dm"].get((idx_d, t))
                                  for idx_d, d in enumerate(demands)
                                 )
                              for t in syst.T
                             ) / len(syst.T)
    

    # print(f"dist_bottle: {dist_bottle}")
    # print(f"solution: {solution['dm']}")
    # print(f"dist_avg: {dist_sum / len(gss)}")

    return (
        solution,
        calculate_key_rate_planning("theoretical", 0, dist_bottle,         sum(haps[0].H) / len(haps[0].H)),
        calculate_key_rate_planning("theoretical", 0, dist_sum / len(gss), sum(haps[0].H) / len(haps[0].H))
    )

    
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
    # demands.append(
    #     demand(
    #         100,
    #         gss[0],
    #         gss[1]
    #     )
    # )
    # demands.append(
    #     demand(
    #         100,
    #         gss[0],
    #         gss[2]
    #     )
    # )

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

    dpts = np.linspace(15 * COORDINATE_SCALE, 2e3 * COORDINATE_SCALE, 10)
    kpts = [calculate_key_rate_planning("theoretical", 0, d/COORDINATE_SCALE, t) * KEY_RATE_SCALE for d in dpts]

    # k (key rate) and d (distance)
    k, di, ke, de = {}, {}, {}, {}
    for idx_l, l in enumerate(links):
        for t in syst.T:
            di[idx_l, t] = m.addVar(name=f"di_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15 * COORDINATE_SCALE, ub=2e3 * COORDINATE_SCALE) # LoS distance (Min height in strat.)
                
            k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

            #print(f"kpts: {kpts}")
            #m.addGenConstrPWL(di[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

            slope = (kpts[9] - kpts[0]) / (dpts[9] - dpts[0])
            print(f"slope: {slope}")

    plt.plot(dpts, kpts)
    plt.show()

    # for idx_d, d in enumerate(demands):
    #     for t in syst.T:
    #         ke[idx_d, t] = m.addVar(name=f"ke_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

    # for t in syst.T:
    #     ke[t] = m.addVar(name=f"ke_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
    for t in syst.T:
        de[t] = m.addVar(name=f"de_{t}", vtype=GRB.CONTINUOUS, lb=15.0, ub=2e3)

    # Coordinate decision variables for each HAP and time
    c1, c2 = {}, {}  # x, y in km
    for idx_h, hnode in enumerate(haps):
        for t in syst.T:
            c1[idx_h, t] = m.addVar(lb=-8e2*COORDINATE_SCALE, ub=8e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
            c2[idx_h, t] = m.addVar(lb=-8e2*COORDINATE_SCALE, ub=8e2*COORDINATE_SCALE, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")

    # # Secondary objective: maximize d
    m.setObjective(gp.quicksum(de[t]
                               for t in syst.T
                              )
                   , GRB.MINIMIZE)

    # m.setObjective(gp.quicksum(ke[t]
    #                            for t in syst.T
    #                           )
    #                , GRB.MAXIMIZE)

    # m.setObjective(gp.quicksum(gp.quicksum(ke[idx_d, t]
    #                                        for t in syst.T
    #                                       )
    #                            for idx_d, d in enumerate(demands)
    #                           )
    #                , GRB.MAXIMIZE)

    # m.setObjective(gp.quicksum(gp.quicksum(di[idx_l, t]
    #                                        for t in syst.T
    #                                       )
    #                            for idx_l, l in enumerate(links)
    #                           )
    #                , GRB.MINIMIZE)

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

    # m.setObjective(gp.quicksum(ke[t]
    #                            for t in syst.T
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

    m.setParam("FeasibilityTol", 1e-6)
    m.setParam("IntFeasTol", 1e-6)
    m.setParam("OptimalityTol", 1e-6)


    # m.setParam("MIPGap", 1e-4)
    # m.setParam("MIPGapAbs", 1e-4)
    # m.setParam("FeasibilityTol", 1e-4)
    # m.setParam("IntFeasTol", 1e-4)
    # m.setParam("OptimalityTol", 1e-4)
    # m.setParam("TimeLimit", 60)  # seconds
    # m.Params.Presolve = 2
    # m.Params.Method = 2
    # m.Params.Cuts = 2
    # m.Params.Heuristics = 0.5
    # m.Params.MIPFocus = 1
    # m.Params.NumericFocus = 1
    # m.Params.Threads = 0
    # m.Params.NodefileStart = 0.5
    # m.Params.NoRelHeurTime = 60
    # m.Params.ConcurrentMIP = 1

    # # One-hop paths
    # m.addConstrs(
    #     (k[idx_l, t] == slope * (di[idx_l, t] - 15) + kpts[0]
    #      for idx_l, l in enumerate(links)
    #      for t in syst.T),
    #     name="keyrate_distance"
    # )

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

    # m.addConstrs(
    #     (
    #         ke[idx_d, t] <= k[idx_l, t] + kpts[0] * (1 - z[idx_l, idx_d, t])
    #         for idx_d, d in enumerate(demands)
    #         for idx_l, l in enumerate(links)
    #         for t in syst.T
    #     ), name="keyrate_active_link_3"
    # )
    # m.addConstrs(
    #     (
    #         ke[t] <= k[idx_l, t] + kpts[0] * (1 - z[idx_l, idx_d, t])
    #         for idx_d, d in enumerate(demands)
    #         for idx_l, l in enumerate(links)
    #         for t in syst.T
    #     ), name="keyrate_active_link_3"
    # )
    m.addConstrs(
        (
            de[t] >= di[idx_l, t] - (1 - z[idx_l, idx_d, t])
            for idx_d, d in enumerate(demands)
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
            # m.addConstr(di[idx_l, t] == dx + dy + haps[hap_idx].H[t]*COORDINATE_SCALE,
            #             name=f"dist_cone_{idx_l}_{t}")
                         
    ## Solve
    m.optimize()

    if m.status in (GRB.OPTIMAL, GRB.TIME_LIMIT) and m.SolCount > 0:
        print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        # solution = {
        #     "z":   {k: v.X for k, v in z.items()},
        #     #"o":   {k: round(v.X, 3) for k, v in o.items()},
        #     "ke":  {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in ke.items()},
        #     "k":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in k.items()},
        #     "di":   {k: round(v.X / COORDINATE_SCALE, 3) for k, v in di.items()},
        #     "c1":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c1.items()},
        #     "c2":  {k: round(v.X / COORDINATE_SCALE, 3) for k, v in c2.items()}
        # }
        solution = {
            "z":   {k: v.X for k, v in z.items()},
            #"o":   {k: round(v.X, 3) for k, v in o.items()},
            "de":  {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in ke.items()},
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
            #print(coord_text)
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

def calculate_key_rate_planning(method, link, d_los, altitude):
    K_MAX = 0
    if method == "plob":
        L_geo = 20 * max(math.log10((R_TX + d_los * 1000 * THETA) / R_RX), 0)
        L_ma  = 0.01 * d_los
        L_t   = L_geo + L_ma
        
        ETA = 10**(-L_t/10)
        
        K_MAX = -ts.ratesources * math.log2(1 - ETA)
    elif method == "theoretical":
        # Compute efficiencies
        eta_theory = ts.channel_theory(direction="downlink", gs_alt=0, hap_alt=altitude,
                                       distance=d_los, n_correction=6, params_in=ts.params)
        # Compute SKRs
        K_MAX = ts.compute_skr(eta_theory)
    elif method == "simulation":
        # Compute efficiencies
        eta_sim = ts.channel_simulation(direction="downlink", gs_alt=0, hap_alt=altitude,
                                       distance=d_los, n_correction=6, params_in=ts.params)
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











