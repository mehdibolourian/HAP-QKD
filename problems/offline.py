from libraries import *

def plot_z_timeline(z, links, demands, figsize=(9,3)):
    """
    Visualize z[idx_l, idx_d, t] binary results as a timeline,
    showing only links that are active (z=1 for any demand/time).

    Parameters
    ----------
    z : 3D array-like (n_links x n_demands x n_time)
        Binary or 0/1 values.
    links : list
        List of link objects with attributes n1.tag and n2.tag
    demands : list
        List of demand objects with attributes n1.tag and n2.tag
    syst : object
        System object containing syst.T (list or range of time steps)
    figsize : tuple, optional
        Figure size.
    """
    n_links = len(links)
    n_demands = len(demands)
    n_time = len(syst.T)

    # Identify links with at least one z=1
    active_links = [
        idx_l for idx_l in range(n_links)
        if any(z[idx_l, idx_d, t] > 0.5 for idx_d in range(n_demands) for t in syst.T)
    ]
    if not active_links:
        print("No active links (all z = 0).")
        return

    # Assign each demand a color
    cmap = plt.cm.get_cmap("tab10", n_demands)
    demand_colors = [cmap(i) for i in range(n_demands)]

    fig, ax = plt.subplots(figsize=figsize)

    for new_idx, idx_l in enumerate(active_links):
        link = links[idx_l]
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
                        y=new_idx,
                        width=t2 - t1 + 1,
                        left=t1,
                        height=0.6,
                        color=demand_colors[idx_d],
                        alpha=0.7,
                        edgecolor='k'
                    )

    # --- Y-axis: only active links ---
    y_labels = [
        f"{links[idx_l].n1.tag}→{links[idx_l].n2.tag}"
        if hasattr(links[idx_l].n1, "tag") and hasattr(links[idx_l].n2, "tag")
        else f"L{idx_l}"
        for idx_l in active_links
    ]
    ax.set_yticks(np.arange(len(active_links)))
    ax.set_yticklabels(y_labels, fontsize=9)

    # --- X-axis: time steps ---
    ax.set_xticks(syst.T)
    ax.set_xlabel("Time Step", fontsize=11)
    ax.set_ylabel("Active Links", fontsize=11)

    # --- Legend: demands ---
    patches = [
        mpatches.Patch(color=demand_colors[i], label=f"Demand {i}: {demands[i].n1.tag}→{demands[i].n2.tag}")
        for i in range(n_demands)
    ]
    ax.legend(handles=patches, loc="upper right", fontsize=9)

    ax.grid(True, axis='x', alpha=0.3)
    plt.title("Active Link Utilization Timeline by Demand (z variables)", fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.show()


def offline(gss, haps, links, demands, f_qkp, problem, demand_active, link_active):
    # Create Optimization Model
    m = gp.Model("hap-qkd")
    
    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    r_1, r_2, r_h, a, z, u, o = {}, {}, {}, {}, {}, {}, {}

    for idx_l, l in enumerate(links):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                # demand_active is np.ones(MAX_DEMANDS), so we check the idx_d
                is_active = (demand_active[idx_d][t] == 1)
                
                # If inactive, we'll force ub to 0.0 later or set it here
                upper_bound = d.K_REQ[t] * KEY_RATE_SCALE if is_active else 0.0
                
                r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
                r_2[idx_l, idx_d, t] = m.addVar(name=f"r_2_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
                z[idx_l, idx_d, t]   = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}",   vtype=GRB.BINARY)
                
    for idx_d, d in enumerate(demands):
        for t in syst.T:
            # demand_active is np.ones(MAX_DEMANDS), so we check the idx_d
            is_active = (demand_active[idx_d][t] == 1)
            
            # If inactive, we'll force ub to 0.0 later or set it here
            upper_bound = d.K_REQ[t] * KEY_RATE_SCALE if is_active else 0.0
            r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
            
    for idx_l, l in enumerate(links):
        for t in syst.T:
            a[idx_l, t] = m.addVar(name=f"a_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=l.K_MAX[t] * KEY_RATE_SCALE * syst.THETA)

    ## Demand satisfaction variable
    if problem == 1:
        for idx_d, d in enumerate(demands):
            u[idx_d] = m.addVar(name=f"u_{idx_d}", vtype=GRB.BINARY)

    nodes = gss + haps
    
    ## Node order variable --> To prevent subtours
    for idx_n, n in enumerate(nodes):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                o[idx_n, idx_d, t] = m.addVar(name=f"o_{idx_n}_{idx_d}_{t}", vtype=GRB.CONTINUOUS)

    # Primary objective: maximize u
    ## Case 1: Maximize demand satisfaction
    if problem == 1:
        m.ModelSense = GRB.MAXIMIZE

        # Primary objective: maximize demand satisfaction
        m.setObjectiveN(gp.quicksum(u[idx_d]
                                    for idx_d, d in enumerate(demands)
                                   ), index=0, priority=3, weight=1.0, abstol=1e-2, reltol=1e-2, name="Primary")
        
        # Secondary objective: maximize keys served
        m.setObjectiveN(gp.quicksum(gp.quicksum(r_h[idx_d, t]
                                                for idx_d, d in enumerate(demands)
                                               )
                                    for t in syst.T
                                   ), index=1, priority=2, weight=1.0, abstol=1e-2, reltol=1e-2, name="Secondary")

    ## Case 2: Maximize total served keys
    elif problem == 2:
        # Objective: maximize keys served
        m.setObjective(sum(sum(r_h[idx_d, t]
                               for idx_d, d in enumerate(demands)
                              )
                           for t in syst.T
                          ) * syst.THETA, GRB.MAXIMIZE
                      )

    m.Params.MIPGap = 0.01      # 1% optimality
    m.Params.MIPFocus = 1       # focus on finding feasible solutions
    m.Params.Heuristics = 0.5   # increase heuristics
    m.Params.Cuts = 1           # reduce cut aggressiveness

    # ### Tuning the accuracy and convergence of the solver
    # m.setParam("MIPGap", 1e-2)
    # m.setParam("MIPGapAbs", 1e-2)
    # m.setParam("FeasibilityTol", 1e-2)
    # m.setParam("IntFeasTol", 1e-2)
    # m.setParam("OptimalityTol", 1e-2)

    # m.Params.Presolve = 2
    # m.Params.Method = 2
    # m.Params.Cuts = 2
    # m.Params.Heuristics = 0.25
    # m.Params.MIPFocus = 1
    # m.Params.NumericFocus = 1
    # m.Params.Threads = 0
    # m.Params.NodefileStart = 0.5
    # m.Params.NoRelHeurTime = 20
    # m.Params.ConcurrentMIP = 1

    ## Constraints
    # --- NEW: Explicitly zero out flows for inactive demands ---
    # This prevents the solver from using links for demands that aren't "there"
    for idx_d, d in enumerate(demands):
        for t in syst.T:
            if demand_active[idx_d][t] == 0:
                m.addConstr(r_h[idx_d,t] == 0, name=f"mask_rh_{idx_d}_{t}")
                for idx_l in range(len(links)):
                    m.addConstr(r_1[idx_l, idx_d, t] == 0, name=f"mask_r1_{idx_l}_{idx_d}_{t}")
                    m.addConstr(r_2[idx_l, idx_d, t] == 0, name=f"mask_r2_{idx_l}_{idx_d}_{t}")
                    m.addConstr(z[idx_l, idx_d, t]   == 0, name=f"mask_z__{idx_l}_{idx_d}_{t}")

    ## --- NEW: Explicitly zero out flows for inactive links ---
    for idx_l, l in enumerate(links):
        for t in syst.T:
            if link_active[idx_l][t] == 0:
                # Force storage growth to zero for this link
                m.addConstr(a[idx_l, t] == 0, name=f"mask_a_link_{idx_l}_{t}")
                for idx_d in range(len(demands)):
                    # Force all flow (direct and from storage) to zero for this link-demand pair
                    m.addConstr(r_1[idx_l, idx_d, t] == 0, name=f"mask_r1_link_{idx_l}_{idx_d}_{t}")
                    m.addConstr(r_2[idx_l, idx_d, t] == 0, name=f"mask_r2_link_{idx_l}_{idx_d}_{t}")
                    m.addConstr(z[idx_l, idx_d, t]   == 0, name=f"mask_z__{idx_l}_{idx_d}_{t}")
    
    # Maximum Link Capacity --> Enforces only r_1
    m.addConstrs(
        (
            gp.quicksum(r_1[idx_l, idx_d, t]
                        for idx_d, d in enumerate(demands)
                       ) <= l.K_MAX[t] * KEY_RATE_SCALE
            for idx_l, l in enumerate(links)
            for t        in syst.T
        ), name="max_link_capacity"
    )

    # Maximum Tx/Rx Connection
    m.addConstrs(
        (
            gp.quicksum(z[idx_l, idx_d, t]
                        for idx_l, l in enumerate(links)
                        for idx_d, d in enumerate(demands)
                        if  l.n1 == h
                       ) <= h.N_TX
            for idx_h, h in enumerate(haps)
            for t        in syst.T
        ), name="max_tx_connections"
    )
    
    m.addConstrs(
        (
            gp.quicksum(z[idx_l, idx_d, t]
                        for idx_l, l in enumerate(links)
                        for idx_d, d in enumerate(demands)
                        if l.n2 == h
                       ) <= h.N_RX
            for idx_h, h in enumerate(haps)
            for t        in syst.T
        ), name="max_rx_connections"
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

    # MTZ subtour elimination --> Eliminates pointless single/multi-hop loops in the flows --> Uses an ordering values for all nodes
    # --> The order values should only increase on the path --> A decrease in order value == a loop (X)
    M = len(nodes)
    m.addConstrs(
        (
            o[nodes.index(l.n2), idx_d, t] >= o[nodes.index(l.n1), idx_d, t] + 1 - M * (1 - z[idx_l, idx_d, t])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="ordering_constraint_1"
    )
    m.addConstrs(
        (
            o[nodes.index(d.n1), idx_d, t] == 1
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="ordering_constraint_2"
    )
    m.addConstrs(
        (
            o[nodes.index(d.n2), idx_d, t] == M
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="ordering_constraint_2"
    )

    # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
    m.addConstrs(
        (
            r_h[idx_d, t] >= r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]# + d.K_REQ[t] * KEY_RATE_SCALE * (1 - z[idx_l, idx_d, t])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="demand_link_coordination_1"
    )

    # # Key rate and routing coordination (1)
    # m.addConstrs(
    #     (
    #         r_1[idx_l, idx_d, t] >= 1e-1 * z[idx_l, idx_d, t]
    #         for idx_l, l in enumerate(links)
    #         for idx_d, d in enumerate(demands)
    #         for t        in syst.T
    #     ), name="key_rate_routing_coordination_1"
    # )
    # Key rate and routing coordination (2)
    m.addConstrs(
        (
            r_1[idx_l, idx_d, t] <= d.K_REQ[t] * KEY_RATE_SCALE * z[idx_l, idx_d, t]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="key_rate_routing_coordination_2"
    )
    
    # Demand satisfaction check
    if problem == 1:
        M = [sum(d.K_REQ[t]
                 for t in syst.T
                )
             for d in demands
            ]
        
        m.addConstrs(
            (sum(r_h[idx_d, t]
                 for t in syst.T
                ) >= sum(d.K_REQ[t]
                         for t in syst.T
                        ) - M[idx_d] * (1 - u[idx_d])
             for idx_d, d in enumerate(demands)
            ), name="satisfaction_lower"
        )
        m.addConstrs(
            (sum(r_h[idx_d, t]
                 for t in syst.T
                ) <= sum(d.K_REQ[t]
                         for t in syst.T
                        ) + M[idx_d] * u[idx_d]
             for idx_d, d in enumerate(demands)
            ), name="satisfaction_upper"
        )
    
    # Whether to deploy QKP or not
    if f_qkp:
        # QKP on HAPs and GSs
        m.addConstrs(
                (r_2[idx_l, idx_d, 0] == 0
                 for idx_l, l in enumerate(links)
                 for idx_d, d in enumerate(demands)),
                name="initial_empty_QKP"
            )

        m.addConstrs(
            (
                gp.quicksum(a[idx_l, tp]
                            for tp in range(t)
                           ) >= syst.THETA * gp.quicksum(r_2[idx_l, idx_d, t]
                                                         for idx_d, d in enumerate(demands)
                                                        ) * STORAGE_SCALE
                for idx_l, l in enumerate(links)
                for t        in syst.T[1:]
            ), name="qkp_min_capacity"
        )
        
        m.addConstrs(
            (
                a[idx_l, t] == syst.THETA * (l.K_MAX[t] * gp.quicksum(z[idx_l, idx_d, t]
                                                                      for idx_d, d in enumerate(demands)
                                                                     ) * KEY_RATE_SCALE - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                                                                                                      for idx_d, d in enumerate(demands)
                                                                                                     )
                                            ) * STORAGE_SCALE
                for idx_l, l in enumerate(links)
                for t        in syst.T
            ), name="qkp_sequence"
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
    else:
        m.addConstrs(
            (
                a[idx_l, t] == 0
                for idx_l, l in enumerate(links)
                for t        in syst.T
            ), name="No_QKP_2"
        )

        m.addConstrs(
                (r_2[idx_l, idx_d, t] == 0
                 for idx_l, l in enumerate(links)
                 for idx_d, d in enumerate(demands)
                 for t in syst.T
                ),
                name="No_QKP_1"
            )

    m.optimize()

    k_srv = 0
    a_lst = 0
    if m.status == GRB.OPTIMAL:
        # print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        # Store solutions as dict of numpy arrays
        solution_all = {
            "o":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in o.items() if abs(v.X) > 0},
            "u":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in u.items()},
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
            "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items()},
            "z":   {k: v.X for k, v in z.items()}
        }
        solution_filtered = {
            "o":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in o.items() if abs(v.X) > 0},
            "u":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in u.items() if abs(v.X) > 0},
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items() if abs(v.X) > 0},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items() if abs(v.X) > 0},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items() if abs(v.X) > 0},
            "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items() if abs(v.X) > 0}
        }

        r_total = {k: solution_all["r_1"].get(k, 0) + solution_all["r_2"].get(k, 0)
           for k in set(solution_all["r_1"]) | set(solution_all["r_2"])}

        # plot_z_timeline(solution_all["z"], links, demands)
        # plot_z_timeline(r_total, links, demands)

        # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
        # pp.pprint(solution_filtered)

        k_srv = sum(solution_all["r_h"][idx_d, t]
                    for idx_d, d in enumerate(demands)
                    for t in syst.T
                   ) * syst.THETA

        k_req = sum(d.K_REQ[t]
                    for idx_d, d in enumerate(demands)
                    for t in syst.T
                   ) * syst.THETA

        a_lst = sum(sum(solution_all["a"][idx_l, t]
                        for t in syst.T
                       )
                    for idx_l, l in enumerate(links)
                   )

        #print(f"k_req: {k_req}, k_srv: {k_srv}")

        # for idx_l, l in enumerate(links):
        #     for t in syst.T:
        #         print(f"K_MAX[{idx_l}][{t}]: {l.K_MAX[t]}, {l.n1.tag}, {l.n2.tag}")
        
        #print(solution)
    else:
        print("No optimal solution found.")
        solution_all = None
        
    return solution_all, k_srv, a_lst

def offline_lp(gss, haps, links, demands, f_qkp, problem, Z, demand_active, link_active):
    # Create Optimization Model
    m = gp.Model("hap-qkd")
    
    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    r_1, r_2, r_h, a, u = {}, {}, {}, {}, {}

    for idx_l, l in enumerate(links):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                # demand_active is np.ones(MAX_DEMANDS), so we check the idx_d
                is_active = (demand_active[idx_d][t] == 1)
                
                # If inactive, we'll force ub to 0.0 later or set it here
                upper_bound = d.K_REQ[t] * KEY_RATE_SCALE if is_active else 0.0
                
                r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
                r_2[idx_l, idx_d, t] = m.addVar(name=f"r_2_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
                
    for idx_d, d in enumerate(demands):
        for t in syst.T:
            # demand_active is np.ones(MAX_DEMANDS), so we check the idx_d
            is_active = (demand_active[idx_d][t] == 1)
            
            # If inactive, we'll force ub to 0.0 later or set it here
            upper_bound = d.K_REQ[t] * KEY_RATE_SCALE if is_active else 0.0
            r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
            
    for idx_l, l in enumerate(links):
        for t in syst.T:
            a[idx_l, t] = m.addVar(name=f"a_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=l.K_MAX[t] * KEY_RATE_SCALE * syst.THETA)

    ## Demand satisfaction variable
    if problem == 1:
        for idx_d, d in enumerate(demands):
            u[idx_d] = m.addVar(name=f"u_{idx_d}", vtype=GRB.BINARY)

    nodes = gss + haps

    # Primary objective: maximize u
    ## Case 1: Maximize demand satisfaction
    if problem == 1:
        m.ModelSense = GRB.MAXIMIZE

        # Primary objective: maximize demand satisfaction
        m.setObjectiveN(gp.quicksum(u[idx_d]
                                    for idx_d, d in enumerate(demands)
                                   ), index=0, priority=3, weight=1.0, abstol=1e-2, reltol=1e-2, name="Primary")
        
        # Secondary objective: maximize keys served
        m.setObjectiveN(gp.quicksum(gp.quicksum(r_h[idx_d, t]
                                                for idx_d, d in enumerate(demands)
                                               )
                                    for t in syst.T
                                   ), index=1, priority=2, weight=1.0, abstol=1e-2, reltol=1e-2, name="Secondary")

    ## Case 2: Maximize total served keys
    elif problem == 2:
        # Objective: maximize keys served
        m.setObjective(sum(sum(r_h[idx_d, t]
                               for idx_d, d in enumerate(demands)
                              )
                           for t in syst.T
                          ) * syst.THETA, GRB.MAXIMIZE
                      )

    # ### Tuning the accuracy and convergence of the solver
    # m.setParam("MIPGap", 1e-2)
    # m.setParam("MIPGapAbs", 1e-2)
    # m.setParam("FeasibilityTol", 1e-2)
    # m.setParam("IntFeasTol", 1e-2)
    # m.setParam("OptimalityTol", 1e-2)

    # m.Params.Presolve = 2
    # m.Params.Method = 2
    # m.Params.Cuts = 2
    # m.Params.Heuristics = 0.25
    # m.Params.MIPFocus = 1
    # m.Params.NumericFocus = 1
    # m.Params.Threads = 0
    # m.Params.NodefileStart = 0.5
    # m.Params.NoRelHeurTime = 20
    # m.Params.ConcurrentMIP = 1

    ## Constraints
    # --- NEW: Explicitly zero out flows for inactive demands ---
    # This prevents the solver from using links for demands that aren't "there"
    for idx_d, d in enumerate(demands):
        for t in syst.T:
            if demand_active[idx_d][t] == 0:
                m.addConstr(r_h[idx_d,t] == 0, name=f"mask_rh_{idx_d}_{t}")
                for idx_l in range(len(links)):
                    m.addConstr(r_1[idx_l, idx_d, t] == 0, name=f"mask_r1_{idx_l}_{idx_d}_{t}")
                    m.addConstr(r_2[idx_l, idx_d, t] == 0, name=f"mask_r2_{idx_l}_{idx_d}_{t}")

    ## --- NEW: Explicitly zero out flows for inactive links ---
    for idx_l, l in enumerate(links):
        for t in syst.T:
            if link_active[idx_l][t] == 0:
                # Force storage growth to zero for this link
                m.addConstr(a[idx_l, t] == 0, name=f"mask_a_link_{idx_l}_{t}")
                for idx_d in range(len(demands)):
                    # Force all flow (direct and from storage) to zero for this link-demand pair
                    m.addConstr(r_1[idx_l, idx_d, t] == 0, name=f"mask_r1_link_{idx_l}_{idx_d}_{t}")
                    m.addConstr(r_2[idx_l, idx_d, t] == 0, name=f"mask_r2_link_{idx_l}_{idx_d}_{t}")
    
    # Maximum Link Capacity --> Enforces only r_1
    m.addConstrs(
        (
            gp.quicksum(r_1[idx_l, idx_d, t]
                        for idx_d, d in enumerate(demands)
                       ) <= l.K_MAX[t] * KEY_RATE_SCALE
            for idx_l, l in enumerate(links)
            for t        in syst.T
        ), name="max_link_capacity"
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
        (
            r_h[idx_d, t] >= r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]# + d.K_REQ[t] * KEY_RATE_SCALE * (1 - z[idx_l, idx_d, t])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="demand_link_coordination_1"
    )

    # # Key rate and routing coordination (1)
    # m.addConstrs(
    #     (
    #         r_1[idx_l, idx_d, t] >= 1e-1 * z[idx_l, idx_d, t]
    #         for idx_l, l in enumerate(links)
    #         for idx_d, d in enumerate(demands)
    #         for t        in syst.T
    #     ), name="key_rate_routing_coordination_1"
    # )
    # Key rate and routing coordination (2)
    m.addConstrs(
        (
            r_1[idx_l, idx_d, t] <= d.K_REQ[t] * KEY_RATE_SCALE * Z[idx_l][idx_d][t]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="key_rate_routing_coordination_2"
    )
    
    # Demand satisfaction check
    if problem == 1:
        M = [sum(d.K_REQ[t]
                 for t in syst.T
                )
             for d in demands
            ]
        
        m.addConstrs(
            (sum(r_h[idx_d, t]
                 for t in syst.T
                ) >= sum(d.K_REQ[t]
                         for t in syst.T
                        ) - M[idx_d] * (1 - u[idx_d])
             for idx_d, d in enumerate(demands)
            ), name="satisfaction_lower"
        )
        m.addConstrs(
            (sum(r_h[idx_d, t]
                 for t in syst.T
                ) <= sum(d.K_REQ[t]
                         for t in syst.T
                        ) + M[idx_d] * u[idx_d]
             for idx_d, d in enumerate(demands)
            ), name="satisfaction_upper"
        )
    
    # Whether to deploy QKP or not
    if f_qkp:
        # QKP on HAPs and GSs
        m.addConstrs(
                (r_2[idx_l, idx_d, 0] == 0
                 for idx_l, l in enumerate(links)
                 for idx_d, d in enumerate(demands)),
                name="initial_empty_QKP"
            )

        m.addConstrs(
            (
                gp.quicksum(a[idx_l, tp]
                            for tp in range(t)
                           ) >= syst.THETA * gp.quicksum(r_2[idx_l, idx_d, t]
                                                         for idx_d, d in enumerate(demands)
                                                        ) * STORAGE_SCALE
                for idx_l, l in enumerate(links)
                for t        in syst.T[1:]
            ), name="qkp_min_capacity"
        )
        
        m.addConstrs(
            (
                a[idx_l, t] == syst.THETA * (l.K_MAX[t] * sum(Z[idx_l][idx_d][t]
                                                              for idx_d, d in enumerate(demands)
                                                             ) * KEY_RATE_SCALE - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                                                                                              for idx_d, d in enumerate(demands)
                                                                                             )
                                            ) * STORAGE_SCALE
                for idx_l, l in enumerate(links)
                for t        in syst.T
            ), name="qkp_sequence"
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
    else:
        m.addConstrs(
            (
                a[idx_l, t] == 0
                for idx_l, l in enumerate(links)
                for t        in syst.T
            ), name="No_QKP_2"
        )

        m.addConstrs(
                (r_2[idx_l, idx_d, t] == 0
                 for idx_l, l in enumerate(links)
                 for idx_d, d in enumerate(demands)
                 for t in syst.T
                ),
                name="No_QKP_1"
            )

    m.optimize()

    k_srv = 0
    a_lst = 0
    if m.status == GRB.OPTIMAL:
        #print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        # Store solutions as dict of numpy arrays
        solution_all = {
            "u":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in u.items()},
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
            "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items()},
        }
        solution_filtered = {
            "u":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in u.items() if abs(v.X) > 0},
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items() if abs(v.X) > 0},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items() if abs(v.X) > 0},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items() if abs(v.X) > 0},
            "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items() if abs(v.X) > 0}
        }

        r_total = {k: solution_all["r_1"].get(k, 0) + solution_all["r_2"].get(k, 0)
           for k in set(solution_all["r_1"]) | set(solution_all["r_2"])}

        #plot_z_timeline(solution_all["z"], links, demands)
        #plot_z_timeline(r_total, links, demands)

        # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
        # pp.pprint(solution_filtered)

        k_srv = sum(solution_all["r_h"][idx_d, t]
                    for idx_d, d in enumerate(demands)
                    for t in syst.T
                   ) * syst.THETA

        k_req = sum(d.K_REQ[t]
                    for idx_d, d in enumerate(demands)
                    for t in syst.T
                   ) * syst.THETA

        a_lst = sum(sum(solution_all["a"][idx_l, t]
                        for t in syst.T
                       )
                    for idx_l, l in enumerate(links)
                   )

        #print(f"k_req: {k_req}, k_srv: {k_srv}")

        # for idx_l, l in enumerate(links):
        #     for t in syst.T:
        #         print(f"K_MAX[{idx_l}][{t}]: {l.K_MAX[t]}, {l.n1.tag}, {l.n2.tag}")
        
        #print(solution)
    else:
        print("No optimal solution found.")
        solution_all = None
        
    return solution_all, k_srv, a_lst

def offline_relaxed(gss, haps, links, demands, f_qkp, problem, demand_active, link_active):
    # Create Optimization Model
    m = gp.Model("hap-qkd")
    
    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    r_1, r_2, r_h, a, z, u, o = {}, {}, {}, {}, {}, {}, {}

    for idx_l, l in enumerate(links):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                # demand_active is np.ones(MAX_DEMANDS), so we check the idx_d
                is_active = (demand_active[idx_d][t] == 1)
                
                # If inactive, we'll force ub to 0.0 later or set it here
                upper_bound = d.K_REQ[t] * KEY_RATE_SCALE if is_active else 0.0
                
                r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
                r_2[idx_l, idx_d, t] = m.addVar(name=f"r_2_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
                z[idx_l, idx_d, t]   = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}",   vtype=GRB.CONTINUOUS, lb=0.0, ub=is_active)
                
    for idx_d, d in enumerate(demands):
        for t in syst.T:
            # demand_active is np.ones(MAX_DEMANDS), so we check the idx_d
            is_active = (demand_active[idx_d][t] == 1)
            
            # If inactive, we'll force ub to 0.0 later or set it here
            upper_bound = d.K_REQ[t] * KEY_RATE_SCALE if is_active else 0.0
            r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
            
    for idx_l, l in enumerate(links):
        for t in syst.T:
            a[idx_l, t] = m.addVar(name=f"a_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=-GRB.INFINITY, ub=l.K_MAX[t] * KEY_RATE_SCALE * syst.THETA)

    ## Demand satisfaction variable
    if problem == 1:
        for idx_d, d in enumerate(demands):
            u[idx_d] = m.addVar(name=f"u_{idx_d}", vtype=GRB.BINARY)

    nodes = gss + haps
    
    ## Node order variable --> To prevent subtours
    for idx_n, n in enumerate(nodes):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                o[idx_n, idx_d, t] = m.addVar(name=f"o_{idx_n}_{idx_d}_{t}", vtype=GRB.CONTINUOUS)

    # Primary objective: maximize u
    ## Case 1: Maximize demand satisfaction
    if problem == 1:
        m.ModelSense = GRB.MAXIMIZE

        # Primary objective: maximize demand satisfaction
        m.setObjectiveN(gp.quicksum(u[idx_d]
                                    for idx_d, d in enumerate(demands)
                                   ), index=0, priority=3, weight=1.0, abstol=1e-2, reltol=1e-2, name="Primary")
        
        # Secondary objective: maximize keys served
        m.setObjectiveN(gp.quicksum(gp.quicksum(r_h[idx_d, t]
                                                for idx_d, d in enumerate(demands)
                                               )
                                    for t in syst.T
                                   ), index=1, priority=2, weight=1.0, abstol=1e-2, reltol=1e-2, name="Secondary")

    ## Case 2: Maximize total served keys
    elif problem == 2:
        # Objective: maximize keys served
        m.setObjective(sum(sum(r_h[idx_d, t]
                               for idx_d, d in enumerate(demands)
                              )
                           for t in syst.T
                          ) * syst.THETA, GRB.MAXIMIZE
                      )

    m.Params.Method = 2        # Barrier
    m.Params.Crossover = 0    # Skip crossover (LP only!)

    m.Params.Presolve = 2
    m.Params.Cuts = 0
    m.Params.Heuristics = 0
    m.Params.MIPFocus = 0

    ### Tuning the accuracy and convergence of the solver
    m.setParam("MIPGap", 1e-4)
    m.setParam("MIPGapAbs", 1e-4)
    m.setParam("FeasibilityTol", 1e-4)
    m.setParam("IntFeasTol", 1e-4)
    m.setParam("OptimalityTol", 1e-4)

    # m.Params.Presolve = 2
    # m.Params.Method = 2
    # m.Params.Cuts = 2
    # m.Params.Heuristics = 0.25
    # m.Params.MIPFocus = 1
    # m.Params.NumericFocus = 1
    # m.Params.Threads = 0
    # m.Params.NodefileStart = 0.5
    # m.Params.NoRelHeurTime = 20
    # m.Params.ConcurrentMIP = 1

    ## Constraints
    # --- NEW: Explicitly zero out flows for inactive demands ---
    # This prevents the solver from using links for demands that aren't "there"
    for idx_d, d in enumerate(demands):
        for t in syst.T:
            if demand_active[idx_d][t] == 0:
                m.addConstr(r_h[idx_d,t] == 0, name=f"mask_rh_{idx_d}_{t}")
                for idx_l in range(len(links)):
                    m.addConstr(r_1[idx_l, idx_d, t] == 0, name=f"mask_r1_{idx_l}_{idx_d}_{t}")
                    m.addConstr(r_2[idx_l, idx_d, t] == 0, name=f"mask_r2_{idx_l}_{idx_d}_{t}")
                    m.addConstr(z[idx_l, idx_d, t]   == 0, name=f"mask_z__{idx_l}_{idx_d}_{t}")

    ## --- NEW: Explicitly zero out flows for inactive links ---
    for idx_l, l in enumerate(links):
        for t in syst.T:
            if link_active[idx_l][t] == 0:
                # Force storage growth to zero for this link
                m.addConstr(a[idx_l, t] == 0, name=f"mask_a_link_{idx_l}_{t}")
                for idx_d in range(len(demands)):
                    # Force all flow (direct and from storage) to zero for this link-demand pair
                    m.addConstr(r_1[idx_l, idx_d, t] == 0, name=f"mask_r1_link_{idx_l}_{idx_d}_{t}")
                    m.addConstr(r_2[idx_l, idx_d, t] == 0, name=f"mask_r2_link_{idx_l}_{idx_d}_{t}")
                    m.addConstr(z[idx_l, idx_d, t]   == 0, name=f"mask_z__{idx_l}_{idx_d}_{t}")

    
    # Maximum Link Capacity --> Enforces only r_1
    m.addConstrs(
        (
            gp.quicksum(r_1[idx_l, idx_d, t]
                        for idx_d, d in enumerate(demands)
                       ) <= l.K_MAX[t] * KEY_RATE_SCALE
            for idx_l, l in enumerate(links)
            for t        in syst.T
        ), name="max_link_capacity"
    )

    # Maximum Tx/Rx Connection
    m.addConstrs(
        (
            gp.quicksum(z[idx_l, idx_d, t]
                        for idx_l, l in enumerate(links)
                        for idx_d, d in enumerate(demands)
                        if  l.n1 == h
                       ) <= h.N_TX
            for idx_h, h in enumerate(haps)
            for t        in syst.T
        ), name="max_tx_connections"
    )
    
    m.addConstrs(
        (
            gp.quicksum(z[idx_l, idx_d, t]
                        for idx_l, l in enumerate(links)
                        for idx_d, d in enumerate(demands)
                        if l.n2 == h
                       ) <= h.N_RX
            for idx_h, h in enumerate(haps)
            for t        in syst.T
        ), name="max_rx_connections"
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

    # MTZ subtour elimination --> Eliminates pointless single/multi-hop loops in the flows --> Uses an ordering values for all nodes
    # --> The order values should only increase on the path --> A decrease in order value == a loop (X)
    M = len(nodes)
    m.addConstrs(
        (
            o[nodes.index(l.n2), idx_d, t] >= o[nodes.index(l.n1), idx_d, t] + 1 - M * (1 - z[idx_l, idx_d, t])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="ordering_constraint_1"
    )
    m.addConstrs(
        (
            o[nodes.index(d.n1), idx_d, t] == 1
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="ordering_constraint_2"
    )
    m.addConstrs(
        (
            o[nodes.index(d.n2), idx_d, t] == M
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="ordering_constraint_2"
    )

    # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
    m.addConstrs(
        (
            r_h[idx_d, t] >= r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]# + d.K_REQ[t] * KEY_RATE_SCALE * (1 - z[idx_l, idx_d, t])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="demand_link_coordination_1"
    )

    # # Key rate and routing coordination (1)
    # m.addConstrs(
    #     (
    #         r_1[idx_l, idx_d, t] >= 1e-1 * z[idx_l, idx_d, t]
    #         for idx_l, l in enumerate(links)
    #         for idx_d, d in enumerate(demands)
    #         for t        in syst.T
    #     ), name="key_rate_routing_coordination_1"
    # )
    # Key rate and routing coordination (2)
    m.addConstrs(
        (
            r_1[idx_l, idx_d, t] <= d.K_REQ[t] * KEY_RATE_SCALE * z[idx_l, idx_d, t]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="key_rate_routing_coordination_2"
    )
    
    # Demand satisfaction check
    if problem == 1:
        M = [sum(d.K_REQ[t]
                 for t in syst.T
                )
             for d in demands
            ]
        
        m.addConstrs(
            (sum(r_h[idx_d, t]
                 for t in syst.T
                ) >= sum(d.K_REQ[t]
                         for t in syst.T
                        ) - M[idx_d] * (1 - u[idx_d])
             for idx_d, d in enumerate(demands)
            ), name="satisfaction_lower"
        )
        m.addConstrs(
            (sum(r_h[idx_d, t]
                 for t in syst.T
                ) <= sum(d.K_REQ[t]
                         for t in syst.T
                        ) + M[idx_d] * u[idx_d]
             for idx_d, d in enumerate(demands)
            ), name="satisfaction_upper"
        )
    
    # Whether to deploy QKP or not
    if f_qkp:
        # QKP on HAPs and GSs
        m.addConstrs(
                (r_2[idx_l, idx_d, 0] == 0
                 for idx_l, l in enumerate(links)
                 for idx_d, d in enumerate(demands)),
                name="initial_empty_QKP"
            )

        m.addConstrs(
            (
                gp.quicksum(a[idx_l, tp]
                            for tp in range(t)
                           ) >= syst.THETA * gp.quicksum(r_2[idx_l, idx_d, t]
                                                         for idx_d, d in enumerate(demands)
                                                        ) * STORAGE_SCALE
                for idx_l, l in enumerate(links)
                for t        in syst.T[1:]
            ), name="qkp_min_capacity"
        )
        
        m.addConstrs(
            (
                a[idx_l, t] == syst.THETA * (l.K_MAX[t] * gp.quicksum(z[idx_l, idx_d, t]
                                                                      for idx_d, d in enumerate(demands)
                                                                     ) * KEY_RATE_SCALE - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                                                                                                      for idx_d, d in enumerate(demands)
                                                                                                     )
                                            ) * STORAGE_SCALE
                for idx_l, l in enumerate(links)
                for t        in syst.T
            ), name="qkp_sequence"
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
    else:
        m.addConstrs(
            (
                a[idx_l, t] == 0
                for idx_l, l in enumerate(links)
                for t        in syst.T
            ), name="No_QKP_2"
        )

        m.addConstrs(
                (r_2[idx_l, idx_d, t] == 0
                 for idx_l, l in enumerate(links)
                 for idx_d, d in enumerate(demands)
                 for t in syst.T
                ),
                name="No_QKP_1"
            )

    m.optimize()

    k_srv = 0
    a_lst = 0
    if m.status == GRB.OPTIMAL:
        # print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        # Store solutions as dict of numpy arrays
        solution_all = {
            "o":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in o.items() if abs(v.X) > 0},
            "u":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in u.items()},
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
            "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items()},
            "z":   {k: v.X for k, v in z.items()}
        }
        solution_filtered = {
            "o":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in o.items() if abs(v.X) > 0},
            "u":   {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in u.items() if abs(v.X) > 0},
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items() if abs(v.X) > 0},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items() if abs(v.X) > 0},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items() if abs(v.X) > 0},
            "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items() if abs(v.X) > 0}
        }

        r_total = {k: solution_all["r_1"].get(k, 0) + solution_all["r_2"].get(k, 0)
           for k in set(solution_all["r_1"]) | set(solution_all["r_2"])}

        # plot_z_timeline(solution_all["z"], links, demands)
        # plot_z_timeline(r_total, links, demands)

        # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
        # pp.pprint(solution_filtered)

        k_srv = sum(solution_all["r_h"][idx_d, t]
                    for idx_d, d in enumerate(demands)
                    for t in syst.T
                   ) * syst.THETA

        k_req = sum(d.K_REQ[t]
                    for idx_d, d in enumerate(demands)
                    for t in syst.T
                   ) * syst.THETA

        a_lst = sum(sum(solution_all["a"][idx_l, t]
                        for t in syst.T
                       )
                    for idx_l, l in enumerate(links)
                   )

        #print(f"k_req: {k_req}, k_srv: {k_srv}")

        # for idx_l, l in enumerate(links):
        #     for t in syst.T:
        #         print(f"K_MAX[{idx_l}][{t}]: {l.K_MAX[t]}, {l.n1.tag}, {l.n2.tag}")
        
        #print(solution)
    else:
        print("No optimal solution found.")
        solution_all = None
        
    return solution_all, k_srv, a_lst