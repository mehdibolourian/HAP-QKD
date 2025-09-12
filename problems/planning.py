from libraries import *

def calculate_key_rate_planning(link, d_los):  
    L_geo = 20 * max(math.log10((R_TX + d_los * 1000 * THETA) / R_RX), 0)
    L_ma  = 0.01 * d_los
    
    L_t   = L_geo + L_ma
    
    ETA = 10**(-L_t/10)
    
    K_MAX = -B * math.log2(1 - ETA)
    
    return K_MAX

def planning(gss, haps, links, demands):
    # Create Optimization Model
    m = gp.Model("hap-qkd")

    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    r_1, r_2, r_h, a, z = {}, {}, {}, {}, {}

    for idx_l, l in enumerate(links):
        for idx_d, d in enumerate(demands):
            for t in sys.T:
                r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
                r_2[idx_l, idx_d, t] = m.addVar(name=f"r_2_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
                z[idx_l, idx_d, t]   = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}",   vtype=GRB.BINARY,     lb=0.0)

                
    for idx_d, d in enumerate(demands):
        for t in sys.T:
            r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)
            
    for idx_l, l in enumerate(links):
        for t in sys.T:
            a[idx_l, t] = m.addVar(name=f"a_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

    # k (key rate) and d (distance)
    k, d = {}, {}
    for idx_l, l in enumerate(links):
        for t in sys.T:
            d[idx_l, t] = m.addVar(name=f"d_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=15.0) # LoS distance (Min height in strat.)
            k[idx_l, t] = m.addVar(name=f"k_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

            dpts = np.linspace(15, 100, 10)
            kpts = [calculate_key_rate_planning(0, v) * KEY_RATE_SCALE for v in dpts]

            m.addGenConstrPWL(d[idx_l, t], k[idx_l, t], dpts, kpts, name=f"pwl_key_rate_{idx_l}_{t}")

    d_min, d_max = 15.0, 100.0
    breakpoint_counts = [5, 10, 20]  # numbers of breakpoints to test
    plt.figure(figsize=(8,5))
    for n_bp in breakpoint_counts:
        # breakpoints for this case
        dpts = np.linspace(d_min, d_max, n_bp)
        kpts = [calculate_key_rate_planning(0, d) * KEY_RATE_SCALE for d in dpts]
    
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
    c1, c2, c3 = {}, {}, {}  # x, y, z in km
    for idx_h, hnode in enumerate(haps):
        for t in sys.T:
            c1[idx_h, t] = m.addVar(lb=-1e3, ub=1e3, vtype=GRB.CONTINUOUS, name=f"c1_{idx_h}_{t}")
            c2[idx_h, t] = m.addVar(lb=-1e3, ub=1e3, vtype=GRB.CONTINUOUS, name=f"c2_{idx_h}_{t}")
            c3[idx_h, t] = m.addVar(lb=15,   ub=25,  vtype=GRB.CONTINUOUS, name=f"c3_{idx_h}_{t}")  # altitude in km

    ## Objective
    m.ModelSense = GRB.MAXIMIZE
    
    # Primary objective: maximize x
    m.setObjectiveN(sum(sum(k[idx_l, t]
                           for idx_l, l in enumerate(links)
                          ) + sum(r_h[idx_d, t]
                                  for idx_d, d in enumerate(demands)
                                 )
                       for t in sys.T
                      ) * sys.THETA, index=0, priority=1, weight=1.0, name="Primary")

    # Secondary objective: maximize y
    m.setObjectiveN(sum(sum(a[idx_l, t]
                           for idx_l, l in enumerate(links)
                          )
                       for t in sys.T
                      ), index=2, priority=1, weight=1.0, name="Secondary")

    ## Constraints
    # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
    # r_h = min_{l:z_l=1}(r_1+r_2)
    m.addConstrs(
        (
            r_h[idx_d, t] <= r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] + d.K_REQ[t] * (1 - z[idx_l, idx_d, t])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in sys.T
        ), name="demand_link_coordination"
    )
    
    EPSILON = 1e-3
    # Key rate and routing coordination (1)
    m.addConstrs(
        (
            r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] >= EPSILON * z[idx_l, idx_d, t]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in sys.T
        ), name="demand_link_coordination_1"
    )
    
    # Key rate and routing coordination (2)
    m.addConstrs(
        (
            r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] <= d.K_REQ[t] * z[idx_l, idx_d, t]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in sys.T
        ), name="key_rate_routing_coordination_2"
    )
    
    # Max key rate
    m.addConstrs(
        (
            sum(r_1[idx_l, idx_d, t]
                for idx_d, d in enumerate(demands)
               ) <= k[idx_l, t]
            for idx_l, l in enumerate(links)
            for t        in sys.T
        ), name="max_key_rate"
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
            for t        in sys.T
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
            for t        in sys.T
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
            for t        in sys.T
        ), name="flow_conservation_3"
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
            for t        in sys.T
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
            for t        in sys.T
        ), name="max_rx_connections"
    )
    
    # QKP on HAPs and GSs
    m.addConstrs(
        (
            a[idx_l, t] >= sys.THETA * sum(r_2[idx_l, idx_d, t]
                                           for idx_d, d in enumerate(demands)
                                          )
            for idx_l, l in enumerate(links)
            for t        in sys.T
        ), name="qkp_min_capacity"
    )
    
    m.addConstrs(
        (
            a[idx_l, t+1] == a[idx_l, t] + l.K_MAX[t] - sys.THETA * sum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                                                                        for idx_d, d in enumerate(demands)
                                                                       )
            for idx_l, l in enumerate(links)
            for t        in sys.T[:-1]
        ), name="qkp_sequence"
    )
    
    m.addConstrs(
        (
            a[idx_l, 0] == 0
            for idx_l, l in enumerate(links)
        ), name="qkp_init"
    )
    
    m.addConstrs(
        (
            a[idx_l, t] <= min(l.n1.A_MAX, l.n2.A_MAX)
            for idx_l, l in enumerate(links)
        ), name="qkp_max_capacity"
    )

    # For each link l = (hap, gs), add the SOCP constraint tying d to (c1,c2,c3)
    for idx_l, l in enumerate(links):
        # identify which endpoint is HAP and which is GS
        if isinstance(l.n1, hap) and isinstance(l.n2, gs):
            hap_idx, gs_node = haps.index(l.n1), l.n2
        elif isinstance(l.n2, hap) and isinstance(l.n1, gs):
            hap_idx, gs_node = haps.index(l.n2), l.n1

        # precompute GS coordinates in same (x,y,z) frame and units (km)
        [cg1, cg2] = lonlat_to_xy(gs_node.lg, gs_node.la)  # ensure these are constants in km
        
        print(cg1, cg2)

        for t in sys.T:
            # Quadratic cone: d^2 >= (c1-cg1)^2 + (c2-cg2)^2 + (c3-cg3)^2
            dx = c1[hap_idx, t] - cg1*1e-4
            dy = c2[hap_idx, t] - cg2*1e-4
            dz = c3[hap_idx, t]
            m.addQConstr(d[idx_l, t]*d[idx_l, t] >= dx*dx + dy*dy + 15*15,
                         name=f"dist_cone_{idx_l}_{t}")

    ## Solve
    m.optimize()

    if m.status == GRB.OPTIMAL:
        print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        solution = {
            "r_1": {k: v.X / KEY_RATE_SCALE for k, v in r_1.items()},
            "r_2": {k: v.X / KEY_RATE_SCALE for k, v in r_2.items()},
            "r_h": {k: v.X / KEY_RATE_SCALE for k, v in r_h.items()},
            "a":   {k: v.X / KEY_RATE_SCALE for k, v in a.items()},
            "z":   {k: v.X for k, v in z.items()},
            "k":   {k: v.X / KEY_RATE_SCALE for k, v in k.items()},
            "d":   {k: v.X * 1e3 for k, v in d.items()},
            "c1":  {k: v.X * 1e4 for k, v in c1.items()},
            "c2":  {k: v.X * 1e4 for k, v in c2.items()},
            "c3":  {k: v.X * 1e4 for k, v in c3.items()}
        }
    else:
        print("No optimal solution found.")
        solution = None

    return solution