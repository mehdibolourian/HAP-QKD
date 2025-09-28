from libraries import *

def offline_wind(gss, haps, links, demands):
    # Create Optimization Model
    m = gp.Model("hap-qkd")
    
    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    r_1, r_2, r_h, a, z = {}, {}, {}, {}, {}

    for idx_l, l in enumerate(links):
        for idx_d, d in enumerate(demands):
            for t in syst.T:
                r_1[idx_l, idx_d, t] = m.addVar(name=f"r_1_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
                r_2[idx_l, idx_d, t] = m.addVar(name=f"r_2_{idx_l}_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
                z[idx_l, idx_d, t]   = m.addVar(name=f"z_{idx_l}_{idx_d}_{t}",   vtype=GRB.BINARY)
                
    for idx_d, d in enumerate(demands):
        for t in syst.T:
            r_h[idx_d, t] = m.addVar(name=f"r_h_{idx_d}_{t}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
            
    for idx_l, l in enumerate(links):
        for t in syst.T:
            a[idx_l, t] = m.addVar(name=f"a_{idx_l}_{t}", vtype=GRB.CONTINUOUS, lb=0.0)

    m.ModelSense = GRB.MAXIMIZE
    
    # Primary objective: maximize r_h
    m.setObjectiveN(sum(sum(r_h[idx_d, t]
                           for idx_d, d in enumerate(demands)
                          )
                       for t in syst.T
                      ) * syst.THETA, index=0, priority=2, weight=1.0, abstol=1e-9, reltol=1e-9, name="Primary")

    # m.setObjectiveN(sum(sum(sum(r_1[idx_l, idx_d, t]
    #                             for idx_l, l in enumerate(links)
    #                            )
    #                         for idx_d, d in enumerate(demands)
    #                        )
    #                     for t in syst.T
    #                    ) * syst.THETA, index=1, priority=1, weight=1.0, abstol=1e-9, reltol=1e-9, name="Secondary")

    # # Secondary objective: maximize a
    # m.setObjectiveN(sum(sum(a[idx_l, t]
    #                        for idx_l, l in enumerate(links)
    #                       )
    #                    for t in syst.T
    #                   ), index=1, priority=1, weight=1.0, abstol=1e-6, reltol=1e-6, name="Secondary")

    m.setParam("MIPGap", 1e-9)          # force very tight gap
    m.setParam("MIPGapAbs", 1e-9)
    m.setParam("FeasibilityTol", 1e-9)
    m.setParam("IntFeasTol", 1e-9)
    m.setParam("OptimalityTol", 1e-9)

    ## Constraints
    # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
    # r_h = min_{l:z_l=1}(r_1+r_2)
    m.addConstrs(
        (
            r_h[idx_d, t] <= r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] + d.K_REQ[t] * KEY_RATE_SCALE * (1 - z[idx_l, idx_d, t])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="demand_link_coordination"
    )

    EPSILON = 1e-8
    # Key rate and routing coordination (1)
    m.addConstrs(
        (
            r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] >= EPSILON * z[idx_l, idx_d, t]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="demand_link_coordination_1"
    )
    # Key rate and routing coordination (2)
    m.addConstrs(
        (
            r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t] <= d.K_REQ[t] * KEY_RATE_SCALE * z[idx_l, idx_d, t]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
            for t        in syst.T
        ), name="key_rate_routing_coordination_2"
    )
    
    # Max Key Rate
    m.addConstrs(
        (
            sum(r_1[idx_l, idx_d, t]
                for idx_d, d in enumerate(demands)
               ) <= l.K_MAX[t] * KEY_RATE_SCALE
            for idx_l, l in enumerate(links)
            for t        in syst.T
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
            for t in syst.T
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
                if  l.n1 == n
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
    
    # QKP on HAPs and GSs
    m.addConstrs(
        (
            sum(a[idx_l, tp]
                for tp in range(t)
               ) >= syst.THETA * sum(r_2[idx_l, idx_d, t]
                                     for idx_d, d in enumerate(demands)
                                    ) * STORAGE_SCALE
            for idx_l, l in enumerate(links)
            for t        in syst.T
        ), name="qkp_min_capacity"
    )
    
    m.addConstrs(
        (
            a[idx_l, t] == syst.THETA * (l.K_MAX[t] * KEY_RATE_SCALE - sum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d, t]
                                                                           for idx_d, d in enumerate(demands)
                                                                          )
                                        ) * STORAGE_SCALE
            for idx_l, l in enumerate(links)
            for t        in syst.T
        ), name="qkp_sequence"
    )
    
    # m.addConstrs(
    #     (
    #         sum(a[idx_l, tp]
    #             for tp in range(t)
    #            ) <= min(l.n1.A_MAX, l.n2.A_MAX) * STORAGE_SCALE
    #         for idx_l, l in enumerate(links)
    #         for t        in syst.T
    #     ), name="qkp_max_capacity"
    # )

    m.optimize()
    
    if m.status == GRB.OPTIMAL:
        print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        # Store solutions as dict of numpy arrays
        solution = {
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
            "a": {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items()},
            "z": {k: v.X for k, v in z.items()}
        }

        pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
        pp.pprint(solution)

        # for idx_l, l in enumerate(links):
        #     for t in syst.T:
        #         print(f"K_MAX[{idx_l}][{t}]: {l.K_MAX[t]}")
        
        #print(solution)
    else:
        print("No optimal solution found.")
        solution = None
        
    return solution