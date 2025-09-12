from libraries import *

def offline_wind(gss, haps, links, demands):
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

    m.ModelSense = GRB.MAXIMIZE
    
    # Primary objective: maximize x
    m.setObjectiveN(sum(sum(r_h[idx_d, t]
                           for idx_d, d in enumerate(demands)
                          )
                       for t in sys.T
                      ) * sys.THETA, index=0, priority=2, weight=1.0, name="Primary")

    # Secondary objective: maximize y
    m.setObjectiveN(sum(sum(a[idx_l, t]
                           for idx_l, l in enumerate(links)
                          )
                       for t in sys.T
                      ), index=1, priority=1, weight=1.0, name="Secondary")

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
    
    # Max Key Rate
    m.addConstrs(
        (
            sum(r_1[idx_l, idx_d, t]
                for idx_d, d in enumerate(demands)
               ) <= l.K_MAX[t]
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
                if  l.n1 == n
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

    m.optimize()
    
    if m.status == GRB.OPTIMAL:
        print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        # Store solutions as dict of numpy arrays
        solution = {
            "r_1": {k: v.X / KEY_RATE_SCALE for k, v in r_1.items()},
            "r_2": {k: v.X / KEY_RATE_SCALE for k, v in r_2.items()},
            "r_h": {k: v.X / KEY_RATE_SCALE for k, v in r_h.items()},
            "a": {k: v.X / KEY_RATE_SCALE for k, v in a.items()},
            "z": {k: v.X for k, v in z.items()}
        }
        
        print(solution)
    else:
        print("No optimal solution found.")
        solution = None
        
    return solution