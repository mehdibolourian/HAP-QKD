from libraries import *

def online_new(gss, haps, links, state, action, t, f_qkp):
    demands = state[t]["demands"]

    A       = state[t]["a"]
    z       = action[t]

    A_next  = A.copy()

    # print(f"z: {z}")
    # print(f"A: {A}")

    reward  = -10
    
    # Create Optimization Model
    m = gp.Model("hap-qkd")
    
    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    r_1, r_2, r_h, a = {}, {}, {}, {}

    for idx_l, l in enumerate(links):
        for idx_d, d in enumerate(demands):
            r_1[idx_l, idx_d] = m.addVar(name=f"r_1_{idx_l}_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
            r_2[idx_l, idx_d] = m.addVar(name=f"r_2_{idx_l}_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
                
    for idx_d, d in enumerate(demands):
        r_h[idx_d] = m.addVar(name=f"r_h_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)

    K_MAX_SUM = sum(l.K_MAX[tp]
                    for tp in range(t+1)
                   )
    for idx_l, l in enumerate(links):
        a[idx_l] = m.addVar(name=f"a_{idx_l}", vtype=GRB.CONTINUOUS, lb=-K_MAX_SUM, ub=l.K_MAX[t] * KEY_RATE_SCALE * syst.THETA)

    nodes = gss + haps

    m.ModelSense = GRB.MAXIMIZE

    # Primary objective: maximize demand satisfaction
    m.setObjectiveN(gp.quicksum(r_h[idx_d]
                                for idx_d, d in enumerate(demands)
                               ), index=0, priority=2, weight=1.0, abstol=1e-2, reltol=1e-2, name="Primary")
    
    # Secondary objective: maximize keys served
    m.setObjectiveN(gp.quicksum(a[idx_l]
                                for idx_l, l in enumerate(links)
                               ), index=1, priority=1, weight=1.0, abstol=1e-2, reltol=1e-2, name="Secondary")

    ### Tuning the accuracy and convergence of the solver
    m.setParam("MIPGap", 1e-2)
    m.setParam("MIPGapAbs", 1e-2)
    m.setParam("FeasibilityTol", 1e-2)
    m.setParam("IntFeasTol", 1e-2)
    m.setParam("OptimalityTol", 1e-2)

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
            gp.quicksum(r_1[idx_l, idx_d]
                        for idx_d, d in enumerate(demands)
                       ) <= l.K_MAX[t] * KEY_RATE_SCALE
            for idx_l, l in enumerate(links)
        ), name="max_link_capacity"
    )

    # Flow conservation
    m.addConstrs(
        (gp.quicksum(r_1[idx_l, idx_d] + r_2[idx_l, idx_d]
                     for idx_l, l in enumerate(links)
                     if l.n1 == d.n1)
         - gp.quicksum(r_1[idx_l, idx_d] + r_2[idx_l, idx_d]
                       for idx_l, l in enumerate(links)
                       if l.n2 == d.n1)
         == r_h[idx_d]
         for idx_d, d in enumerate(demands)
        ),
        name="flow_conservation_1"
    )
    m.addConstrs(
        (gp.quicksum(r_1[idx_l, idx_d] + r_2[idx_l, idx_d]
                     for idx_l, l in enumerate(links)
                     if l.n2 == d.n2)
         - gp.quicksum(r_1[idx_l, idx_d] + r_2[idx_l, idx_d]
                       for idx_l, l in enumerate(links)
                       if l.n1 == d.n2)
         == r_h[idx_d]
         for idx_d, d in enumerate(demands)
        ),
        name="flow_conservation_2"
    )
    m.addConstrs(
        (gp.quicksum((r_1[idx_l, idx_d] + r_2[idx_l, idx_d])
                     for idx_l, l in enumerate(links)
                     if l.n1 == n
                    )
         - gp.quicksum((r_1[idx_l, idx_d] + r_2[idx_l, idx_d])
                       for idx_l, l in enumerate(links)
                       if l.n2 == n
                      )
         == 0
         for idx_d, d in enumerate(demands)
         for n in gss + haps
         if n != d.n1 and n != d.n2
        ),
        name="flow_conservation_3"
    )

    # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
    m.addConstrs(
        (
            r_h[idx_d] >= r_1[idx_l, idx_d] + r_2[idx_l, idx_d]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
        ), name="demand_link_coordination_1"
    )

    # # Key rate and routing coordination (1)
    # m.addConstrs(
    #     (
    #         r_1[idx_l, idx_d] >= 1e-1 * z[idx_l, idx_d]
    #         for idx_l, l in enumerate(links)
    #         for idx_d, d in enumerate(demands)
    #     ), name="key_rate_routing_coordination_1"
    # )
    # Key rate and routing coordination (2)
    m.addConstrs(
        (
            r_1[idx_l, idx_d] <= d.K_REQ[t] * KEY_RATE_SCALE * z[idx_l]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
        ), name="key_rate_routing_coordination_2"
    )
    
    # Whether to deploy QKP or not
    if f_qkp:
        if t == 0:
            # QKP on HAPs and GSs
            m.addConstrs(
                    (r_2[idx_l, idx_d] == 0
                     for idx_l, l in enumerate(links)
                     for idx_d, d in enumerate(demands)),
                    name="initial_empty_QKP"
                )

        if t >= 1:
            m.addConstrs(
                (
                    A[idx_l] >= syst.THETA * gp.quicksum(r_2[idx_l, idx_d]
                                                     for idx_d, d in enumerate(demands)
                                                    ) * STORAGE_SCALE
                    for idx_l, l in enumerate(links)
                ), name="qkp_min_capacity"
            )
        
        m.addConstrs(
            (
                a[idx_l] == syst.THETA * (l.K_MAX[t] * z[idx_l] * KEY_RATE_SCALE - gp.quicksum(r_1[idx_l, idx_d] + r_2[idx_l, idx_d]
                                                                                           for idx_d, d in enumerate(demands)
                                                                                          )
                                            ) * STORAGE_SCALE
                for idx_l, l in enumerate(links)
            ), name="qkp_sequence"
        )

        m.addConstrs(
            (
                A[idx_l] + a[idx_l] >= 0
                for idx_l, l in enumerate(links)
            ), name="positive_storage"
        )
    else:
        m.addConstrs(
            (
                a[idx_l] == 0
                for idx_l, l in enumerate(links)
            ), name="No_QKP_2"
        )

        m.addConstrs(
                (r_2[idx_l, idx_d] == 0
                 for idx_l, l in enumerate(links)
                 for idx_d, d in enumerate(demands)
                ),
                name="No_QKP_1"
            )

    m.optimize()
    
    if m.status == GRB.OPTIMAL:
        #print("\n=========== OPTIMAL SOLUTION FOUND ===========")

        # Store solutions as dict of numpy arrays
        solution_all = {
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
            "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items()},
        }
        solution_filtered = {
            "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items() if abs(v.X) > 0},
            "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items() if abs(v.X) > 0},
            "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items() if abs(v.X) > 0},
            "a":   {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items() if abs(v.X) > 0}
        }

        r_total = {k: solution_all["r_1"].get(k, 0) + solution_all["r_2"].get(k, 0)
           for k in set(solution_all["r_1"]) | set(solution_all["r_2"])}

        # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
        # pp.pprint(solution_filtered)

        k_srv = sum(solution_all["r_h"][idx_d]
                    for idx_d, d in enumerate(demands)
                   )

        k_req = sum(d.K_REQ[t]
                    for idx_d, d in enumerate(demands)
                   )

        #print(f"k_req: {k_req}, k_srv: {k_srv}")

        Obj_REQ = sum(d.K_REQ[idx_d]
                     for idx_d, d in enumerate(demands)
                    )
        if m.ObjVal < 0.99 * Obj_REQ:
            reward = m.ObjVal / Obj_REQ
        else:
            reward = 10

        A_next = {k: A.get(k, 0) + solution_all["a"].get(k, 0)
           for k in set(solution_all["a"])}
    else:
        #print("No optimal solution found.")
        solution_all = {
            "a": state[t]["a"].copy()
        }
        
    return solution_all, reward, A_next









# def online(history, gss, haps, links, demand, t, delta, state):
#     # Create Optimization Model
#     m = gp.Model("hap-qkd")
    
#     ## Decision Variables
#     # Dictionaries of decision variables instead of MVar arrays
#     r_1, r_2, r_h, a, z = {}, {}, {}, {}, {}

#     demands = [demand]

#     for idx_l, l in enumerate(links):
#         for idx_d, d in enumerate(demands):
#             r_1[idx_l] = m.addVar(name=f"r_1_{idx_l}_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
#             r_2[idx_l] = m.addVar(name=f"r_2_{idx_l}_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)
#             z[idx_l]   = m.addVar(name=f"z_{idx_l}_{idx_d}",   vtype=GRB.BINARY)
                
#     for idx_d, d in enumerate(demands):
#         r_h = m.addVar(name=f"r_h_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=d.K_REQ[t] * KEY_RATE_SCALE)

#     m.ModelSense = GRB.MAXIMIZE
    
#     # Primary objective: maximize r_h
#     m.setObjectiveN(sum(sum(r_h[idx_d]
#                             for idx_d, d in enumerate(demands)
#                            )
#                        ) * syst.THETA, index=0, priority=2, weight=1.0, abstol=1e-9, reltol=1e-9, name="Primary")

#     m.setParam("MIPGap", 1e-9)          # force very tight gap
#     m.setParam("MIPGapAbs", 1e-9)
#     m.setParam("FeasibilityTol", 1e-9)
#     m.setParam("IntFeasTol", 1e-9)
#     m.setParam("OptimalityTol", 1e-9)

#     ## Constraints
#     # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
#     # r_h = min_{l:z_l=1}(r_1+r_2)
#     m.addConstrs(
#         (
#             r_h[idx_d] <= r_1[idx_l, idx_d] + r_2[idx_l, idx_d] + d.K_REQ[t] * KEY_RATE_SCALE * (1 - z[idx_l, idx_d])
#             for idx_l, l in enumerate(links)
#             for idx_d, d in enumerate(demands)
#         ), name="demand_link_coordination"
#     )

#     EPSILON = 1e-8
#     # Key rate and routing coordination (1)
#     m.addConstrs(
#         (
#             r_1[idx_l, idx_d] + r_2[idx_l, idx_d] >= EPSILON * z[idx_l, idx_d]
#             for idx_l, l in enumerate(links)
#             for idx_d, d in enumerate(demands)
#         ), name="demand_link_coordination_1"
#     )
#     # Key rate and routing coordination (2)
#     m.addConstrs(
#         (
#             r_1[idx_l, idx_d] + r_2[idx_l, idx_d, t] <= d.K_REQ * KEY_RATE_SCALE * z[idx_l, idx_d]
#             for idx_l, l in enumerate(links)
#             for idx_d, d in enumerate(demands)
#         ), name="key_rate_routing_coordination_2"
#     )
    
#     # Max Key Rate
#     m.addConstrs(
#         (
#             gp.quicksum(r_1[idx_d]
#                         for idx_d, d in enumerate(demands)
#                        ) <= l.K_MAX * KEY_RATE_SCALE
#             for idx_l, l in enumerate(links)
#         ), name="max_key_rate"
#     )
    
#     # Flow conservation
#     m.addConstrs(
#         (
#             gp.quicksum(z[idx_l, idx_d]
#                         for idx_l, l in enumerate(links)
#                         if isinstance(l.n1, gs)
#                         if gss.index(l.n1) == gss.index(d.n1)
#                        ) - gp.quicksum(z[idx_l, idx_d]
#                                        for idx_l, l in enumerate(links)
#                                        if isinstance(l.n2, gs)
#                                        if gss.index(l.n2) == gss.index(d.n1)
#                                       ) == 1
#             for idx_d, d in enumerate(demands)
#         ), name="flow_conservation_1"
#     )
#     m.addConstrs(
#         (
#             gp.quicksum(z[idx_l, idx_d]
#                         for idx_l, l in enumerate(links)
#                         if isinstance(l.n2, gs)
#                         if gss.index(l.n2) == gss.index(d.n2)
#                        ) - gp.quicksum(z[idx_l, idx_d]
#                                        for idx_l, l in enumerate(links)
#                                        if isinstance(l.n1, gs)
#                                        if gss.index(l.n1) == gss.index(d.n2)
#                                       ) == 1
#             for idx_d, d in enumerate(demands)
#         ), name="flow_conservation_2"
#     )
#     m.addConstrs(
#         (
#             gp.quicksum(z[idx_l, idx_d]
#                         for idx_l, l in enumerate(links)
#                         if  l.n1 == n
#                        ) - gp.quicksum(z[idx_l, idx_d]
#                                        for idx_l, l in enumerate(links)
#                                        if  l.n2 == n
#                                       ) == 0
#             for idx_d, d in enumerate(demands)
#             for n in gss + haps
#             if  n != d.n1 and n != d.n2
#         ), name="flow_conservation_3"
#     )
#     m.addConstrs(
#         (
#             z[idx_l_1, idx_d] + z[idx_l_2, idx_d] <= 1
#             for idx_l_1, l_1 in enumerate(links)
#             for idx_l_2, l_2 in enumerate(links)
#             if  idx_l_1 < idx_l_2 and (l_1.n1 == l_2.n2) and (l_1.n2 == l_2.n1)
#             for idx_d, d     in enumerate(demands)
#         ), name="loop_prevention"
#     )
    
#     # Maximum Tx/Rx Connection
#     m.addConstrs(
#         (
#             gp.quicksum(z[idx_l, idx_d]
#                         for idx_l, l in enumerate(links)
#                         for idx_d, d in enumerate(demands)
#                         if  l.n1 == n
#                        ) <= n.N_TX
#             for idx_n, n in enumerate(gss + haps)
#         ), name="max_tx_connections"
#     )
    
#     m.addConstrs(
#         (
#             gp.quicksum(z[idx_l, idx_d]
#                         for idx_l, l in enumerate(links)
#                         for idx_d, d in enumerate(demands)
#                         if l.n2 == n
#                        ) <= n.N_RX
#             for idx_n, n in enumerate(gss + haps)
#         ), name="max_rx_connections"
#     )
    
#     # QKP on HAPs and GSs
#     m.addConstrs(
#         (
#             a[idx_l] >= delta * gp.quicksum(r_2[idx_l, idx_d]
#                                             for idx_d, d in enumerate(demands)
#                                            ) * STORAGE_SCALE
#             for idx_l, l in enumerate(links)
#             for t        in syst.T
#         ), name="qkp_min_capacity"
#     )
    
#     m.addConstrs(
#         (
#             a[idx_l] == syst.THETA * (l.K_MAX[t] * KEY_RATE_SCALE - gp.quicksum(r_1[idx_l, idx_d, t] + r_2[idx_l, idx_d]
#                                                                                 for idx_d, d in enumerate(demands)
#                                                                                )
#                                         ) * STORAGE_SCALE
#             for idx_l, l in enumerate(links)
#             for t        in syst.T
#         ), name="qkp_sequence"
#     )
    
#     # m.addConstrs(
#     #     (
#     #         sum(a[idx_l, tp]
#     #             for tp in range(t)
#     #            ) <= min(l.n1.A_MAX, l.n2.A_MAX) * STORAGE_SCALE
#     #         for idx_l, l in enumerate(links)
#     #         for t        in syst.T
#     #     ), name="qkp_max_capacity"
#     # )

#     m.optimize()
    
#     if m.status == GRB.OPTIMAL:
#         print("\n=========== OPTIMAL SOLUTION FOUND ===========")

#         # Store solutions as dict of numpy arrays
#         solution = {
#             "r_1": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_1.items()},
#             "r_2": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_2.items()},
#             "r_h": {k: round(v.X / KEY_RATE_SCALE, 3) for k, v in r_h.items()},
#             "a": {k: round(v.X / KEY_RATE_SCALE / STORAGE_SCALE, 3) for k, v in a.items()},
#             "z": {k: v.X for k, v in z.items()}
#         }

#         # Append to history
#         history[t] = solution

#         # pp = pprint.PrettyPrinter(indent=2, width=120, sort_dicts=False)
#         # pp.pprint(solution)

#         # for idx_l, l in enumerate(links):
#         #     for t in syst.T:
#         #         print(f"K_MAX[{idx_l}][{t}]: {l.K_MAX[t]}")
        
#         #print(solution)
#     else:
#         print("No optimal solution found.")
#         solution = None
        
#     return solution, history