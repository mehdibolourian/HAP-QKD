from libraries import *

def online_new(gss, haps, links, state, action, t, f_qkp, k_srv_heur1, k_srv_heur2, demand_active, link_active):
    demands = state[t]["demands"]
    z       = action

    #print(z)

    if f_qkp:
        A       = state[t]["a"]
        A_next  = A.copy()

    # print(f"z: {z}")
    # print(f"A: {A}")

    reward  = -100
    
    # Create Optimization Model
    m = gp.Model("hap-qkd")
    
    ## Decision Variables
    # Dictionaries of decision variables instead of MVar arrays
    r_1, r_2, r_h, a, o = {}, {}, {}, {}, {}

    for idx_d, d in enumerate(demands):
        # NEW: Check if this specific demand is active
        # demand_active is np.ones(MAX_DEMANDS), so we check the idx_d
        is_active = (demand_active[idx_d][t] == 1)
        
        # If inactive, we'll force ub to 0.0 later or set it here
        upper_bound = d.K_REQ[t] * KEY_RATE_SCALE if is_active else 0.0
        
        r_h[idx_d] = m.addVar(name=f"r_h_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
        
        for idx_l, l in enumerate(links):
            r_1[idx_l, idx_d] = m.addVar(name=f"r_1_{idx_l}_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)
            r_2[idx_l, idx_d] = m.addVar(name=f"r_2_{idx_l}_{idx_d}", vtype=GRB.CONTINUOUS, lb=0.0, ub=upper_bound)

    K_MAX_SUM = sum(l.K_MAX[tp]
                    for tp in range(t+1)
                   )
    for idx_l, l in enumerate(links):
        a[idx_l] = m.addVar(name=f"a_{idx_l}", vtype=GRB.CONTINUOUS, lb=-K_MAX_SUM, ub=l.K_MAX[t] * KEY_RATE_SCALE * syst.THETA)

    nodes = gss + haps

    ## Node order variable --> To prevent subtours
    for idx_n, n in enumerate(nodes):
        for idx_d, d in enumerate(demands):
            o[idx_n, idx_d] = m.addVar(name=f"o_{idx_n}_{idx_d}", vtype=GRB.CONTINUOUS)

    m.setObjective(sum(r_h[idx_d]
                       for idx_d, d in enumerate(demands)
                       if demand_active[idx_d][t] == 1
                      ) * syst.THETA, GRB.MAXIMIZE
                  )

    m.Params.MIPGap = 0.001     # 1% optimality
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
        if demand_active[idx_d][t] == 0:
            m.addConstr(r_h[idx_d] == 0, name=f"mask_rh_{idx_d}")
            for idx_l in range(len(links)):
                m.addConstr(r_1[idx_l, idx_d] == 0, name=f"mask_r1_{idx_l}_{idx_d}")
                m.addConstr(r_2[idx_l, idx_d] == 0, name=f"mask_r2_{idx_l}_{idx_d}")

    ## --- NEW: Explicitly zero out flows for inactive links ---
    for idx_l, l in enumerate(links):
        if link_active[idx_l][t] == 0:
            # Force storage growth to zero for this link
            m.addConstr(a[idx_l] == 0, name=f"mask_a_link_{idx_l}")
            
            for idx_d in range(len(demands)):
                # Force all flow (direct and from storage) to zero for this link-demand pair
                m.addConstr(r_1[idx_l, idx_d] == 0, name=f"mask_r1_link_{idx_l}_d_{idx_d}")
                m.addConstr(r_2[idx_l, idx_d] == 0, name=f"mask_r2_link_{idx_l}_d_{idx_d}")

    
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

    # MTZ subtour elimination --> Eliminates pointless single/multi-hop loops in the flows --> Uses an ordering values for all nodes
    # --> The order values should only increase on the path --> A decrease in order value == a loop (X)
    M = len(nodes)
    m.addConstrs(
        (
            o[nodes.index(l.n2), idx_d] >= o[nodes.index(l.n1), idx_d] + 1 - M * (1 - z[idx_d][idx_l])
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
        ), name="ordering_constraint_1"
    )
    m.addConstrs(
        (
            o[nodes.index(d.n1), idx_d] == 1
            for idx_d, d in enumerate(demands)
        ), name="ordering_constraint_2"
    )
    m.addConstrs(
        (
            o[nodes.index(d.n2), idx_d] == M
            for idx_d, d in enumerate(demands)
        ), name="ordering_constraint_2"
    )

    # Demand-level and link-level key rate coordination (Note that r_h is a part of the maximization objective)
    m.addConstrs(
        (
            r_h[idx_d] >= r_1[idx_l, idx_d] + r_2[idx_l, idx_d]
            for idx_l, l in enumerate(links)
            for idx_d, d in enumerate(demands)
        ), name="demand_link_coordination_1"
    )

    # Key rate and routing coordination (2)
    m.addConstrs(
        (
            r_1[idx_l, idx_d] <= d.K_REQ[t] * KEY_RATE_SCALE * z[idx_d][idx_l] #z[idx_l]
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
                a[idx_l] == syst.THETA * (l.K_MAX[t] * sum(z[idx_d][idx_l]
                                                           for idx_d, d in enumerate(demands)
                                                          ) * KEY_RATE_SCALE - gp.quicksum(r_1[idx_l, idx_d] + r_2[idx_l, idx_d]
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

        k_srv = sum(solution_all["r_h"][idx_d] * syst.THETA
                    for idx_d, d in enumerate(demands)
                    if demand_active[idx_d][t] == 1
                   )

        k_req = sum(d.K_REQ[t] * syst.THETA
                    for idx_d, d in enumerate(demands)
                    if demand_active[idx_d][t] == 1
                   )

        # if k_srv_heur != 0:
        #     print(f"k_req: {k_req}, k_srv: {k_srv}, k_srv_heur: {k_srv_heur}")

        Obj_REQ = sum(d.K_REQ[t]
                     for idx_d, d in enumerate(demands)
                    ) * syst.THETA
        # if m.ObjVal < 0.99 * Obj_REQ:
        #reward = m.ObjVal / Obj_REQ
        # else:
        #     reward = 10
        ######################

        # reward = 0
        # for idx_d, d in enumerate(demands):
        #     if solution_all["r_h"][idx_d] >= 0.95 * d.K_REQ[t]:
        #         reward += 1
        #     else:
        #         reward += -1 + solution_all["r_h"][idx_d] / d.K_REQ[t]

        # reward = reward / len(demands)
        ###################

        # reward = 0
        # if k_srv >= 0.95 * k_req:
        #     reward = 1
        # else:
        #     reward = -1 + k_srv / k_req

        #print(f"k_srv: {k_srv}, k_srv_heur: {k_srv_heur}, k_req: {k_req}")

        #######################
        # ratio_agent = k_srv      / (k_req + 1e-8)
        # ratio_heur  = k_srv_heur / (k_req + 1e-8)
    
        # if ratio_agent >= 0.99 or k_req == 0:
        #     reward = 1.0
        # elif ratio_agent <= 1e-6:
        #     reward = -1.0
        # else:
        #     reward = (ratio_agent - ratio_heur)
        #     reward = max(-1.0, min(1.0, reward))

        ratio_agent = k_srv / (k_req + 1e-8)

        ratio_h1 = k_srv_heur1 / (k_req + 1e-8)
        ratio_h2 = k_srv_heur2 / (k_req + 1e-8)
        
        h_low  = min(ratio_h1, ratio_h2)
        h_high = max(ratio_h1, ratio_h2)
        
        if k_req == 0:
            reward = 0.0 * 6
        
        elif ratio_agent >= 0.95:
            reward = 1.0 * 6
        
        elif ratio_agent >= h_high:
            reward = (0.6 + 0.4 * (ratio_agent - h_high)) * 6
        
        elif ratio_agent >= h_low:
            reward = (0.2 + 0.4 * (ratio_agent - h_low) / (h_high - h_low + 1e-8)) * 6
        
        else:
            reward = (-1 * (h_low - ratio_agent)) * 6

        if f_qkp:
            A_next = {k: A.get(k, 0) + solution_all["a"].get(k, 0)
               for k in set(solution_all["a"])}
        else:
            A_next = {k: solution_all["a"].get(k, 0)
               for k in set(solution_all["a"])}
    else:
        #print("No optimal solution found.")
        solution_all = {
            "a": state[t]["a"].copy()
        }
        
    return solution_all, reward, A_next

def calculate_ppo_reward(m, demands, t, syst):
    """
    m: The solved Gurobi model for the current timestep
    demands: List of demand objects
    t: Current timestep
    syst: System parameters
    """
    # 1. Calculate the Target (The "Goal")
    # We use a per-demand satisfaction check to give the agent granular feedback
    total_req = sum(d.K_REQ[t] for d in demands) * syst.THETA
    total_achieved = m.ObjVal
    
    # 2. Individual Satisfaction Check (Granular Feedback)
    # Note: This requires you to have the solved values for each demand from Gurobi
    # Assuming z_vars or similar are accessible to check flow per demand
    satisfied_count = 0
    for d in demands:
        # Check if this specific demand met its threshold
        # You'll need to pull the specific flow value for demand 'd' from your model
        achieved_d = get_flow_for_demand(m, d) 
        if achieved_d >= d.K_REQ[t] * syst.THETA:
            satisfied_count += 1
            
    # 3. Construct the Reward Signal
    # Component A: Throughput Ratio (0.0 to 1.0)
    throughput_ratio = total_achieved / total_req if total_req > 0 else 1.0
    
    # Component B: Success Bonus 
    # Encourages satisfying as many GS-GS pairs as possible
    success_bonus = satisfied_count / len(demands)
    
    # Component C: The Penalty
    # We only penalize if we are significantly below the target
    penalty = 0
    if total_achieved < 0.95 * total_req:
        penalty = -0.5 * (1.0 - throughput_ratio) # Scales with how much we missed
        
    # Final Reward Calculation
    # We weigh success heavily to drive the agent toward satisfying demands
    reward = (0.3 * throughput_ratio) + (0.7 * success_bonus) + penalty
    
    return reward







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