    if tag is "UID":
        mutations = np.zeros(1)
    else:
        mutations = np.zeros(2)

    n_molecule = len(nsamples_per_molecule)
    cut_off = 0
    ns_max =  np.max(nsamples_per_molecule)
    control_mean = ns_max * (cycles-np.floor(np.log2(ns_max))+2) * bases_per_amplicon * error_rate
    while n_molecule * poisson.pmf(cut_off, control_mean) > 1e-4:
        cut_off += 1
    cut_off = max(cut_off,2)
    print("cut_off:",cut_off, "iteration:",n_molecule, "control_mean:",control_mean)

    for i, ns in enumerate(nsamples_per_molecule):
        if i % int(n_molecule/20) == 0:
            print("Processing:", i * 100 / n_molecule, "%")
        # zero_control_mean = ns * (cycles-np.floor(np.log2(ns))+2) * bases_per_amplicon * error_rate
        # if np.random.poisson(zero_control_mean) == 0:
        #     if tag is "UID":
        #         UIDs = translate(sampleUID(cycles, ns), cycles)
        #         mutation_count = np.zeros_like(UIDs, dtype=np.int16)
        #         mutations = np.vstack((mutations, np.hstack( [UIDs, mutation_count] )))
        #     else:
        #         ids = sampleUID(cycles, ns)
        #         is_left_side = ids < 2**(cycles-1)
        #         is_left_side = is_left_side.reshape((ns, 1))
        #         WBCs = translate_well(ids, cycles, tag[i])
        #         mutation_count = np.zeros_like(WBCs, dtype=np.int16)
        #         mutations = np.vstack(( mutations, np.hstack( [WBCs, is_left_side, mutation_count] )))                   
        
        # else:
        #     ids = sampleUID(cycles, ns)
        #     ancestors = [np.array(list(set(np.floor(ids / (2**c)))), dtype=np.int32) for c in range(cycles)]
        #     counts = [len(x) for x in ancestors]
        #     true_mean = np.sum(counts) * bases_per_amplicon * error_rate
        #     control_pval = np.array([np.exp(-zero_control_mean), 1.0-np.exp(-zero_control_mean)], dtype=np.float64)
        #     poisson_pval = np.array([poisson.pmf(i, true_mean) for i in range(cut_off)], dtype=np.float64)
        #     poisson_pval[-1] = 1.0 - np.sum(poisson_pval[:-1])
        #     adjust_pval = np.array([(poisson_pval[0] - control_pval[0]) / control_pval[1]] +\
        #                            [poisson_pval[i] / control_pval[1] for i in range(cut_off) if i > 0], dtype=np.float64)

        #     rand = np.random.choice(cut_off, p=adjust_pval)
        #     if rand == 0:
        #         if tag is "UID":
        #             UIDs = translate(ids, cycles)
        #             mutation_count = np.zeros_like(UIDs, dtype=np.int16)
        #             mutations = np.vstack((mutations, np.hstack( [UIDs, mutation_count] )))
        #         else:
        #             is_left_side = ids < 2**(cycles-1)
        #             is_left_side = is_left_side.reshape((ns, 1))
        #             WBCs = translate_well(ids, cycles, tag[i])
        #             mutation_count = np.zeros_like(WBCs, dtype=np.int16)
        #             mutations = np.vstack(( mutations, np.hstack( [WBCs, is_left_side, mutation_count] )))
        #     else:
        #         # print(i, rand)
        #         locs = np.random.choice(len(counts), size=(rand,), p=np.array(counts,dtype=np.float32) / np.sum(counts))
        #         aids = [np.random.choice(ancestors[loc]) for loc in locs]
        #         mutation_loc = [[aid, loc, 1, int(aid % 2 == 0)] for loc, aid in zip(locs,aids)]
        #         for i in range(rand):
        #             while True:
        #                 update_list = mutation_loc[i]
        #                 if update_list[1] == 0:
        #                     if not isinstance(update_list[0], list):
        #                         update_list[0] = [update_list[0]]
        #                     break
        #                 if update_list[2] == 1:
        #                     loc = update_list[1]
        #                     naid = update_list[0] * 2 + update_list[3]
        #                     if naid not in ancestors[loc-1]:
        #                         mutation_loc[i] = []
        #                         break
        #                     else:
        #                         update_list[0] = [naid]
        #                         update_list[1] -= 1
        #                         update_list[2] = 2
        #                 else:
        #                     loc = update_list[1]
        #                     naid = [aid * 2 for aid in update_list[0] if aid*2 in ancestors[loc-1]] + \
        #                         [aid*2 + 1 for aid in update_list[0] if aid*2+1 in ancestors[loc-1]]
        #                     if len(naid) == 0:
        #                         mutation_loc[i] = []
        #                         break
        #                     update_list[0] = naid
        #                     update_list[1] -= 1
                
        #         mutation_count = np.zeros((ns, 2), dtype=np.int16)
        #         for ml in mutation_loc:
        #             if len(ml) == 0:
        #                 continue
        #             else:
        #                 if ml[2] == 1:
        #                     for aid in ml[0]:
        #                         mutation_count[np.where(ids == aid)[0][0], ml[3]] += 1
        #                 else:
        #                     for aid in ml[0]:
        #                         mutation_count[np.where(ids == aid)[0][0], :] += 1
                
        #         if tag is "UID":
        #             UIDs = translate(ids, cycles)
        #             mutations = np.vstack((mutations, np.hstack( [UIDs, mutation_count] )))
        #         else:
        #             is_left_side = ids < 2**(cycles-1)
        #             is_left_side = is_left_side.reshape((ns, 1))
        #             WBCs = translate_well(ids, cycles, tag[i])
        #             mutations = np.vstack(( mutations, np.hstack( [WBCs, is_left_side, mutation_count] )))

        if tag is "UID":
            UIDs = sampleUID(cycles, ns).reshape((ns, 1))
            mutations = np.vstack(( mutations, UIDs ))
        else:
            ids = sampleUID(cycles, ns).reshape((ns, 1))
            is_left_side = ids < 2**(cycles-1)
            is_left_side = is_left_side.reshape((ns, 1))
            mutations = np.vstack(( mutations, np.hstack( [ids, is_left_side] )))
    
    # if tag is "UID":
    #     mutations = pd.DataFrame( mutations[1:,], columns= ["leftUID","rightUID","UID_left_mutation","UID_right_mutation"] )
    # else:
    #     mutations = pd.DataFrame( mutations[1:,], columns= ["leftWBC","rightWBC","is_left_side","WBC_left_mutation","WBC_right_mutation"] )
    if tag is "UID":
        mutations = pd.DataFrame( mutations[1:,], columns= ["UID"] )
    else:
        mutations = pd.DataFrame( mutations[1:,], columns= ["WBC","is_left_side"] )

    mutations = mutations.reset_index(drop=True)
    return mutations
