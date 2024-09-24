
import numpy as np
import csv

def state_indices_to_merge(means, n_states_keep):
    print(means)
    means_sorted_indices = np.argsort(means)
    state_indices_to_merge = np.full(means.shape, False)
    if n_states_keep > len(means) - 1:
        return state_indices_to_merge
    for i in range(len(means) - n_states_keep + 1):
        state_indices_to_merge[means_sorted_indices[i]] = True
    return state_indices_to_merge

def merge_states(means, vars, transition_matrix, stat_probs, n_states_keep, betas_list):
    merge_ind = state_indices_to_merge(means, n_states_keep)
    print(merge_ind)
    new_means = np.zeros(n_states_keep)
    new_vars = np.zeros(n_states_keep)
    new_stat_probs = np.zeros(n_states_keep)
    alphas = np.zeros(n_states_keep)
    new_transition_matrix = np.zeros((n_states_keep, n_states_keep))
    tm_merged_cols = np.zeros((len(means), n_states_keep))
    new_betas_list = []
    for elem in betas_list:
        new_betas_list.append(np.zeros(n_states_keep))
    ind_new = 0
    merge_mean = 0
    merge_var = 0
    merge_alpha = 0
    merge_stat_prob = 0
    merge_state_ind = n_states_keep - 1
    for i in range(len(means)):
        if merge_ind[i]:
            merge_mean = max(merge_mean, means[i])
            merge_var = max(merge_var, vars[i])
            merge_alpha = max(merge_alpha, means[i])
            merge_stat_prob += stat_probs[i]
            for b_i in range(len(betas_list)):
                new_betas_list[b_i][merge_state_ind] += betas_list[b_i][i]
        else:
            new_means[ind_new] = means[i]
            new_vars[ind_new] = vars[i]
            new_stat_probs[ind_new] = stat_probs[i]
            for b_i in range(len(betas_list)):
                new_betas_list[b_i][ind_new] = betas_list[b_i][i]
            ind_new += 1
    merge_state_ind = n_states_keep - 1
    alphas[merge_state_ind] = merge_alpha
    new_means[merge_state_ind] = merge_mean
    new_vars[merge_state_ind] = merge_var
    new_stat_probs[merge_state_ind] = merge_stat_prob
    # transition matrix, first step merge columns in each row - just add upp
    for r in range(len(means)):
        ind_new_c = 0
        for c in range(len(means)):
            if merge_ind[c]:
                tm_merged_cols[r][merge_state_ind] += transition_matrix[r][c]
            else:
                tm_merged_cols[r][ind_new_c] += transition_matrix[r][c]
                ind_new_c += 1
    # transition matrix, second step, combine rows, scale with normalized stat probs
    ind_new_r = 0
    for r in range(len(means)):
        if merge_ind[r]:
            for c in range(n_states_keep):
                new_transition_matrix[merge_state_ind][c] += stat_probs[r]*tm_merged_cols[r][c]/merge_stat_prob
        else:
            for c in range(n_states_keep):
                new_transition_matrix[ind_new_r][c] += tm_merged_cols[r][c]
            ind_new_r += 1
    return (new_means, new_vars, alphas, new_stat_probs, new_transition_matrix, new_betas_list)



if __name__=='__main__':
    nStates = 8
    transitionMatrix = np.zeros((nStates,nStates))
    with open('../../data/models/cont/transitionMatrix.txt', 'r') as csvfile:
        reader = csv.reader(csvfile)
        r = 0
        for row in reader:
            rowElems = row[0].split(sep = ' ')
            for c in range(nStates):
                transitionMatrix[r][c] = float(rowElems[c])
            r += 1

    statProbs = np.zeros(nStates)
    with open('../../data/models/cont/stationaryDistr.txt', 'r') as csvfile:
        reader = csv.reader(csvfile)
        r = 0
        for row in reader:
            rowElems = row[0].split(sep = ' ')
            statProbs[r] = float(rowElems[0])
            r+= 1

    means = np.zeros(nStates)
    variances = np.zeros(nStates)
    with open('../../data/models/cont/normalParams.txt', 'r') as csvfile:
        reader = csv.reader(csvfile)
        r = 0
        for row in reader:
            rowElems = row[0].split(sep = ' ')
            means[r] = float(rowElems[0])
            variances[r] = pow(float(rowElems[1]), 2)
            r+= 1

    merged_means, merged_vars, alphas, merged_stat_probs, merged_transition_mat, betas_list = \
        merge_states(means, variances, transitionMatrix, statProbs, 2, [])
    print(merged_means)
    print(merged_vars)
    print(alphas)
    print(merged_stat_probs)
    print(merged_transition_mat)
 