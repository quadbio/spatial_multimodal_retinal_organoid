import scvelo as scv
import cellrank as cr
import numpy as np
import pandas as pd

adata = scv.read("data/scRNAseq/h5ad/onlyRetina_timecourse.h5ad")
scv.tl.recover_dynamics(adata, n_jobs=50)

# velocity kernel
vk = cr.tl.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()
# pseudotime kernel
pk = cr.tl.kernels.PseudotimeKernel(adata, time_key = 'velocity_pseudotime')
pk.compute_transition_matrix(threshold_scheme='hard')
# velocity-pseudotime hybrid kernel
combined_kernel = 0.5 * vk + 0.5 * pk

# obsorption probability with the hybrid kernel
g = cr.tl.estimators.GPCCA(combined_kernel, write_to_adata=False)
terminal_states = adata.obs['annot_ct_refine']
terminal_states[terminal_states.isin(["Immature rods","Immature cones","Immature BC","PR/BC-PC","Late PR/BC-IPC","Early PR/BC-IPC","AC/HC","Immature IN","RPC","Prolif. RPC","PC-MG"])] = np.NAN
terminal_states = terminal_states.cat.add_categories('BC')
terminal_states[terminal_states.isin(["BC-ON","BC-OFF"])] = "BC"
terminal_states = terminal_states.cat.remove_unused_categories()
g.set_terminal_states(terminal_states)
g.compute_absorption_probabilities(n_jobs=60)

# output
pd.DataFrame(adata.obsm['to_terminal_states'], columns = adata.obsm['to_terminal_states'].names, index = adata.obs.index).to_csv("data/scRNAseq/cellrank_lineages_absorption_prob.tsv", sep="\t")
adata.obs['terminal_states_probs'] = np.NAN
adata.write("data/scRNAseq/h5ad/onlyRetina_timecourse.h5ad")
