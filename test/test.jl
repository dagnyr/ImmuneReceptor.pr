# load data
# ensure that the column w/ the cdr3 is entitled "cdr3" and the TCR chain column is entitled "chain"

df = load_cdr3s(/path/to/your/data.csv)
list_of_cdrs = df.cdr3s

# create list of motifs

## first, find every single motif based on your sample cdr3s.
## you define a minimum length and maximum length interval for motifs
all_motifs = get_motifs(list_of_cdrs, min_length, max_length)
sig_motifs = find_significant_motifs(all_motifs, sample_cdrs, reference_cdrs)

# make graph

# make a graph, defining your significant motifs as identified above.
# define using isglobal and islocal if you want edges based on hamming distance (global edges) or just shared significant motifs (local edges)
graph = make_edges(df, sig_motifs, isglobal = true, islocal = true)

# scoring graphs

## provide the graph generated and then the reference CDR3s to compare length distributions
## define a simulation depth for randomly drawing CDRs from the reference dataset for empirical p-val generation

length_pvals = find_length_pvals(graph, reference_cdrs, sim_depth)

## define the graph previously created and a simulation depth (e.g. 1000)
## hypergeometric test is used when there are < 200 members in a cluster, empirical p-val based on sims is done when > 200

vgene_pvals = score_vgene(graph, sim_depth)

## provide the graph and a dataframe
## ensure that your original dataset has a "sample" variable for your samples/donors
## the hla_df file should have your samples as the first column and the alleles as the subsequent columns

hla_pvals = score_hla(graph, hla_df)

# maybe make function to make heat maps from the outputs of these?

# generate summary

# get the graph, define:
# members per cluster
# donors per cluster
# length mode
# length hla_pvals
# vgene pval and matched vgene
# all significant HLAs for the cluster as a list
