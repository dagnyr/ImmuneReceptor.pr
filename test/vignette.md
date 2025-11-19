# Vignette

Below is a short vignette/tutorial on using TCR for motif and hamming distance based clustering of TCRs based on their CDR3s. Before using the tool, make sure to look at the file format section to make sure all of your data is formatted and prepared correctly in advance.

# Introduction

## Formatting TCR data

**Your CDRs should be in a .csv format. Ensure the following:**

 - The column containing the **TCR chain** (e.g. TRA, TRB, TRG, TRD, etc.) should be called **"chain"**.
 -  The columns containing the vgenes, jgenes, and dgenes if present (e.g. TRGC1, TRDJ3) should be called **"v_gene"**, **"j_gene"**, and **"d_gene"** respectively.
 - The column containing the cdr3 amino acid sequences (e.g. "CALGDNLAFLRIGGAYTDKLIF", "CALGEVGRGIRAKLIF") should be called **"cdr3"**.
 - The column containing the donor or sample information (e.g. donor1, sample1, etc.) should be called **"sample"** - this is used for the HLA enrichment function.
 
 Example  of a suitable format:

| barcode | chain | v_gene | j_gene | d_gene | cdr3 |  ... |  sample |
|--|--|--|--|--|--|--|--|
| CATT | TRB | TRBV1 | TRBJ2 | TRBD3 | CATTAGF | ... | sample_1 |
| CGTA | TRB | TRBV4 | TRBJ5 | TRBD6 | CGCATAGCF | ... | sample_2 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |

As of now, we suggest looking at one type of chain at a time, ie. filter your data frame to include only TRG for example to analyze the TRGs first. Then redo for TRD and so on. Currently we do not support looking at paired TRA/TRBB and TRG/TRD chains, etc.

## Formatting HLA typing data

HLA typing data should be provided as a separate csv. For donors and samples with donor information available, provide a .csv file with the donors/samples as the first columns and the genes as the subsequent rows, as shown below with mock data:
| donor/sample | HLADPA1 | HLADPB1 | ... |
|--|--|--|--|
| sample_1 | DPA1_*01:03 | DPB1_*02:05 | ... |
| sample_2 | DPA1_*01:06 | DPB1_*02:05 | ... |
| ... | ... | ... | ... |

The titles and convention you use for naming will not impact the HLA enrichment function, as long as you are **1)** consistent with the naming convention and **2)** the donor is the first column and **3)** the subsequent columns are the HLA genes.

## Time considerations
The functions to find motifs and finding significant motifs will be the most time consuming and computationally consuming tasks - these may take longer, so just keep that in mind when you are planning your analysis. There is a progress bar to keep track of progress for these functions.


# What it actually does

## Simple Overview
Overall, this is what the package is doing:

 1. Loads your data, double checks it, and removes and counts duplicates.
 2. Find motifs present in your sample dataset based on motif lengths you define.
 3. Identifies from those motifs, motifs that are enriched in your sample compared to a reference dataset.
 4. Builds a graph based on the hamming distances and presence of shared significant motifs of your CDR sequences (either or both).
 5. Scores the clusters of TCRs identified from the graph to look at enrichment of most common CDR3 length, v-gene utilization, and optionally, HLA type enrichment.

## Comparing to a reference
Comparison to our reference dataset is used again and again in the below functions and is very important to answer the question **"is X the result of antigen specific T-cell responses, or is it occurring just by chance?"**. One way to answer this question is by having a naive T-cell reference dataset of CDR3s that are unselected and represented in a sort of unbiased way, what a truly random population of CDR3s would look like. So below, in many of the functions, *we are essentially just checking if a motif is present more so than we would expect by chance, is the modal CDR3 length longer than we would expect just by chance*, and so on - it is a key aspect of the package's function.

# How to use it

Firstly, you will need your sample dataset (the CSV containing the TCR sequencing information, like the CDR3s and TCR chains, etc.), your own reference dataset (or you can use ours), and **optionally** a file with HLA typing information - then you can follow along as shown below:

## Loading the data

    # load the dataframe with your sample data using load_cdr3s()
    # this function double checks everything is looking good + loads data
    df = load_cdr3s(/path/to/your/data.csv)
    
    # optionally, if you have it, load your HLA typing here too:
    hla_df = CSV.read(/path/to/your/hla_data.csv, DataFrame)

## Finding motifs

    ## get_motifs() takes your cdrs, chops them into motifs + counts them
    ## you will define a minimum length and maximum length of motifs
    all_motifs = get_motifs(list_of_cdrs, min_length, max_length)

    # find_significant_motifs() finds motifs enriched in your dataset compared to the reference dataset
    # you can use your own reference or the ones we used
    sig_motifs = find_significant_motifs(all_motifs, sample_cdrs, reference_cdrs)

## Making a graph

**Global Edges** -> Edges defined based on hamming distance (e.g. draw an edge between CDR3s with very similar global sequence).
**Local Edges** -> Edges based on presence of shared significant motifs (e.g. draw an edge between CDR3s with shared motifs).

    # using Graphs.jl and MetaGraphsNext.jl, we make a graph
    # it is made using the cdr3 dataframe and significant motifs
    # set "isglobal" as true if you want global edges
    # set "islocal" as true if you want local edges
    graph = make_edges(df, sig_motifs, isglobal = true, islocal = true)


# Scoring the clusters

After the graphing, we will have clusters of TCRs/CDRs based on the global and local edges combined. Using several functions, we can then score the enrichment of the **most common CDR3 length** in the cluster as compared to the reference dataset, the **most enriched v-genes** in the cluster compared to the overall sample dataset, and optionally, the **enriched HLA genes** in each cluster as compared to the overall sample dataset.

Note that for many of these functions you will be defining a **sim_depth**. In these functions, at some point, there is repeated resampling of either your sample dataset or the reference dataset to create fake "clusters" to compare to the real clusters (found earlier), to see if length, v-genes, etc. are enriched in your cluster as compared to those reference clusters.

## Length P-Values
**find_length_pvals()** finds the most common length per cluster and then compares this to the frequency in the reference dataset, then providing a p-value.

    # finds the most common length per cluster + compares to reference
    # you define the graph, reference, and simulation depth
    length_pvals = find_length_pvals(graph, reference_cdrs, sim_depth)

## V-Gene P-Vals
**score_vgene()** compares your clusters to the overall v-gene utilization in your sample, providing adjusted p-vals for each v-gene present in each cluster.

    # define a sim depth and provide your graph
    vgene_pvals = score_vgene(graph, sim_depth)

## HLA P-Vals
**score_hla()** compares how enriched each HLA allele is for every HLA gene in your clusters compared to the overall dataset using a hypergeometric test, providing a list of HLA genes and matching p-values for each cluster as an output.

    # make sure you look at how to prepare the hla_df above
    hla_pvals = score_hla(graph, hla_df)
