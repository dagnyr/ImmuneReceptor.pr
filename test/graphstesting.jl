# testing for graphs

g = MetaGraph(
    Graph();
    label_type=Symbol,
    vertex_data_type=Vector{String},  # test storing multiple labels?
    edge_data_type=Tuple{String,Int},
)

# plans?
# labels will be cdr3s
# vertices will have vector data storing v-gene, j-gene, d-gene (None if no gene), patient,
# edges will store label (global -> distance, local -> motif) and p-val
# maybe store vector of dictionary with String, Any so I can store label followed by thing?
# vertices will store

vertex_data_type = Dict{Symbol,Any}
