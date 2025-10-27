# testing for graphs

g = MetaGraph(
    Graph();
    label_type=Symbol,
    vertex_data_type=Vector{String},  # test storing multiple labels?
    edge_data_type=Vector{String},
)

# plans?
# labels will be cdr3s
# vertices will have vector data storing v-gene, j-gene, d-gene (None if no gene), patient,
# edges will store label (global -> distance, local -> motif) and p-val
# maybe store vector of dictionary with String, Any so I can store label followed by thing?
# vertices will store

#vertex_data_type = Dict{Symbol,Any}

add_vertex!(g, :Owo, ["hey", "hello", "hi"])
add_vertex!(g, :Uwu, ["goodbye", "seeya", "bye"])

add_edge!(g, :Uwu, :Owo, ["tidings", "bonjour", "nihao"])

edge_data = g[:Uwu, :Owo]
push!(edge_data, "shalom")
g[:Uwu, :Owo] = edge_data

g[:Uwu, :Owo]


# now trying with dictionary
# Dict{Symbol,Any}

g = MetaGraph(
    Graph();
    label_type=Symbol,
    vertex_data_type=Dict{Symbol,Any},  # test storing multiple labels?
    edge_data_type=Dict{Symbol,Any},
)

add_vertex!(g, :Owo, Dict(:label => "ACGATAG", :vgene => "V13B")) # cdr3 and gene is fake for testing
add_vertex!(g, :Uwu, Dict(:label => "CCGGTAT", :vgene => "V12A"))

add_edge!(g, :Uwu, :Owo, Dict(:type => "local", :pval => 0.02))
add_edge!(g, :Uwu, :Owo, Dict(:type => "global", :distance => 2)) # replaces prior edge

# solution? add a vector of dicionaries to edges and have a global distance, global p-val, motif and motif p-val entry.
# to determine if theres a global or local edge use bool to check for each edge?
