# The format of simpleweightedgraph files is as follows: for each graph,
# a one line header: "LightGraphs.SimpleWeightedGraph", <num_vertices>, <num_edges>, {"d" | "u"}, <name>[, <ver>, <vdatatype>, <wdatatype>, <graphcode>]
#   - "LightGraphs.SimpleWeightedGraph" is a fixed string
#   - num_vertices is an integer
#   - num_edges is an integer
#   - "d" for directed graph, "u" for undirected. Note that this
#       option does not perform any additional edge construction; it's
#       merely used to return the correct type of graph.
#   - name is a string
#   - ver is an int
#   - vdatatype is a string ("UInt8", etc.)
#   - wdatatype is a string describing the data type of the weights
#   - graphcode is a string.
# header followed by a list of (comma-delimited) edges - src,dst,weight
# Multiple graphs may be present in one file.

const FIXEDSTR = "LightGraphs.SimpleWeightedGraph"

struct SWGFormat <: AbstractGraphFormat end

struct SWGHeader
    nv::Int
    ne::Int
    is_directed::Bool
    name::String
    ver::Int
    vdtype::DataType    # vertex data type
    wdtype::DataType    # weight data type
    code::String
end
function show(io::IO, h::SWGHeader)
    isdir = h.is_directed ? "d" : "u"
    print(io, "$FIXEDSTR,$(h.nv),$(h.ne),$isdir,$(h.name),$(h.ver),$(h.vdtype),$(h.wdtype),$(h.code)")
end

function _swg_read_one_graph(f::IO, header::SWGHeader)
    T = header.vdtype
    U = header.wdtype
    if header.is_directed
        g = SimpleWeightedDiGraph{T, U}(header.nv)
    else
        g = SimpleWeightedGraph{T, U}(header.nv)
    end
    for i = 1:header.ne
        line = chomp(readline(f))
        if length(line) > 0
            src_s, dst_s, weight_s = split(line, r"\s*,\s*")
            src = parse(T, src_s)
            dst = parse(T, dst_s)
            weight = parse(U, weight_s)
            add_edge!(g, src, dst, weight)
        end
    end
    return g
end

function _swg_skip_one_graph(f::IO, n_e::Integer)
    for _ in 1:n_e
        readline(f)
    end
end

function _parse_header(s::AbstractString)
    addl_info = false
    fixedstr, nvstr, nestr, dirundir, graphname, _ver, _vdtype, _wdtype, graphcode  = split(s, r"s*,s*", limit=9)
    fixedstr != FIXEDSTR && error("Error parsing header.")

    n_v = parse(Int, nvstr)
    n_e = parse(Int, nestr)
    dirundir = strip(dirundir)
    directed = !(dirundir == "u")
    graphname = strip(graphname)
    ver = parse(Int, _ver)
    vdtype = eval(Symbol(_vdtype))
    wdtype = eval(Symbol(_wdtype))

    return SWGHeader(n_v, n_e, directed, graphname, ver, vdtype, wdtype, graphcode)
end

"""
    loadswg_mult(io)

Return a dictionary of (name=>graph) loaded from IO stream `io`.
"""
function loadswg_mult(io::IO)
    graphs = Dict{String,AbstractGraph}()
    while !eof(io)
        line = strip(chomp(readline(io)))
        if !(startswith(line, "#") || line == "")
            header = _parse_header(line)
            g = _swg_read_one_graph(io, header)
            graphs[header.name] = g
        end
    end
    return graphs
end

function loadswg(io::IO, gname::String)
    while !eof(io)
        line = strip(chomp(readline(io)))
        (startswith(line, "#") || line == "") && continue
        header = _parse_header(line)
        if gname == header.name
            return _swg_read_one_graph(io, header)
        else
            _swg_skip_one_graph(io, header.ne)
        end
    end
    error("Graph $gname not found")
end

"""
    saveswg(io, g, gname)

Write a graph `g` with name `gname` in a proprietary format
to the IO stream designated by `io`. Return 1 (number of graphs written).
"""
function saveswg(io::IO, g::AbstractGraph, gname::String)
    header = SWGHeader(nv(g), ne(g), is_directed(g), gname, 1, eltype(g), weighttype(g), "simpleweightedgraph")
    # write header line
    line = string(header)
    write(io, "$line\n")
    # write edges
    for e in edges(g)
        write(io, "$(src(e)),$(dst(e)),$(weight(e))\n")
    end
    return 1
end

"""
    saveswg_mult(io, graphs)

Write a dictionary of (name=>graph) to an IO stream `io`,
with default `GZip` compression. Return number of graphs written.
"""
function saveswg_mult(io::IO, graphs::Dict)
    ng = 0
    for (gname, g) in graphs
        ng += saveswg(io, g, gname)
    end
    return ng
end


loadgraph(io::IO, gname::String, ::SWGFormat) = loadswg(io, gname)
loadgraphs(io::IO, ::SWGFormat) = loadswg_mult(io)
savegraph(io::IO, g::AbstractGraph, gname::String, ::SWGFormat) = saveswg(io, g, gname)
savegraph(io::IO, g::AbstractGraph, ::SWGFormat) = saveswg(io, g, "graph")
savegraph(io::IO, d::Dict, ::SWGFormat) = saveswg_mult(io, d)