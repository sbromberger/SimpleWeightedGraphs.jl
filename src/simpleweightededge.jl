import Base: Pair, Tuple, show, ==
import LightGraphs: AbstractEdge, src, dst, reverse

abstract type AbstractSimpleWeightedEdge{T} <: AbstractEdge{T} end

struct SimpleWeightedEdge{T<:Integer, U<:Real} <: AbstractSimpleWeightedEdge{T}
    src::T
    dst::T
    weight::U
end

SimpleWeightedEdge(t::NTuple{2}) = SimpleWeightedEdge(t[1], t[2], one(Float64))
SimpleWeightedEdge(t::NTuple{3}) = SimpleWeightedEdge(t[1], t[2], t[3])
SimpleWeightedEdge(p::Pair) = SimpleWeightedEdge(p.first, p.second, one(Float64))
SimpleWeightedEdge{T, U}(p::Pair) where T<:Integer where U <: Real = SimpleWeightedEdge(T(p.first), T(p.second), one(U))
SimpleWeightedEdge{T, U}(t::NTuple{3}) where T<:Integer where U <: Real = SimpleWeightedEdge(T(t[1]), T(t[2]), U(t[3]))
SimpleWeightedEdge{T, U}(t::NTuple{2}) where T<:Integer where U <: Real = SimpleWeightedEdge(T(t[1]), T(t[2]), one(U))
SimpleWeightedEdge(x, y) = SimpleWeightedEdge(x, y, one(Float64))
eltype(e::T) where T<:AbstractSimpleWeightedEdge= eltype(src(e))

# Accessors
src(e::AbstractSimpleWeightedEdge) = e.src
dst(e::AbstractSimpleWeightedEdge) = e.dst
weight(e::AbstractSimpleWeightedEdge) = e.weight

# I/O
show(io::IO, e::AbstractSimpleWeightedEdge) = print(io, "Edge $(e.src) => $(e.dst) with weight $(e.weight)")

# Conversions
Tuple(e::AbstractSimpleWeightedEdge) = (src(e), dst(e), weight(e))

(::Type{SimpleWeightedEdge{T, U}})(e::AbstractSimpleWeightedEdge) where {T <: Integer, U <: Real} = SimpleWeightedEdge{T, U}(T(e.src), T(e.dst), U(e.weight))

# Convenience functions - note that these do not use weight.
reverse(e::T) where T<:AbstractSimpleWeightedEdge = T(dst(e), src(e), weight(e))
==(e1::AbstractSimpleWeightedEdge, e2::AbstractSimpleWeightedEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2))
==(e1::AbstractSimpleWeightedEdge, e2::AbstractEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2))
==(e1::AbstractEdge, e2::AbstractSimpleWeightedEdge) = (src(e1) == src(e2) && dst(e1) == dst(e2))