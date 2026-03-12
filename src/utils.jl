#
# Copyright (c) 2023 Tobias Thummerer, Lars Mikelsons
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

function isZeroTangent(d)
    return false
end
function isZeroTangent(d::ZeroTangent)
    return true
end
function isZeroTangent(d::AbstractArray{<:ZeroTangent})
    return true
end

# Zero-copy type conversion helpers (shared by forwarddiff.jl and reversediff.jl)
@inline _as_f64_vec(v::Vector{Float64}) = v
@inline _as_f64_vec(v::AbstractVector{Float64}) = v
@inline _as_f64_vec(v::AbstractVector{<:Real}) = Float64.(v)

@inline _as_u32_vec(v::Vector{UInt32}) = v
@inline _as_u32_vec(v::AbstractVector{<:UInt32}) = convert(Array{UInt32,1}, v)
@inline _as_u32_vec(v::AbstractVector) = convert(Array{UInt32,1}, unsense(v))

# In-place value extraction from ReverseDiff TrackedReal vectors
@inline function _rd_extract_values(v::AbstractVector{<:ReverseDiff.TrackedReal})
    out = Vector{Float64}(undef, length(v))
    @inbounds for i in eachindex(v)
        out[i] = ReverseDiff.value(v[i])
    end
    out
end
