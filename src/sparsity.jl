using SparseArrays, FMIBase, SparseDiffTools

function fmi2dependencyKindToDependencyIndex(kind::fmi2DependencyKind)
    if kind == fmi2DependencyKindDependent
        return UInt32(5)
    end
    return kind
end

function fmi3dependencyKindToDependencyIndex(kind::fmi3DependencyKind)
    return kind
end

struct DependencyMatrix
    matrix::SparseMatrixCSC
    vr_idx_dict::Dict
end

# From FMI3-Standard:  
# If dependencies is not present, it must be assumed that the unknown depends on all knowns. If dependencies is present as empty list, the unknown depends on none of the knowns.

DependencyMatrix(md::fmi2ModelDescription) = begin
    vrs = [mV.valueReference for mV in md.modelVariables]
    sort!(vrs); unique!(vrs)
    vr_idx_dict = Dict(zip(vrs, 1:length(vrs)))
    dep_mtx = spzeros(UInt32, length(vrs), length(vrs))
    @info "Constructing Dependency Matrix"
    for dependency_category in [md.modelStructure.derivatives, md.modelStructure.outputs, md.modelStructure.initialUnknowns]
        for dep_info in dependency_category
            dependent_vR = md.modelVariables[dep_info.index].valueReference
            if !isnothing(dep_info.dependencies)
                for (idx, dependency) in enumerate(dep_info.dependencies)
                    dependency_vR = md.modelVariables[dependency].valueReference
                    dependency_kind = dep_info.dependenciesKind[idx]
                    dep_mtx[vr_idx_dict[dependent_vR], vr_idx_dict[dependency_vR]] = fmi2dependencyKindToDependencyIndex(dependency_kind)
                end
            else
                dep_mtx[vr_idx_dict[dependent_vR], :] .= 5
            end
        end
    end
    DependencyMatrix(dep_mtx, vr_idx_dict)
end

DependencyMatrix(md::fmi3ModelDescription) = begin
    vrs = [mV.valueReference for mV in md.modelVariables]
    sort!(vrs); unique!(vrs)
    vr_idx_dict = Dict(zip(vrs, 1:length(vrs)))
    dep_mtx = spzeros(UInt32, length(vrs), length(vrs))
    @info "Constructing Dependency Matrix"
    for dependency_category in [md.modelStructure.continuousStateDerivatives, md.modelStructure.outputs, md.modelStructure.initialUnknowns]
        for dep_info in dependency_category
            dependent_vR = dep_info.index
            if !isnothing(dep_info.dependencies)
                for (idx, dependency) in enumerate(dep_info.dependencies)
                    dependency_vR = dependency
                    dependency_kind = dep_info.dependenciesKind[idx]
                    dep_mtx[vr_idx_dict[dependent_vR], vr_idx_dict[dependency_vR]] = fmi3dependencyKindToDependencyIndex(dependency_kind)
                end
            else
                dep_mtx[vr_idx_dict[dependent_vR], :] .= 5
            end
        end
    end
    DependencyMatrix(dep_mtx, vr_idx_dict)
end

function Base.getindex(D::DependencyMatrix, i::UInt32, j::UInt32)
    return D.matrix[D.vr_idx_dict[i], D.vr_idx_dict[j]]
end
function Base.getindex(D::DependencyMatrix, i::Vector{UInt32}, j::Vector{UInt32})
    return D.matrix[getindex.(Ref(D.vr_idx_dict), i), getindex.(Ref(D.vr_idx_dict), j)]
end

function get_coloring(D::DependencyMatrix, f_refs::Vector{UInt32}, x_refs::Vector{UInt32}, row_based::Bool = false)
    mat = D[f_refs, x_refs]
    mat[mat.<5] .= 0
    return matrix_colors(mat, partition_by_rows = row_based), mat
end

## Plan für Decrompression:
# Funktion decompress, die aus Sparse Compressed und Compression Map Volle Matrix macht
# Funktion decompression_map, die aus sparsity und coloring eine Decompression Map macht

## Decompression Map:
# Es reicht eine Matrix, die die Information enthält, in welche Spalte welcher Eintrag der  Matrix komprimierten kommt, 
# weil Zeile ja immer beibehalten wird

#[1 0 0 1
# 1 0 1 0
# 0 1 1 0
# 0 1 0 1] bekommt das Coloring [1 1 2 2] und die Compression Map
#[1 4
# 1 3
# 2 3
# 2 4]
#


# function decompression_map(colorvec::Vector, sparsity::Matrix)
#     map = zeros(Int, size(sparsity)[1], maximum(colorvec))
#     nonzeros = findall(x-> x!=0, sparsity)
#     for i in 1:maximum(colorvec)
#         current_cols = findall(x->x==i, colorvec)
#         for nonzero in nonzeros[findall(x->x[2] in current_cols, nonzeros)]
#             map[nonzero[1],i] = nonzero[2] 
#         end
#     end
#     map
# end

# function decompress_sparse(compressed, decompression_map::Matrix{Int})
#     decompressed = zeros(size(decompression_map)[1], maximum(decompression_map))
#     for idx in eachindex(IndexCartesian(), decompression_map)
#         if decompression_map[idx] == 0
#             decompressed[idx] = 0.0
#         else
#             decompressed[idx[1], decompression_map[idx]] = compressed[idx]
#         end
#     end
#     return decompressed
# end

# sparsity = [1 0 0 1
#             1 0 1 0
#             0 1 1 0
#             0 0 0 1]
            
# colorvec = [1, 1, 2, 2]

# matrix_colors(sparsity, partition_by_rows = false)
# map = decompression_map(colorvec, sparsity)
# decom = decompress_sparse(ones(4,2), map)



function compute_coo_from_coloring(sparsity_pattern, coloring)
    # Initialize lists to store COO format data
    # Row-Indices
    I = Int[]
    # Col-Indices
    J = Int[]
    # Loop through each color
    for color in 1:maximum(coloring)
        color_cols = findall(coloring .== color)
        # Loop through each row for each color, and identify the nonzero column
        for row in 1:size(sparsity_pattern, 1)
            nzidxs = findall(sparsity_pattern[row, color_cols] .!= 0)
            println(nzidxs)
            if length(nzidxs) > 1
                error("Got more than one nonzero entry in Row $row at colums $nzidxs. This indicates that coloring and sparsity_pattern do not match.")
            elseif length(nzidxs) == 0
                # This is needed in Case there is a Zero in the Compressed Jacobian. In this case we want to keep its coordinates, so the Value-Vector keeps its shape n_rows x n_colors
                push!(I, row)
                push!(J, color_cols[1])
            else
                push!(I, row)
                push!(J, color_cols[nzidxs[1]])                
            end
        end
    end
    return I, J
end


I, J = compute_coo_from_coloring(sparsity, [1,1,1,1])
values = ones(8)

values = [1,2,3,4,5,6,7,8]
A = sparse(I, J, values)

# https://stackoverflow.com/questions/77510614/is-there-an-efficient-way-to-reuse-the-sparsity-pattern-of-a-sparse-matrix-in-ju
function sparse_sched(I,J,n = maximum(I), m = maximum(J))
    D1 = collect(1:length(I))
    Z = sparse(I,J,1:length(I),n, m,(x,y)->(D1[y] = x; x))
    D2 = similar(D1)
    for i in 1:nnz(Z)
        D2[Z.nzval[i]] = i
    end
    for i in eachindex(D2)
        D1[i] == 0 || (D2[i] = D2[D1[i]])
    end
    return D2
end

sched = sparse_sched(I, J)


function apply_sched!(A,S,V,op = +)
    A.nzval .= zero(eltype(A))
    for i in eachindex(V)
        v = A.nzval[S[i]]
        A.nzval[S[i]] = iszero(v) ? V[i] : op(v, V[i])
    end
    return A
end

apply_sched!(A, sched, [8, 7, 6, 5, 4, 3, 2, 1])