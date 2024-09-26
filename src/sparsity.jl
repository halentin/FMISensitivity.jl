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
    return matrix_colors(D[f_refs, x_refs], partition_by_rows = row_based)
end