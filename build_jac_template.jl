using FMIZoo, FMIImport
ENV["EXPORTINGTOOL"] = "Dymola"
ENV["EXPORTINGVERSION"] = "2023x"
ENV["FMIVERSION"] = 2.0
ENV["FMUSTRUCT"] = "FMUCOMPONENT"
function getFMUStruct(modelname, mode, tool=ENV["EXPORTINGTOOL"], version=ENV["EXPORTINGVERSION"], fmiversion=ENV["FMIVERSION"], fmustruct=ENV["FMUSTRUCT"]; kwargs...)
    
    # choose FMU or FMUComponent
    if endswith(modelname, ".fmu")
        fmu = loadFMU(modelname; kwargs...)
    else
        fmu = loadFMU(modelname, tool, version, fmiversion; kwargs...) 
    end

    if fmustruct == "FMU"
        return fmu, fmu

    elseif fmustruct == "FMUCOMPONENT"
        inst, _ = FMIImport.prepareSolveFMU(fmu, nothing, mode; loggingOn=true)
        return inst, fmu

    else
        @assert false "Unknown fmuStruct, variable `FMUSTRUCT` = `$(fmustruct)`"
    end
end
c, fmu = getFMUStruct("SpringFrictionPendulumExtForce1D", :ME)
stateVRs = fmu.modelDescription.stateValueReferences
derVRs = [fmu.modelDescription.modelVariables[der.index].valueReference for der in fmu.modelDescription.modelStructure.derivatives]

modelVariablesForValueReference.(Ref(fmu), stateVRs)
modelVariablesForValueReference.(Ref(fmu), derVRs)

deps = fmu.modelDescription.modelStructure.derivatives
md = fmu.modelDescription
# find all state references in derivative dependencies
vr_dep_1 = [md.modelVariables[dep].valueReference for dep in fmu.modelDescription.modelStructure.derivatives[1].dependencies]
filtered_vr = [vr for vr in vr_dep_1 if vr in md.stateValueReferences]
deps_vr = [dep for dep in fmu.modelDescription.modelStructure.derivatives]
test = filtered_vr[1]
using SparseArrays

vrs = [mV.valueReference for mV in md.modelVariables]
sort!(vrs); unique!(vrs)
vr_idx_dict = Dict(zip(vrs, 1:length(vrs)))
vr_idx_dict[md.modelVariables[6].valueReference]
dep_mtx = spzeros(length(vrs), length(vrs))
dep_mtx[vr_idx_dict[md.modelVariables[1].valueReference], vr_idx_dict[md.modelVariables[1].valueReference]] = 1.0
dep_mtx[vr_idx_dict[md.modelVariables[1].valueReference], vr_idx_dict[md.modelVariables[2].valueReference]] = 1.0
# das bedeutet jetzt quasi vR[1] hängt ab von vR[1] und vR[1] hängt ab von vR[2]

# vR[12] hängt ab von vR[16]:
dep_mtx[vr_idx_dict[md.modelVariables[12].valueReference], vr_idx_dict[md.modelVariables[16].valueReference]] = 1.0

dep_mtx

dep_mtx = spzeros(UInt32, length(vrs), length(vrs))
for dependency_category in [md.modelStructure.derivatives, md.modelStructure.outputs, md.modelStructure.initialUnknowns]
    for dep_info in dependency_category
        dependent_vR = md.modelVariables[dep_info.index].valueReference
        for (idx, dependency) in enumerate(dep_info.dependencies)
            dependency_vR = md.modelVariables[dependency].valueReference
            dependency_kind = dep_info.dependenciesKind[idx]
            dep_mtx[vr_idx_dict[dependent_vR], vr_idx_dict[dependency_vR]] = dependencyKindToDependencyIndex(dependency_kind)
        end
    end
end


function dependencyKindToDependencyIndex(kind::fmi2DependencyKind)
    if kind == fmi2DependencyKindDependent
        return UInt32(5)
    end
    return kind
end

md.stateValueReferences
md.derivativeValueReferences
dep_mtx[[vr_idx_dict[vR] for vR in md.derivativeValueReferences], [vr_idx_dict[vR] for vR in md.stateValueReferences]]
dep_mtx

struct DependencyMatrix
    matrix::SparseMatrixCSC
    vr_idx_dict::Dict
end

function Base.getindex(D::DependencyMatrix, i::UInt32, j::UInt32)
    return D.matrix[D.vr_idx_dict[i], D.vr_idx_dict[j]]
end
function Base.getindex(D::DependencyMatrix, i::Vector{UInt32}, j::Vector{UInt32})
    return D.matrix[getindex.(Ref(D.vr_idx_dict), i), getindex.(Ref(D.vr_idx_dict), j)]
end

DependencyMatrix(md::fmi2ModelDescription) = begin
    vrs = [mV.valueReference for mV in md.modelVariables]
    sort!(vrs); unique!(vrs)
    vr_idx_dict = Dict(zip(vrs, 1:length(vrs)))
    dep_mtx = spzeros(UInt32, length(vrs), length(vrs))
    @info "Constructing Dependency Matrix"
    for dependency_category in [md.modelStructure.derivatives, md.modelStructure.outputs, md.modelStructure.initialUnknowns]
        for dep_info in dependency_category
            dependent_vR = md.modelVariables[dep_info.index].valueReference
            for (idx, dependency) in enumerate(dep_info.dependencies)
                dependency_vR = md.modelVariables[dependency].valueReference
                dependency_kind = dep_info.dependenciesKind[idx]
                dep_mtx[vr_idx_dict[dependent_vR], vr_idx_dict[dependency_vR]] = dependencyKindToDependencyIndex(dependency_kind)
            end
        end
    end
    DependencyMatrix(dep_mtx, vr_idx_dict)

end
dep_mtx = DependencyMatrix(md)
typeof(md.modelVariables[1].valueReference)

dep_mtx[md.derivativeValueReferences[2], md.stateValueReferences[1]]
dep_mtx[md.derivativeValueReferences, md.stateValueReferences]

using SparseDiffTools
colorvec = matrix_colors(dep_mtx[md.derivativeValueReferences, md.stateValueReferences])
colorvec = [1,1,1,2,3,1,1]
for i in 1:maximum(colorvec)
    n_hot = zeros(length(colorvec))
    n_hot[colorvec .== i] .= 1.0
    println(n_hot)
end

function n_hot(len::Int, i::Vector{Int})
    ret = zeros(len)
    ret[i] .= 1.0
    return ret
end

n_hot(5, [1,4])

test = [1, 2, 3, 4, 2,4, 1]
findall(x->x==3,test)