using LinearAlgebra
using SparseDiffTools, SparseArrays

A = I
A = Matrix(I, 5, 5)
A[2,2] = 0
A[2,1] = 1
A = ones(5,5)
A = sparse(A)
c = matrix_colors(A)
maximum(c)
B
A



fmu.modelDescription.derivativeValueReferences

length(fmu.modelDescription.modelVariables)
Int(ans)

# Konstruktion der Adjacency Matrix n x n: 
n = length(fmu.modelDescription.modelVariables)
adj_mtx = zeros(Int, n, n)
for (var_idx, output) in enumerate(fmu.modelDescription.modelStructure.outputs)
    print("Output: ")
    println(output.index)
    for dependency in fmu.modelDescription.modelStructure.outputs[var_idx].dependencies
        println(dependency)
        adj_mtx[output.index,dependency] = 1
    end
end

for (var_idx, output) in enumerate(fmu.modelDescription.modelStructure.derivatives)
    print("Output: ")
    println(output.index)
    for dependency in fmu.modelDescription.modelStructure.derivatives[var_idx].dependencies
        println(dependency)
        adj_mtx[output.index,dependency] = 1
    end
end

for (var_idx, output) in enumerate(fmu.modelDescription.modelStructure.initialUnknowns)
    print("Output: ")
    println(output.index)
    for dependency in fmu.modelDescription.modelStructure.initialUnknowns[var_idx].dependencies
        println(dependency)
        adj_mtx[output.index,dependency] = 1
    end
end

fmu.modelDescription.modelVariables
fmu.modelDescription.stateValueReferences
fmu.modelDescription.modelVariables[1].valueReference
state_idxs = findall(x->x.valueReference in fmu.modelDescription.stateValueReferences, fmu.modelDescription.modelVariables)
adj_mtx[state_idxs, state_idxs]
matrix_colors(adj_mtx', partition_by_rows=false)
maximum(ans)

