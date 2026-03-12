#
# Copyright (c) 2023 Tobias Thummerer, Lars Mikelsons
# Licensed under the MIT license. See LICENSE file in the project root for details.
#

module FMISensitivity

# load modules for reusability in other frameworks
import SciMLSensitivity
import SciMLSensitivity: ForwardDiff
import SciMLSensitivity: FiniteDiff
import SciMLSensitivity: ReverseDiff
import SciMLSensitivity: Zygote

import FMIBase.ChainRulesCore
using FMIBase.ChainRulesCore: ZeroTangent, NoTangent, @thunk

using SciMLSensitivity.LinearAlgebra

using FMIBase
using FMIBase.FMICore
using FMIBase: undual, unsense, untrack, FMUEvaluationOutput

include("utils.jl")
include("types.jl")
include("primitives.jl")
include("chainrules.jl")
include("forwarddiff.jl")
include("reversediff.jl")

end # module
