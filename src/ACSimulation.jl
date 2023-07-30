#############################################################################
#############################################################################
#
# This is the main include file for the project. It exports all 
#   necessary functionality for running experiments.
#
# Work is based of this paper:
#   "An exact solution of the Grad-Shafranov-Helmhotz equation with 
#       central current density reversal"
#   doi: 10.1063/1.1924554
#                                                                               
#############################################################################
#############################################################################

using Bessels, Roots

import Base: -

include("System.jl")
include("Data.jl")
include("MagneticField.jl")
include("Current.jl")
include("Pressure.jl")
include("OptimisationUtil.jl")

nothing