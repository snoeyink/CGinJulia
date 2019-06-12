# run CGinJulia unit tests   Jack Snoeyink Jun 2016
using CGinJulia
using Base.Test

tests = ["GG_2dTypes"]   # the test file names are stored as strings...
for t in tests
 include("$(t).jl")                         # ... so that they can be evaluated in a loop
end
