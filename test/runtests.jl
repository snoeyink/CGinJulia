# run CGinJulia unit tests   Jack Snoeyink Jun 2016
using CGinJulia, Test

tests = ["CG_2dTypes","CG_3dTypes"]   # the test file names are stored as strings...
for t in tests
 include("$(t).jl")                         # ... so that they can be evaluated in a loop
end
