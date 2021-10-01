#=
Indexer struct for array type entities

An indexer has a symbol, range and current value
The range can be numbers or other indexers.
It can be a min,max pair, or an array of values.
=#

mutable struct Indexer
    symbol::Symbol   # symbol used in expressions
    range::Array     # Range can be a min/max pair, or an array of values.
    value::Int       # Current value (may not be needed?)
end
