#=
Indexer struct for array type entities

An indexer has a symbol, range and current value
The range is an array of integers or any other Array subclass that can be used like "for i in I.range"
=#

mutable struct Indexer
    symbol::Symbol   # symbol used in expressions
    range::Array     # Range of integers.
    value::Int       # Current value
end
