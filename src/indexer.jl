#=
Indexer struct for array type entities

An indexer has a symbol, range and current value
The range is an array of integers or any other Array subclass that can be used like "for i in I.range"
=#

"""
    Indexer

An entity representing an index to be applied to a variable or coefficient. It has a 
symbol that can be used in expressions, a range of integer values, and a current 
value that can be accessed with `I.value` for an index labeled I.

This should be built with the `index` function.
"""
mutable struct Indexer
    symbol::Symbol   # symbol used in expressions
    range::Array     # Range of integers.
    value::Int       # Current value
    tag::Int         # This indexer's position in the array of indexers?
end
