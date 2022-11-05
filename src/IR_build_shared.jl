# These are pieces that are shared by FEM and FVM and anything else.

# A useful pattern: (a-b)*c + d
function generate_the_one_pattern(a, b, c, d)
    IRtypes = IR_entry_types();
    result = IR_operation_node(IRtypes.math_op, [:+, 
                d,
                IR_operation_node(IRtypes.math_op, [:*,
                    c, 
                    IR_operation_node(IRtypes.math_op, [:-, a, b])
                    ])
                ]);
    
    return result;
end

# TDM(A,b,C) = transpose(A) * diagm(b) * C
# P_ij =  A'_ik * b_k * C_kj = A_ki * b_k * C_kj
# return the IR for A[k,i] * b[k] * C[k,j]
# A,b,C are symbols or IR_data_node
# i,j,k are symbols, numbers, IR_part
function generate_linalg_TDM_product(A, b, C, i, j, k)
    IRtypes = IR_entry_types();
    if typeof(A) == IR_data_node
        A_part = IR_data_node(IRtypes.float_data, A.label, [:?,:?], [k,i]);
    elseif typeof(A) <: IR_part
        A_part = apply_indexed_access(A, [k,i], IRtypes);
    else
        A_part = IR_data_node(IRtypes.float_data, A, [:?,:?], [k,i]);
    end
    if typeof(b) == IR_data_node
        b_part = IR_data_node(IRtypes.float_data, b.label, [:?], [k]);
    elseif typeof(b) <: IR_part
        b_part = apply_indexed_access(b, [k], IRtypes);
    else
        b_part = IR_data_node(IRtypes.float_data, b, [:?], [k]);
    end
    if typeof(C) == IR_data_node
        C_part = IR_data_node(IRtypes.float_data, C.label, [:?,:?], [k,j]);
    elseif typeof(C) <: IR_part
        C_part = apply_indexed_access(C, [k,j], IRtypes);
    else
        C_part = IR_data_node(IRtypes.float_data, C, [:?,:?], [k,j]);
    end
    
    return IR_operation_node(IRtypes.math_op, [:*, A_part, b_part, C_part]);
end

# Tv(A,b) = transpose(A) * b
# P_i =  A'_ij * b_j = A_ji * b_j
# return the IR for A[j,i] * b[k]
# A,b are symbols or IR_data_node
# i,j are symbols, numbers, IR_parts
function generate_linalg_Tv_product(A, b, i, j)
    IRtypes = IR_entry_types();
    if typeof(A) == IR_data_node
        A_part = IR_data_node(IRtypes.float_data, A.label, [:?,:?], [j,i]);
    elseif typeof(A) <: IR_part
        A_part = apply_indexed_access(A, [j,i], IRtypes);
    else
        A_part = IR_data_node(IRtypes.float_data, A, [:?,:?], [j,i]);
    end
    if typeof(b) == IR_data_node
        b_part = IR_data_node(IRtypes.float_data, b.label, [:?], [j]);
    elseif typeof(b) <: IR_part
        b_part = apply_indexed_access(b, [j], IRtypes);
    else
        b_part = IR_data_node(IRtypes.float_data, b, [:?], [j]);
    end
    
    return IR_operation_node(IRtypes.math_op, [:*, A_part, b_part]);
end
