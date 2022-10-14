#=
# Uses the elemental order to determine the nodal order.
# Order within elements is lexicographic. Other options may be added.
=#
export reorder_grid_element_first!

function reorder_grid_element_first!(grid, porder)
    dim = size(grid.allnodes,1);
    new_node_order= element_first_order_map(grid, porder, dim, true);
    reorder_grid_nodes!(grid, new_node_order);
    return grid;
end

function element_first_order_map(grid, porder, dim, invert = true)
    nnodes = size(grid.allnodes,2);
    eforder = zeros(Int, nnodes); # The ordering
    node_done = zeros(Int, nnodes);
    
    # Loop over the elements in elemental_order, filling in the eforder as we go
    l2g = grid.loc2glb;
    efind = 0;
    for ei=grid.elemental_order
        for ni=1:size(grid.loc2glb,1)
            # node's old index
            old_ind = grid.loc2glb[ni,ei];
            if node_done[old_ind] == 0 # not done yet
                # node's new index
                efind = efind+1;
                
                if invert
                    eforder[old_ind] = efind;
                else
                    eforder[efind] = old_ind;
                end
                node_done[old_ind] = efind;
                
            else # already mapped 
                
            end
            
        end
        
        
        
        
        # # Add interior first
        # for k=1:porder+1
        #     if dim > 1
        #         for j = 1:porder+1
        #             if dim > 2
        #                 for i = 1:porder+1
        #                     # 3D
        #                     lind = i + (porder+1)*((j-1) + (porder+1)*(k-1)); # local index
        #                     gind = l2g[lind, ei]; # global index
        #                     if !node_done[gind]
        #                         efind = efind+1;
        #                         if invert
        #                             eforder[gind] = efind;
        #                         else
        #                             eforder[efind] = gind;
        #                         end
        #                         node_done[gind] = true;
        #                     end
        #                 end # i loop
                        
        #             else # 2D
        #                 lind = j + (porder+1)*(k-1); # local index
        #                 gind = l2g[lind, ei]; # global index
        #                 if !node_done[gind]
        #                     efind = efind+1;
        #                     if invert
        #                         eforder[gind] = efind;
        #                     else
        #                         eforder[efind] = gind;
        #                     end
        #                     node_done[gind] = true;
        #                 end
        #             end
                    
        #         end # j loop
                
        #     else # 1D
        #         lind = k; # local index
        #         gind = l2g[lind, ei]; # global index
        #         if !node_done[gind]
        #             efind = efind+1;
        #             if invert
        #                 eforder[gind] = efind;
        #             else
        #                 eforder[efind] = gind;
        #             end
        #             node_done[gind] = true;
        #         end
        #     end
            
        # end # k loop
    end # elements
    
    return eforder;
end
