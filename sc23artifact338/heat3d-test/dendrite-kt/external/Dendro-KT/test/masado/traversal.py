## Hopefully easier than the C++ version I was in the middle of.
##
## Masado Ishii  --  UofU SoC, 2018-11-15

## 2018-11-16
## ----------
##   It is relatively straightforward to enumerate the possible paths for
## a single level of refinement. It is possible to count the minimum
## number of orientations contained in each of these paths. My program
## currently does both of these tasks. However, this is not the same as
## finding the minimum number of orientations needed overall,
## since it only captures the complexity of a single level of refinement.
##
##   What we are after is the smallest orientation table that is closed
## under multiplication.
##
##   Note that by unique 'orientation' I mean a pair:
## (region starting node, curve net displacement).
## Such a pair determines the connectivity of the curve segment and 
## the volume filled. Different choices in the hyperplane are left flexible.
##
##   With that clarification, our problem can be formulated as a graph problem. 
## The graph to consider is a bipartite graph. On the left is the set
## of orientations. On the right is a set of groupings of orientations.
## An edge from an orientation (left) to a grouping (right) means
## one of the ways to subdivide that orientation is through the grouping on the
## right. Edges from a grouping (right) to orientations (left) mean
## the orientations are contained in that grouping. In my program thus far,
## an initial orientation on the left is chosen, namely (0, +highest_dim),
## and all groupings which receive an edge from this orientation are explored,
## and we count the number of orientations contained in each grouping. The
## larger problem is to find a minimal orientation table that is closed under
## multiplication. In terms of our graph, we want to select a subset of
## groupings on the right, subject to two constraints. The closure constraint
## is: Every orientation in the range of the subset of groupings must have an
## edge back into the same subset of groupings. 'Range' here means all the
## orientations which receive an edge from any one in the subset of groupings.
## The minimization constraint is: The range of the subset of groupings must
## be minimal.
##

import math
import collections
import itertools

__DEBUG__=False;

def __main__():
    ##dim = 2;   ## The answer it gives is exactly 1 possible traversal (correct).
    dim = 3;   ## The answer it gives is 5 possible traversals.
    ##dim = 4;   ## The answer it gives is 5733 possible traversals.
    ##dim = 5;   ## An infeasible number of traversals.
    (num_solutions, min_num_repl, min_lines, max_num_repl, max_lines, ss_repl) = traverse(dim);
    print("Dimension==%d" % dim);
    print("Total # of valid displacement sequences: %d." % num_solutions);
    print("Minimum # of orientations (%d) on lines" % min_num_repl, min_lines);
    print("Maximum # of orientations (%d) on lines" % max_num_repl, max_lines);
    print("Number of sets of SAT orientations: %d" % len(ss_repl));

    ### for s_repl in ss_repl:
    ###     print(s_repl);

    orig_task = (0, dim-1);
    refinements = ssrepl2refinements(ss_repl);
    (min_size, solutions) = solve(dim, orig_task, refinements);
    print("Victory! min_size==%d, solutions:" % min_size);
    print(*solutions, sep='\n');
    

    ## Results so far:
    ## 3D -- My program claims there are two possible orientation tables with
    ##       the minimum length of 12 rows:
    ##
    ##           Victory! min_size==12, solutions:
    ##           ((0, 0), (0, 1), (0, 2), (3, 0), (3, 1), (3, 2), (5, 0), (5, 1), (5, 2), (6, 0), (6, 1), (6, 2))
    ##           ((1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (4, 0), (4, 1), (4, 2), (7, 0), (7, 1), (7, 2))


class Progress:
    def __init__(self, total_count, frac_freq):
        self.total_count = total_count;
        self.frac_freq = frac_freq;
        self.progress_count = 0;
        self.last_emit = 0;

    def begin(self):
        print("0%", end='');

    def emit(self):
        print("....%d%%" % round(self.progress_count / self.total_count * 100), end='', flush=True);
        self.last_emit = self.progress_count;

    def update_inc(self, inc=1):
        self.progress_count += inc;
        if ((self.progress_count - self.last_emit) / self.total_count >= self.frac_freq):
            self.emit();

    def update_prog(self, prog):
        self.progress_count = prog;
        if ((self.progress_count - self.last_emit) / self.total_count >= self.frac_freq):
            self.emit();

    def finish(self):
        print();
        


def generate_base_path(dim):
    if (dim == 0):
        return [];
    else:
      lower_level = generate_base_path(dim-1);
      return lower_level + [1 << (dim-1)] + lower_level;

def generate_moves(inregion_pos, outregion_pos, region_disp, dim):
    if ((inregion_pos & region_disp) == (outregion_pos & region_disp)):
        return [region_disp];
    else:
        return (1 << d for d in range(dim-1,-1,-1) if (1 << d) != region_disp);

def traverse(dim):
    my_base_path = generate_base_path(dim);
    my_base_order = (
            [0] + list(itertools.accumulate(my_base_path, lambda x,y : x ^ y)));

    ### print("Displacement sequences:");

    num_regions = 1 << dim;
    assert num_regions == len(my_base_order), "num_regions does not match region visitation list!";

    inregion_start = 0;
    inregion_end = (1 << (dim-1));  ## Net disp is major dim, both inner and outer.
    num_solutions = 0;

    move_hist = [];
    moves = collections.deque();

    ## Count how many orientations we will need.
    ## NOTE: This does not take into account further iterations, which
    ##       fill out the multiplication table of orientations.
    inregion_hist = [];
    min_num_repl = None;
    min_lines = [];
    max_num_repl = None;
    max_lines = [];

    ## Enumerate all sets of orientations that can be used
    ## to refine this permutation of this orientation.
    ss_repl = set();

    distance = 0;
    inregion_pos = inregion_start;

    save_stuff = dict();

    while (True):
        if __DEBUG__: print("      State: d=%d, [out,in] == [%d,%d]" % (distance, my_base_order[distance], inregion_pos), end=' ');
        ## Is region terminal?
        if (distance == num_regions-1):
            ## If yes, and we reached the goal, then record a solution.
            final_move = my_base_path[-1];
            if (inregion_pos == inregion_end ^ final_move):
                if __DEBUG__: print("Applying the final move (0,%d)." % final_move);
                move_hist.append(final_move);
                inregion_hist.append(inregion_pos);
                ### print("%5d  " % num_solutions, move_hist);   ## Comment out to save on i/o

                ## Keep track of min and max num of orientations,
                ## which are determined by net displacement and the
                ## starting corner.
                s_repl = frozenset(zip(inregion_hist, move_hist));
                num_repl = len(s_repl);
                if (min_num_repl is None or num_repl < min_num_repl):
                    min_num_repl = num_repl;
                    min_lines = [num_solutions];
                elif (num_repl == min_num_repl):
                    min_lines.append(num_solutions);
                if (max_num_repl is None or num_repl > max_num_repl):
                    max_num_repl = num_repl;
                    max_lines = [num_solutions];
                elif (num_repl == max_num_repl):
                    max_lines.append(num_solutions);

                ## Record the exact set of orientations.
                ss_repl.add(s_repl);
                if (dim == 4 and num_solutions == 796):
                    save_stuff[796] = list(zip(inregion_hist, move_hist));

                num_solutions += 1;
            else:
                if __DEBUG__: print("Stuck at %d" % inregion_pos);

        else:
            ## Otherwise, we generate new moves.
            new_moves = list(generate_moves(
                    inregion_pos, my_base_order[distance],
                    my_base_path[distance], dim));
            for m in new_moves:
                moves.append((inregion_pos, m, distance));

        ## Apply the next move.
        if (not moves):
            break;
        (inregion_pos, next_move, distance) = moves.pop();
        if __DEBUG__: print("Applying the next move (%d,%d)." % (my_base_path[distance], next_move));
        move_hist[distance:] = [next_move];
        inregion_hist[distance:] = [inregion_pos];
        inregion_pos ^= next_move;
        inregion_pos ^= my_base_path[distance];
        distance += 1;

    print("Base path:", my_base_path, sep='\n');

    print(save_stuff);

    return (num_solutions, min_num_repl, min_lines, max_num_repl, max_lines, ss_repl);

## The second component of a task (i, highest_dim) can be expressed
## either as a power of two or as the index of the dimension.
## This function transforms from the former representation to the latter.
def ssrepl2refinements(ss_repl):
    return set(map(lambda r: frozenset(map(lambda i_hdp2: (i_hdp2[0], i_hdp2[1].bit_length() - 1), r)), ss_repl));

def binary_decompose(i):
    while (i):
        yield i & 1;
        i >>= 1;

def permute_bits(i, p):
    return sum(map(lambda b,d: b << d, binary_decompose(i), p));


## Tasks are represented as tuples (i, highest_dim). The index 'i' can be
## considered as a compressed bit-vector of the starting point in D coordinates.
##
## This function is a generator of lambdas that extend the permute/reflect
## transformations to the task tuples. There are (2^D)(D!) transformations.
def get_permute_reflect(dim):
    perms = itertools.permutations(list(range(0,dim)));
    flips = list(range(0, 1<<dim));
    ## The group of permute/reflect transformations is a semidirect product
    ## of the permutations and the flips. This implies that each transformation
    ## can be expressed as a permutation followed by a flip.
    for (p, f) in itertools.product(perms, flips):
        if __DEBUG__: print(p, f);
        yield lambda i__hd: (f ^ permute_bits(i__hd[0],p), p[i__hd[1]]);


def solve(dim, orig_task, refinements):
    ## refinements is a set of satisfying sets of orientations for the
    ## task (0,+highest_dim) with the default under-permutation.
    ##
    ## Not easy to isolate the under-permutations for a single task.
    ## But we can rotate through all full-permutations to produce
    ## all tasks with all under-permutations. That's what we ultimately
    ## need to work with anyway.
    ##
    ## This dictionary stores answers to the question:
    ##   Using a given set of subtasks, what parent tasks can be implemented
    ##   with a configuration that is built from this exact set of subtasks?
    ##
    progress_counter = Progress((1 << dim)*math.factorial(dim), 0.01);
    progress_counter.begin();
    satisfied_parents = dict();
    for phi in get_permute_reflect(dim):
        rot_parent = phi(orig_task);
        rot_refinements = map(lambda r: tuple(sorted(set(map(phi, r)))), refinements);
        for r in rot_refinements:
            try:
                satisfied_parents[r].add(rot_parent);
            except:
                satisfied_parents[r] = set([rot_parent]);
        progress_counter.update_inc(1);
    progress_counter.finish();
    ## Note: Changed the type of the dictionary keys: previously frozensets of
    ## tuples, now monotone tuples of tuples.
    ## The reason is that next we will operate on monotone tuples of tuples.
    ## Would like to not have to cast to frozenset for each one, now that the
    ## keys are computed.

    ## For 3D there is exactly one satisfied parent per set of refining tasks.
    ## For 4D there are 1 or 2 satisfied parents per set of refining tasks.
    ## So we have to keep the dictionary on sets, even though it's slow.

    ## DEBUG
    print("Maximum # satisfied parents per refinement: %d"
            % max(map(len, satisfied_parents.values())));
    print("Total # of distinct values in the table: %d"
            % len(set().union(*satisfied_parents.values())));
    ### for (k,v) in satisfied_parents.items():
    ###     print(k,v, sep=': ', end='\n\n');

    ## Since there are so many refinements, it becomes infeasible to 
    ## loop through all (lower) subsets of refinements. While this approach
    ## never tests sets of tasks that would form incomplete refinements, the
    ## approach actually does redundant work. There are fewer tasks than
    ## possible refinements (in 4D, it is 384 versus a million). So a better
    ## approach is to loop through all (lower) subsets of tasks rather than
    ## subsets of refinements.
    ##
    ## Accumulate results with a dynamic-programming algorithm.
    ## Induction on task subset size, so we encounter minimal solutions first.
    all_tasks = sorted(set((phi(orig_task) for phi in get_permute_reflect(dim))));
    sat_parents_tasks = {tuple(): set()};  ## Base case is null.
    cand_size = 0;
    solutions = [];
    while(not solutions):
        sat_parents_tasks_lower = sat_parents_tasks;
        sat_parents_tasks = {};
        cand_size += 1;
        level_combos = itertools.combinations(all_tasks, cand_size);
        for combo in level_combos:
            if (combo in satisfied_parents):
                new_sat_p = satisfied_parents[combo];
            else:
                new_sat_p = set();
            sat_parents_tasks[combo] = new_sat_p.union(
                    *map(lambda subcombo: sat_parents_tasks_lower[subcombo],
                        itertools.combinations(combo, cand_size-1)));       ## Recursion
            ## The CLOSURE condition:
            if (sat_parents_tasks[combo]):
                print(sat_parents_tasks[combo]);
            if (set(combo).issubset(sat_parents_tasks[combo])):
                ## Declare victory.
                solutions.append(combo);

            ## Add all solutions in the same (lowest) level.
        print("Finished level %d.\n" % cand_size);

    ##TODO when we build the satisfied_parents table, we should include
    ##     a tag in the parents about which permutation of the parent was
    #      satisfied. We need to be able to recover a full construction.

    return (cand_size, solutions);


__main__();
