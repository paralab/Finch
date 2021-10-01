if !@isdefined(Jpycachesim)
    include("Jpycachesim.jl")
    using .Jpycachesim
end

#pcs_get_cachesim_from_file("cachesim/cachedef");

# This is my laptop
# l3 = Cache("L3", 12288, 16, 64, "LRU");
# l2 = Cache("L2", 256, 4, 64, "LRU");
# l1 = Cache("L1", 32, 8, 64, "LRU");

al3 = Cache("L3", 32, 1, 64, "LRU");
al2 = Cache("L2", 16, 1, 64, "LRU");
al1 = Cache("L1", 8, 1, 64, "LRU");

bl3 = Cache("L3", 32, 1, 64, "LRU");
bl2 = Cache("L2", 16, 1, 64, "LRU");
bl1 = Cache("L1", 8, 1, 64, "LRU");

c1 = pcs_build_cache([al1,al2,al3]);
c2 = pcs_build_cache([bl1,bl2,bl3]);

pcs_load(0,1,c1);
(m1,h1) = pcs_get_stats(c1);

pcs_load(0,130,c2);
(m2,h2) = pcs_get_stats(c2);

pcs_load(0,1,c1);
(m3,h3) = pcs_get_stats(c1);

println("test stats: should be miss=1, ?, 0, hit = 0, ?, 1");
println("miss= "*string(m1[1])*", "*string(m2[1])*", "*string(m3[1]));
println("hit= "*string(h1[1])*", "*string(h2[1])*", "*string(h3[1]));

pcs_load(0,8,c1);
pcs_store(0,8,c1);
pcs_load(28, 32,c1);

pcs_load(0,128,c2);
pcs_store(0,16,c2);
pcs_load(28, 1024,c2);

pcs_print_stats(c1);
pcs_print_stats(c2);

cobj = unsafe_load(c1);
lev2 = unsafe_load(cobj.load_from);
lev3 = unsafe_load(lev2.load_from);

println("miss: "*string(Int(cobj.MISS.count)));
println("hit: "*string(Int(cobj.HIT.count)));
println("load: "*string(Int(cobj.LOAD.count)));
println("store: "*string(Int(cobj.STORE.count)));