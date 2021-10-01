from cachesim import CacheSimulator, Cache, MainMemory

#########################
verbose = False
model = 1
#########################

mem = MainMemory()
if model == 1:
    l3 = Cache("L3", 20480, 16, 64, "LRU")  # 20MB: 20480 sets, 16-ways with cacheline size of 64 bytes
    mem.load_to(l3)
    mem.store_from(l3)
    l2 = Cache("L2", 512, 8, 64, "LRU", store_to=l3, load_from=l3)  # 256KB
    l1 = Cache("L1", 64, 8, 64, "LRU", store_to=l2, load_from=l2)  # 32KB
elif model == 2:
    l3 = Cache("L3", 12288, 16, 64, "LRU")
    mem.load_to(l3)
    mem.store_from(l3)
    l2 = Cache("L2", 256, 4, 64, "LRU", store_to=l3, load_from=l3)
    l1 = Cache("L1", 32, 8, 64, "LRU", store_to=l2, load_from=l2)
elif model == 3:
    l3 = Cache("L3", 20, 4, 64, "LRU")
    mem.load_to(l3)
    mem.store_from(l3)
    l2 = Cache("L2", 2, 1, 64, "LRU", store_to=l3, load_from=l3)
    l1 = Cache("L1", 1, 1, 64, "LRU", store_to=l2, load_from=l2)

cs = CacheSimulator(l1, mem)


infile = open("cachesim_output.out", "rb")

while True:
    chunk = infile.read(9)
    if len(chunk) < 9:
        break
    
    t = chunk[0]
    a = int.from_bytes(chunk[1:8], byteorder='little', signed=False)
    
    if t==0:
        if verbose==True:
            print('load  '+str(a))
        cs.load(a, length=8)
    
    if t==1:
        if verbose==True:
            print('store '+str(a))
        cs.store(a, length=8)

cs.force_write_back()
cs.print_stats()