Periodic Domain Decomposition Test
----------------------------------

This test checks that the mapping of periodic partners is correct across domains. It does this by forming a 2D box domain with a known periodic configuration and directly checking that the periodic partners are correct. It supports multiple variables, which may be set to be periodic independently of each other. It also supports mixed periodic/non-periodic boundaries.

The node layout is for this test is:
```
12--13--14--15
 |   |   |   |
 8-- 9--10--11
 |   |   |   |
 4-- 5-- 6-- 7
 |   |   |   |
 0-- 1-- 2-- 3
 
 
^ y direction
|
|
|
----> x direction
```

With periodic bounds in the _x_ direction:

* 3 is mapped to 0
* 7 is mapped to 4
* 11 is mapped to 8
* 15 is mapped to 12

With periodic bounds in the _y_ direction:

* 12 is mapped to 0
* 13 is mapped to 1
* 14 is mapped to 2
* 15 is mapped to 3

With periodic bounds in _both_ directions:

* 12 is mapped to 0
* 13 is mapped to 1
* 14 is mapped to 2
* 15 is mapped to 0
* 3 is mapped to 0
* 7 is mapped to 4
* 11 is mapped to 8

With domain decomposition, the "index_from" array will be set to the global
solution index that corresponds to the index of the unknown in the the solution array. If these values are correctly set, the mapping will automatically apply periodic bounds during the copying of the results to the node data arrays.

With the known grid configuaration above, the periodic pairs can be determined. Unfortunately, the domain decomposition process uses random values, so the exact mapping may vary across systems. However, the values for the from array must be the same for any nodes that are periodic partners. This fact is used for the test. The mapping in constructed using the library as normal. An expected mapping is then created based on the solution mapping and the periodic information. Those mappings are compared. If there are no discrepancies, the test passes.

This test assumes the solution_map is correct, and may give unreliable results if that is not the case.

The test specifically confirms:

* The solution_map and from mapping on all processes is sane (i.e. all values
  are in the proper range)

* At least one process has each unknown

* The values are consistant across all processes (i.e. all processes that have
  an unknown will have the same entry in the mappings)

* The global "from" array mapping (i.e. the values from all the nodes combined) match the expected results for this system.

* The "from" array mapping on each process is correct for the unknowns on that process.
