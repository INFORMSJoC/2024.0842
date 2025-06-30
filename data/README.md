The name of data file is prefixed with `rrap_` followed by several (optional) parts separated by '_':
- ns<number>: the number of subsystems is <number>
- nh<number>: the number of component types is <number>
- m<number>: the number of resource types is <number>
- g<number>: the value of $\gamma_{jh}$ based on which component reliabilities are randomly drawn.
- seed<number>: the index of instance associated with the given combination of configurations 


The structure of data file is: 
- The first line: the number of resource types, the number of subsystems, and the number of component types.
- The second line: the amount of each resource.
- The next <number of subsystems> lines, each with <number of component types> values representing the reliability of each component type (r_jh).
- The next <number of subsystems $\times$ number of resource types> lines: the units of resource i consumed by each type-h component in subsystem j (a_
  {ijh})
