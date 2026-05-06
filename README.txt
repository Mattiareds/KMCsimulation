================================================================================
  KINETIC MONTE CARLO SIMULATION OF ATOMIC DIFFUSION ON A NANOPYRAMID SURFACE
  README - Code Logic and Architecture
================================================================================

--------------------------------------------------------------------------------
1. OVERVIEW
--------------------------------------------------------------------------------

This program implements a Kinetic Monte Carlo (kMC) simulation of atomic
diffusion and deposition on the shell of a truncated octahedron nanoparticle,
but can be adapted to others geometry under request. 
The core of the nanoparticle is fixed; only the shell sites are active. Atoms are
deposited stochastically on the shell surface and then diffuse according to
thermally activated hopping rates between nearest-neighbour (NN) sites.
The activation barriers reported are calculated via Density Functional Theory
relax calculations, combined with Drag Method.

The simulation follows the standard Residence Time Algorithm (BKL algorithm):
at each step, all possible processes and their rates are known, a process is
chosen with probability proportional to its rate, and the simulation clock is
advanced by a random time drawn from an exponential distribution.

The codebase is split into three main classes:
  - coordinates  : reads and stores atomic positions (xyz files)
  - geometry     : builds the geometric structure (planes, NN tables, site types)
  - metropolis   : runs the kMC simulation itself


--------------------------------------------------------------------------------
2. INPUT FILES
--------------------------------------------------------------------------------

piramide.xyz
    XYZ file for the CORE nanoparticle. Defines the fixed inner structure.
    Used to determine which shell sites lie on the core surface (border sites).

siti_jmol.xyz
    XYZ file for the SHELL diffusion sites. These are the sites on which
    atoms can be deposited and between which atoms can hop.

coordinates_settings.in
    Settings for the CORE geometry:
      Line 1: d_pv (nearest-neighbour distance, Angstrom)  spigolo (edge length)
      Line 2: n_cut[0..3] (number of planes cut from each face family)

shell_settings.in
    Same format as above, but for the SHELL geometry.

mc_settings.in
    Monte Carlo parameters (in order):
      n_barriers       : number of distinct energy barriers
      barriers[0..N]   : barrier values in eV (see Section 5 for mapping)
      E_b              : additional energy per nearest neighbour (eV)
      F                : deposition rate (atoms per site per unit time)
      filling          : target fractional coverage (stops deposition when reached)
      type             : simulation type (3 = kMC)
      temperature      : temperature in Kelvin
      nu_0             : attempt frequency prefactor (Hz), typically 1e12
      n_class          : total number of process classes
      total_time       : total simulation time


--------------------------------------------------------------------------------
3. CLASS: coordinates
--------------------------------------------------------------------------------

Stores the raw atomic positions read from an xyz file.

Key data:
  siti             : Nx3 vector of (x, y, z) coordinates
  specie_chimica   : chemical species label per site
  N                : number of atoms (from the xyz header)
  d_pv             : nearest-neighbour distance (from settings file)
  spigolo          : edge length of the pyramid without cuts
  n_cut[4]         : cuts applied to top(0), bottom(1), lateral(2), other(3)

Key methods:
  initialize()     : reads the xyz file
  settings_reader(): reads d_pv, spigolo, n_cut from the settings file
  output_writer()  : writes current configuration to output (xyz format)


--------------------------------------------------------------------------------
4. CLASS: geometry
--------------------------------------------------------------------------------

Builds and stores the geometric relationships between shell sites.

--- 4.1 Plane definitions (make_planes) ---

The nanopyramid has faces belonging to two crystallographic families:
  {100} faces : top face (plane 0), bottom face (plane 1), lateral faces (2-5)
  {111} faces : eight diagonal faces (planes 6-13)
  Second layer: plane 14, a virtual plane one layer below the bottom face

Each plane is defined by a normal vector (nx, ny, nz) and an offset d such
that a site i lies ON the plane if:
  |nx*xi + ny*yi + nz*zi - d| / ||(nx,ny,nz)|| < epsilon  (epsilon = 1e-3)

The offsets are derived from d_pv, spigolo, and n_cut using:
  factor = d_pv / sqrt(2)
  offset_top    = factor * (spigolo - 1 - n_cut[0])
  offset_bottom = factor * (spigolo - 1 - n_cut[1])
  offset_lateral= factor * (spigolo - 1 - n_cut[2])
  offset_111    = factor * (spigolo - 1)

--- 4.2 Site classification (site_characterisation) ---

Each shell site is tested against all 15 planes. The number of planes it
belongs to determines its type (stored in type_of_plane):

  0  : internal {100} site (on exactly one {100} plane of the shell)
  1  : internal {111} site (on exactly one {111} plane of the shell)
  2  : edge between two {111} faces
  3  : edge between one {111} and one {100} face
  4  : edge between two {100} faces
  5  : corner between two {111} and one {100} face
  6  : corner between two {111} and two {100} faces
  14 : second-layer site (on the virtual second-layer plane)

After test_border() is called (see below), border sites get additional types:
  7  : {100} border site (on a shell {100} plane AND on a core surface plane)
  8  : {111} border site (on a shell {111} plane AND on a core surface plane)
  15 : second-layer border site

Sites of type 14 are stored in upper_sites and initially deactivated.

--- 4.3 Border detection (test_border) ---

After classifying shell sites against shell planes, each site on a single
shell plane (type 0 or 1) is tested against the core planes (planes 0-13,
NOT plane 14 which is internal). If the site also lies on a core surface
plane, it is a border site and its type is incremented by 7:
  0 -> 7  (100 border)
  1 -> 8  (111 border)

This identifies the ring of shell sites that sit directly on the core facets.

--- 4.4 Nearest-neighbour table (initialize_nn) ---

For each pair of sites (i, j), if their distance is within 0.08 Angstrom of
d_pv, site j is added to the NN list of site i. This gives table_of_nn[i],
the list of geometric nearest neighbours of each site.

--- 4.5 Edge substitution (pv_substitution) ---

Edge sites (type > 1 and type < 7) are geometrically unstable and are not
used as diffusion endpoints. For each internal site i (type < 2), if one of
its NNs is an edge site (1 < type < 7), that edge NN is replaced in the NN
list of i by the nearest non-edge site on the ADJACENT face (not the same
face as i). This is recorded in edge_map[i][j], which stores the original
edge site index if a substitution was made, or -1 if no substitution occurred.

The substitution ensures that diffusion rates between two face-internal sites
correctly account for the geometry of the edge between them, using the
adjacent-face site as the effective diffusion partner.


--------------------------------------------------------------------------------
5. CLASS: metropolis - PROCESS CLASSIFICATION
--------------------------------------------------------------------------------

The kMC simulation uses a classification scheme where each possible hop is
assigned to a class index. The class encodes both the TYPE of hop (which
determines the base energy barrier) and the NUMBER OF NEAREST NEIGHBOURS of
the hopping atom (which adds E_b per neighbour to the barrier).

The class index is:
  class = base_process_index + nnn_atoms[site]

where nnn_atoms[site] is the number of occupied NN sites of the hopping atom
(counted only among same-face neighbours, i.e. edge_map entry < 0).

Base process indices and their barrier mapping (barrier index in parentheses):

  Base  {111}c -> {111}c, same facet             (barrier[0])
  Base  {111}b -> {111}c, same facet              (barrier[6])
  Base  {111}b -> {111}b, different facets        (barrier[2])
  Base  {100}b -> {111}b  (escape from 100)       (barrier[4])
  Base  {111}b -> {100}b  (enter 100)             (barrier[3])
  Base  {100}c -> {100}c, same facet              (barrier[1])
  Base  {100}b -> {100}c                          (barrier[5])
  Base  {100}c -> {100}b                          (barrier[7])
  Base  push ({111}b pushes atom on {100}b)       (barrier[8])
  Base  mixed site -> second-layer site           (barrier[9])
  Base  kink (second-layer -> first-layer         (barrier[10])
             even if destination is occupied)

The hop rate for class c is:
  k(c) = nu_0 * exp( -(barrier[b] + nnn * E_b) / (kB * T) )

where b is the barrier index, nnn is the neighbour count, kB is Boltzmann's
constant, and T is temperature.

The barrier file (mc_settings.in line 3) lists barriers in this order:
  [0]=111C-111C  [1]=100C-100C  [2]=111<->111  [3]=111->100  [4]=100->111
  [5]=100B-100C  [6]=111B-111C  [7]=100C-100B  [8]=push      [9]=111mix->100
  [10]=kink


--------------------------------------------------------------------------------
6. CLASS: metropolis - DATA STRUCTURES
--------------------------------------------------------------------------------

atoms[i]
    Boolean vector. True if site i is currently occupied.

nnn_atoms[i]
    Integer vector. Number of occupied NN sites of site i, counting only
    same-face neighbours (edge_map entry < 0 for that NN slot).

table_of_processes[i]
    For each site i, the list of BASE process indices geometrically available
    from i (regardless of current occupation). This is filled once at startup
    by table_of_processes_filler() and updated only when the geometry changes
    (e.g. second-layer activation). Push and kink processes are NOT stored
    here; they are managed separately as dynamic processes.

table_of_end_pos[i]
    Parallel to table_of_processes[i]. For each process of site i, the
    destination site index.

table_of_initial_pos[j]
    Inverse map: for each site j, the list of sites from which j can be
    reached. Used by map_of_class_next_eraser to remove processes whose
    destination has become occupied.

map_class_processes[c][k][0,1]
    3D array. For class c, the k-th available process is the hop from site
    map_class_processes[c][k][0] to site map_class_processes[c][k][1].
    Entries with value -1 are empty slots. Size: 100 classes x 1000 processes.

movements_per_class[c]
    Integer. Number of currently available processes in class c.

P[c]
    Double. Total rate contribution of class c:
      P[c] = movements_per_class[c] * k(c)

dynamic_processes
    Temporary list (cleared each MC step) of push and kink processes that
    are available given the current configuration. Each entry is a tuple
    (site, destination, class). These are filled in classification() and
    erased at the start of the next classification() call.

deactivated_sites
    List of second-layer site indices not yet activated. Initially populated
    with all upper_sites. A site is removed from this list when activated.


--------------------------------------------------------------------------------
7. CLASS: metropolis - SIMULATION FLOW
--------------------------------------------------------------------------------

start_of_the_sim()
  1. Opens output files (MC_processes.txt, debug.out).
  2. Calls table_of_processes_filler() to build the permanent process table.
  3. Initialises atoms[] to false, deactivated_sites from upper_sites.
  4. Deposits the first atom (deposition_func).
  5. Computes interested_sites and nnn_atoms for the first configuration.

algorithm()  [main loop, runs while time < total_time]
  For each MC step:
  1. second_layer_activation() : checks if any second-layer site should now
                                  be activated (>= 4 occupied NN sites).
  2. classification()          : rebuilds the available process list for all
                                  interested sites (see Section 8).
  3. probability_filler()      : computes P[c] for all classes.
  4. time_prob_calc()          : advances time, selects deposition or diffusion,
                                  picks a specific process (see Section 9).
  5. Execute the chosen process:
       - Deposition: atoms[pos] = true, write xyz snapshot.
       - Diffusion:  atoms[pos] = false, atoms[next] = true.
  6. interested_sites_calc()   : computes the set of sites whose process list
                                  may have changed (NN of pos and NN of next).
  7. nn_updater()              : recomputes nnn_atoms for all affected sites.
  8. pos = next.


--------------------------------------------------------------------------------
8. classification() - INCREMENTAL PROCESS TABLE UPDATE
--------------------------------------------------------------------------------

This function updates map_class_processes for every site in interested_sites.
It does NOT recompute geometry; it only updates which processes are currently
available given the occupation state.

For each occupied site in interested_sites:
  1. Erase all current entries for this site from map_class_processes:
       for t = 0 to nnn_atoms[site]: map_of_class_position_eraser(site, t)
     This removes the site's processes from ALL classes it could have been in
     (since nnn_atoms may have changed, the old class is uncertain).
  2. Re-insert each geometrically possible process (from table_of_processes)
     into the correct class:
       if destination pv is empty:
         map_of_class_filler(base_process + nnn_atoms[site], site, pv)
  3. Check for push and kink opportunities:
       if destination pv is occupied and conditions are met, add a dynamic
       process to dynamic_processes and call map_of_class_filler.

At the START of classification(), before the loop:
  All dynamic_processes from the previous step are erased from map_class_processes
  and the dynamic_processes list is cleared.

map_of_class_position_eraser(site, n)
  Iterates over table_of_processes[site] and calls map_of_class_eraser for
  each process with class = base_process + n.

map_of_class_next_eraser(destination)
  When a destination site becomes occupied, removes from map_class_processes
  all processes that point TO that destination. Iterates over
  table_of_initial_pos[destination], and for each source site erases the
  process at class = base_process + nnn_atoms[source].
  Called in time_prob_calc() immediately after selecting a diffusion move,
  before classification() updates things.

map_of_class_eraser(c, i, j)
  Finds the entry (i, j) in map_class_processes[c], sets it to (-1, -1),
  decrements movements_per_class[c], and shifts subsequent entries left
  (via shifter) to maintain a compact array.

map_of_class_filler(c, i, j)
  Finds the first empty slot in map_class_processes[c], writes (i, j),
  and increments movements_per_class[c].


--------------------------------------------------------------------------------
9. time_prob_calc() - TIME ADVANCE AND PROCESS SELECTION
--------------------------------------------------------------------------------

1. Compute total deposition rate:
     P_dep = n_sites * F
2. Compute total diffusion rate:
     P_diff = sum of P[c] for all classes c
3. Total rate:
     P_tot = P_dep + P_diff
4. Advance time:
     dt = -log(U1) / P_tot    where U1 ~ Uniform(0,1)
5. Draw a second uniform random number U2.
   If U2 * P_tot > P_diff : deposition event
   Else                   : diffusion event
6. For diffusion: select class c with probability P[c] / P_diff,
   then select a uniformly random process within that class from
   map_class_processes[c].
7. Before executing the move, immediately remove the chosen process from
   map_class_processes (map_of_class_position_eraser for pos,
   map_of_class_next_eraser for next).


--------------------------------------------------------------------------------
10. SECOND LAYER ACTIVATION
--------------------------------------------------------------------------------

Sites of type 14 (second-layer {100}) start deactivated. Each MC step,
second_layer_activation() checks every deactivated site: if 4 or more of its
NN sites are occupied, the site is activated.

second_layer_updates(upper_site) then:
  1. Adds hop processes between upper_site and its NN second-layer sites
     (type 14->14 as {100}c, type 14->15 as {100}b).
  2. Adds kink processes (upper_site to occupied first-layer sites below).
  3. For each NN site of type 7 ({100} border of the first layer), calls
     activate_edges() to convert the old edge geometry into a passable
     transition site.

activate_edges(upper_site, pv)
  When the second layer is active, some previously excluded edge sites between
  the {111} and the new {100} become valid transition sites (type 10, "mixed").
  For each such edge in edge_map[pv]:
    - Restore the original NN (before pv_substitution) as the new NN of pv.
    - Set the edge site's type to 10 (mixed).
    - Append upper_site as a new NN of the edge site.
    - Append the edge site as a new NN of the old {111}b site.
    - Call table_of_processes_updater() for both affected sites.
    - Reset the edge_map entry to -1 (substitution no longer active).

IMPORTANT: whenever append_nn() adds a new NN to a site, the corresponding
edge_map entry for that site must also be extended with -1 to keep the two
vectors in sync. This is critical for nn_updater() to count second-layer
neighbours correctly.


--------------------------------------------------------------------------------
11. DEPOSITION
--------------------------------------------------------------------------------

deposition_func() picks a random site index, rejects it if:
  - Its type is >= 2 (edge or deactivated site)
  - It is already occupied
and retries until a valid empty surface site is found. The atom is placed,
n_deposited is incremented, and the event is logged.

Deposition stops (F is set to 0) when n_deposited >= filling * n_sites.


--------------------------------------------------------------------------------
12. OUTPUT FILES
--------------------------------------------------------------------------------

MC_processes.out
    One line per MC event:
      - For diffusion: "process_name  from  X  to  Y"
      - For deposition: "New atom: X  at time: T"
        followed by the full xyz snapshot of all shell sites.
    Useful if you are interested to follow each process, but can be very large, 
    its writing can be deactivated from mc_settings file.


occ*.out
  .xyz-like file with the coordinates of the configuration each time a new atom
  is deposited

debug.out
    Verbose log of internal state: process table at startup, site types,
    classification and erasure calls, neighbour counts, second-layer events.
    Used for debugging; can be very large.

--------------------------------------------------------------------------------
AUTHORS
--------------------------------------------------------------------------------
This code was written entirely by me, Mattia Piero Rossi, as part of my Master's 
degree thesis in Physics. This work was supervised by Prof. Riccardo Ferrando

--------------------------------------------------------------------------------
REFERENCES
--------------------------------------------------------------------------------
Here some interesting papers about this algorithm, that was fitted for this 
specific problem. 

https://doi.org/10.1016/0021-9991(75)90060-1

https://arxiv.org/abs/0904.2556 


================================================================================
END OF README
================================================================================