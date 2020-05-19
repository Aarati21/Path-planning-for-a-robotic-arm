# Path-planning-for-a-robotic-arm-
An effective planner for a planar arm to move from its start joint angles to the goal joint angles. It returns a plan that is collision-free and supports upto 7 degrees of freedom. 4 algorithms have been implemented.

(a) RRT <br />
(b) RRT-Connect <br />
(c) RRT∗ <br />
(d) PRM <br />

For 20 randomly generated start and goal pairs on map2.txt, and for each of the four planners, a table of results is reported showingthe comparison of:

(a) Average planning times <br />
(b) Average path qualities <br />
(c) Average number of vertices generated (in constructed graph/tree) <br />
(d) Success rates for generating solutions within 5 seconds <br />

