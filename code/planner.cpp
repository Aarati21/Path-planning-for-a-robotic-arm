
#include <iostream> 
#include <math.h>
#include "mex.h"
#include <cstdio>
#include <chrono>
#include <stdio.h>
#include <vector>
#include <ctime>
#include <queue>
#include <algorithm>

 /* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10
using namespace std;




int x_size, y_size;
double* map;
int numofDOFs;
double* armstart_anglesV_rad;
double* armgoal_anglesV_rad;
int planner_id;

typedef struct {
    int X1, Y1;
    int X2, Y2;
    int Increment;
    int UsingYIndex;
    int DeltaX, DeltaY;
    int DTerm;
    int IncrE, IncrNE;
    int XIndex, YIndex;
    int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int* pY, int x_size, int y_size)
{
    double cellsize = 1.0;
    //take the nearest cell
    *pX = (int)(x / (double)(cellsize));
    if (x < 0) *pX = 0;
    if (*pX >= x_size) *pX = x_size - 1;

    *pY = (int)(y / (double)(cellsize));
    if (y < 0) *pY = 0;
    if (*pY >= y_size) *pY = y_size - 1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t* params)
{
    params->UsingYIndex = 0;

    if (fabs((double)(p2y - p1y) / (double)(p2x - p1x)) > 1)
        (params->UsingYIndex)++;

    if (params->UsingYIndex)
    {
        params->Y1 = p1x;
        params->X1 = p1y;
        params->Y2 = p2x;
        params->X2 = p2y;
    }
    else
    {
        params->X1 = p1x;
        params->Y1 = p1y;
        params->X2 = p2x;
        params->Y2 = p2y;
    }

    if ((p2x - p1x) * (p2y - p1y) < 0)
    {
        params->Flipped = 1;
        params->Y1 = -params->Y1;
        params->Y2 = -params->Y2;
    }
    else
        params->Flipped = 0;

    if (params->X2 > params->X1)
        params->Increment = 1;
    else
        params->Increment = -1;

    params->DeltaX = params->X2 - params->X1;
    params->DeltaY = params->Y2 - params->Y1;

    params->IncrE = 2 * params->DeltaY * params->Increment;
    params->IncrNE = 2 * (params->DeltaY - params->DeltaX) * params->Increment;
    params->DTerm = (2 * params->DeltaY - params->DeltaX) * params->Increment;

    params->XIndex = params->X1;
    params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t* params, int* x, int* y)
{
    if (params->UsingYIndex)
    {
        *y = params->XIndex;
        *x = params->YIndex;
        if (params->Flipped)
            *x = -*x;
    }
    else
    {
        *x = params->XIndex;
        *y = params->YIndex;
        if (params->Flipped)
            *y = -*y;
    }
}

int get_next_point(bresenham_param_t* params)
{
    if (params->XIndex == params->X2)
    {
        return 0;
    }
    params->XIndex += params->Increment;
    if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
        params->DTerm += params->IncrE;
    else
    {
        params->DTerm += params->IncrNE;
        params->YIndex += params->Increment;
    }
    return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double* map,
    int x_size,
    int y_size)

{
    bresenham_param_t params;
    int nX, nY;
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);

    //make sure the line segment is inside the environment
    if (x0 < 0 || x0 >= x_size ||
        x1 < 0 || x1 >= x_size ||
        y0 < 0 || y0 >= y_size ||
        y1 < 0 || y1 >= y_size)
        return 0;

    ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
    ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

    //iterate through the points on the segment
    get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
    do {
        get_current_point(&params, &nX, &nY);
        if (map[GETMAPINDEX(nX, nY, x_size, y_size)] == 1)
            return 0;
    } while (get_next_point(&params));

    return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double* map,
    int x_size, int y_size)
{
    double x0, y0, x1, y1;
    int i;

    //iterate through all the links starting with the base
    x1 = ((double)x_size) / 2.0;
    y1 = 0;
    for (i = 0; i < numofDOFs; i++)
    {
        //compute the corresponding line segment
        x0 = x1;
        y0 = y1;
        x1 = x0 + LINKLENGTH_CELLS * cos(2 * PI - angles[i]);
        y1 = y0 - LINKLENGTH_CELLS * sin(2 * PI - angles[i]);

        //check the validity of the corresponding line segment
        if (!IsValidLineSegment(x0, y0, x1, y1, map, x_size, y_size))
            return 0;
    }
    return 1;
}


struct state {
    state* prev_node;
    int num_node;
    double* angles;
    double cost;
};

double comp_distance(double* p, double* q) {
    double dist = 0.0;
    for (int i = 0; i < numofDOFs; i++) {
        dist += pow(p[i] - q[i], 2);
    }
    return sqrt(dist);
}


double planning_quality(double*** plan, int* planlength) {
    double dist = 0;

    for (int i = 0; i < *planlength - 1; i++) {
        double* curr_step = (*plan)[i];
        double* next_step = (*plan)[i + 1];
        dist += comp_distance(curr_step, next_step);
    }
    return dist;
}



pair<int, double*> random_configuration(double sample_goal) {

    double* angles = new double[numofDOFs];
    double random_value = ((double)rand() / (double)RAND_MAX);

    if (random_value >= sample_goal) {
        
            return make_pair(1, armgoal_anglesV_rad);
    }

    while (1) {
        
        for (int i = 0; i < numofDOFs; i++)
            angles[i] = (double)rand() * 2.0 * PI / (double)RAND_MAX;

        if (IsValidArmConfiguration(angles, numofDOFs, map, x_size, y_size)) {
            return make_pair(0, angles);

        }
    }
    
}


pair<double, state*> closest_node(double* q_rand, vector<state*> tree_graph) {
    state* q_node;
    double min_dist = 100000;

for (int it = 0; it < tree_graph.size(); it++) {

    double rms = 0;

    for (int j = 0; j < numofDOFs; j++) {
        rms += pow((((tree_graph)[it])->angles)[j] - (q_rand)[j], 2);
    }
    rms = sqrt(rms);
    if (rms < min_dist) {
        min_dist = rms;

        q_node = (tree_graph)[it];
    }

}

return make_pair(min_dist, q_node);

}








pair<int, int> valid_move(double closest_distance, state* closest_state, double* unit_vector, int steps) {
    double* tempJoint = new double[numofDOFs];
    int last_good_config;

    for (int i = 0; i <= steps; i++) {
        last_good_config = i;
        for (int j = 0; j < numofDOFs; j++) {
            tempJoint[j] = closest_state->angles[j] + (closest_distance / steps * i * unit_vector[j]);
        }
        if (!IsValidArmConfiguration(tempJoint, numofDOFs, map, x_size, y_size)) {

            return make_pair(0, last_good_config);
        }
    }

    return make_pair(1, last_good_config);
}



namespace rrt {
    int steps = 100;
    double sample_goal = 0.85;
    double epsilon = 0.5;
    vector<state*> tree_graph;


    void planner(double* map,
        int x_size,
        int y_size,
        double* armstart_anglesV_rad,
        double* armgoal_anglesV_rad,
        int numofDOFs,
        double*** plan,
        int* planlength)
    {



        clock_t start = clock();
        
        state* start_state = (state*)malloc(sizeof(state));
        start_state->angles = armstart_anglesV_rad;
        start_state->prev_node = 0;
        start_state->num_node = 1;
        tree_graph.push_back(start_state);

        int goal_reached = 0;
        
        int status = 0;
        int exceed = 0;

        state* new_node_main = (state*)malloc(sizeof(state));

        while (!goal_reached) {
            

            state* random_state = (state*)malloc(sizeof(state));
            pair<int, double*> random_node = random_configuration(sample_goal);
            random_state->angles = random_node.second;
            

            pair<double, state*> closest = closest_node(random_state->angles, tree_graph);
            double closest_distance = closest.first;
            state* closest_state = closest.second;


            double* unit_vector = new double[numofDOFs];
            for (int j = 0; j < numofDOFs; j++) {
                unit_vector[j] = ((random_state->angles[j] - closest_state->angles[j]) / closest_distance);
            }
            
            double duration1 = (clock() - start) / (double)CLOCKS_PER_SEC;
            if (duration1 > 100){
                printf("Took too long.");
                exceed = 1;
                return;
            }
            
            if (closest_distance < epsilon)
                status = 0;
            else {
                status = 1;
                closest_distance = epsilon;
            }
            pair<int, int> validity = (valid_move(closest_distance, closest_state, unit_vector, steps));
            int valid = validity.first;
            int ls = validity.second;
            if (!valid) {
                continue;
                
            }
            else {
                if (status == 0) {
                    random_state->prev_node = closest_state;
                    random_state->num_node = closest_state->num_node + 1;
                    tree_graph.push_back(random_state);

                    if (comp_distance(random_state->angles, armgoal_anglesV_rad) <= 0.0001) {
                        goal_reached = 1;
                        new_node_main = closest_state;
                    }
                    delete[] unit_vector;
                    continue;
                }
                else {
                    double* new_vec = new double[numofDOFs];
                    for (int i = 0; i < numofDOFs; i++) {
                        new_vec[i] = (epsilon / steps * ls* unit_vector[i])  + closest_state->angles[i];
                    }

                    state* new_node = (state*)malloc(sizeof(state));
                    new_node->angles = new_vec;
                    new_node->prev_node = closest_state;
                    new_node->num_node = closest_state->num_node + 1;
                    tree_graph.push_back(new_node);

                    if (comp_distance(new_node->angles, armgoal_anglesV_rad) <= 0.0001) {
                        goal_reached = 1;
                        new_node_main = new_node;
                    }
                    delete[] unit_vector;
                    delete[] new_vec;
                    continue;
                }
            }

        }

        if (goal_reached) {
            
            *plan = NULL;
            *planlength = 0;

            *planlength = new_node_main->num_node;
            *plan = (double**)malloc(*planlength * sizeof(double*));
            for (int i = *planlength - 1; i >= 0; i--) {
                (*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
                for (int j = 0; j < numofDOFs; j++) {
                    (*plan)[i][j] = new_node_main->angles[j];
                }
                new_node_main = new_node_main->prev_node;

            }

            double quality = planning_quality(plan, planlength);
            int num_nodes = tree_graph.size();
            double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
            free(new_node_main);
            free(start_state);
            //printf("Time elapsed: %f   |   Node number  %d    |    Path Quality  %f   |    Plan Length  %d \n", duration, num_nodes, quality, *planlength);
           
            return;
        }

    }

}
















namespace rrtconnect {

    int steps = 100;
    int tree_status = 0;
    int trap_status = 0;
    int connect_status = 0;
    double sample_goal = 0.85;
    double epsilon = 0.5;
    int exceed = 0;
    int goal_reached = 0;
    int status = 0;


    
    
    state* start_state = (state*)malloc(sizeof(state));
    state* goal_state = (state*)malloc(sizeof(state));
    vector<state*>tree_graphA ;
    vector<state*>tree_graphB ;
    vector<state*>tree_graph ;
    state* new_node_main = (state*)malloc(sizeof(state));

    state* treeA_endnode= (state*)malloc(sizeof(state));
    state* treeB_endnode = (state*)malloc(sizeof(state));
    vector<state*>path;
      
    void initialize() {
        start_state->angles = armstart_anglesV_rad;
        start_state->prev_node = 0;
        start_state->num_node = 1;
        tree_graphA.push_back(start_state);
        goal_state->angles = armgoal_anglesV_rad;
        goal_state->prev_node = 0;
        goal_state->num_node = 1;
        tree_graphB.push_back(goal_state);

    }



    void swap() {
        if (tree_status == 0)
           tree_status = 1;
        else
            tree_status = 0;
    }




    void extend(clock_t start) {

        while (1) {
            state* random_state = (state*)malloc(sizeof(state));
            pair<int, double*> random_node = random_configuration(sample_goal);
            random_state->angles = random_node.second;

            tree_graph.clear();
            if (tree_status == 0)
                tree_graph = tree_graphA;
            else
                tree_graph = tree_graphB;
            pair<double, state*> closest = closest_node(random_state->angles, tree_graph);
            double closest_distance = closest.first;
            state* closest_state = closest.second;
            

            double* unit_vector = new double[numofDOFs];
            for (int j = 0; j < numofDOFs; j++) {
                unit_vector[j] = ((random_state->angles[j] - closest_state->angles[j]) / closest_distance);
            }

            
           
            
            double duration1 = (clock() - start) / (double)CLOCKS_PER_SEC;
            if (duration1 > 5) {
                
                exceed = 1;
                return;
            }

            
            if (closest_distance < epsilon)
                status = 0;
            else {
                status = 1;
                closest_distance = epsilon;
            }
            pair<int, int> validity = (valid_move(closest_distance, closest_state, unit_vector, steps));
            int valid = validity.first;
            int ls = validity.second;
            if (!valid) {
                if (ls <= 1) {
                    swap();
                    
                }
                delete[] unit_vector;
                continue;
               
            }
            else {
                if (status == 0) {
                    random_state->prev_node = closest_state;
                    random_state->num_node = closest_state->num_node + 1;
                    new_node_main = random_state;
                    if (tree_status == 0) {
                        tree_graphA.push_back(random_state);
                    }
                    else {
                        tree_graphB.push_back(random_state);
                    }

                    
                    delete[] unit_vector;
                    return;
                }
                else {
                    double* new_vec = new double[numofDOFs];
                    for (int i = 0; i < numofDOFs; i++) {
                        new_vec[i] = unit_vector[i] * epsilon + closest_state->angles[i];
                    }

                    state* new_node = (state*)malloc(sizeof(state));;
                    new_node->angles = new_vec;
                    new_node->prev_node = closest_state;
                    new_node->num_node = closest_state->num_node + 1;
                    new_node_main = new_node;
                    if (tree_status == 0) {
                        tree_graphA.push_back(new_node);
                    }
                    else {
                        tree_graphB.push_back(new_node);
                    }

                    delete[] unit_vector;
                    delete[] new_vec;
                    return;
                }
            }

        }


    }

    void connect() {
        
        double* unit_vector = new double[numofDOFs];
        double* q_rand_angles = new double[numofDOFs];
        double* new_angles = new double[numofDOFs];
        
        
        if (tree_status == 0)
            tree_graph = tree_graphA;
        else
            tree_graph = tree_graphB;
        pair<double, state*> closest = closest_node(new_node_main->angles, tree_graph);
        double closest_distance = closest.first;
        state* closest_state = closest.second;
        q_rand_angles = new_node_main->angles;


        for (int j = 0; j < numofDOFs; j++) {
            unit_vector[j] = ((q_rand_angles[j] - closest_state->angles[j]) / closest_distance);
        }

        
        pair<int, int> validity = (valid_move(closest_distance, closest_state, unit_vector, 5000));
        int valid = validity.first;
        trap_status = validity.second;

        
        if (!valid) {
            
                return;
            
        }
        else {
            connect_status = 1;
            if (tree_status == 0) {
                treeA_endnode = closest_state;
                treeB_endnode = new_node_main;
            }
            else {
                treeB_endnode = closest_state;
                treeA_endnode = new_node_main;
            }
            
            return;

        }

    }


    void backtrack(double*** plan, int* planlength, clock_t start) {
        int treeA_size = treeA_endnode->num_node;
        
        *planlength = treeA_endnode->num_node + treeB_endnode->num_node;
        *plan = (double**)malloc(*planlength * sizeof(double*));

        for (int i = treeA_size; i < *planlength; i++) {
            (*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
            for (int j = 0; j < numofDOFs; j++) {
                (*plan)[i][j] = treeB_endnode->angles[j];
            }
            treeB_endnode = treeB_endnode->prev_node;
        }
        for (int i = treeA_size - 1; i >= 0; i--) {
            (*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
            for (int j = 0; j < numofDOFs; j++) {
                (*plan)[i][j] = treeA_endnode->angles[j];
            }
            treeA_endnode = treeA_endnode->prev_node;

        }
        double quality = planning_quality(plan, planlength);
        int num_nodes = tree_graphA.size() + tree_graphB.size();
        double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
        tree_graphA.clear();
        tree_graphB.clear();
        tree_graph.clear();
        free(new_node_main);
        free(start_state);
        free(goal_state);
        free(treeA_endnode);
        free(treeB_endnode);
        //printf("RRT CONNECT \n");
        //printf("Time elapsed: %f   |   Node number  %d    |    Path Quality  %f    | Plan Length %d \n", duration, num_nodes, quality, *planlength);
        
      
    }

    void planner(
        double* map,
        int x_size,
        int y_size,
        double* armstart_anglesV_rad,
        double* armgoal_anglesV_rad,
        int numofDOFs,
        double*** plan,
        int* planlength)
    {
        clock_t start = clock();
        initialize();
        while ((!connect_status) ) {      
            extend(start);
            if (exceed == 1)
                return;
            swap();
            connect();
        }
        
        if (connect_status) {
            backtrack(plan, planlength, start);
        }
       

        double duration1 = (clock() - start) / (double)CLOCKS_PER_SEC;
        if (duration1 > 5) {
            printf("Planner too too long");
            return;
        }
    }
}






namespace rrtstar {
    
    pair<int, double*> random_configuration1(double sample_goal, int goal_reach) {

        double* angles = new double[numofDOFs];
        double random_value = ((double)rand() / (double)RAND_MAX);

        if ((random_value >= sample_goal) && (!goal_reach)) {

            return make_pair(1, armgoal_anglesV_rad);
        }

        while (1) {

            for (int i = 0; i < numofDOFs; i++)
                angles[i] = (double)rand() * 2.0 * PI / (double)RAND_MAX;

            if (IsValidArmConfiguration(angles, numofDOFs, map, x_size, y_size)) {
                return make_pair(0, angles);

            }
        }

    }
      
    
    
    

   
    void planner(double* map, int x_size, int y_size,
        double* armstart_anglesV_rad, double* armgoal_anglesV_rad, int numofDOFs,
        double*** plan, int* planlength) {

        clock_t start = clock();

        
        *plan = NULL;
        *planlength = 0;

        
        double epsilon = 0.5;
        int steps = 500;
        int goal_reached = 0;
        int goal_once = 0;
        int num_samples = 0;
        
        double sample_goal = 0.85;


        state* start_state = (state*)malloc(sizeof(state));
        state* goal_state = (state*)malloc(sizeof(state));
        vector<state*>tree_graph;
        start_state->angles = armstart_anglesV_rad;
        start_state->prev_node = 0;
        start_state->num_node = 1;
        start_state->cost = 0;
        tree_graph.push_back(start_state);
        

        double* random_angle;
        state* closest_state;
       
        
        while (1) {


            if (goal_once) {
                num_samples++;
                if (num_samples>50)
                    break;
            }


            vector<state*> neigh_pack;
            vector<double> corr_dist;
            random_angle = (double*)malloc(numofDOFs * sizeof(double));
            
            pair<int, double*> random_node = random_configuration1(sample_goal, goal_once);
            for (int i = 0; i < numofDOFs; i++) {
                random_angle[i] = random_node.second[i];
            }
            goal_reached = random_node.first;
            
            double duration1 = (clock() - start) / (double)CLOCKS_PER_SEC;
            if (duration1 > 100) {
                printf("Planner too too long");
                return;
            }
            

           
            double radius = min(pow(( log(tree_graph.size()) / tree_graph.size())*500, (1.0 / numofDOFs)), epsilon);
            double closest_distance = (pow(2 * PI, 2) * numofDOFs);

            for (auto it = tree_graph.begin(); it != tree_graph.end(); it++) {


                double currNeighborDistance = comp_distance((*it)->angles, random_angle);
                if (currNeighborDistance <= radius) {
                    neigh_pack.push_back(*it);
                    corr_dist.push_back(currNeighborDistance);
                }

                if (currNeighborDistance < closest_distance) {
                    closest_state = (*it);
                    closest_distance = currNeighborDistance;
                }
                
            }


                                                                              
            double* unit_vector = new double[numofDOFs];
            for (int j = 0; j < numofDOFs; j++) {
                unit_vector[j] = ((random_angle[j] - closest_state->angles[j]) / closest_distance);
            }
            if (closest_distance > epsilon) {
                goal_reached = 0;
                for (int j = 0; j < numofDOFs; j++) {
                    random_angle[j] = closest_state->angles[j] + epsilon * unit_vector[j];
                }
                closest_distance = epsilon;
            }


            int valid = (valid_move(closest_distance, closest_state, unit_vector, steps)).first;

            if (valid) {

               

                state* cheap_state = closest_state;
                double cheapest = closest_state->cost + closest_distance;
                for (int i = 0; i < neigh_pack.size(); i++) {
                    double* unit_vector1 = new double[numofDOFs];
                    for (int j = 0; j < numofDOFs; j++) {
                        unit_vector1[j] = ((random_angle[j] - ((neigh_pack)[i]->angles)[j]) / corr_dist[i]);
                    }
                    int valid1 = (valid_move(corr_dist[i], (neigh_pack)[i], unit_vector, steps)).first;
                    

                    if (valid) {
                        double currCost = (neigh_pack)[i]->cost + corr_dist[i];
                        if (currCost < cheapest) {
                            cheap_state = (neigh_pack)[i];
                            cheapest = currCost;
                        }
                    }
                }

                state* new_node = (state*)malloc(sizeof(state));
                new_node->cost = cheapest;
                new_node->angles = random_angle;
                new_node->prev_node = cheap_state;
                new_node->num_node = cheap_state->num_node + 1;
                tree_graph.push_back(new_node);

                for (int i = 0; i < neigh_pack.size(); i++) {
                    
                    double currCost = new_node->cost + corr_dist[i];
                    double* unit_vector2 = new double[numofDOFs];
                    for (int j = 0; j < numofDOFs; j++) {
                        unit_vector2[j] = ((random_angle[j] - ((neigh_pack)[i]->angles)[j]) / corr_dist[i]);
                    }
                    int valid2 = (valid_move(corr_dist[i], (neigh_pack)[i], unit_vector2, steps)).first;
                    if ((valid2) && currCost < (neigh_pack)[i]->cost) {
                        (neigh_pack)[i]->cost = currCost;
                        (neigh_pack)[i]->prev_node = new_node;
                        (neigh_pack)[i]->num_node = new_node->num_node + 1;
                    }
                }

                

                if (goal_reached) {
                    goal_state = new_node;
                    goal_once = 1; 
                   
                }
            }
        }
        
        
         
    
        vector<state*> critical_path;
        state* new_node_main = goal_state;
        while (new_node_main->angles != start_state->angles) {
            critical_path.push_back(new_node_main);
            new_node_main = new_node_main->prev_node;
        }
        critical_path.push_back(start_state);
       
        reverse(critical_path.begin(), critical_path.end());
        *planlength = critical_path.size();
        *plan = (double**)malloc(*planlength * sizeof(double*));
        for (int i = 0; i <= *planlength - 1; i++) {
            (*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
            for (int j = 0; j < numofDOFs; j++) {
                (*plan)[i][j] = (critical_path[i]->angles)[j];
            }
        }


        double quality = planning_quality(plan, planlength);
        int num_nodes = tree_graph.size();
        double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
                
        printf("Time elapsed: %f   |   Node number  %d    |    Path Quality  %f   |    Plan Length  %d \n", duration, num_nodes, quality, *planlength);
        
        return;
    }

    
}
















namespace prm {
    int steps = 100;
    int goal_reached = 0;
    int radius = 2;
    int neighbors = 3;
    int num_nodes = 5000;

    class prm_state {
    private:
        prm_state* prev_node;
        vector<prm_state*> knn;
        int num_node;
        double cost;
        double* angles;

    public:
        prm_state();
        prm_state(double* a);
        void make_edge(prm_state* node);
        int num_edges();
        double get_cost();
        double* get_angle();
        prm_state* get_prev_node();
        int get_num_node();
        void set_cost(double c);
        void set_prev_node(prm_state* p);
        void set_num_node(int n);
        vector<prm_state*> get_edge();
    };


    prm_state::prm_state() {}
    prm_state::prm_state(double* a) {
        prev_node = 0;
        cost = DBL_MAX;
        angles = (double*)malloc(numofDOFs * sizeof(double));
        for (int i = 0; i < numofDOFs; i++) {
            angles[i] = a[i];
        }
    }
    vector<prm_state*> prm_state::get_edge() { return knn; }
    void prm_state::make_edge(prm_state* node) { knn.push_back(node); }
    int prm_state::num_edges() { return knn.size(); }

    double prm_state::get_cost() { return cost; }
    double* prm_state::get_angle() { return angles; }
    int prm_state::get_num_node() { return num_node; }
    prm_state* prm_state::get_prev_node() { return prev_node; }

    void prm_state::set_cost(double c) { cost = c; }
    void prm_state::set_prev_node(prm_state* p) { prev_node = p; }
    void prm_state::set_num_node(int i) { num_node = i; }





    struct distance_heap {
        template <class T>
        bool operator()(T s1, T s2)
        {
            return s1.first < s2.first;
        }
    };
    struct cost_heap {
        template <class T>
        bool operator()(T s1, T s2)
        {
            return s1->get_cost() > s2->get_cost();
        }
    };




    vector<prm_state*>graph;
    prm_state* start_state = (prm_state*)malloc(sizeof(prm_state));
    prm_state* goal_state = (prm_state*)malloc(sizeof(prm_state));

    pair<int, prm_state*> random_configuration() {
        double* angles = new double[numofDOFs];

        while (1) {

            for (int i = 0; i < numofDOFs; i++)
                angles[i] = (double)rand() * 2.0 * PI / (double)RAND_MAX;


            if (IsValidArmConfiguration(angles, numofDOFs, map, x_size, y_size)) {
                prm_state* q_new = new prm_state(angles);
                return make_pair(0, q_new);

            }
        }
    }




    pair<int, int> valid_move(double closest_distance, prm_state* closest_state, double* unit_vector) {
        double* tempJoint = new double[numofDOFs];

        int ls;
       

        for (int i = 0; i <= steps; i++) {
            ls = i;
            for (int j = 0; j < numofDOFs; j++) {
                tempJoint[j] = (closest_state->get_angle())[j] + (closest_distance / steps * i * (unit_vector[j]));
            }
            if (!IsValidArmConfiguration(tempJoint, numofDOFs, map, x_size, y_size)) {
                return make_pair(0, ls);
            }
        }
        return make_pair(1, ls);
    }


    void extension(prm_state* q_rand) {

        priority_queue<pair<double, prm_state*>, vector<pair<double, prm_state*>>, distance_heap> neighbor_priority;
        double closest_distance;
        q_rand->set_num_node(graph.size());
        graph.push_back(q_rand);

        for (auto it = graph.begin(); it != graph.end(); it++) {
            closest_distance = comp_distance((*it)->get_angle(), q_rand->get_angle());
            double* unit_vector = new double[numofDOFs];
            if (closest_distance < radius) {

                for (int j = 0; j < numofDOFs; j++) {
                    unit_vector[j] = (((q_rand->get_angle())[j] - ((*it)->get_angle())[j]) / closest_distance);
                }

                int valid = (valid_move(closest_distance, *it, unit_vector)).first;

                if (!valid) {
                    continue;
                }
                if (valid) {
                    pair<double, prm_state*> neighor_cell = make_pair(closest_distance, (*it));
                    neighbor_priority.push(neighor_cell);

                }
            }
        }


        while (!neighbor_priority.empty()) {
            while (neighbor_priority.size() > neighbors)
                neighbor_priority.pop();
            prm_state* chosen = neighbor_priority.top().second;
            neighbor_priority.pop();
            q_rand->make_edge(chosen);
            chosen->make_edge(q_rand);
        }

        return;

    }



    void search_path(prm_state* start_state, prm_state* goal_state, double*** plan, int* planlength, clock_t start) {

        priority_queue<prm_state*, vector<prm_state*>, cost_heap> OPEN;
        vector<bool> EXPANDED(graph.size(), false);
        goal_reached = 0;
        start_state->set_cost(0);
        OPEN.push(start_state);
        prm_state* s;

        while ((OPEN.empty() == false) && (!EXPANDED[goal_state->get_num_node()])) {
            s = OPEN.top();
            OPEN.pop();


            if (!EXPANDED[s->get_num_node()]) {
                EXPANDED[s->get_num_node()] = 1;
                vector<prm_state*> connected_vertex = (s->get_edge());
                while (!connected_vertex.empty()) {
                    prm_state* i = connected_vertex.back();
                    connected_vertex.pop_back();
                    if (!EXPANDED[(i)->get_num_node()]) {
                        double cost = comp_distance((i)->get_angle(), s->get_angle()) + s->get_cost();
                        if (cost < (i)->get_cost()) {
                            (i)->set_cost(cost);
                            (i)->set_prev_node(s);
                        }
                        OPEN.push(i);
                    }
                }


            }
        }
        if (EXPANDED[goal_state->get_num_node()])
            goal_reached = 1;
        if (goal_reached) {
            vector<prm_state*> critical_path;
            prm_state* new_node_main = goal_state;
            while ((new_node_main->get_angle() != start_state->get_angle())) {
                critical_path.push_back(new_node_main);
                new_node_main = new_node_main->get_prev_node();
            }
            critical_path.push_back(start_state);

            reverse(critical_path.begin(), critical_path.end());
            *planlength = critical_path.size();
            *plan = (double**)malloc((*planlength) * sizeof(double*));

            for (int i = 0; i <= *planlength - 1; i++) {
                (*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
                for (int j = 0; j < numofDOFs; j++) {
                    (*plan)[i][j] = (critical_path[i]->get_angle())[j];
                }
            }


            double quality = planning_quality(plan, planlength);

            double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
            
            //printf(" PRM Time elapsed: %f   |   PlanLength  %d \n", duration,  *planlength);
           
        }
        return;


    }



    void planner(
        double* map,
        int x_size,
        int y_size,
        double* armstart_anglesV_rad,
        double* armgoal_anglesV_rad,
        int numofDOFs,
        double*** plan,
        int* planlength)
    {
        
        clock_t start1= clock();
        prm_state* start_state = new prm_state(armstart_anglesV_rad);
        prm_state* goal_state = new prm_state(armgoal_anglesV_rad);
        start_state->set_num_node(0);
        start_state->set_cost(0);

        prm_state* q_rand = random_configuration().second;
        q_rand->set_num_node(0);
        graph.push_back(q_rand);
        for (int i = 1; i < num_nodes; i++) {
            prm_state* q_rand = random_configuration().second;
            extension(q_rand);
        }
        
        
        
        extension(start_state);
        extension(goal_state);
        
        if ((start_state->num_edges() > 0) && (goal_state->num_edges() > 0)) {
            clock_t start = clock();
            
            search_path(start_state, goal_state, plan, planlength, start);
        }
        else {
            printf("THe graph built does not have desired path");
            return;
        }

    }






















}
















static void starter_planner(
    double* map,
    int x_size,
    int y_size,
    double* armstart_anglesV_rad,
    double* armgoal_anglesV_rad,
    int numofDOFs,
    double*** plan,
    int* planlength)
{
    //no plan by default
    *plan = NULL;
    *planlength = 0;
    srand(time(NULL));
    
    if (planner_id == 0) {
        
        if (!IsValidArmConfiguration(armgoal_anglesV_rad, numofDOFs, map, x_size, y_size)) {
            printf("Goal nor reached because Infeasible  goal state. No plan");
            return;
        }
        if (!IsValidArmConfiguration(armstart_anglesV_rad, numofDOFs, map, x_size, y_size)) {
            printf("Goal nor reached because Infeasible  goal state. No plan");
            return;
        }
        rrt::planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
            
            
        
    }

    if (planner_id == 1) {

        if (!IsValidArmConfiguration(armgoal_anglesV_rad, numofDOFs, map, x_size, y_size)) {
            printf("Goal nor reached because Infeasible  goal state. No plan");
            return;
        }
        if (!IsValidArmConfiguration(armstart_anglesV_rad, numofDOFs, map, x_size, y_size)) {
            printf("Goal nor reached because Infeasible  goal state. No plan");
            return;
        }
        rrtconnect::planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
    }
    
    if (planner_id == 2) {
        if (!IsValidArmConfiguration(armgoal_anglesV_rad, numofDOFs, map, x_size, y_size)) {
            printf("Goal nor reached because Infeasible  goal state. No plan");
            return;
        }
        if (!IsValidArmConfiguration(armstart_anglesV_rad, numofDOFs, map, x_size, y_size)) {
            printf("Goal nor reached because Infeasible  goal state. No plan");
            return;
        }
        rrtstar::planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
    }
    
    if (planner_id == 3) {
        if (!IsValidArmConfiguration(armgoal_anglesV_rad, numofDOFs, map, x_size, y_size)) {
            printf("Goal nor reached because Infeasible  goal state. No plan");
            return;
        }
        if (!IsValidArmConfiguration(armstart_anglesV_rad, numofDOFs, map, x_size, y_size)) {
            printf("Goal nor reached because Infeasible  goal state. No plan");
            return;
        }
        prm::planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, plan, planlength);
    }

    return;
}

//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])

{

    /* Check for proper number of arguments */
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:planner:invalidNumInputs",
            "Four input arguments required.");
    }
    else if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:planner:maxlhs",
            "One output argument required.");
    }

    /* get the dimensions of the map and the map matrix itself*/
    x_size = (int)mxGetM(MAP_IN);
    y_size = (int)mxGetN(MAP_IN);
    map = mxGetPr(MAP_IN);

    /* get the start and goal angles*/
    numofDOFs = (int)(MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if (numofDOFs <= 1) {
        mexErrMsgIdAndTxt("MATLAB:planner:invalidnumofdofs",
            "it should be at least 2");
    }
    armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))) {
        mexErrMsgIdAndTxt("MATLAB:planner:invalidnumofdofs",
            "numofDOFs in startangles is different from goalangles");
    }
    armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);

    //get the planner id
    planner_id = (int)*mxGetPr(PLANNER_ID_IN);
    if (planner_id < 0 || planner_id > 4) {
        mexErrMsgIdAndTxt("MATLAB:planner:invalidplanner_id",
            "planner id should be between 0 and 3 inclusive");
    }

    
    double** plan = NULL;
    int planlength = 0;
    
    
    starter_planner(map, x_size, y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);

   
    if (planlength > 0)
    {
        PLAN_OUT = mxCreateNumericMatrix((mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL);
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int i, j;
        for (i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j * planlength + i] = plan[i][j];
            }
        }
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix((mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL);
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for (j = 0; j < numofDOFs; j++)
        {
            plan_out[j] = armstart_anglesV_rad[j];
        }
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix((mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL);
    int* planlength_out = (int*)mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;


    return;

}





