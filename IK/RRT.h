#ifndef RRT_H
#define RRT_H
#endif
#include <math.h> 
#include "stdafx.h"
#include "IK_group.h"
#include <iostream>
#include <Eigen/Dense>
#include<ctime>
#include<vector>
#include"interpolation_quater.h"
using namespace Eigen;
using namespace std;


vector<joint> motion_planning_RRT(struct joint initial, struct joint goal);