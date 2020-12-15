#include <math.h>
#include "stdafx.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <Eigen/Dense>
#include <array>
#include<cmath>
#include"interpolation_quater.h"
#include<Eigen/Core> 
#include<Eigen/SVD>
using namespace Eigen;
using namespace std;



struct joint         //六个关节角度
{
	float joint1;
	float joint2;
	float joint3;
	float joint4;
	float joint5;
	float joint6;
	float idx;
};


bool isobs_free(struct joint pose1,struct joint pose2)
{





	return true;

}

void random_jointangle(double theta[6],int num)
{
	srand(time(NULL)+num);
	for (int i=0;i<6;i++)
	theta[0]=(rand()%200-100);
	theta[1] = (rand() % 200-100);
	theta[2] = (rand() % 200-100);
	theta[3] = (rand() % 200-100);
	theta[4] = (rand() % 200-100);
	theta[5] = (rand() % 200-100);
	
}

float distance_joint(struct joint point1, struct joint point2)
{
	float distance;
	distance = sqrt(pow((point1.joint1 - point2.joint1), 2) + pow((point1.joint2 - point2.joint2), 2) + pow((point1.joint3 - point2.joint3), 2) + pow((point1.joint4 - point2.joint4), 2) + pow((point1.joint5 - point2.joint5), 2) + pow((point1.joint6 - point2.joint6), 2));

	return distance;

}

vector<joint> motion_planning_RRT(struct joint initial, struct joint goal)
{
	vector<joint>tree;
	joint point_rand;
	joint point;
	double theta[6];
	int idx=0;
	float distance;
	float min_distance=1000;
	point.joint1 =  initial.joint1;
	point.joint2 =  initial.joint2;
	point.joint3 =  initial.joint3;
	point.joint4 =  initial.joint4;
	point.joint5 =  initial.joint5;
	point.joint6 =  initial.joint6;
	tree.push_back(point);

	for (int i=0; i < 5000; i++) //迭代次数，需要设置
	{
		random_jointangle(theta,i);
		//cout << theta[1]<<"**"<<endl;
		point_rand.joint1 = theta[0]/100 * 3.14;   //步长设置参数
		point_rand.joint2 = theta[1]/100 * 3.14;
		point_rand.joint3 = theta[2] / 100 * 3.14;
		point_rand.joint4 = theta[3]/100 * 3.14;
		point_rand.joint5 = theta[4]/100 * 3.14;
		point_rand.joint6 = theta[5]/100 * 3.14;
		for ( int j= 0; j < tree.size(); j++)
		{
			distance=distance_joint(point, tree[j]);
			if (distance < min_distance)
				min_distance = distance;
			idx = j;
		}
		//cout << distance << endl;
		//cout << point_rand.joint1 << endl;
		point.joint1 = 0.01*(point_rand.joint1 - tree[idx].joint1) / distance + tree[idx].joint1;
		point.joint2 = 0.01*(point_rand.joint2 - tree[idx].joint2) / distance + tree[idx].joint2;
		point.joint3 = 0.01*(point_rand.joint3 - tree[idx].joint3) / distance + tree[idx].joint3;
		point.joint4 = 0.01*(point_rand.joint4 - tree[idx].joint4) / distance + tree[idx].joint4;
		point.joint5 = 0.01*(point_rand.joint5 - tree[idx].joint5) / distance + tree[idx].joint5;
		point.joint6 = 0.01*(point_rand.joint6 - tree[idx].joint6) / distance + tree[idx].joint6;


		if (isobs_free(point,tree[idx]))
		{
			point.idx = idx; 
			tree.push_back(point);
		
		if (distance_joint(goal, point) < 0.1)//目标参数，需要设置
		{
			break;
		}
         }



	}
	return tree;
}

