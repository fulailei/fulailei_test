
#ifndef new_IK_group_H
#define new_IK_group_H
#endif
#include <math.h>
#include "stdafx.h"
#include <iostream>
#include <Eigen/Dense>
#include <array>
#include<cmath>
#include<Eigen/Core> 
#include<Eigen/SVD>
#include<vector>

struct DH
{

	float alpha1;
	float a1;
	float theta1;
	float d1;
	float alpha2;
	float a2;
	float theta2;
	float d2;
	float alpha3;
	float a3;
	float theta3;
	float d3;
	float alpha4;
	float a4;
	float theta4;
	float d4;
	float alpha5;
	float a5;
	float theta5;
	float d5;
	float alpha6;
	float a6;
	float theta6;
	float d6;
};
struct Position
{
	float x;
	float y;
	float z;
	float Euler_z;
	float Euler_y;
	float Euler_x;
	float Quaternion[4];
};
struct joint         //六个关节角度
{
	float joint1;
	float joint2;
	float joint3;
	float joint4;
	float joint5;
	float joint6;
};

enum RotationStatus1
{
	RXYZn, RXZYn, RYXZn, RYZXn, RZXYn, RZYXn
};
Eigen::Matrix<float, 4, 4> pose2rot_euler(struct Position pose, struct DH dh1, enum rotationstatus1 choose_r);
Eigen::Matrix<float, 4, 4> pose2rot_quaternion(struct Position pose, struct DH dh1);
bool theta1_3(float pose_x, float pose_y, float pose_z, std::vector<float>& vec_angles, struct DH dh1, int joint1_d, int joint3_d);
void theta4_5(float theta1, float theta2, float theta3, Eigen::Matrix4f pose_rot, std::vector<float>& vec_angles, struct DH dh1, int joint4_d);
bool getIKPose_UR(Position pose, joint & vec_angles, float ME[3], DH dh1, float theta[6][8]);
float getikpose(struct Position pose, std::vector<float>& vec_angles, int me[3], struct DH dh1, int choose_type);

void AxisAngToQuaternion1(float omg[3], float theta, float q[4]);

//Matrix<double, 4, 4> QuaternionToRot(double q[4]);
