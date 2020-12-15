#pragma once
#ifndef interpolation_quater_H
#define interpolation_quater_H
#endif
#include <math.h> 
#include "stdafx.h"
//#include "IK_group.h"
#include <iostream>
#include <Eigen/Dense>
#include<ctime>
#include<vector>
using namespace Eigen;
using namespace std;


void RPYToRot(double roll, double pitch, double yaw, double R[3][3], enum RotationStatus choose_R); // 欧拉角转换为旋转矩阵
void RotToQuaternion(double R[3][3], double q[4]);     // 旋转矩阵转换四元数
void QuaternionToRot(double q[4], double R[3][3]);     //四元数转换旋转矩阵
void AxisAngToQuaternion(double omg[3], double theta, double q[4]);//旋转轴和角，转换为四元数
void RotToAxisAng(double R[3][3], double omghat[3], double *theta);//旋转矩阵转换为旋转轴和旋转角
//vector<Position> pose_initial_goal(struct Position initial, struct Position goal, int num, int ME[3], struct DH dh1, enum RotationStatus choose_R);//初始值，目标值（位置+欧拉角），插值个数
//vector<joint> joint_initial_goal(struct joint initial, struct joint goal, int num);//初始值，目标值（关节角度），插值个数
void quat2euler(double q[4], double p[3]);
