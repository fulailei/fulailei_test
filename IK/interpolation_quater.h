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


void RPYToRot(double roll, double pitch, double yaw, double R[3][3], enum RotationStatus choose_R); // ŷ����ת��Ϊ��ת����
void RotToQuaternion(double R[3][3], double q[4]);     // ��ת����ת����Ԫ��
void QuaternionToRot(double q[4], double R[3][3]);     //��Ԫ��ת����ת����
void AxisAngToQuaternion(double omg[3], double theta, double q[4]);//��ת��ͽǣ�ת��Ϊ��Ԫ��
void RotToAxisAng(double R[3][3], double omghat[3], double *theta);//��ת����ת��Ϊ��ת�����ת��
//vector<Position> pose_initial_goal(struct Position initial, struct Position goal, int num, int ME[3], struct DH dh1, enum RotationStatus choose_R);//��ʼֵ��Ŀ��ֵ��λ��+ŷ���ǣ�����ֵ����
//vector<joint> joint_initial_goal(struct joint initial, struct joint goal, int num);//��ʼֵ��Ŀ��ֵ���ؽڽǶȣ�����ֵ����
void quat2euler(double q[4], double p[3]);
