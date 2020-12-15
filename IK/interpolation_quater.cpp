
#include <math.h>
#include "stdafx.h"
#include <iostream>
#include <Eigen/Dense>
#include <array>
#include<vector>
#include "IK_group.h"
using namespace Eigen;
using namespace std;
#define			PI					3.14159265358979323846
//if the norm of vector is near zero(< 1.0E-6),regard as zero.
#define			ZERO_VECTOR			1.0E-6	
#define			ZERO_ELEMENT		1.0E-6		
#define			ZERO_ANGLE			1.0E-6
#define			ZERO_DISTANCE		1.0E-6
#define			ZERO_VALUE			1.0E-6
/**
* @brief 			Description: Computes the unit vector of Euler axis and rotation angle corresponding to rotation matrix.
* @param[in]		R				A rotation matrix.
* @param[out]		omghat			the unit vector of Euler axis .
* @param[out]		theta			the rotation angle.
* @return			No return value.
* @retval			0
* @note:			if  theta is zero ,the unit axis is undefined and set it as a zero vector [0;0;0].
*@warning:
*/
struct DH      //DH表
{

	float alpha1; float a1; float theta1; float d1;
	float alpha2; float a2; float theta2; float d2;
	float alpha3; float a3; float theta3; float d3;
	float alpha4; float a4; float theta4; float d4;
	float alpha5; float a5; float theta5; float d5;
	float alpha6; float a6; float theta6; float d6;
};
struct Position
{
	float x;
	float y;
	float z;
	float Euler_z;
	float Euler_y;
	float Euler_x;
};

struct joint
{
	float joint1;
	float joint2;
	float joint3;
	float joint4;
	float joint5;
	float joint6;
};

void RotToAxisAng(double R[3][3], double omghat[3], double *theta)//旋转矩阵转换为旋转轴和旋转角
{
	double tmp;
	double omg[3] = { 0 };
	double acosinput = (R[0][0] + R[1][1] + R[2][2] - 1.0) / 2.0;
	if (fabs(acosinput - 1.0)<ZERO_VALUE)
	{
		memset(omghat, 0, 3 * sizeof(double));
		*theta = 0.0;
	}
	else if (acosinput <= -1.0)
	{
		if ((1.0 + R[2][2]) >= ZERO_VALUE)
		{
			omg[0] = 1.0 / sqrt(2 * (1.0 + R[2][2]))*R[0][2];
			omg[1] = 1.0 / sqrt(2 * (1.0 + R[2][2]))*R[1][2];
			omg[2] = 1.0 / sqrt(2 * (1.0 + R[2][2]))*(1.0 + R[2][2]);
		}
		else if ((1.0 + R[1][1] >= ZERO_VALUE))
		{
			omg[0] = 1.0 / sqrt(2 * (1.0 + R[1][1]))*R[0][1];
			omg[1] = 1.0 / sqrt(2 * (1.0 + R[1][1]))*(1.0 + R[1][1]);
			omg[2] = 1.0 / sqrt(2 * (1.0 + R[1][1]))*R[2][1];
		}
		else
		{
			omg[0] = 1.0 / sqrt(2 * (1.0 + R[0][0]))*(1.0 + R[0][0]);
			omg[1] = 1.0 / sqrt(2 * (1.0 + R[0][0]))*R[1][0];
			omg[2] = 1.0 / sqrt(2 * (1.0 + R[0][0]))*R[2][0];
		}
		omghat[0] = omg[0];
		omghat[1] = omg[1];
		omghat[2] = omg[2];
		*theta = PI;
	}
	else
	{
		*theta = acos(acosinput);
		tmp = 2.0*sin(*theta);
		omghat[0] = (R[2][1] - R[1][2]) / tmp;
		omghat[1] = (R[0][2] - R[2][0]) / tmp;
		omghat[2] = (R[1][0] - R[0][1]) / tmp;

	}

	return;
}


/**
* @brief 			Description: Computes the unit quaternion corresponding to the Euler axis and rotation angle.
* @param[in]		omg				Unit vector of Euler axis.
* @param[in]		theta			Rotation angle.
* @param[in]		q				The unit quaternion
* @return			No return value.
* @note:
*@warning:
*/
void AxisAngToQuaternion(double omg[3], double theta, double q[4])//旋转轴和角，转换为四元数
{
	q[0] = cos(theta / 2.0);
	q[1] = omg[0] * sin(theta / 2.0);
	q[2] = omg[1] * sin(theta / 2.0);
	q[3] = omg[2] * sin(theta / 2.0);
	return;
}


/**
* @brief 			Description:Computes the unit quaternion corresponding to a rotation matrix.
* @param[in]		q				Unit quaternion.
* @param[out]		R				Rotation matrix.
* @return			No return value.
* @note:
* @warning:
*/
void QuaternionToRot(double q[4], double R[3][3])       //四元数转换旋转矩阵
{
	R[0][0] = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
	R[0][1] = 2.0*(q[1] * q[2] - q[0] * q[3]);
	R[0][2] = 2.0*(q[0] * q[2] + q[1] * q[3]);
	R[1][0] = 2.0*(q[0] * q[3] + q[1] * q[2]);
	R[1][1] = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
	R[1][2] = 2.0*(q[2] * q[3] - q[0] * q[1]);
	R[2][0] = 2.0*(q[1] * q[3] - q[0] * q[2]);
	R[2][1] = 2.0*(q[0] * q[1] + q[2] * q[3]);
	R[2][2] = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
	return;
}

/**
* @brief 			Description: Computes the unit quaternion corresponding to the rotation matrix.
* @param[in]		R				The rotation matrix.
* @param[out]		q				The unit quaternion.
* @return			No return value.
* @note:
* @warning:
*/
void RotToQuaternion(double R[3][3], double q[4])      // 旋转矩阵转换四元数
{
	double omghat[3];
	double theta;
	RotToAxisAng(R, omghat, &theta);
	AxisAngToQuaternion(omghat, theta, q);
	return;
}
enum RotationStatus
{
	RXYZ, RXZY, RYXZ, RYZX, RZXY, RZYX
};
void RPYToRot(double roll, double pitch, double yaw, double R[3][3],enum RotationStatus choose_R) // 欧拉角转换为旋转矩阵
{
	double z = yaw;
	double y = pitch;
	double x = roll;

	MatrixXd R_z(4, 4);
	MatrixXd R_y(4, 4);
	MatrixXd R_x(4, 4);
	MatrixXd R_T(4, 4);
	R_z << cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	R_y << cos(y), 0, sin(y), 0, 0, 1, 0, 0, -sin(y), 0, cos(y), 0, 0, 0, 0, 1;
	R_x << 1, 0, 0, 0, 0, cos(x), -sin(x), 0, 0, sin(x), cos(x), 0, 0, 0, 0, 1;
	switch (choose_R)
	{
	case(RXYZ):
		R_T = R_x*R_y*R_z;
		break;
	case(RXZY):
		R_T = R_x*R_z*R_y;
		break;
	case(RYXZ):
		R_T = R_y*R_x*R_z;
		break;
	case(RYZX):
		R_T = R_y*R_z*R_x;
		break;
	case(RZXY):
		R_T = R_z*R_x*R_y;
		break;
	case(RZYX):
		R_T = R_z*R_y*R_x;
		break;
	}
	R[0][0] = R_T(0, 0);
	R[0][1] = R_T(0, 1);
	R[0][2] = R_T(0, 2);
	R[1][0] = R_T(1, 0);
	R[1][1] = R_T(1, 1);
	R[1][2] = R_T(1, 2);
	R[2][0] = R_T(2, 0);
	R[2][1] = R_T(2, 1);
	R[2][2] = R_T(2, 2);
	return;
}
void rotationMatrixToEulerAngles(double ROT[3][3],double p[3]) //旋转矩阵转欧拉角
{
	float sy = sqrt(ROT[0][0] * ROT[0][0] + ROT[1][0] * ROT[1][0]);   
	bool singular = sy < 1e-6;
	// If    
	float x, y, z;   
	if (!singular) 
	{       
	x = atan2(ROT[2][1], ROT[2][2]);   
	y = atan2(-ROT[2][0], sy);       
	z = atan2(ROT[1][0], ROT[0][0]); 
	} 
	else 
	{       
	x = atan2(-ROT[1][2], ROT[1][1]);  
	y = atan2(-ROT[2][0], sy);     
	z = 0;   
	} 
	p[0] = x;
	p[1] = y;
	p[2] = z;
	return ;  
}

void quat2euler(double q[4], double p[3])
{
	double R[3][3];

	QuaternionToRot(q,R);

	rotationMatrixToEulerAngles(R,p);


}




vector<Position> initial_goal_pose(struct Position initial, struct Position goal,int num, int ME[3],struct DH dh1,enum RotationStatus choose_R)//初始值，目标值（位置+欧拉角），插值个数
{
	Position pose;
	vector<Position> path;
	double roll_initial,roll_goal;
	double pitch_initial,pitch_goal;
	double yaw_initial,yaw_goal;
	double rot_initial[3][3];
	double rot_goal[3][3];
	double delta_met[3][3];
	roll_initial = initial.Euler_z;
	pitch_initial = initial.Euler_y;
	yaw_initial = initial.Euler_x;
	roll_goal = goal.Euler_z;
	pitch_goal = goal.Euler_y;
	yaw_goal = goal.Euler_x;
	RPYToRot(roll_initial, pitch_initial, yaw_initial, rot_initial,choose_R);
	RPYToRot(roll_goal, pitch_goal, yaw_goal, rot_goal,choose_R);
	MatrixXf eig_initial(3,3);
	MatrixXf eig_goal(3,3);
	eig_initial << rot_initial[0][0], rot_initial[0][1], rot_initial[0][2], rot_initial[1][0], rot_initial[1][1], rot_initial[1][2], rot_initial[2][0], rot_initial[2][1], rot_initial[2][2];
	eig_goal<< rot_goal[0][0], rot_goal[0][1], rot_goal[0][2], rot_goal[1][0], rot_goal[1][1], rot_goal[1][2], rot_goal[2][0], rot_goal[2][1], rot_goal[2][2];
	MatrixXf delta;
	delta = eig_initial.inverse()*eig_goal;
	delta_met[0][0] = delta(0, 0);
	delta_met[0][1] = delta(0, 1);
	delta_met[0][2] = delta(0, 2);
	delta_met[1][0] = delta(1, 0);
	delta_met[1][1] = delta(1, 1);
	delta_met[1][2] = delta(1, 2);
	delta_met[2][0] = delta(2, 0);
	delta_met[2][1] = delta(2, 1);
	delta_met[2][2] = delta(2, 2);
	double omghat[3];
	double theta;
	RotToAxisAng(delta_met,omghat, &theta);//求出theta值，也就是初始位置和目标位置的转换矩阵在四元数中的旋转角。
	//位置直线插值//
	double x_initial, y_initial, z_initial;
	double x_goal, y_goal, z_goal;
	double x_delta, y_delta, z_delta;
	x_initial = initial.x;
	y_initial = initial.y;
	z_initial = initial.z;
	x_goal = goal.x;
	y_goal = goal.y;
	z_goal = goal.z;
	x_delta = (x_goal - x_initial)/num;
	y_delta = (y_goal - y_initial)/num;
	z_delta = (z_goal - z_initial)/num;
    double theta_delta;
	theta_delta = theta / num;
	double theta_num;
	double qua[4];
	double ROT_num[3][3];
	double ROT_point[3][3];
	double eul[3];
	for (int i = 0; i <= num; i++)
	{
		pose.x = x_initial + i*x_delta;
		pose.y = y_initial + i*y_delta;
		pose.z = z_initial + i*z_delta;
	
	//位置直线插值//
	
	//姿态插值//
		theta_num = i*theta_delta;
		AxisAngToQuaternion(omghat, theta_num, qua);
		QuaternionToRot(qua, ROT_num);
		MatrixXf eig_ROT_num(3, 3);
		MatrixXf ROT_p;
		eig_ROT_num << ROT_num[0][0], ROT_num[0][1], ROT_num[0][2], ROT_num[1][0], ROT_num[1][1], ROT_num[1][2], ROT_num[2][0], ROT_num[2][1], ROT_num[2][2];
		ROT_p = eig_initial*eig_ROT_num;
		ROT_point[0][0] = ROT_p(0, 0);
		ROT_point[0][1] = ROT_p(0, 1);
		ROT_point[0][2] = ROT_p(0, 2);
		ROT_point[1][0] = ROT_p(1, 0);
		ROT_point[1][1] = ROT_p(1, 1);
		ROT_point[1][2] = ROT_p(1, 2);
		ROT_point[2][0] = ROT_p(2, 0);
		ROT_point[2][1] = ROT_p(2, 1);
		ROT_point[2][2] = ROT_p(2, 2);
		rotationMatrixToEulerAngles(ROT_point, eul);
		pose.Euler_z = eul[0];
		pose.Euler_y = eul[1];
		pose.Euler_x = eul[2];
		path.push_back(pose);
	}
	return path;

}

vector<joint> initial_goal_joint_time(struct joint initial, struct joint goal, int num,int time)//初始值，目标值（关节角度），插值个数(三项插值)
{
	float t = time;
	float t1;
	joint angle;
	vector<joint>path;
	for (int i = 0; i <= num; i++) {
		t1 = i*t / num;
		angle.joint1 = initial.joint1 + (initial.joint1 - goal.joint1) * 3 / (-t*t)*t1*t1 + (initial.joint1 - goal.joint1) * 2 / (t*t*t)*t1*t1*t1;
		angle.joint2 = initial.joint2 + (initial.joint2 - goal.joint2) * 3 / (-t*t)*t1*t1 + (initial.joint2 - goal.joint2) * 2 / (t*t*t)*t1*t1*t1;
		angle.joint3 = initial.joint3 + (initial.joint3 - goal.joint3) * 3 / (-t*t)*t1*t1 + (initial.joint3 - goal.joint3) * 2 / (t*t*t)*t1*t1*t1;
		angle.joint4 = initial.joint4 + (initial.joint4 - goal.joint4) * 3 / (-t*t)*t1*t1 + (initial.joint4 - goal.joint4) * 2 / (t*t*t)*t1*t1*t1;
		angle.joint5 = initial.joint5 + (initial.joint5 - goal.joint5) * 3 / (-t*t)*t1*t1 + (initial.joint5 - goal.joint5) * 2 / (t*t*t)*t1*t1*t1;
		angle.joint6 = initial.joint6 + (initial.joint6 - goal.joint6) * 3 / (-t*t)*t1*t1 + (initial.joint6 - goal.joint6) * 2 / (t*t*t)*t1*t1*t1;
		path.push_back(angle);
	}
	return path;
}

vector<Position> initial_goal_line(struct Position initial, struct Position goal, int num, int ME[3], struct DH dh1, enum RotationStatus choose_R)//初始值，目标值（位置+欧拉角），插值个数
{
	Position pose;
	vector<Position> path;
	double theta;
	//位置直线插值//
	double x_initial, y_initial, z_initial;
	double x_goal, y_goal, z_goal;
	double x_delta, y_delta, z_delta;
	x_initial = initial.x;
	y_initial = initial.y;
	z_initial = initial.z;
	x_goal = goal.x;
	y_goal = goal.y;
	z_goal = goal.z;
	x_delta = (x_goal - x_initial) / num;
	y_delta = (y_goal - y_initial) / num;
	z_delta = (z_goal - z_initial) / num;
	for (int i = 0; i <= num; i++)
	{
		pose.x = x_initial + i*x_delta;
		pose.y = y_initial + i*y_delta;
		pose.z = z_initial + i*z_delta;
		pose.Euler_x = initial.Euler_x;
		pose.Euler_y = initial.Euler_y;
		pose.Euler_z = initial.Euler_z;
		path.push_back(pose);
	}
	return path;
}
		//位置直线插值//
vector<joint> initial_goal_joint_line(struct joint initial, struct joint goal, int num)//初始值，目标值（关节角度），插值个数(三项插值)
{
	joint angle;
	vector<joint>path;
	for (int i = 0; i <= num; i++) {
		angle.joint1 = initial.joint1 + i*(goal.joint1 - initial.joint1) / num;
		angle.joint2 = initial.joint2 + i*(goal.joint2 - initial.joint2) / num;
		angle.joint3 = initial.joint3 + i*(goal.joint3 - initial.joint3) / num;
		angle.joint4 = initial.joint4 + i*(goal.joint4 - initial.joint4) / num;
		angle.joint5 = initial.joint5 + i*(goal.joint5 - initial.joint5) / num;
		angle.joint6 = initial.joint6 + i*(goal.joint6 - initial.joint6) / num;
		path.push_back(angle);
	}
	return path;
}
