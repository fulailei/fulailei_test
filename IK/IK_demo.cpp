// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
// inv_link6逆解函数       pose_initial_goal 插补函数

#include <math.h> 
#include "stdafx.h"
#include <iostream>
#include <Eigen/Dense>
#include<ctime>
#include<vector>
#include"interpolation_quater.h"
#include "IK_group.h"
//#include"RRT.h"
using namespace Eigen;
using namespace std;
Eigen::MatrixXd pinv_eigen_based(Eigen::MatrixXd & origin, const float er = 0) {
	// 进行svd分解
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_holder(origin,
		Eigen::ComputeThinU |
		Eigen::ComputeThinV);
	// 构建SVD分解结果
	Eigen::MatrixXd U = svd_holder.matrixU();
	Eigen::MatrixXd V = svd_holder.matrixV();
	Eigen::MatrixXd D = svd_holder.singularValues();

	// 构建S矩阵
	Eigen::MatrixXd S(V.cols(), U.cols());
	S.setZero();

	for (unsigned int i = 0; i < D.size(); ++i) {

		if (D(i, 0) > er) {
			S(i, i) = 1 / D(i, 0);
		}
		else {
			S(i, i) = 0;
		}
	}

	// pinv_matrix = V * S * U^T
	return V * S * U.transpose();
}
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
struct Position      //位置+欧拉角
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
Matrix<double, 4, 4> QuaternionToRot(double q[4])       //四元数转换旋转矩阵
{
	MatrixXd R(4, 4);
	R(0, 0) = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
	R(0, 1) = 2.0*(q[1] * q[2] - q[0] * q[3]);
	R(0, 2) = 2.0*(q[0] * q[2] + q[1] * q[3]);
	R(0, 3) = 0;
	R(1, 0) = 2.0*(q[0] * q[3] + q[1] * q[2]);
	R(1, 1) = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
	R(1, 2) = 2.0*(q[2] * q[3] - q[0] * q[1]);
	R(1, 3) = 0;
	R(2, 0) = 2.0*(q[1] * q[3] - q[0] * q[2]);
	R(2, 1) = 2.0*(q[0] * q[1] + q[2] * q[3]);
	R(2, 2) = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
	R(2, 3) = 0;
	R(3, 0) = 0;
	R(3, 1) = 0;
	R(3, 2) = 0;
	R(3, 3) = 1;
	return R;
}
struct Quaternion1
{
	float x;
	float y;
	float z;
	float w;
};

struct rotMatrix
{
	float m11, m12, m13;
	float m21, m22, m23;
	float m31, m32, m33;
};


//fuction:旋转矩阵到四元数的转换
void rotMatrixToQuternion(Quaternion1 &q, rotMatrix &r)
{
	float tr = r.m11 + r.m22 + r.m33;
	float temp = 0.0;
	if (tr > 0.0)
	{
		temp = 0.5f / sqrtf(tr + 1);
		q.w = 0.25f / temp;
		q.x = (r.m23 - r.m32) * temp;
		q.y = (r.m31 - r.m13) * temp;
		q.z = (r.m12 - r.m21) * temp;
	}
	else
	{
		if (r.m11 > r.m22 && r.m11 > r.m33)
		{
			temp = 2.0f * sqrtf(1.0f + r.m11 - r.m22 - r.m33);
			q.w = (r.m32 - r.m23) / temp;
			q.x = 0.25f * temp;
			q.y = (r.m12 + r.m21) / temp;
			q.z = (r.m13 + r.m31) / temp;
		}
		else if (r.m22 > r.m33)
		{
			temp = 2.0f * sqrtf(1.0f + r.m22 - r.m11 - r.m33);
			q.w = (r.m13 - r.m31) / temp;
			q.x = (r.m12 + r.m21) / temp;
			q.y = 0.25f * temp;
			q.z = (r.m23 + r.m32) / temp;
		}
		else
		{
			temp = 2.0f * sqrtf(1.0f + r.m33 - r.m11 - r.m22);
			q.w = (r.m21 - r.m12) / temp;
			q.x = (r.m13 + r.m31) / temp;
			q.y = (r.m23 + r.m32) / temp;
			q.z = 0.25f * temp;
		}
	}
}
void rotToQuaternion1(double R[3][3], double q[4])
{
	float tr = R[0][0] + R[1][1] + R[2][2];
	float temp = 0.0;
	if (tr > 0.0)
	{
		temp = 0.5f / sqrtf(tr + 1);
		q[0] = 0.25f / temp;
		q[1] = (R[2][1] - R[1][2]) * temp;
		q[2] = (R[0][2] - R[2][0]) * temp;
		q[3] = (R[1][0] - R[0][1]) * temp;
	}
	else
	{
		if (R[0][0] > R[1][1] && R[0][0] > R[2][2])
		{
			temp = 2.0f * sqrtf(1.0f + R[0][0] - R[1][1] - R[2][2]);
			q[0] = (R[2][1] - R[1][2]) / temp;
			q[1] = 0.25f * temp;
			q[2] = (R[0][1] + R[1][0]) / temp;
			q[3] = (R[0][2] + R[2][0]) / temp;
		}
		else if (R[1][1] >R[2][2])
		{
			temp = 2.0f * sqrtf(1.0f + R[1][1] - R[0][0] - R[2][2]);
			q[0] = (R[0][2] - R[2][0]) / temp;
			q[1] = (R[0][1] + R[1][0]) / temp;
			q[2] = 0.25f * temp;
			q[3] = (R[1][2] + R[2][1]) / temp;
		}
		else
		{
			temp = 2.0f * sqrtf(1.0f + R[2][2] - R[0][0] - R[1][1]);
			q[0] = (R[1][0] - R[0][1]) / temp;
			q[1] = (R[0][2] + R[2][0]) / temp;
			q[2] = (R[1][2] + R[2][1]) / temp;
			q[3] = 0.25f * temp;
		}
	}
}
int main()
{
	struct DH dh2;//定义工业机器人（带球形腕关节）的dh表参数，d-h表的参数定义并不一定相同，实际使用中，需要调整d-h参数匹配该逆解。关节的旋转参考如下矩阵
	dh2.theta1 = 0;
	dh2.d1 = 495;
	dh2.a1 = 175; 
	dh2.alpha1 = 3.1415 / 2;
	dh2.theta2 = 0;
	dh2.d2 = 0; 
	dh2.a2 = 900; 
	dh2.alpha2 = 0;
	dh2.theta3 = 0;
	dh2.d3 = 175;
	dh2.a3 = 960;
	dh2.alpha3 = 3.1415 / 2;
	dh2.theta4 = 0; 
	dh2.d4 = 0; 
	dh2.a4 = 0; 
	dh2.alpha4 = -3.1415 / 2;
	dh2.theta5 = 0;
	dh2.d5 = 0; 
	dh2.a5 = 0;
	dh2.alpha5 = 3.1415 / 2;
	dh2.theta6 = 0;
	dh2.d6 = 135; 
	dh2.a6 = 0; 
	dh2.alpha6 = 0;
	////dh表参数定义

	Position pose_initial;
	//position pose_goal;
	float *vec_angles;
	int me[3] = { 0,0,0 };
	Position pose1;
    pose_initial.x = 806.29;
	pose_initial.y = 0;
	pose_initial.z = 1154;
	pose_initial.Euler_z = 3.1415;
	pose_initial.Euler_y = 1.0472;
	pose_initial.Euler_x = 3.1415;
	

	////起始位置//位置+欧拉角（弧度）
	
	////起始位置//
	//float rx;
	//float ry;
	//double q[4];
	//q[0] = 0.5;
	//q[1] = 0;
	//q[2] = 0.866;
	//q[3] = 0;
	//double p[3];
	////rx= 2.0*(q[1] * q[3] - q[0] * q[2]);
	//ry= 2.0*(q[0] * q[1] + q[2] * q[3]);


	////quat2euler(q, p);
	////cout << p[0] << endl;
	////cout << p[1] << endl;
	////cout << p[2] << endl;
	////cout << rx << endl;
	////cout << ry << endl;


	////目标位置//
	//pose_goal.x =1019.61;
	//pose_goal.y = 0;
	//pose_goal.z = 1217.5;
	//pose_goal.euler_z =3.1415;
	//pose_goal.euler_y = 1.0472;
	//pose_goal.euler_x = 3.1415;
	//pose_goal.quaternion[0] = q[0];
	//pose_goal.quaternion[1] = q[1];
	//pose_goal.quaternion[2] = q[2];
	//pose_goal.quaternion[3] = q[3];
	////目标位置//

	//joint joint_6link;//路径关节点
	//joint joint_initial;//初始关节点
	//joint joint_goal;//目标关节点
	//vector<joint> joint_path;//测试直线插值后的求逆的解
	//vector<joint> joint_path2;//测试关节插值
	//vector<joint> joint_path3;//测试直线插值
 //   int me[3] = {0,0,0};//轴1，3，4的正反向，机器人形态的控制参数
 //   int issolved;
	////rotationstatus choose;
	//double omg[3] = {0,1,0};
	//double theta = 0;
	double q1[4]= { 0.3705 ,   -0.0118076  , -0.917115930   , 0.050726};
	//AxisAngToQuaternion1(omg, theta, q1);
	double R[3][3];
	MatrixXd R1(4, 4);
	R1=QuaternionToRot(q1);
	cout << R1 << endl;
	//std::cout << QuaternionToRot(q1) << endl;
	/*R[0][0] = R1(0, 0);
	R[0][1] = R1(0, 1);
	R[0][2] = R1(0, 2);
	R[1][0] = R1(1, 0);
	R[1][1] = R1(1, 1);
	R[1][2] = R1(1, 2);
	R[2][0] = R1(2, 0);
	R[2][1] = R1(2, 1);
	R[2][2] = R1(2, 2);*/
	R[0][0] =0.567;
	R[0][1] = 0.8157;
	R[0][2] = 0.115;
	R[1][0] = 0.8015;
	R[1][1] = -0.5785;
	R[1][2] = 0.1513;
	R[2][0] = 0.1899;
	R[2][1] = 0.0064;
	R[2][2] = -0.9818;


	rotMatrix r1;
	Quaternion1 q2;
	r1.m11 = 0.7358;
	r1.m12 = 0.5651;
	r1.m13 = -0.373;
	r1.m21 = 0.3738;
	r1.m22 = 0.12;
	r1.m23 = 0.9197;
	r1.m31 = 0.5646;
	r1.m32 = -0.8162;
	r1.m33 = -0.1227;
	rotMatrixToQuternion(q2, r1);
	RotToQuaternion(R, q1);
	//cout << q1[0] << "::" << q1[1] << "::" << q1[2] << "::" << q1[3] << endl;
	cout << q2.w << "::" << q2.x << "::" << q2.y << "::" << q2.z << endl;

	rotToQuaternion1(R, q1);
	cout << q1[0] << "::R::" << q1[1] << "::" << q1[2] << "::" << q1[3] << endl;
	//q1[0] = q2.w;
	//q1[1] = q2.x;
	//q1[2] = q2.y;
	//q1[3] = q2.z;
	q1[1] = -q1[1];
	q1[2] = -q1[2];
	q1[3] = -q1[3];
	R1 = QuaternionToRot(q1);
	cout << R1 << endl;



	
    //joint_path=pose_initial_goal(pose_initial, pose_goal, 1000,ME,dh1,RZYX);//姿态插补程序，输入初始姿态结构体，输入目标姿态结构体，输入插补个数   ME=姿态选择，dh1=DH表，RZYX=欧拉角旋转顺序
	
	



	
    clock_t startTime = clock();
	

   clock_t endTime = clock();
   Eigen::MatrixXd M(4, 3) ;
   M << 1, 0, 1,
	   0, 1, 2,
	   2, 3, 4,
	   5, 6, 7;
   cout << M << endl;
   cout << "M::" << pinv_eigen_based(M) << endl;

//cout << "整个程序用时：" << double(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	system("pause");
	return 0;
}

