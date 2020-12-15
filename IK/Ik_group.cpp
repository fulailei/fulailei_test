#include <math.h>
#include "stdafx.h"
#include <iostream>
#include <Eigen/Dense>
#include <array>
#include<cmath>
#include"interpolation_quater.h"
#include<Eigen/Core> 
#include<Eigen/SVD>
using namespace Eigen;
using namespace std;
struct DH
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
struct joint         //六个关节角度
{
	float joint1;
	float joint2;
	float joint3;
	float joint4;
	float joint5;
	float joint6;
};

enum RotationStatus
{
	RXYZ, RXZY, RYXZ, RYZX, RZXY, RZYX
};

float *pose_ouler(struct Position pose,float *const buf,struct DH dh1,enum RotationStatus choose_R )     //  将欧拉角和位置转换为旋转矩阵
{



	float x, y, z,x1,y1,z1;
	x = pose.Euler_x;
	y = pose.Euler_y;
	z = pose.Euler_z;
	x1 = pose.x;
	y1 = pose.y;
	z1 = pose.z;
	MatrixXd R_z(4, 4);
	MatrixXd R_y(4, 4);
	MatrixXd R_x(4, 4);
	MatrixXd R_T(4, 4);
	R_z << cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	R_y << cos(y), 0, sin(y), 0, 0, 1, 0, 0, -sin(y), 0, cos(y), 0, 0, 0, 0, 1;
	R_x << 1, 0, 0, 0, 0, cos(x), -sin(x), 0, 0, sin(x), cos(x), 0, 0, 0, 0, 1;
	switch(choose_R)
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
	buf[0] = R_T(0,0);
	buf[1] = R_T(0,1);
	buf[2] = R_T(0,2);
	buf[3] = x1- R_T(0,2)*dh1.d6;
	buf[4] = R_T(1,0);
	buf[5] = R_T(1,1);
	buf[6] = R_T(1,2);
	buf[7] = y1- R_T(1,2)*dh1.d6;
	buf[8] = R_T(2,0);
	buf[9] = R_T(2,1);
	buf[10] = R_T(2,2);
	buf[11] = z1- R_T(2,2)*dh1.d6;
	buf[12] = 0;
	buf[13] = 0;
	buf[14] = 0;
	buf[15] = 1;
	return buf;
}


float *theta1_3(float x0, float y0, float z0, float *const buf,struct DH dh1,int joint1_d = 0,int joint3_d=0)       //根据位置求前三个关节的角度，1-3轴角度	float l1 = 0.7;float l2 = 0.6;h=0.4865
{
	float l1 ; float l2 ; float h ; float d ;
	l1 = dh1.a2; l2 = dh1.a3; h = dh1.d1; d = dh1.a1;
	float x1, y1;
	h = dh1.d1; d = dh1.a1;
	l1 = dh1.a2;
	l2 = sqrt(dh1.a3*dh1.a3 + dh1.d3*dh1.d3);
	x1 = dh1.d2 * (-y0 / (sqrt(x0*x0 + y0*y0)))+x0;
	y1 = dh1.d2 * (x0 / (sqrt(x0*x0 + y0*y0)))+y0;
	float x,y,zlink, ylink, j1;
	zlink = z0 - h;//0.4865
	j1 = atan2(y1, x1);//
	ylink = sqrt(x1*x1 + y1*y1) - d;//根据一轴的位置确定这个0.15部分的正负性
	if (joint1_d > 0)
	{
		x1 = -dh1.d2 * (-y0 / (sqrt(x0*x0 + y0*y0))) + x0;
		y1 = -dh1.d2 * (x0 / (sqrt(x0*x0 + y0*y0))) + y0;
		j1 = atan2(y1, x1)+3.1415;//
		ylink = sqrt(x1*x1 + y1*y1) +d;
	}
	x = ylink; 
	y = zlink;

	
	float theta1, theta2, theta3, theta4;//根据2-3的位置关系，确定2-3轴的一个组合
	theta1 = acos((x*x + y*y + l1*l1 - l2*l2) / (2 * l1*sqrt(x*x + y*y))) + atan(y / x);
	theta2 = acos((x*x + y*y + l2*l2 - l1*l1) / (2 * l2*sqrt(x*x + y*y))) + acos((x*x + y*y + l1*l1 - l2*l2) / (2 * l1*sqrt(x*x + y*y))) ;
	theta3 = -acos((x*x + y*y + l1*l1 - l2*l2) / (2 * l1*sqrt(x*x + y*y))) + atan(y / x);
	theta4 = -theta2;
	if (joint1_d > 0)
	{
		theta1 = -theta1+3.1415926;
		theta2 = -theta2;;
		theta3 = -theta3+3.1415926;
		theta4 = -theta4 ;
	}
	buf[0] = j1;
	buf[1] = 3.1415/2-theta1;
	buf[2] = theta2-3.1415/2;
	if (joint3_d > 0)
	{
		buf[1] = 3.1415 / 2 - theta3;
		buf[2] = theta4 - 3.1415 / 2;
	}
	buf[2] = buf[2] + atan2(dh1.d3, dh1.a3);
	return buf;//,theta2, theta3, theta4;//
}


float *theta4_5(float theta1, float theta2, float theta3,float *A,float *const buf,struct DH dh1,int joint4_d=0)   //求出前三个关节后，根据前三个关节的的末端位姿和整体全五轴的末端位姿的关系，求4-5轴的角度
{
	float j4, j5,j6;
	MatrixXf C1(4, 4);
	MatrixXf C2(4, 4);
	MatrixXf C3(4, 4);
	MatrixXf C4(4, 4);
	MatrixXf C5(4, 4);
	float z,z1,z2,al1,al2;
	MatrixXf B1(4, 4);
	MatrixXf B2(4, 4);
	MatrixXf B3(4, 4);
	MatrixXf B4(4, 4);
	MatrixXf thans_y(4, 4);
	thans_y << 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, dh1.d3, 0, 0, 0, 1;
	z = theta1 + 3.14159;
	float d ; 
	d = dh1.d1;
	float a ; 
	a = -dh1.a1;
	float al ; 
	al = dh1.alpha1;
	B1 << cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B2 << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, d, 0, 0, 0, 1;
	B3 << 1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B4 << 1, 0, 0, 0, 0, cos(al), -sin(al), 0, 0, sin(al), cos(al), 0, 0, 0, 0, 1;
	MatrixXf M0;
	M0 = B1*B2*B3*B4;
	z1 = theta2 + 3.1415 / 2;
	d = dh1.d2;
	a = dh1.a2;
	al1 = dh1.alpha2;
	B1 << cos(z1), -sin(z1), 0, 0, sin(z1), cos(z1), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B2 << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, d, 0, 0, 0, 1;
	B3 << 1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B4 << 1, 0, 0, 0, 0, cos(al1), -sin(al1), 0, 0, sin(al1), cos(al1), 0, 0, 0, 0, 1;
	MatrixXf M1;
	M1 = B1*B2*B3*B4;
	z2 = theta3 + 3.1415 / 2;
	d = dh1.d3;
	a = dh1.a3;
	al2 = dh1.alpha3;
	B1 << cos(z2), -sin(z2), 0, 0, sin(z2), cos(z2), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B2 << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, d, 0, 0, 0, 1;
	B3 << 1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B4 << 1, 0, 0, 0, 0, cos(al2), -sin(al2), 0, 0, sin(al2), cos(al2), 0, 0, 0, 0, 1;
	MatrixXf M2;
	M2 = B1*B2*B3*B4;
	MatrixXf M3;
	M3 = M0*M1*M2*thans_y;
	C1 = M3;
	C2 << *A, *(A + 1), *(A + 2), *(A + 3), *(A + 4), *(A + 5), *(A + 6),*(A + 7),*(A + 8),*(A + 9),*(A + 10),*(A + 11),*(A+12),*(A+13),*(A+14),*(A+15);
	C3 = C1.inverse()*C2;
	int j4_int;
	j4 = atan2(C3(1,2) , C3(0,2));//四轴+180度，五轴取反
	j5 = atan2(C3(0,2)*cos(j4) + C3(1,2)*sin(j4), C3(2,2));
	buf[0] = theta1;
	buf[1] = theta2;
	buf[2] = theta3;
	if (joint4_d > 0)
	{
		j4 = j4 -3.1415;
		j5 = -j5;
	}
	z2 = j4;
	B1 << cos(z2), -sin(z2), 0, 0, sin(z2), cos(z2), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B4 << cos(j5), 0, sin(j5), 0, 0, 1, 0, 0, -sin(j5), 0, cos(j5), 0, 0, 0, 0, 1;
	C4 = C1*B1*B4;
	C5 = C4.inverse()*C2;
	j6 = atan2(C5(1, 0),C5(1,1));
	if (joint4_d > 0)
	{
		j6 = j6+3.14;
	}
	j4_int = j4 / 3.1415;
	j4 = j4 - j4_int*3.1415926 * 2;
	//if (j6 < ((-135) / 180 * 3.14))
	//{
		//j6 = j6 + 6.28;
	//}
	buf[3] = j4+dh1.theta4;
	buf[4] = j5+dh1.theta5;
	buf[5] = j6+dh1.theta6;
	return buf;
}
float inv_link6(struct Position pose, float *const theta, int ME[3], struct DH dh1, enum RotationStatus choose_R)   //输入位置和欧拉角数组【位置+欧拉角】，输出六个轴关节角度theta[0]--theta[5]
{
	int m1, m2, m3;
	m1 = ME[0];
	m2 = ME[1];
	m3 = ME[2];
	float *first_step = new float[16];
	float *second_step = new float[5];
	float *third_step = new float[5];
	float *p1 = pose_ouler(pose, first_step, dh1,choose_R);
	float *p2 = theta1_3((*(first_step + 3)), (*(first_step + 7)), (*(first_step + 11)), second_step, dh1, m1, m2);
	float *p3 = theta4_5(*second_step, *(second_step + 1), *(second_step + 2), first_step, third_step, dh1, m3);

	theta[0] = *third_step+0;
	theta[1] = *(third_step + 1);
	theta[2] = *(third_step + 2);
	theta[3] = *(third_step + 3);
	theta[4] = *(third_step + 4);
	theta[5] = *(third_step + 5);
	delete first_step, second_step, third_step, p1, p2, p3;
	int A=0;
	if (isnan(theta[1])||isnan(theta[2])||isnan(theta[3])||isnan(theta[4])||isnan(theta[5]))
		{
			A = 1;
        }
	
	return A;
}
float obs_cost(float obs[3], float robot_points[3],float rad)//球形障碍物绝对值度量函数
{
	float cost;
	float distance;
	distance = sqrt((obs[0] - robot_points[0])*(obs[0] - robot_points[0]) + (obs[1] - robot_points[1])*(obs[1] - robot_points[1]) + (obs[2] - robot_points[2])*(obs[2] - robot_points[2]));
	cost = exp(rad - distance);
	return cost;
}
MatrixXd liman_mat(double num,MatrixXd a)    //光滑度黎曼度量函数
{
	MatrixXd C;
	MatrixXd A;
	C.setIdentity(num, num);
	for (int i = 0; i < num; i++)
	{
		C(i + 1, i) = -1;
	}
	A = C.transpose()*C;
	return A;
}
Matrix<double, 4, 4> DH_Matrix(float z, float a, float al, float d)//每两个轴的DH值转换的矩阵
{
	MatrixXd B1(4, 4);
	MatrixXd B2(4, 4);
	MatrixXd B3(4, 4);
	MatrixXd B4(4, 4);
	B1 << cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B2 << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, d, 0, 0, 0, 1;
	B3 << 1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	B4 << 1, 0, 0, 0, 0, cos(al), -sin(al), 0, 0, sin(al), cos(al), 0, 0, 0, 0, 1;
	MatrixXd M0;
	M0 = B1*B2*B3*B4;
	return M0;
}

Matrix<double, Dynamic, Dynamic> fK_link6(struct DH dh1, struct joint j1)//正运动学
{
	MatrixXd A0_1;
	MatrixXd A1_2;
	MatrixXd A2_3;
	MatrixXd A3_4;
	MatrixXd A4_5;
	MatrixXd A5_6;
	MatrixXd A0_6;
	MatrixXd thans_y(4, 4);
	thans_y << 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, dh1.d3, 0, 0, 0, 1;
	A0_1 = DH_Matrix(j1.joint1+3.14, -dh1.a1, dh1.alpha1, dh1.d1);
	A1_2 = DH_Matrix(j1.joint2+3.14/2, dh1.a2, dh1.alpha2, dh1.d2);
	A2_3 = DH_Matrix(j1.joint3+3.14/2, dh1.a3, dh1.alpha3, 0);
	A3_4 = DH_Matrix(j1.joint4, dh1.a4, dh1.alpha4, dh1.d4);
	A4_5 = DH_Matrix(j1.joint5, dh1.a5, dh1.alpha5, dh1.d5);
	A5_6 = DH_Matrix(j1.joint6, dh1.a6, dh1.alpha6, dh1.d6);
	A0_6 = A0_1*A1_2*A2_3*thans_y*A3_4*A4_5*A5_6;
	return A0_6;
}
Matrix<double,Dynamic,Dynamic> jacobi_mat(struct DH dh2,struct joint j2)//6关节的雅克比矩阵，3*6矩阵
{
	MatrixXd jaco;
	MatrixXd jaco_delta1;
	MatrixXd jaco_delta2;
	MatrixXd jaco_delta3;
	MatrixXd jaco_delta4;
	MatrixXd jaco_delta5;
	MatrixXd jaco_delta6;
	MatrixXd jaco_fin(3,6);
	jaco = fK_link6(dh2, j2);
	joint j2_joint1;
	joint j2_joint2;
	joint j2_joint3;
	joint j2_joint4;
	joint j2_joint5;
	joint j2_joint6;
	j2_joint1 = j2;
	j2_joint2 = j2;
	j2_joint3 = j2;
	j2_joint4 = j2;
	j2_joint5 = j2;
	j2_joint6 = j2;
	j2_joint1.joint1 = j2.joint1 + 0.0001;
	j2_joint2.joint2 = j2.joint2 + 0.0001;
	j2_joint3.joint3 = j2.joint3 + 0.0001;
	j2_joint4.joint4 = j2.joint4 + 0.0001;
	j2_joint5.joint5 = j2.joint5 + 0.0001;
	j2_joint6.joint6 = j2.joint6 + 0.0001;
	jaco_delta1 = fK_link6(dh2, j2_joint1);
	jaco_delta2 = fK_link6(dh2, j2_joint2);
	jaco_delta3 = fK_link6(dh2, j2_joint3);
	jaco_delta4 = fK_link6(dh2, j2_joint4);
	jaco_delta5 = fK_link6(dh2, j2_joint5);
	jaco_delta6 = fK_link6(dh2, j2_joint6);
	jaco_fin << jaco_delta1.block(0, 3, 3, 1) - jaco.block(0, 3, 3, 1), jaco_delta2.block(0, 3, 3, 1) - jaco.block(0, 3, 3, 1), jaco_delta3.block(0, 3, 3, 1) - jaco.block(0, 3, 3, 1),
		jaco_delta4.block(0, 3, 3, 1) - jaco.block(0, 3, 3, 1), jaco_delta5.block(0, 3, 3, 1) - jaco.block(0, 3, 3, 1), jaco_delta6.block(0, 3, 3, 1) - jaco.block(0, 3, 3, 1);
	return jaco_fin/0.0001;

}
float mat2(struct DH dh3, struct joint j3)//测试雅克比矩阵和正运动学
{
	MatrixXd jaco;
	MatrixXd jaco1;
	jaco = jacobi_mat( dh3,j3);
	jaco1 = fK_link6(dh3, j3);
	cout << jaco <<endl;
	cout << jaco1 << endl;
	return jaco(0, 3);
}
template<typename _Matrix_Type_>                  //求伪逆矩阵
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon =
	std::numeric_limits<double>::epsilon())
{
	Eigen::JacobiSVD< _Matrix_Type_ > svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}

struct Vec3d      //三坐标的结构体，用来表示点和向量
{
	double x, y, z;

	Vec3d()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}
	Vec3d(double dx, double dy, double dz)
	{
		x = dx;
		y = dy;
		z = dz;
	}
	void Set(double dx, double dy, double dz)
	{
		x = dx;
		y = dy;
		z = dz;
	}
};


//计算三点成面的法向量
void Cal_Normal_3D(const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, Vec3d &vn)
{
	//v1(n1,n2,n3);
	//平面方程: na * (x – n1) + nb * (y – n2) + nc * (z – n3) = 0 ;
	double na = (v2.y - v1.y)*(v3.z - v1.z) - (v2.z - v1.z)*(v3.y - v1.y);
	double nb = (v2.z - v1.z)*(v3.x - v1.x) - (v2.x - v1.x)*(v3.z - v1.z);
	double nc = (v2.x - v1.x)*(v3.y - v1.y) - (v2.y - v1.y)*(v3.x - v1.x);

	//平面法向量
	vn.Set(na, nb, nc);
}

void solveCenterPointOfCircle(const Vec3d& r1, const Vec3d& r2, const Vec3d& r3, Vec3d &r_n)//已知三个点求空间圆心坐标
{ 
double a1, b1, c1, d1;
double a2, b2, c2, d2;
double a3, b3, c3, d3; 
double na, nb, nc;
double x1 = r1.x, y1 = r1.y, z1 = r1.z; 
double x2 = r2.x, y2 = r2.y, z2 = r2.z; 
double x3 = r3.x, y3 = r3.y, z3 = r3.z;
a1 = (y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2); 
b1 = -(x1*z2 - x2*z1 - x1*z3 + x3*z1 + x2*z3 - x3*z2);
c1 = (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
d1 = -(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
a2 = 2 * (x2 - x1); 
b2 = 2 * (y2 - y1); 
c2 = 2 * (z2 - z1);
d2 = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2;
a3 = 2 * (x3 - x1);
b3 = 2 * (y3 - y1);
c3 = 2 * (z3 - z1); 
d3 = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3;
na = -(b1*c2*d3 - b1*c3*d2 - b2*c1*d3 + b2*c3*d1 + b3*c1*d2 - b3*c2*d1) / (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1);
nb = (a1*c2*d3 - a1*c3*d2 - a2*c1*d3 + a2*c3*d1 + a3*c1*d2 - a3*c2*d1) / (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1); 
nc = -(a1*b2*d3 - a1*b3*d2 - a2*b1*d3 + a2*b3*d1 + a3*b1*d2 - a3*b2*d1) / (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1); 
r_n.Set(na, nb, nc);
}
void unit_vector(const Vec3d& u1, Vec3d &u_n)   //普通矢量转换单位矢量
{
	double a, b, c,m;
	double ua, ub, uc;
	a = u1.x;
	b = u1.y;
	c = u1.z;
	m = sqrt(a*a + b*b + c*c);
	ua = a / m;
	ub = b / m;
	uc = c / m;
	u_n.Set(ua, ub, uc);
}
MatrixXd Rot_7(Vec3d &n1,Vec3d &n2,Vec3d &n3)  //取n1为零位，即第七轴的关节角度为0，求出7轴连接处的位姿。
{
	MatrixXd ROT7(4, 4);
	Vec3d sn;
	Vec3d rn;
	Vec3d xn;
	Cal_Normal_3D(n1, n2, n3, sn);
	solveCenterPointOfCircle(n1, n2, n3, rn);
	Vec3d pn;
	pn.x = sn.x + rn.x;
	pn.y = sn.y + rn.y;
	pn.z = sn.z + rn.z;
	Cal_Normal_3D(n1, rn, pn, xn);
	Vec3d x1;
	Vec3d y1;
	Vec3d z1;
	Vec3d k1;
	k1.x = rn.x - n1.x;
	k1.y = rn.y - n1.y;
	k1.z = rn.z - n1.z;
	unit_vector(xn, x1);
	unit_vector(k1, y1);
	unit_vector(sn, z1);

	ROT7 << x1.x, y1.x, z1.x, rn.x, x1.y, y1.y, z1.y, rn.y, x1.z, y1.z, z1.z, rn.z, 0, 0, 0, 1;
	return ROT7;
}
MatrixXd ROT_7_2(MatrixXd ROT6,Vec3d &p1,Vec3d &p2,Vec3d &p3)  //7轴连接处到末端抓手的转换矩阵
{
	MatrixXd ROT7_2;
	MatrixXd ROT7_1;
	MatrixXd B4(4,4);
      //B4 << 1, 0, 0, 0, 0, cos(1.57), -sin(1.57), 0, 0, sin(1.57), cos(1.57), 0, 0, 0, 0, 1;
		ROT7_2 = ROT6;
		ROT7_2(0, 3) = p1.x;
		ROT7_2(1, 3) = p1.y;
		ROT7_2(2, 3) = p1.z;

		ROT7_1=Rot_7(p1, p2, p3);

		


		//cout << "这个是"<<ROT7_1 << endl;
		return ROT7_1.inverse()*ROT7_2;//ROT7_1.inverse()


}
MatrixXd ROT_7_1(MatrixXd ROT6, Vec3d &p1, Vec3d &p2, Vec3d &p3)//法兰盘到7轴连接处的位姿
{
	MatrixXd ROT7_2;
	MatrixXd ROT7_1;

	ROT7_2 = ROT6;
	ROT7_2(0, 3) = p1.x;
	ROT7_2(1, 3) = p1.y;
	ROT7_2(2, 3) = p1.z;
	ROT7_1 = Rot_7(p1, p2, p3);
	return ROT6.inverse()*ROT7_1;


}
MatrixXd ROT6_7(double joint7, MatrixXd ROT7_1, MatrixXd ROT7_2)//通过三个点标定以后，获得了两个转换矩阵，再根据7轴的旋转角度，可以得出法兰盘到末端抓手的转换矩阵
{
	MatrixXd ROT(4,4);
	ROT << cos(joint7), -sin(joint7), 0, 0, sin(joint7), cos(joint7), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	//ROT << 1, 0, 0, 0, 0, cos(joint7), -sin(joint7), 0, 0, sin(joint7), cos(joint7), 0, 0, 0, 0, 1;
	return ROT7_1*ROT*ROT7_2;
}
MatrixXd RPYToROT(struct Position p1) // 欧拉角转换为旋转矩阵
{
	MatrixXd R(4, 4);
	double alpha = p1.Euler_z;
	double beta = p1.Euler_y;
	double gamma = p1.Euler_x;
	R(0,0) = cos(alpha)*cos(beta);
	R(0,1) = cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma);
	R(0,2) = cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma);
	R(1,0) = sin(alpha)*cos(beta);
	R(1,1) = sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma);
	R(1,2) = sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma);
	R(2,0) = -sin(beta);
	R(2,1) = cos(beta)*sin(gamma);
	R(2,2) = cos(beta)*cos(gamma);
	R(0, 3) = p1.x;
	R(1, 3) = p1.y;
	R(2, 3) = p1.z;
	R(3, 3) = 1;
	R(3, 0) = 0;
	R(3, 1) = 0;
	R(3, 2) = 0;
	return R;
}



int Singular_test(struct joint value)//判断五轴会不会卡死，五轴不能为零。//常用的方案可以在0位附近走关节。
{
	int i = 0;
	if (abs(value.joint5)<0.001)
	{
		i = 1;
	}
	return i;
}
