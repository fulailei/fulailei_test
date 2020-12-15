
#include "new_IK_group.h"
#define PI 3.1415926
using namespace Eigen;
using namespace std;
//if the norm of vector is near zero(< 1.0E-6),regard as zero.
#define			ZERO_VECTOR			1.0E-6	
#define			ZERO_ELEMENT		1.0E-6		
#define			ZERO_ANGLE			1.0E-6
#define			ZERO_DISTANCE		1.0E-6
#define			ZERO_VALUE			1.0E-6


void RotToAxisAng1(float R[3][3], float omghat[3], float *theta)//旋转矩阵转换为旋转轴和旋转角
{
	float tmp;
	float omg[3] = { 0 };
	float acosinput = (R[0][0] + R[1][1] + R[2][2] - 1.0) / 2.0;
	if (fabs(acosinput - 1.0)<ZERO_VALUE)
	{
		memset(omghat, 0, 3 * sizeof(float));
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
void AxisAngToQuaternion1(float omg[3], float theta, float q[4])//旋转轴和角，转换为四元数
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
Matrix<float, 4, 4> QuaternionToRot(float q[4])       //四元数转换旋转矩阵
{
	MatrixXf R(4, 4);
	R(0,0) = q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
	R(0,1) = 2.0*(q[1] * q[2] - q[0] * q[3]);
	R(0,2) = 2.0*(q[0] * q[2] + q[1] * q[3]);
	R(0, 3) = 0;
	R(1,0) = 2.0*(q[0] * q[3] + q[1] * q[2]);
	R(1,1) = q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3];
	R(1,2) = 2.0*(q[2] * q[3] - q[0] * q[1]);
	R(1, 3) = 0;
	R(2,0) = 2.0*(q[1] * q[3] - q[0] * q[2]);
	R(2,1) = 2.0*(q[0] * q[1] + q[2] * q[3]);
	R(2,2) = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3];
	R(2, 3) = 0;
	R(3, 0) = 0;
	R(3, 1) = 0;
	R(3, 2) = 0;
	R(3, 3) = 1;
	return R;
}

/**
* @brief 			Description: Computes the unit quaternion corresponding to the rotation matrix.
* @param[in]		R				The rotation matrix.
* @param[out]		q				The unit quaternion.
* @return			No return value.
* @note:
* @warning:
*/
void RotToQuaternion1(float R[3][3], float q[4])      // 旋转矩阵转换四元数
{
	float omghat[3];
	float theta;
	RotToAxisAng1(R, omghat, &theta);
	AxisAngToQuaternion1(omghat, theta, q);
	return;
}
Matrix<float, 4, 4> DH_Matrix1(float z, float a, float al, float d)//每两个轴的DH值转换的矩阵
{
	MatrixXf Rot_Z_thetai(4, 4);
	MatrixXf Trans_Z_di(4, 4);
	MatrixXf Trans_X_ai(4, 4);
	MatrixXf Rot_X_alphai(4, 4);
	Rot_Z_thetai << cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	Trans_Z_di << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, d, 0, 0, 0, 1;
	Trans_X_ai << 1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	Rot_X_alphai << 1, 0, 0, 0, 0, cos(al), -sin(al), 0, 0, sin(al), cos(al), 0, 0, 0, 0, 1;
	MatrixXf single_;
	single_ = Rot_Z_thetai*Trans_Z_di*Trans_X_ai*Rot_X_alphai;
	return single_;
}

Matrix<float, Dynamic, Dynamic> fK_link(struct DH dh1, std::vector<float>& vec_angles)//正运动学
{
	MatrixXf Rot0_1;
	MatrixXf Rot1_2;
	MatrixXf Rot2_3;
	MatrixXf Rot3_4;
	MatrixXf Rot4_5;
	MatrixXf Rot5_6;
	MatrixXf Rot0_6;
	MatrixXf Rot_y(4, 4);
	Rot_y << 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, dh1.d3, 0, 0, 0, 1;
	Rot0_1 = DH_Matrix1(vec_angles[0] + PI, -dh1.a1, dh1.alpha1, dh1.d1);
	Rot1_2 = DH_Matrix1(vec_angles[1] + PI / 2, dh1.a2, dh1.alpha2, dh1.d2);
	Rot2_3 = DH_Matrix1(vec_angles[2] + PI / 2, dh1.a3, dh1.alpha3, 0);
	Rot3_4 = DH_Matrix1(vec_angles[3], dh1.a4, dh1.alpha4, dh1.d4);
	Rot4_5 = DH_Matrix1(vec_angles[4], dh1.a5, dh1.alpha5, dh1.d5);
	Rot5_6 = DH_Matrix1(vec_angles[5], dh1.a6, dh1.alpha6, dh1.d6);
	Rot0_6 = Rot0_1*Rot1_2*Rot2_3*Rot_y*Rot3_4*Rot4_5*Rot5_6;
	return Rot0_6;
}



Matrix<float, 4, 4> pose2Rot_Euler(struct Position pose, struct DH dh1, enum RotationStatus1 choose_R)     //  将欧拉角和位置转换为旋转矩阵
{
	float x, y, z, x1, y1, z1;
	x = pose.Euler_x;
	y = pose.Euler_y;
	z = pose.Euler_z;
	x1 = pose.x;
	y1 = pose.y;
	z1 = pose.z;
	Eigen::MatrixXf R_z(4, 4);
	Eigen::MatrixXf R_y(4, 4);
	Eigen::MatrixXf R_x(4, 4);
	Eigen::MatrixXf R_T(4, 4);
	R_z << cos(z), -sin(z), 0, 0, sin(z), cos(z), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	R_y << cos(y), 0, sin(y), 0, 0, 1, 0, 0, -sin(y), 0, cos(y), 0, 0, 0, 0, 1;
	R_x << 1, 0, 0, 0, 0, cos(x), -sin(x), 0, 0, sin(x), cos(x), 0, 0, 0, 0, 1;
	switch (choose_R)
	{
	case(RXYZn):
		R_T = R_x*R_y*R_z;
		break;
	case(RXZYn):
		R_T = R_x*R_z*R_y;
		break;
	case(RYXZn):
		R_T = R_y*R_x*R_z;
		break;
	case(RYZXn):
		R_T = R_y*R_z*R_x;
		break;
	case(RZXYn):
		R_T = R_z*R_x*R_y;
		break;
	case(RZYXn):
		R_T = R_z*R_y*R_x;
		break;
	}
	cout << R_T << endl;

	R_T(0, 3) = x1 -R_T(0, 2)*dh1.d6;
	R_T(1,3)= y1 - R_T(1, 2)*dh1.d6;
	R_T(2,3) = z1 - R_T(2, 2)*dh1.d6;
	
	return R_T;
}
Matrix<float, 4, 4> pose2Rot_Quaternion(struct Position pose, struct DH dh1)     //  将欧拉角和位置转换为旋转矩阵
{
	float x, y, z;
	float Quaternion[4];
	Eigen::MatrixXf R_z(4, 4);
	Eigen::MatrixXf R_y(4, 4);
	Eigen::MatrixXf R_x(4, 4);
	Eigen::MatrixXf R_T(4, 4);
	x = pose.x;
	y = pose.y;
	z = pose.z;
	Quaternion[0] = pose.Quaternion[0];
	Quaternion[1] = pose.Quaternion[1];
	Quaternion[2] = pose.Quaternion[2];
	Quaternion[3] = pose.Quaternion[3];
	R_T=QuaternionToRot(Quaternion);       //四元数转换旋转矩阵
	cout << R_T << endl;
	R_T(0, 3) = x - R_T(0, 2)*dh1.d6;
	R_T(1, 3) = y - R_T(1, 2)*dh1.d6;
	R_T(2, 3) = z - R_T(2, 2)*dh1.d6;
	return R_T;
}






bool theta1_3(float pose_x, float pose_y, float pose_z, std::vector<float>& vec_angles, struct DH dh1, int joint1_d, int joint3_d)       //根据位置求前三个关节的角度，1-3轴角度	
{
	float Array_angles[3];
	float length1; float length2; float hight; float offset;
	float x, y, zlink, ylink, j1;
	float x1, y1;
	length1 = dh1.a2; length2 = dh1.a3; hight = dh1.d1; offset = dh1.a1;
	hight = dh1.d1; offset = dh1.a1;
	length1 = dh1.a2;
	length2 = sqrt(dh1.a3*dh1.a3 + dh1.d3*dh1.d3);
	x1 = dh1.d2 * (-pose_y / (sqrt(pose_x*pose_x + pose_y*pose_y))) + pose_x;
	y1 = dh1.d2 * (pose_x / (sqrt(pose_x*pose_x + pose_y*pose_y))) + pose_y;
	zlink = pose_z - hight;//
	j1 = atan2(y1, x1);//
	ylink = sqrt(x1*x1 + y1*y1) - offset;//
	if (joint1_d > 0)
	{
		x1 = -dh1.d2 * (-pose_y / (sqrt(pose_x*pose_x + pose_y*pose_y))) + pose_x;
		y1 = -dh1.d2 * (pose_x / (sqrt(pose_x*pose_x + pose_y*pose_y))) + pose_y;
		j1 = atan2(y1, x1) + PI;//
		ylink = sqrt(x1*x1 + y1*y1) + offset;
	}
	x = ylink;
	y = zlink;
	float theta1, theta2, theta3, theta4;//根据2-3的位置关系，确定2-3轴的一个组合
	theta1 = acos((x*x + y*y + length1*length1 - length2*length2) / (2 * length1*sqrt(x*x + y*y))) + atan(y / x);
	theta2 = acos((x*x + y*y + length2*length2 - length1*length1) / (2 * length2*sqrt(x*x + y*y))) + acos((x*x + y*y + length1*length1 - length2*length2) / (2 * length1*sqrt(x*x + y*y)));
	theta3 = -acos((x*x + y*y + length1*length1 - length2*length2) / (2 * length1*sqrt(x*x + y*y))) + atan(y / x);
	theta4 = -theta2;
	if (joint1_d > 0)
	{
		theta1 = -theta1 + PI;
		theta2 = -theta2;;
		theta3 = -theta3 + PI;
		theta4 = -theta4;
	}
	Array_angles[0] = j1;
	Array_angles[1] = 3.1415 / 2 - theta1;
	Array_angles[2] = theta2 - PI / 2;
	if (joint3_d > 0)
	{
		Array_angles[1] = PI / 2 - theta3;
		Array_angles[2] = theta4 - PI / 2;
	}
	Array_angles[2] = Array_angles[2] + atan2(dh1.d3, dh1.a3);
	vec_angles.push_back(Array_angles[0]);
	vec_angles.push_back(Array_angles[1]);
	vec_angles.push_back(Array_angles[2]);
	return true;
}


void theta4_5(float theta1, float theta2, float theta3, Eigen::Matrix4f pose_ROT, std::vector<float>& vec_angles, struct DH dh1, int joint4_d)   //求出前三个关节后，根据前三个关节的的末端位姿和整体全五轴的末端位姿的关系，求4-5轴的角度
{
	float joint_4, joint_5, joint_6;
	float Array_angles[6];
	float z, z1, z2;
	Eigen::MatrixXf Rot4_6(4, 4);
	Eigen::MatrixXf Rot0_5(4, 4);
	Eigen::MatrixXf Rot5_6(4, 4);
	Eigen::MatrixXf Rot3_4(4, 4);
	Eigen::MatrixXf Rot4_5(4, 4);
	Eigen::MatrixXf Rot0_3(4, 4);
	Eigen::MatrixXf Rot_y(4, 4);
	Rot_y << 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0, dh1.d3, 0, 0, 0, 1;
	z = theta1 + PI;
	z1 = theta2 + PI / 2;
	z2 = theta3 + PI / 2;
	Rot0_3 = DH_Matrix1(z, -dh1.a1, dh1.alpha1, dh1.d1)*DH_Matrix1(z1, dh1.a2, dh1.alpha2, dh1.d2)*DH_Matrix1(z2, dh1.a3, dh1.alpha3, dh1.d3)*Rot_y;
	Rot4_6 = Rot0_3.inverse()* pose_ROT;
	int j4_int;
	joint_4 = atan2(Rot4_6(1, 2), Rot4_6(0, 2));//四轴+180度，五轴取反
	joint_5 = atan2(Rot4_6(0, 2)*cos(joint_4) + Rot4_6(1, 2)*sin(joint_4), Rot4_6(2, 2));
	Array_angles[0] = theta1;
	Array_angles[1] = theta2;
	Array_angles[2] = theta3;
	if (joint4_d > 0)
	{
		joint_4 = joint_4 - PI;
		joint_5 = -joint_5;
	}
	z2 = joint_4;
	Rot3_4 << cos(z2), -sin(z2), 0, 0, sin(z2), cos(z2), 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
	Rot4_5 << cos(joint_5), 0, sin(joint_5), 0, 0, 1, 0, 0, -sin(joint_5), 0, cos(joint_5), 0, 0, 0, 0, 1;
	Rot0_5 = Rot0_3*Rot3_4*Rot4_5;
	Rot5_6 = Rot0_5.inverse()*pose_ROT;
	joint_6 = atan2(Rot5_6(1, 0), Rot5_6(1, 1));
	if (joint4_d > 0)
	{
		joint_6 = joint_6 + PI;
	}
	j4_int = joint_4 / PI;
	joint_4 = joint_4 - j4_int*PI * 2;
	Array_angles[3] = joint_4 + dh1.theta4;
	Array_angles[4] = joint_5 + dh1.theta5;
	Array_angles[5] = joint_6 + dh1.theta6;
	std::vector<float> vec_angle(Array_angles, Array_angles + 6);
	vec_angles = vec_angle;
}
float getIKPose(struct Position pose, std::vector<float>& vec_angles, int ME[3], struct DH dh1, int choose_type)   //输入位置和欧拉角数组【位置+欧拉角】，输出六个轴关节角度theta[0]--theta[5]
{
	int joint1_d, joint3_d, joint4_d;
	joint1_d = ME[0];
	joint3_d = ME[1];
	joint4_d = ME[2];
	Matrix<float, 4, 4> first_step;
	std::vector<float> second_step;
	std::vector<float>  third_step;
	if (choose_type < 0.5)
	{
		first_step = pose2Rot_Quaternion(pose, dh1);
	}
	else 
	{
		first_step = pose2Rot_Euler(pose, dh1, RZYXn);
	}
	theta1_3(first_step(0,3), first_step(1,3), first_step(2,3), second_step, dh1, joint1_d, joint3_d);
	theta4_5(second_step[0], (second_step[1]), (second_step[2]), first_step, third_step, dh1, joint4_d);
	vec_angles = third_step;
	int A = 0;
	if (isnan(vec_angles[1]) || isnan(vec_angles[2]) || isnan(vec_angles[3]) || isnan(vec_angles[4]) || isnan(vec_angles[5]))
	{
		A = 1;
	}

	return A;
}
bool getIKPose_UR(Position pose, joint& vec_angles, float ME[3], DH dh1, float theta[6][8])
{
	Matrix4f first_step;
	Matrix4f pose_mat;
	float psi;
	float phi;
	float joint1;
	float joint2;
	float joint3;
	float joint4;
	float joint5;
	float joint6;
	first_step = pose2Rot_Quaternion(pose, dh1);
	pose_mat = first_step;
	pose_mat(0,3) = pose.x;
	pose_mat(1,3) = pose.y;
	pose_mat(2,3) = pose.z;
	//joint 1//
	psi = atan2(first_step(1, 3), first_step(0, 3));
	phi = acos(dh1.d4 / sqrt(first_step(1, 3)*first_step(1, 3) + first_step(0, 3)*first_step(0, 3)));
	joint1 = PI / 2 + psi + phi;
	joint1 = PI / 2 + psi - phi;
	for (int i = 0; i < 4; i++)
	{
		theta[0][i] = PI / 2 + psi + phi;
		theta[0][i + 4] = PI / 2 + psi - phi;
	}
	//joint 5//
	Matrix4f T1_0;
	Matrix4f T1_6;
	float p16z;
	T1_0 = (DH_Matrix1(joint1, dh1.a1, dh1.alpha1, dh1.d1)).inverse();
	T1_6 = T1_0*pose_mat;
	p16z = T1_6(2, 3);
	for (int i = 0; i < 2; i++)
	{
		joint5 = acos((p16z - dh1.d4) / dh1.d6);
		theta[4][i * 4] = joint5;
		theta[4][i * 4 + 1] = joint5;
		joint5 = -joint5;
		theta[4][i * 4 + 2] = joint5;
		theta[4][i * 4 + 3] = joint5;
	}
	//joint6//
	Matrix4f T0_1;
	Matrix4f T6_1;
	float T61zy;
	float T61zx;
	float t5;
	T0_1 = DH_Matrix1(joint1, dh1.a1, dh1.alpha1, dh1.d1);
	T6_1 = pose_mat.inverse()*T0_1;
	T61zy = T6_1(1, 2);
	T61zx = T6_1(0, 2);
	for (int i = 0; i < 4; i++)
	{
		t5 = theta[5][2 * i];
		theta[5][2 * i + 1] = atan2(-T61zy / sin(t5), T61zx / sin(t5));
		theta[5][2 * i + 1] = atan2(-T61zy / sin(t5), T61zx / sin(t5));
	}
	//joint3//
	Matrix4f T6_5;
	Matrix4f T5_4;
	Matrix4f T1_4;
	T6_5 = (DH_Matrix1(joint6, dh1.a6, dh1.alpha6, dh1.d6)).inverse();
	T5_4 = (DH_Matrix1(joint5, dh1.a5, dh1.alpha5, dh1.d5)).inverse();
	T1_4 = T1_0*pose_mat*T6_5*T5_4;
	float T1_3x;
	float T1_3y;
	float T1_3z;
	float p13norm;
	T1_3x = T1_4(0, 3) - T1_4(0, 2)*dh1.d4;
	T1_3y = T1_4(1, 3) - T1_4(1, 2)*dh1.d4;
	T1_3z = T1_4(2, 3) - T1_4(2, 2)*dh1.d4;
	p13norm = sqrt(T1_3x*T1_3x + T1_3y*T1_3y + T1_3z*T1_3z);
	joint3 = acos((p13norm - dh1.a2*dh1.a2 - dh1.a3*dh1.a3) / (2 * dh1.a2*dh1.a3));
	for (int i = 0; i < 4; i++)
	{
		theta[2][2 * i] = joint3;
		theta[2][2 * i + 1] = -joint3;
	}
	//joint2 ,joint4//
	Matrix4f T3_2;
	Matrix4f T2_1;
	Matrix4f T3_4;
	for (int i = 0; i < 8; i++)
	{
		theta[1][i] = -atan2(T1_3y, -T1_3x) + asin(dh1.a3*sin(theta[3][i]) / p13norm);
		T3_2 = (DH_Matrix1(theta[3][i], dh1.a3, dh1.alpha3, dh1.d3)).inverse();
		T2_1 = (DH_Matrix1(theta[2][i], dh1.a2, dh1.alpha2, dh1.d2)).inverse();
		T3_4 = T3_2*T2_1*T1_4;
		theta[3][i] = atan2(T3_2(2, 1), T3_4(1, 1));
	}
	for (int i = 0; i < 8; i++)
	{
		theta[0][i] = theta[0][i] - PI;
	}
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			if (theta[i][j] <= -PI)
			{
				theta[i][j] = theta[i][j] + 2 * PI;
			}
			else if (theta[i][j] > PI)
			{
				theta[i][j] = theta[i][j] - 2 * PI;
			}
		}
	}
}