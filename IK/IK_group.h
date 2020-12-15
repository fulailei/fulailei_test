#ifndef IK_group_H
#define IK_group_H
#endif

float *pose_ouler(struct Position pose, float *const buf, struct DH dh1, enum RotationStatus choose_R);
float *theta4_5(float theta1, float theta2, float theta3, float *A, float *const buf, struct DH dh1, int joint4_d);
float *theta1_3(float x1, float y1, float z1, float *const buf, struct DH dh1, int joint1_d, int joint3_d );
float inv_link6(struct Position pose, float *const theta, int ME[3], struct DH dh1, enum RotationStatus choose_R);   //输入矩阵【位置+欧拉角】，指针指向六个关节的矩阵，ME【3】，表示轴1，3，4的正负关系，默认状态为【0，0，0】，当取1时，比如【1，0，0】表示轴一取负，其他轴发生相应改变，【0，1，0】表示轴三取负，其他轴发生相应改变，同理一共八种情况。
float mat2(struct DH dh3, struct joint j3);

