#ifndef IK_group_H
#define IK_group_H
#endif

float *pose_ouler(struct Position pose, float *const buf, struct DH dh1, enum RotationStatus choose_R);
float *theta4_5(float theta1, float theta2, float theta3, float *A, float *const buf, struct DH dh1, int joint4_d);
float *theta1_3(float x1, float y1, float z1, float *const buf, struct DH dh1, int joint1_d, int joint3_d );
float inv_link6(struct Position pose, float *const theta, int ME[3], struct DH dh1, enum RotationStatus choose_R);   //�������λ��+ŷ���ǡ���ָ��ָ�������ؽڵľ���ME��3������ʾ��1��3��4��������ϵ��Ĭ��״̬Ϊ��0��0��0������ȡ1ʱ�����硾1��0��0����ʾ��һȡ���������ᷢ����Ӧ�ı䣬��0��1��0����ʾ����ȡ���������ᷢ����Ӧ�ı䣬ͬ��һ�����������
float mat2(struct DH dh3, struct joint j3);

