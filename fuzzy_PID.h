#ifndef FUZZY_PID_H_
#define FUZZY_PID_H_
#include<iostream>
#include<string>
using std::string;
using std::cout;
using std::cin;
using std::endl;
class FuzzyPID
{
public:
	const static int N=7;
private:
	float target;  //系统的控制目标
	float actual;  //采样获得的实际值
	float e;       //误差
	float e_pre_1; //上一次的误差
	float e_pre_2; //上上次的误差
	float de;      //误差的变化率
	float emax;    //误差基本论域上限
	float demax;   //误差辩化率基本论域的上限
	float delta_Kp_max;   //delta_kp输出的上限
	float delta_Ki_max;   //delta_ki输出上限
	float delta_Kd_max;   //delta_kd输出上限
	float Ke;      //Ke=n/emax,量化论域为[-3,-2,-1,0,1,2,3]
	float Kde;     //Kde=n/demax,量化论域为[-3,-2,-1,0,1,2,3]
	float Ku_p;    //Ku_p=Kpmax/n,量化论域为[-3,-2,-1,0,1,2,3]
	float Ku_i;    //Ku_i=Kimax/n,量化论域为[-3,-2,-1,0,1,2,3]
	float Ku_d;    //Ku_d=Kdmax/n,量化论域为[-3,-2,-1,0,1,2,3]
	int Kp_rule_matrix[N][N];//Kp模糊规则矩阵
	int Ki_rule_matrix[N][N];//Ki模糊规则矩阵
	int Kd_rule_matrix[N][N];//Kd模糊规则矩阵
	string mf_t_e;       //e的隶属度函数类型
	string mf_t_de;      //de的隶属度函数类型
	string mf_t_Kp;      //kp的隶属度函数类型
	string mf_t_Ki;      //ki的隶属度函数类型
	string mf_t_Kd;      //kd的隶属度函数类型
	float *e_mf_paras; //误差的隶属度函数的参数
	float *de_mf_paras;//误差的偏差隶属度函数的参数
	float *Kp_mf_paras; //kp的隶属度函数的参数
	float *Ki_mf_paras; //ki的隶属度函数的参数
	float *Kd_mf_paras; //kd的隶属度函数的参数
	float Kp;
	float Ki;
	float Kd;
	float A;
	float B;
	float C;
	void showMf(const string & type,float *mf_paras);      //显示隶属度函数的信息
	void setMf_sub(const string & type,float *paras,int n);//设置模糊隶属度函数的子函数
public:
	FuzzyPID(float e_max,float de_max,float kp_max,float ki_max,float kd_max,float Kp0,float Ki0,float Kd0);
	FuzzyPID(float *fuzzyLimit,float *pidInitVal);
	~FuzzyPID();
	float trimf(float x,float a,float b,float c);          //三角隶属度函数
	float gaussmf(float x,float ave,float sigma);          //正态隶属度函数
	float trapmf(float x,float a,float b,float c,float d); //梯形隶属度函数
	void setMf(const string & mf_type_e,float *e_mf,
			   const string & mf_type_de,float *de_mf,
			   const string & mf_type_Kp,float *Kp_mf,
		       const string & mf_type_Ki,float *Ki_mf,
			   const string & mf_type_Kd,float *Kd_mf);	//设置模糊隶属度函数的参数
	void setRuleMatrix(int kp_m[N][N],int ki_m[N][N],int kd_m[N][N]);  //设置模糊规则
	float realize(float t,float a);              //实现模糊控制
	void showInfo();                                      //显示该模糊控制器的信息


};

#endif