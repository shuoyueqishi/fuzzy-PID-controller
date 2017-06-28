#include<iostream>
#include"fuzzy_PID.h"

#define NB -3
#define NM -2
#define NS -1
#define ZO 0
#define PS 1
#define PM 2
#define PB 3

int main()
{
	float target=600;
	float actual=0;
	float u=0;
	int deltaKpMatrix[7][7]={{PB,PB,PM,PM,PS,ZO,ZO},
	                         {PB,PB,PM,PS,PS,ZO,NS},
						     {PM,PM,PM,PS,ZO,NS,NS},
	                         {PM,PM,PS,ZO,NS,NM,NM},
	                         {PS,PS,ZO,NS,NS,NM,NM},
	                         {PS,ZO,NS,NM,NM,NM,NB},
	                         {ZO,ZO,NM,NM,NM,NB,NB}};
	int deltaKiMatrix[7][7]={{NB,NB,NM,NM,NS,ZO,ZO},
	                         {NB,NB,NM,NS,NS,ZO,ZO},
						     {NB,NM,NS,NS,ZO,PS,PS},
	                         {NM,NM,NS,ZO,PS,PM,PM},
	                         {NM,NS,ZO,PS,PS,PM,PB},
	                         {ZO,ZO,PS,PS,PM,PB,PB},
	                         {ZO,ZO,PS,PM,PM,PB,PB}};
	int deltaKdMatrix[7][7]={{PS,NS,NB,NB,NB,NM,PS},
	                         {PS,NS,NB,NM,NM,NS,ZO},
						     {ZO,NS,NM,NM,NS,NS,ZO},
	                         {ZO,NS,NS,NS,NS,NS,ZO},
	                         {ZO,ZO,ZO,ZO,ZO,ZO,ZO},
	                         {PB,NS,PS,PS,PS,PS,PB},
	                         {PB,PM,PM,PM,PS,PS,PB}};
	float e_mf_paras[]={-3,-3,-2,-3,-2,-1,-2,-1,0,-1,0,1,0,1,2,1,2,3,2,3,3};
	float de_mf_paras[]={-3,-3,-2,-3,-2,-1,-2,-1,0,-1,0,1,0,1,2,1,2,3,2,3,3};
	float Kp_mf_paras[]={-3,-3,-2,-3,-2,-1,-2,-1,0,-1,0,1,0,1,2,1,2,3,2,3,3};
	float Ki_mf_paras[]={-3,-3,-2,-3,-2,-1,-2,-1,0,-1,0,1,0,1,2,1,2,3,2,3,3};
	float Kd_mf_paras[]={-3,-3,-2,-3,-2,-1,-2,-1,0,-1,0,1,0,1,2,1,2,3,2,3,3};
    FuzzyPID fuzzypid(1500,650,0.3,0.4,0.2,0.02,0.65,0.005);
	fuzzypid.setMf("trimf",e_mf_paras,"trimf",de_mf_paras,"trimf",Kp_mf_paras,"trimf",Ki_mf_paras,"trimf",Kd_mf_paras);
	fuzzypid.setRuleMatrix(deltaKpMatrix,deltaKiMatrix,deltaKdMatrix);
	cout<<"num target    actual"<<endl;
	/*fuzzy.showInfo();*/
	for(int i=0;i<200;i++)
	{
		u=fuzzypid.realize(target,actual);
		actual+=u;
		cout<<i<<"   "<<target<<"    "<<actual<<endl;
	}
	fuzzypid.showInfo();
	system("pause");
	return 0;
}