// Author: Juergen Gall 


#include "NRBM.h"

namespace NRBM {
  
  // RBM2Twist
  void RBM2Twist(CVector<float> &T, CMatrix<float>& RBM) {
    CVector<double> omega;
    double theta;
    CVector<double> v;
    CMatrix<double> Id(3,3);
    CMatrix<double> R(3,3);
    Id.fill(0.0);
    for (int i = 0; i <3; i++)
      Id(i,i)=1.0;
    for (int i = 0; i <3; i++)
      for (int j = 0; j <3; j++)
        R(i,j)=double(RBM(i,j));
    // if (R==Id) 
    if((fabs(R(0,0)-1.0)<0.00001)&&(fabs(R(1,1)-1.0)<0.00001)&&(fabs(R(2,2)-1.0)<0.00001)       )
      {
	T.setSize(6);
	T(0)=double(RBM(3,0)); T(1)=double(RBM(3,1)); T(2)=double(RBM(3,2));
	T(3)=0.0; T(4)=0.0; T(5)=0.0;
      }
    else {
      float diag = (R(0,0)+R(1,1)+R(2,2)-1.0)*0.5;
      if (diag < -1.0) diag = -1.0;
      else if (diag > 1.0) diag = 1.0;
      theta = acos(diag);
      if (sin(theta)==0) theta += 0.001;
      omega.setSize(3);
      omega(0)=(R(1,2)-R(2,1));
      omega(1)=(R(2,0)-R(0,2));
      omega(2)=(R(0,1)-R(1,0));
      omega*=(1.0/(2.0*sin(theta)));
      CMatrix<double> omegaHat(3,3);
      omegaHat.data()[0] = 0.0;            omegaHat.data()[1] = -omega(2); omegaHat.data()[2] = omega(1);
      omegaHat.data()[3] = omega(2);  omegaHat.data()[4] = 0.0;            omegaHat.data()[5] = -omega(0);
      omegaHat.data()[6] = -omega(1); omegaHat.data()[7] = omega(0);  omegaHat.data()[8] = 0.0;
      CMatrix<double> omegaT(3,3);
      for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++)
          omegaT(i,j) = omega(i)*omega(j);
      R = (omegaHat*(double)sin(theta))+((omegaHat*omegaHat)*(double)(1.0-cos(theta)));
      R(0,0) += 1.0; R(1,1) += 1.0; R(2,2) += 1.0;
      CMatrix<double> A(3,3);
      A.fill(0.0);
      A(0,0)=1.0; A(1,1)=1.0; A(2,2)=1.0;
      A -= R;  A*=omegaHat;  A+=omegaT*theta;
      CVector<double> p(3);
      p(0)=RBM(3,0);
      p(1)=RBM(3,1);
      p(2)=RBM(3,2);
      A.inv();
      v=A*p;
      T.setSize(6);
      T(0)=float(v(0)*theta);
      T(1)=float(v(1)*theta);
      T(2)=float(v(2)*theta);
      T(3)=float(theta*omega(0));
      T(4)=float(theta*omega(1));
      T(5)=float(theta*omega(2));
    }
  }

  // Twist2RBM
  void Twist2RBM(CVector<float>& aTwist, CMatrix<float>& A) {
    // Determine main body motion
    CVector<float> moment(3);
    moment(0) = aTwist(0); moment(1) = aTwist(1); moment(2) = aTwist(2);
    CVector<float> omega(3);
    omega(0) = aTwist(3); omega(1) = aTwist(4); omega(2) = aTwist(5);
    A.setSize(4,4);
    if ((omega(0)==0)&&(omega(1)==0)&&(omega(2)==0)) {
      A.data()[0] = 1.0; A.data()[1] = 0.0; A.data()[2] = 0.0; A.data()[3] = moment(0);
      A.data()[4] = 0.0; A.data()[5] = 1.0; A.data()[6] = 0.0; A.data()[7] = moment(1);
      A.data()[8] = 0.0; A.data()[9] = 0.0; A.data()[10]= 1.0; A.data()[11]= moment(2);
      A.data()[12]= 0.0; A.data()[13]= 0.0; A.data()[14]= 0.0; A.data()[15]= 1.0;
    }
    else {
      float theta = sqrt(omega(0)*omega(0)+omega(1)*omega(1)+omega(2)*omega(2));
      float invTheta = 0.0;
      if (theta != 0) invTheta = 1.0/theta;
      omega *= invTheta; moment *= invTheta;
      CMatrix<float> omegaHat(3,3);
      omegaHat.data()[0] = 0.0;       omegaHat.data()[1] = -omega(2); omegaHat.data()[2] = omega(1);
      omegaHat.data()[3] = omega(2);  omegaHat.data()[4] = 0.0;       omegaHat.data()[5] = -omega(0);
      omegaHat.data()[6] = -omega(1); omegaHat.data()[7] = omega(0);  omegaHat.data()[8] = 0.0;
      CMatrix<float> omegaT(3,3);
      for (int j = 0; j < 3; j++)
	for (int i = 0; i < 3; i++)
	  omegaT(i,j) = omega(i)*omega(j);
      CMatrix<float> R(3,3);
      R = (omegaHat*(float)sin(theta))+((omegaHat*omegaHat)*(float)(1.0-cos(theta)));
      R(0,0) += 1.0; R(1,1) += 1.0; R(2,2) += 1.0;
      CMatrix<float> T(R); T *= -1.0;
      T(0,0) += 1.0; T(1,1) += 1.0; T(2,2) += 1.0;
      CVector<float> t(3);
      t = T*(omega/moment)+((omegaT*moment)*theta);
      A.data()[0] = R.data()[0]; A.data()[1] = R.data()[1]; A.data()[2] = R.data()[2]; A.data()[3] = t(0);
      A.data()[4] = R.data()[3]; A.data()[5] = R.data()[4]; A.data()[6] = R.data()[5]; A.data()[7] = t(1);
      A.data()[8] = R.data()[6]; A.data()[9] = R.data()[7]; A.data()[10]= R.data()[8]; A.data()[11]= t(2);
      A.data()[12]= 0.0;         A.data()[13]= 0.0;         A.data()[14]= 0.0;         A.data()[15]= 1.0;
    }
  }

  // Tranbsform Twist with Rigid body motion ...
  void AdjointTwist(CVector<float> &T, CMatrix<float>& RBM) {
    CMatrix<float> XiHat(4,4);
    CMatrix<float> XiNew(4,4);
    
    CMatrix<float> RInv;
    RInv=RBM;
    RInv.inv();
    
    for(int i=0;i<4;i++)
      for(int j=0;j<4;j++)
	XiHat(i,j)=0.0;
    
    XiHat(3,0)=T(0);
    XiHat(3,1)=T(1);
    XiHat(3,2)=T(2);
    XiHat(0,1)=T(5);XiHat(1,0)=-T(5);
    XiHat(0,2)=-T(4);XiHat(2,0)=T(4);
    XiHat(1,2)=T(3);XiHat(2,1)=-T(3);
    
    //XiNew=RInv*XiHat*RBM;
    XiNew=RBM*XiHat*RInv;

    T(0)=XiNew(3,0);
    T(1)=XiNew(3,1);
    T(2)=XiNew(3,2);
    T(3)=XiNew(1,2);
    T(4)=XiNew(2,0);
    T(5)=XiNew(0,1);  

    
    
  }

  double Pitch(CVector<float> &TV)
  {
    double h;
    double n;
    
    h=TV(0)*TV(3)+TV(1)*TV(4)+TV(2)*TV(5);
    
    n=sqrt(TV(3)*TV(3)+TV(4)*TV(4)+TV(5)*TV(5));
    if(n!=0)
      return h/n;
    else
      return 0;
  }
  
  // Converts rigid body motion into rotation vector and translation vector
  void RBM2RVT(const CMatrix<float>* RBM, CVector<double>& rvt) {
    RM2RV(RBM, rvt);
    //std::cout << *RBM << std::endl;
    rvt(3) = (*RBM)(3,0);
    rvt(4) = (*RBM)(3,1);
    rvt(5) = (*RBM)(3,2);
    //std::cout << *RBM << std::endl;
  }
  
  void RVT2RBM(const CVector<double>* rvt, CMatrix<float>& RBM) {
    RV2RM(rvt, RBM);
//     RBM(3,0) = (float)(*rvt)(3);
//     RBM(3,1) = (float)(*rvt)(4);
//     RBM(3,2) = (float)(*rvt)(5);
    RBM.data()[3] = (float)(*rvt)(3);
    RBM.data()[7] = (float)(*rvt)(4);
    RBM.data()[11] = (float)(*rvt)(5);
    RBM.data()[12] = 0.0f;
    RBM.data()[13] = 0.0f;
    RBM.data()[14] = 0.0f;
    RBM.data()[15] = 1.0f;
  }

  // Converts rotation matrix into rotation vector representation
  void RM2RV(const CMatrix<float>* R, CVector<double>& rv) { 
    if((fabs((*R)(0,0)-1.0)<0.00001)&&(fabs((*R)(1,1)-1.0)<0.00001)&&(fabs((*R)(2,2)-1.0)<0.00001)) {
      rv(0)=0.0; rv(1)=0.0; rv(2)=0.0;
    } else {
      double diag = ((double)(*R)(0,0)+(double)(*R)(1,1)+(double)(*R)(2,2)-1.0)*0.5;
      if (diag < -1.0) diag = -1.0;
      else if (diag > 1.0) diag = 1.0;
      
      double theta = acos(diag);
      if (sin(theta)==0) theta += 0.001;
      
      theta /= (2.0*sin(theta));
      rv(0)=((double)(*R)(1,2)-(double)(*R)(2,1))*theta;
      rv(1)=((double)(*R)(2,0)-(double)(*R)(0,2))*theta;
      rv(2)=((double)(*R)(0,1)-(double)(*R)(1,0))*theta;
      //rv*= theta/(2.0*sin(theta));
    }
  }
  
  void RV2RM(const CVector<double>* rv, CMatrix<float>& R) {
    if (((*rv)(0)==0)&&((*rv)(1)==0)&&((*rv)(2)==0)) {    
      R(0,0) = 1.0; R(1,0) = 0.0; R(2,0) = 0.0;
      R(0,1) = 0.0; R(1,1) = 1.0; R(2,1) = 0.0;
      R(0,2) = 0.0; R(1,2) = 0.0; R(2,2) = 1.0;
    } else {
      double theta = sqrt((*rv)(0)*(*rv)(0)+(*rv)(1)*(*rv)(1)+(*rv)(2)*(*rv)(2));
      double invTheta = 0.0;
      if (theta != 0) invTheta = 1.0/theta;
      CVector<double> omega = (*rv)*invTheta;
      CMatrix<double> omegaHat(3,3);
      omegaHat.data()[0] = 0.0;       omegaHat.data()[1] = -omega(2); omegaHat.data()[2] = omega(1);
      omegaHat.data()[3] = omega(2);  omegaHat.data()[4] = 0.0;       omegaHat.data()[5] = -omega(0);
      omegaHat.data()[6] = -omega(1); omegaHat.data()[7] = omega(0);  omegaHat.data()[8] = 0.0;
      CMatrix<double> omegaT(3,3);
      for (int j = 0; j < 3; j++)
	for (int i = 0; i < 3; i++)
	  omegaT(i,j) = omega(i)*omega(j);
      omegaHat = (omegaHat*(double)sin(theta))+((omegaHat*omegaHat)*(double)(1.0-cos(theta)));
      R(0,0) = omegaHat(0,0) + 1.0; R(1,0) = omegaHat(1,0); R(2,0) = omegaHat(2,0);
      R(0,1) = omegaHat(0,1); R(1,1) = omegaHat(1,1) + 1.0; R(2,1) = omegaHat(2,1);
      R(0,2) = omegaHat(0,2); R(1,2) = omegaHat(1,2); R(2,2) = omegaHat(2,2) + 1.0;
    }  
  }

  // mean of rotations
  void meanRotation(const std::vector<CMatrix<float> >* vR, const std::vector<double>* vW, CMatrix<float>& meanR) {
    CMatrix<float> tmpRM(3,3);
    CVector<double> tmpRV(3);
    CVector<double> wRV(3);
    for(int t=0; t<5; ++t) {
      wRV.fill(0);
      for(unsigned int i = 0; i<vR->size(); ++i) {
	tmpRM = trans(meanR) * (*vR)[i];
	RM2RV(&tmpRM, tmpRV);
	wRV += (*vW)[i]*tmpRV;
      }
      RV2RM(&wRV, tmpRM);
      meanR = meanR * tmpRM;
    }
  }

  void invRBM(CMatrix<float>& RBM) {
    CMatrix<float> tmp(RBM); 
    
    //RBM.data()[ 0] = tmp.data()[ 0];
    RBM.data()[ 1] = tmp.data()[ 4];
    RBM.data()[ 2] = tmp.data()[ 8];
    RBM.data()[ 3] = - tmp.data()[0] * tmp.data()[3] - tmp.data()[4] * tmp.data()[7] - tmp.data()[8] * tmp.data()[11];
    
    RBM.data()[ 4] = tmp.data()[ 1];
    //RBM.data()[ 5] = tmp.data()[ 5];
    RBM.data()[ 6] = tmp.data()[ 9];
    RBM.data()[ 7] = - tmp.data()[1] * tmp.data()[3] - tmp.data()[5] * tmp.data()[7] - tmp.data()[9] * tmp.data()[11];
    
    RBM.data()[ 8] = tmp.data()[ 2];
    RBM.data()[ 9] = tmp.data()[ 6];
    //RBM.data()[10] = tmp.data()[10];
    RBM.data()[11] = - tmp.data()[2] * tmp.data()[3] - tmp.data()[6] * tmp.data()[7] - tmp.data()[10] * tmp.data()[11];
    
    RBM.data()[12] = 0.0f;
    RBM.data()[13] = 0.0f;
    RBM.data()[14] = 0.0f;
    RBM.data()[15] = 1.0f;
    
  }

  void invRM(CMatrix<float>& R) {
    float tmp; 
    
    tmp = R.data()[ 1];
    R.data()[ 1] = R.data()[ 3];
    R.data()[ 3] = tmp;

    tmp = R.data()[ 2];
    R.data()[ 2] = R.data()[ 6];
    R.data()[ 6] = tmp;

    tmp = R.data()[ 5];
    R.data()[ 5] = R.data()[ 7];
    R.data()[ 7] = tmp;

  }

  //
  // turns this quaternion into the dual part of the passed quaternion
  //
  void QT2DQ(const CVector<float>& q0, const CVector<float>& t, CVector<float>& dq)
  {
    dq[0] = -0.5*(t[0]*q0[1] + t[1]*q0[2] + t[2]*q0[3]);
    dq[1] = 0.5*( t[0]*q0[0] + t[1]*q0[3] - t[2]*q0[2]);
    dq[2] = 0.5*(-t[0]*q0[3] + t[1]*q0[0] + t[2]*q0[1]);
    dq[3] = 0.5*( t[0]*q0[2] - t[1]*q0[1] + t[2]*q0[0]);
  }

  void QDQ2RBM(const CVector<float>& q0, const CVector<float>& dq, CMatrix<float>& RBM) {
    Quat2RM(q0,RBM);
    RBM.data()[3] =  2*(-dq[0]*q0[1] + dq[1]*q0[0] - dq[2]*q0[3] + dq[3]*q0[2]);
    RBM.data()[7] =  2*(-dq[0]*q0[2] + dq[1]*q0[3] + dq[2]*q0[0] - dq[3]*q0[1]);
    RBM.data()[11] = 2*(-dq[0]*q0[3] - dq[1]*q0[2] + dq[2]*q0[1] + dq[3]*q0[0]);
    RBM.data()[12] = 0.0f;
    RBM.data()[13] = 0.0f;
    RBM.data()[14] = 0.0f;
    RBM.data()[15] = 1.0f;
  }

  void RBM2QDQ(const CMatrix<float>& RBM, CVector<float>& q0, CVector<float>& dq) {
    CVector<float> t(3);
    t[0] = RBM.data()[3];
    t[1] = RBM.data()[7];
    t[2] = RBM.data()[11];
    RM2Quat(RBM,q0);
    QT2DQ(q0,t,dq);
  }


  void RM2Quat(const CMatrix<float>& R, CVector<float>& q) {
    q[0] = sqrt( std::max( 0.0f, 1 + R(0,0) + R(1,1) + R(2,2) ) ) / 2.0; 
    q[1] = sqrt( std::max( 0.0f, 1 + R(0,0) - R(1,1) - R(2,2) ) ) / 2.0; 
    q[2] = sqrt( std::max( 0.0f, 1 - R(0,0) + R(1,1) - R(2,2) ) ) / 2.0; 
    q[3] = sqrt( std::max( 0.0f, 1 - R(0,0) - R(1,1) + R(2,2) ) ) / 2.0; 
    // liuyebin
	if(R(1,2) - R(2,1)<0)
		q[1] = -q[1];
	if(R(2,0) - R(0,2)<0)
		q[2] = -q[2];
	if(R(0,1) - R(1,0)<0)
		q[3] = -q[3];
	//q[1] = copysign( q[1], R(1,2) - R(2,1) ); 
 //   q[2] = copysign( q[2], R(2,0) - R(0,2) ); 
 //   q[3] = copysign( q[3], R(0,1) - R(1,0) ); 
  }
  
  void Quat2RM(const CVector<float>& q, CMatrix<float>& R) {
    float xx = 2.0*q[1]*q[1];
    float yy = 2.0*q[2]*q[2];
    float zz = 2.0*q[3]*q[3];
    float wx = 2.0*q[0]*q[1];
    float wy = 2.0*q[0]*q[2];
    float wz = 2.0*q[0]*q[3];
    float xy = 2.0*q[1]*q[2];
    float xz = 2.0*q[1]*q[3];
    float yz = 2.0*q[2]*q[3];
    R(0,0) = 1-yy-zz; R(1,0) = xy-wz; R(2,0) = xz+wy;
    R(0,1) = xy+wz; R(1,1) = 1-xx-zz; R(2,1) = yz-wx; 
    R(0,2) = xz-wy; R(1,2) = yz+wx; R(2,2) = 1-xx-yy; 
  }

  void toOpenGL(const CMatrix<float>& RBM, float* mat) {
    mat[0] = RBM.data()[0];
    mat[1] = RBM.data()[4];
    mat[2] = RBM.data()[8];
    mat[3] = 0;
    mat[4] = RBM.data()[1];
    mat[5] = RBM.data()[5];
    mat[6] = RBM.data()[9];
    mat[7] = 0;
    mat[8] = RBM.data()[2];
    mat[9] = RBM.data()[6];
    mat[10] = RBM.data()[10];
    mat[11] = 0;
    mat[12] = RBM.data()[3];
    mat[13] = RBM.data()[7];
    mat[14] = RBM.data()[11];
    mat[15] = 1;
  }

}



