
#ifndef CONTROL_H
#define CONTROL_H

#include <math.h>

struct {
  float kp_theta = .0762 ;
  float kd_theta = .0735;
  float ki_theta = 0.0;
  float kp_phi = .0165;
  float kd_phi = .004;
  float kp_psi = .0819;
  float kd_psi = .156;
  float ki_psi = 0.0;
  float km = .274; //0.338; //on my other file
} gains;

#include "tuning_utilities.h"

// physical parameters of the system
static struct {
  float m1=0.108862;
  float ell1=0.247;
  float m2=0.4717;
  float ell2=-0.039;
  float m3=.1905;
  float g = 9.81;
  float ellT = 0.29;
  float ell3x=-.007;
  float ell3y=.018;
  float J1x = 0.000189;
  float J1y = 0.001953;
  float J1z = 0.001894;
  float J2x = 0.00231;
  float J2y = 0.003274;
  float J2z = 0.003416;
  float J3x = 0.0002222;
  float J3y = 0.0001956;
  float J3z = 0.000027;
  float d = 0.12;
  float fe = (m1*ell1+m2*ell2)*g/ellT;  
  float force_max = 0.7; //pwm limit is 0.7 I believe
} P;

// reference structure the reference signals for psi and theta
struct Reference {
  float theta = 0.0;
  float psi = 0.0;
  float phi = 0.0;  
};

// Lateral controller for hummingbird
class CtrlLatLonPID {
  public:
    float theta_d2;
    float theta_d3;
    float theta_dot_d2;
    float theta_dot_d3;
    float phi_d1;
    float phi_dot_d2;
    float phi_dot_d3;
    float psi_d1;
    float psi_dot_d2;
    float psi_dot_d3;    
    float integrator_psi;
    float error_psi_d1;
    float integrator_theta;
    float error_theta_d1;
    
    CtrlLatLonPID() {  
    }

    void init() {
      // persistent variables
      phi_d1 = 0.0;
      phi_dot_d2 = 0.0;
      phi_dot_d3 = 0.0;
      psi_d1 = 0.0;
      psi_dot_d2 = 0.0;
      psi_dot_d3 = 0.0; 
      theta_d2 = 0.0;
      theta_d3 = 0.0;
      theta_dot_d2 = 0.0;
      theta_dot_d3 = 0.0;
      integrator_psi = 0.0;
      integrator_theta = 0.0;
      error_theta_d1 = 0.0;     
      error_psi_d1 = 0.0;      
    }

    void update(float psi_ref, float theta_ref,
                SensorUtilities &sensors, 
                MotorUtilities &rotors, 
                float Ts) {

      // tune gains
      tuneGains();

      // compute phi_dot/psi_dot (with quadratic prediction)
      float phi = sensors.roll;
      float theta_d1 = sensors.pitch;
      float theta = 3*theta_d1 - 3*theta_d2 + theta_d3;
      float theta_dot_d1 = (sensors.pitch - theta_d2) / Ts;
      float theta_dot = 3*theta_dot_d1 - 3*theta_dot_d2 + theta_dot_d3;
      float phi_dot_d1 = (phi - phi_d1) / Ts;
      float phi_dot = 3*phi_dot_d1 - 3*phi_dot_d2 + phi_dot_d3;      
      float psi = sensors.yaw;
      float psi_dot_d1 = (psi - psi_d1) / Ts;
      float psi_dot = 3*psi_dot_d1 - 3*psi_dot_d2 + psi_dot_d3;  

      //compute errors
      float error_psi = psi_ref - psi;
      float error_theta = theta_ref - theta;

      // feedback linearized force
      //(P.m1*P.ell1 + P.m2*P.ell2)*P.g*np.cos(theta) / P.ellT
      float force = (P.m1*P.ell1 + P.m2*P.ell2)*P.g*cos(theta)/P.ellT;//P.fe;
      
      //pitch integrator
      if (abs(theta_dot) < 0.1) {
        integrator_theta += (error_theta + error_theta_d1) * Ts/2.0;
      }
      //pitch control
      float force_unsat = force + gains.kp_theta * error_theta + gains.ki_theta * integrator_theta - gains.kd_theta * theta_dot;

      //yaw integrator
      if (abs(psi_dot) < 0.1) {
        integrator_psi += (error_psi + error_psi_d1) * Ts/2.0;
      }
      
      //yaw control
      float phi_ref = gains.kp_psi * error_psi + integrator_psi * gains.ki_psi - gains.kd_psi * psi_dot;
      float error_phi = phi_ref - phi;
      //T_phi = self.kp_phi * error_phi - self.kd_phi * self.phi_dot
      float torque =gains.kp_phi * error_phi - gains.kd_phi * phi_dot;

      // convert force and torque to pwm and send to motors
      float left_pwm = (force_unsat+torque/P.d)/(2.0*gains.km);
      float right_pwm = (force_unsat-torque/P.d)/(2.0*gains.km);
      left_pwm = saturate(left_pwm, 0.0, 0.7);
      right_pwm = saturate(right_pwm, 0.0, 0.7);
      rotors.update(left_pwm, right_pwm); 

      // update all delayed variables
      phi_d1 = phi;
      phi_dot_d3 = phi_dot_d2;
      phi_dot_d2 = phi_dot_d1;
      psi_d1 = psi;
      psi_dot_d3 = psi_dot_d2;
      psi_dot_d2 = psi_dot_d1;      
      error_psi_d1 = error_psi;
      theta_d3 = theta_d2;
      theta_d2 = theta_d1;
      theta_dot_d3 = theta_dot_d2;
      theta_dot_d2 = theta_dot_d1;
      error_theta_d1 = error_theta;
      // print commanded values
      Serial.print("Psi_ref here:");
      Serial.print(psi_ref*180/PI);
      Serial.print(",");
      Serial.print("Psi:");
      Serial.print(psi*180/PI);
      Serial.print(",");      
      Serial.print("Phi_ref:");
      Serial.println(phi_ref*180/PI);
      Serial.print(",");
      Serial.print("Phi:");
      Serial.print(phi*180/PI); 
      Serial.print("Theta_ref:");
      Serial.print(theta_ref*180/PI);
      Serial.print(",");
      Serial.print("Theta:");
      Serial.println(theta*180/PI);     
    }

    float saturate(float value, float min_value, float max_value) {
      // Implements the saturation function
      return min(max(min_value, value), max_value);
    }  
};

#endif 
