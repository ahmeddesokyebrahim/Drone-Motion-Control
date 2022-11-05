#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"
#include "Math/Angles.h"
#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;

#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();

  // Load parameters (default to 0)
  kpPosXY = config->Get(_config + ".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);

  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
  
  #else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
  #endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to
  //   individual motor thrust commands
  // INPUTS:
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS:
  // - you can access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  /*
  c= f1+f2+f3+f4
  𝑥==(𝐹1+𝐹4−𝐹2−𝐹3)
  𝜏𝑦=(𝐹1+𝐹2−𝐹3−𝐹4)𝑙
  𝜏𝑧=𝜏1+𝜏2+𝜏3+𝜏4
  
  mx = l * (f1,f2,f3,f4)  // l*(f_1 + f_4 - f_2 - f_3)
  my = l * (f1,f2,f3,f4)  // l*(f_1 + f_2 - f_3 - f_4)
  mz -> f, kappa
  mz = -t1 + t2 -t3 + t4
  
  kappa = t/f
  t = kappa *f
  t1 = kappa * f1
  t2 = kappa * f2
  t3 = kappa * f3
  t4 = kappa * f4

  mz = -kappa * (f1-f2+f3-f4)

  c= f1+f2+f3+f4  (1)
  mx = l*(f1 + f4 - f2 - f3) (2)
  my = l*(f1 + f2 - f3 - f4) (3)
  mz = -kappa * (f1-f2+f3-f4) (4)
  
  (0.25 * C)    +    (0.25/L  * mx) +  (0.25/L  * my) +  (-0.25/K   * mz) = u1
  (0.25 * C)    +    (-0.25/L * mx) +  (0.25/L  * my) +  (0.25/K   * mz)  = u2
  (0.25 * C)    +    (-0.25/L * mx) +  (-0.25/L * my) +  (-0.25/K   * mz) = u3
  (0.25 * C)    +    (0.25/L  * mx) +  (-0.25/L * my) +  (0.25/K    * mz) = u4
  */

  float K = -kappa;
  float C = collThrustCmd;
  float mx = momentCmd.x;
  float my = momentCmd.y;
  float mz = momentCmd.z;
  float l = L / sqrt(2.0f);

  cmd.desiredThrustsN[0] = (0.25 * C) + ((0.25 / l) * mx) + ((0.25 / l) * my) + ((0.25 / K) * mz);
  cmd.desiredThrustsN[1] = (0.25 * C) + ((-0.25 / l) * mx) + ((0.25 / l) * my) + ((-0.25 / K) * mz);
  cmd.desiredThrustsN[2] = (0.25 * C) + ((0.25 / l) * mx) + ((-0.25 / l) * my) + ((-0.25 / K) * mz);
  cmd.desiredThrustsN[3] = (0.25 * C) + ((-0.25 / l) * mx) + ((-0.25 / l) * my) + ((0.25 / K) * mz);

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS:
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS:
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)

  V3F momentCmd;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  /*
  p_err= p_c - p_actual 
  u_bar_p = self.k_p_p * p_err
  q_err= q_c - q_actual 
  
  u_bar_q = self.k_p_q * q_err
  
  r_err= r_c - r_actual
  
  u_bar_r = self.k_p_r * r_err
  
  For roll, pitch and yaw  
  𝑢¯p=𝜏𝑥/𝐼𝑥
  𝑢¯𝑞=𝜏𝑦/𝐼𝑦
  𝑢¯𝑟=𝜏𝑧/𝐼𝑧 
  */

  float p_err = pqrCmd.x - pqr.x;
  float u_bar_p = kpPQR.x * p_err;

  momentCmd.x = u_bar_p * Ixx;

  float q_err = pqrCmd.y - pqr.y;
  float u_bar_q = kpPQR.y * q_err;
  momentCmd.y = u_bar_q * Iyy;

  float r_err = pqrCmd.z - pqr.z;
  float u_bar_r = kpPQR.z * r_err;
  momentCmd.z = u_bar_r * Izz;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// returns a desired roll and pitch rate
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS:
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //    return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS:
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first
  /*
    b_x_a = rot_mat[0,2]
    b_x_c_dot = self.k_p_roll * (b_x_c_target - b_x_a)
    
    b_y_a = rot_mat[1,2]
    b_y_c_dot = self.k_p_pitch * (b_y_c_target - b_y_a)
    
    rot_mat1 = np.array([[rot_mat[1,0], -rot_mat[0,0]], [rot_mat[1,1], -rot_mat[0,1]]])/rot_mat[2,2]
    
    ang_vel_bf = np.matmul(rot_mat1, np.array([b_x_c_dot, b_y_c_dot]).transpose())
    
    p_c = ang_vel_bf[0]
    q_c = ang_vel_bf[1]
    
    return p_c, q_c

*/
  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  /*
        rot_mat1 = np.array([[rot_mat[1, 0], -rot_mat[0, 0]], [rot_mat[1, 1], -rot_mat[0, 1]]] ) / rot_mat[2, 2];

        b_x = rot_mat[0,2]
        b_x_err = b_x_c - b_x
        b_x_p_term = self.k_p_roll * b_x_err

        b_y = rot_mat[1,2]
        b_y_err = b_y_c - b_y  
        b_y_p_term = self.k_p_pitch * b_y_err
        
        b_x_commanded_dot = b_x_p_term
        b_y_commanded_dot = b_y_p_term
        
        rot_mat1=np.array([[rot_mat[1,0],-rot_mat[0,0]],[rot_mat[1,1],-rot_mat[0,1]]])/rot_mat[2,2]
        
        rot_rate = np.matmul(rot_mat1,np.array([b_x_commanded_dot,b_y_commanded_dot]).T)
        p_c = rot_rate[0]
        q_c = rot_rate[1]
        */

  ///f = m.a
  // a = f/m
  float a = -1 * (collThrustCmd / mass);

  float b_x_c = accelCmd.x / a;

  float b_x = R(0, 2);
  float b_x_err = b_x_c - b_x;
  float b_x_p_term = kpBank * b_x_err;

  float b_y_c = accelCmd.y / a;
  float b_y = R(1, 2);
  float b_y_err = b_y_c - b_y;
  float b_y_p_term = kpBank * b_y_err;

  // p_c = (r21*b_x_p_term -r11* b_y_p_term) / r33
  // q_c = (r22*b_x_p_term - r12 * b_y_p_term) / r33

  float p_c = (R(1, 0) * b_x_p_term - R(0, 0) * b_y_p_term) / R(2, 2);
  float q_c = (R(1, 1) * b_x_p_term - R(0, 1) * b_y_p_term) / R(2, 2);

  pqrCmd.x = p_c;
  pqrCmd.y = q_c;
  pqrCmd.z = 0.0f;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ,
                                   Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical
  //   acceleration feed-forward command
  // INPUTS:
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS:
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  //u_1_bar = self.z_k_p * (z_target - z_actual) + self.z_k_d * (z_dot_target - z_dot_actual) + z_dot_dot_target
  // c = (u_1_bar - self.g) / rot_mat[2,2]

  float u_1_bar = kpPosZ * (posZCmd - posZ) + kpVelZ * (velZCmd - velZ) + accelZCmd;
  float acceleration = (u_1_bar - CONST_GRAVITY) / R(2, 2);
  CONSTRAIN(acceleration,  maxDescentRate/ dt, maxAscentRate / dt);
  thrust = -mass * acceleration;


  // float z_err = posZCmd - posZ;
  // float p_term = kpPosZ * z_err;

  // float z_dot_err = velZCmd - velZ;
  // integratedAltitudeError += z_err * dt;

  // float d_term = kpVelZ * z_dot_err + velZ;
  // float i_term = KiPosZ * integratedAltitudeError;
  // float b_z = R(2, 2);

  // float u_1_bar = p_term + d_term + i_term + accelZCmd;

  // float acc = (u_1_bar - CONST_GRAVITY) / b_z;

  // thrust = -mass * CONSTRAIN(acc, -maxAscentRate / dt, maxAscentRate / dt);

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS:
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations.
  //     the Z component should be 0
  // HINTS:
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  //  x_dot_dot_command = self.x_k_p * (x_target - x_actual) + self.x_k_d * (x_dot_target - x_dot_actual) + x_dot_dot_target

  //   y_dot_dot_command = self.y_k_p * (y_target - y_actual) + self.y_k_d * (y_dot_target - y_dot_actual) + y_dot_dot_target
  
  V3F accelCmd = accelCmdFF;
  if (velCmd.mag() > maxSpeedXY)
  {
    velCmd = maxSpeedXY * velCmd.norm();
  }

  accelCmd.x = kpPosXY * (posCmd.x - pos.x) + kpVelXY * (velCmd.x - vel.x) + accelCmdFF.x;

  accelCmd.y = kpPosXY * (posCmd.y - pos.y) + kpVelXY * (velCmd.y - vel.y) + accelCmdFF.y;

  if (accelCmd.mag() > maxAccelXY)
  {
    accelCmd = maxAccelXY * accelCmd.norm();
  }
  /////////////////////////////// END STUDENT CODE ////////////////////////////


   return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS:
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS:
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b].
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd = 0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  // psi_error = psi_target - psi_actual
  // r_c = self.k_p_yaw * psi_error

  float yawError = yawCmd - yaw;
  yawError = AngleNormF(yawError);
  yawRateCmd = kpYaw * yawError;
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;
}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f * (maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust + thrustMargin) * 4.f, (maxMotorThrust - thrustMargin) * 4.f);

  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);

  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
