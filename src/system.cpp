#include <cmath>
#include "system.h"

/***************************************************/
/*                                                                              */
/*   Functions defined in class "Initialization"              */
/*                                                                              */
/***************************************************/
//--------------------------------------------------
// Initialization
//--------------------------------------------------
Initialization::Initialization(communicator &world){

/**********************/
  start_time = 0.0;
  end_time   = 20000;
  output_dt  = 100;
/**********************/
  num_output = int (ceil((end_time - start_time)/ output_dt));

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< System parameters are set\n";
    cout<<"<<<<< start time:             "<<start_time<<"\n";
    cout<<"<<<<< end time:               "<<end_time<<"\n";
    cout<<"<<<<< output dt:              "<<output_dt<<"\n";
    cout<<"<<<<< number of output file:  "<<num_output<<"\n";
    cout<<"***********************************************************\n";
  }
}
//--------------------------------------------------
// Get_num_timestep
//--------------------------------------------------
int Initialization::Get_num_timestep(){
  return num_output;
}
//--------------------------------------------------
// Get_start_time
//--------------------------------------------------
Real Initialization::Get_start_time(){
  return start_time;
}
//--------------------------------------------------
// Get_time_interval
//--------------------------------------------------
Real Initialization::Get_dt(){
  return output_dt;
}
