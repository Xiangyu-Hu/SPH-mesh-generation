#include "glbfunc.h"
#include "visualize.h"

/*********************************************************/
/*                                                       */
/*         Functions defined in class "Visualize"        */
/*                                                       */
/*********************************************************/

//-------------------------------------------------------
// set the scene
//-------------------------------------------------------
void Visualize::Set_scene(char *filename, int flag)
{
  char ch;
  char    filename1[256];
  if (flag == 1)
    sprintf(filename1,"./scene/template.pov");
  else if (flag == 2)
    sprintf(filename1,"./scene/template2.pov");

  ifstream in(filename1); // read template
  ofstream out(filename, ios::trunc);  // open outfile

  while (in.get(ch)){
    out<<ch;
  }
  in.close();
//  out.close();
}
//-------------------------------------------------------
// draw boundary
//-------------------------------------------------------
void Visualize::Draw_boundary
(char *filename, Real x_min, Real x_max, Real y_min, Real y_max, Real z_min, Real z_max)
{
  ofstream out(filename, ios::app);
  
  // Draw the outline of the domain
  out<<"union{\n";
    out<<"\tcylinder{<"<<x_min<<","<<y_min<<","<<z_min<<">,<"<<x_max<<","<<y_min<<","<<z_min<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_min<<","<<y_max<<","<<z_min<<">,<"<<x_max<<","<<y_max<<","<<z_min<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_min<<","<<y_max<<","<<z_max<<">,<"<<x_max<<","<<y_max<<","<<z_max<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_min<<","<<y_min<<","<<z_max<<">,<"<<x_max<<","<<y_min<<","<<z_max<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_min<<","<<y_min<<","<<z_min<<">,<"<<x_min<<","<<y_max<<","<<z_min<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_max<<","<<y_min<<","<<z_min<<">,<"<<x_max<<","<<y_max<<","<<z_min<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_max<<","<<y_min<<","<<z_max<<">,<"<<x_max<<","<<y_max<<","<<z_max<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_min<<","<<y_min<<","<<z_max<<">,<"<<x_min<<","<<y_max<<","<<z_max<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_min<<","<<y_min<<","<<z_min<<">,<"<<x_min<<","<<y_min<<","<<z_max<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_max<<","<<y_min<<","<<z_min<<">,<"<<x_max<<","<<y_min<<","<<z_max<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_max<<","<<y_max<<","<<z_min<<">,<"<<x_max<<","<<y_max<<","<<z_max<<">,1.5*r}\n";
    out<<"\tcylinder{<"<<x_min<<","<<y_max<<","<<z_min<<">,<"<<x_min<<","<<y_max<<","<<z_max<<">,1.5*r}\n";
    out<<"\tsphere{<"<<x_min<<","<<y_min<<","<<z_min<<">,1.5*r}\n";
    out<<"\tsphere{<"<<x_max<<","<<y_min<<","<<z_min<<">,1.5*r}\n";
    out<<"\tsphere{<"<<x_min<<","<<y_max<<","<<z_min<<">,1.5*r}\n";
    out<<"\tsphere{<"<<x_max<<","<<y_max<<","<<z_min<<">,1.5*r}\n";
    out<<"\tsphere{<"<<x_min<<","<<y_min<<","<<z_max<<">,1.5*r}\n";
    out<<"\tsphere{<"<<x_max<<","<<y_min<<","<<z_max<<">,1.5*r}\n";
    out<<"\tsphere{<"<<x_min<<","<<y_max<<","<<z_max<<">,1.5*r}\n";
    out<<"\tsphere{<"<<x_max<<","<<y_max<<","<<z_max<<">,1.5*r}\n";
  out<<"\ttexture{\n\t  pigment{\n\t\t    color rgbf <0.9255,0.9412,0.9451>}}}\n";

//  out.close();
}
//-------------------------------------------------------
// draw particle
//-------------------------------------------------------
void Visualize::Draw_particle
(char *filename, my_real pos, int id, Real val, Real min, Real max, Real tr)
{
  ofstream out(filename, ios::app);

  double r = 0.;
  double g = 0.;
  double b = 0.;
  get_rgb_value (r,g,b,Real(id),min,max);

  out<<"// particle id "<<id<<"\n";
  out<<"sphere{<"<<pos.i<<","<<pos.j<<","<<pos.k<<">,"<<val<<" texture{ pigment{ color rgbf <"<<r<<","<<g<<","<<b<<","<<tr<<">}}"<<"}\n\n";

//  out.close();
}
//-------------------------------------------------------
// draw polygon
//-------------------------------------------------------
void Visualize::Draw_polygon
(char *filename,vector<int> &f_vert,vector<double> &v,
 int j, Real val, Real min, Real max, int flag)
{
  // flag 1: output all the info. of the polygon
  // flag 2: output the outline of the polygon

  char s[10][256];
  int k,l,n=f_vert[j];
  double r = 0.;
  double g = 0.;
  double b = 0.;
  double tr  = 0.5; // translucent factor
  get_rgb_value (r,g,b,val,min,max);

  ofstream out(filename, ios::app);

  for(k=0;k<n;k++) {
    l=3*f_vert[j+k+1];
    sprintf(s[k],"<%g,%g,%g>",v[l],v[l+1],v[l+2]);
  }

  if (flag == 1){
    // Draw the interior of the polygon
    out<<"union{\n";
    for(k=2;k<n;k++) out<<"\ttriangle{"<<s[0]<<","<<s[k-1]<<","<<s[k]<<"}\n";
    out<<"\ttexture{\n\t  pigment{\n\t\t    color rgbf <"<<r<<","<<g<<","<<b<<","<<tr<<">}}}\n";
  }
  // Draw the outline of the polygon
  out<<"union{\n";
  for(k=0;k<n;k++) {
    l=(k+1)%n;
    out<<"\tcylinder{"<<s[k]<<","<<s[l]<<",1.5*r}\n\tsphere{"<<s[l]<<",1.5*r}\n";
  }
  out<<"\ttexture{\n\t  pigment{\n\t\t    color rgbf <0.9255,0.9412,0.9451>}}}\n";

 // out.close();

}
