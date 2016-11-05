#include "colors.inc"


#declare lookFrom = <0,0,1>;
#declare lookAt = 0;

//background { CHSL2RGB( <220, 1, .15> ) }
 
light_source { lookFrom White } 
global_settings 
{ 
	assumed_gamma 1.5
	max_trace_level 10
}
          
camera 
{  
	//spherical
	//fisheye
	//ultra_wide_angle                
	//angle 180
	//angle 270                  	
	//angle 360 

	location lookFrom
	sky<0,0,1>                   // z up
	right <-4/3,0,0>             // right-handed coordinate system
	look_at lookAt  
}  

#declare fin = finish
{
   specular .5
}

#declare tex = texture
{
    pigment { White }
    //pigment { Black }
    finish { fin }
}
  


