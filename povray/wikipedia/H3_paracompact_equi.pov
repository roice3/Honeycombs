#include "colors.inc"


#declare lookFrom = 0;
#declare lookAt = <0,0,-1>;

//background { CHSL2RGB( <220, 1, .15> ) }
 
light_source { lookFrom White } 
global_settings 
{ 
	assumed_gamma 1.5
	max_trace_level 10
}

camera
{ 
	spherical
    angle 360  // horizontal
          180  // vertical(optional)
    location lookFrom
    look_at  lookAt
    sky<0,0,-1> 
    right<-2,0,0>
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
  


