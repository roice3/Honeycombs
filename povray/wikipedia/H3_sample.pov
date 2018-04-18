#include "colors.inc"


#declare lookFrom = 0;
#declare lookAt = -z;   // Use this for loop graphs
#declare lookAt = -x;   // Use this when defining the honeycomb with 6 angles


//background { CHSL2RGB( <220, 1, .15> ) } 
background { White } 

 
light_source { lookFrom White } 
global_settings 
{ 
	assumed_gamma 1.5
	max_trace_level 10
}

camera
{   
/*  // Spherical images can be made by uncommenting this.
	spherical	
    angle 360  // horizontal
          180  // vertical(optional)
*/      
    location lookFrom
    look_at  lookAt
    sky<0,0,-1> 
    right<-4/3,0,0>		// Change to chosen aspect ratio.
    //right<-2,0,0>     // For spherical images.
} 

#declare fin = finish
{
   specular .5
}

#declare tex = texture
{
    //pigment { White }
    pigment { Gray40 }
    finish { fin }
}

// If you generate a different honeycomb, change this include!  
#include "2-4-3-2-3-3_0001.pov"

