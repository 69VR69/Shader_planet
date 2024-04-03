#iChannel0 "file://bufferA.glsl"

#define numboids 20.		// Number of boids (must be integer value represented as float)
#define speed 10.0f			// Boid speed
#define a1 0.001f			// Collision factor
#define a2 0.01f			// Cohesion factor	
#define a3 0.05f			// Alignment factor

vec4 getBoid(float id) { return texture(iChannel0, vec2(id + .5f,.5f)/iResolution.xy);}

vec2 hash(float n) { return fract(sin(vec2(n,n*7.))*43758.5f); }
                        
// returns true if fragment coordinate maps to boid data
//bool isBoid(vec2 uv) { return (floor(uv.x) == uv.x) && (uv.x < numboids) && (uv.y == 0.0f);}	

// Very Strange implementation - We use the buffer as a texture that stores and mantains the boid information at each time step.
// Then the main image program queries the boid data from this buffer.
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    fragColor = vec4(0.0f,0.0f,0.0f,0.0f);			// default data
    //float id = floor(fragCoord.x)/iResolution.x;
    
    if (fragCoord.y > 0.5 || fragCoord.x > numboids) discard;
    
    float id = floor(fragCoord.x);
    
    // initialize random boid positions and velocities
    if(iFrame < 5)
    {	vec2 init_p = hash(id)*iResolution.xy;
        vec2 init_v = vec2(sin(id),cos(id));
     	init_v /= speed * length(init_v);
        fragColor = vec4(init_p, init_v);
    }
    else {        
        
        // boid of interest
        vec4 boid = getBoid(id);
        vec2 p = boid.xy;		// boid position
        vec2 v = boid.zw;		// boid velocity
        
        
        // initialize velocity update vectors
        vec2 v_collision,v_cohesion,v_alignment = vec2(0.0f,0.0f);

        // loop over neighboring boids and build update vectors accordingly
        for (float nid = 0.0f; nid < numboids; nid++) {
            
            if (nid != id) {
                vec4 neighbor = getBoid(nid);
            	vec2 p_neighbor = neighbor.xy;	// nieghboring boid position
            	vec2 v_neighbor = neighbor.zw;	// neighboring boid velocity
            	
                vec2 sep = p_neighbor - p;
                float ls = length(sep);
                float r = 50.f;				// separation radius
          
            	// Collision Update
                if (ls < r) {
                	v_collision -= sep;  
                }
            	
                // Cohesion Update
                v_cohesion += sep;
                
                // Alignment Update
                v_alignment += v_neighbor;
            }
            		
            // Perceived values
            v_cohesion = v_cohesion/(numboids - 1.0f);
            v_alignment = v_alignment/(numboids - 1.0f) - v;


            // Apply position and velocity updates
            v = v + a1*v_collision + a2*v_cohesion + a3*v_alignment;
            v = v / length(v) * speed;		// constrain velocity
            p = p + iTimeDelta * v;			// update position
            
            
            
            // Impose wrapped boundary conditions
            if (p.x < 0.) {
            	p.x = iResolution.x + p.x;        
            }
            if (p.x > iResolution.x) {
            	p.x = p.x - iResolution.x;       
            }
             if (p.y < 0.) {
            	p.y = iResolution.y + p.y;       
            }
            if (p.y > iResolution.y) {
            	p.y = p.y - iResolution.y;       
            }
 

        	fragColor = vec4(p,v);
        }
        
    }
}