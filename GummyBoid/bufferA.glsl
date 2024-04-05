#iChannel0 "file://bufferA.glsl"
#iChannel1 "file://bufferB.glsl"

#define numboids 50.		        // Number of boids (must be integer value represented as float)
#define speed 5.		            // Boid speed
#define a1 numboids*0.0000066667    // Collision factor
#define a2 numboids*0.0000666667    // Cohesion factor	
#define a3 numboids*0.0003333333    // Alignment factor
#define a4 numboids*0.0000016667    // Center attraction factor

vec3 getBoidPosition(float id)
{
    return texture(iChannel1, vec2(id + .5, .5) / iResolution.xy).xyz;
}

vec3 getBoidVelocity(float id)
{
    return texture(iChannel0, vec2(id + .5, .5) / iResolution.xy).xyz;
}

// Then the main image program queries the boid data from this buffer.
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{

    fragColor = vec4(0.0, 0.0, 0.0, 1.0);			// default data
    if(fragCoord.y > 1. || fragCoord.x > numboids)
        discard;

    float id = floor(fragCoord.x);

    // initialize random boid positions and velocities
    if(iFrame <= 1)
    {
        vec3 init_v = vec3(sin(id), cos(id), -cos(id));
        init_v /= speed * length(init_v);
        fragColor = vec4(init_v, 1.);
    }
    else
    {        

        // boid of interest
        vec3 v = getBoidVelocity(id);
        vec3 p = getBoidPosition(id);

        // Boundary conditions
        vec3 maxBounds = vec3(1.5, 1.0, 1.0)*6.;
        vec3 minBounds = (-maxBounds) + vec3(0.0, 4.5, 0.0);

        // initialize velocity update vectors
        vec3 v_collision, v_cohesion, v_alignment, v_center_attraction = vec3(0.);

        // loop over neighboring boids and build update vectors accordingly
        for(float nid = 0.0; nid < numboids; nid++)
        {

            if(nid != id)
            {
                vec3 p_neighbor = getBoidPosition(id);	// nieghboring boid position
                vec3 v_neighbor = getBoidVelocity(id);	// neighboring boid velocity

                vec3 sep = p_neighbor - p;
                float ls = length(sep);
                float r = 2.;				// separation radius

            	// Boid Collision Update
                if(ls < r)
                {
                    v_collision -= sep;
                }

                // Boundary Collision Update
                for(int i = 0; i < 3; i++)
                {
                    if(p[i] > maxBounds[i])
                        v_collision[i] -= 1.;
                    if(p[i] < minBounds[i])
                        v_collision[i] += 1.;
                }

                // Cohesion Update
                v_cohesion += sep;

                // Alignment Update
                v_alignment += v_neighbor;
            }

            // Perceived values
            v_cohesion = v_cohesion / (numboids - 1.0);
            v_alignment = v_alignment / (numboids - 1.0) - v;

            // Add an attractive force to the center of the screen
            v_center_attraction = vec3(0., 2., 0.) - p;

            // Apply position and velocity updates
            v += a1 * v_collision + a2 * v_cohesion + a3 * v_alignment + a4 * v_center_attraction;
            v /= length(v) * speed;		// constrain velocity
        }

        // Add a bit of randomness to the boid movement
        v += 0.01 * vec3(sin(id), cos(id), -cos(id));

        fragColor = vec4(v.xyz, 1.);

    }
}