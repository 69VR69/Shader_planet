#iChannel0 "file://bufferA.glsl"

#define numboids 20.		// number of boids (must be integer value represented as float)

float sdBox(vec3 p, vec3 b)
{
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float sdBox( in vec2 p, in vec2 b )
{
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}

// A function that defines whether the screen coordinate uv is inside a boid.
// boid.xyz = boid screenspace coordinate
// boid.w = boid id
bool drawBoid(vec3 uv, vec4 boid)
{
    // Box
    vec3 bounds = vec3(10.0, 10.0, 10.0);
    return sdBox(uv.xy - boid.xy, bounds.xy) < 0.0;
}

bool drawBoid(vec2 uv, vec4 boid)
{
    // Box
    vec2 bounds = vec2(10.0, 10.0);
    return sdBox(uv - boid.xy, bounds) < 0.0;
}

vec4 getBoid(float id)
{
    return texture(iChannel0, vec2(id + .5, .5) / iResolution.xy);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
  
    
    //Background Color
    fragColor = vec4(0.1f, 0.1f, 0.1f, 1.f);
    
    
    for (float i = 0.f; i < numboids; i++) {
    	
        vec4 boid = getBoid(i);
        
        //Boid color
    	vec4 col = vec4(boid.x/iResolution.x, 0., boid.y/iResolution.y,1.f);
    	
        if (drawBoid(fragCoord, boid)) fragColor = col;           
    
    }
    

}