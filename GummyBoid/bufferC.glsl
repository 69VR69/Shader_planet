#iChannel0 "file://bufferC.glsl"
/*
#define GRAVITY -5.0
#define ATTRACTION 60.0


// Random function from https://www.shadertoy.com/view/4ssXRX
float nrand(vec2 n) { 
    return fract(sin(dot(n, vec2(12.9898, 4.1414))) * 43758.5453);
}
vec4 resetParticules(vec2 id)
{
    return vec4(nrand(id),3.*nrand(id.yx),1.,1.);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec2 pixelSize = 1.0 / iResolution.xy;
    
    // Initialize the texture on the first frame with random positions and null velocity.
    if (iFrame == 1)
    {
        fragColor = resetParticules(uv);
        return;
    }
    
    vec4 previousFrameValues = texture(iChannel0, uv);
    vec3 position = previousFrameValues.xyz;
    float velocity = previousFrameValues.w;
 
    // Gravity.
    velocity += GRAVITY;

    float randValue = nrand(uv.yx * iTime) * 0.5;

    // Update position.
    position.y += velocity * iTimeDelta * 0.001 * noise2(position);
    position.x += iTimeDelta * 0.001;

    if(position.y < -1.0)
    {
        position.y = resetParticules(uv).y;
    }
    fragColor = vec4(position.xyz, velocity);
}*/