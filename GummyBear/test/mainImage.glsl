#iChannel0 "file://bufferB.glsl"

float sdBox( in vec2 p, in vec2 b )
{
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}

float opUnion( float d1, float d2 )
{
    return min(d1,d2);
}

// MainImage
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;

    fragColor = vec4(0.0, 0.0, 0.0, 1.0);

    vec4 data = vec4(texture(iChannel0, uv));

    if(data.x > 0.001)
    {
        vec2 pos = data.zw;

        // Set the color depending of the part of the screen using the uv coordinates
        vec3 color = vec3(0.0);
        
        vec2 diff = uv-pos;

        if(diff.x > 0. && diff.y > 0.)
        {
            color = vec3(0.82, 0.07, 0.75);
        }
        else if(diff.x > 0. && diff.y < 0.)
        {
            color = vec3(0.0, 1.0, 0.0);
        }
        else if(diff.x < 0. && diff.y > 0.)
        {
            color = vec3(0.0, 0.0, 1.0);
        }
        else if(diff.x < 0. && diff.y < 0.)
        {
            color = vec3(1.0, 0.0, 0.0);
        }

        // Use sdBox to draw the particles
        float d = sdBox(uv - pos, vec2(0.4));

        float t = opUnion(d,length(uv-pos));
        if(t < 0.0)
        {
            color = vec3(1.0);
        }
        // Draw the particles< size)
         fragColor = vec4(color, 1.0);
    }
}