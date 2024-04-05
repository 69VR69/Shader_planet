
vec4 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return vec4(c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y),1.);
}


void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    // Create a rainbow gradient slowly moving from left to right following an arc
    vec2 uv = (fragCoord.xy / iResolution.xy) * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y;
    float angle = atan(uv.y, uv.x);
    float dist = length(uv);
    float hue = (angle + 3.14159) / 6.28318;
    float sat = dist;
    float val = 1.0;

    // Add the movement of the rainbow
    hue += iTime * 0.1;

    // Add a little bit of noise to the saturation
    sat += 0.2 * cos(iTime * 0.9 + dist * 12.0);

    // Add a little bit of noise to the value
    val += 0.2 * cos(iTime * 0.5 + dist * 10.0);

    // Add a little bit of noise to the hue
    hue += 0.2 * sin(iTime * 0.5 + dist * 6.0);

    // Clamp the values
    hue = mod(hue, 1.0);
    sat = clamp(sat, 0.0, 1.0);
    val = clamp(val, 0.0, 1.0);

    // Convert the HSV color to RGB
    fragColor = hsv2rgb(vec3(hue, sat, val));
}