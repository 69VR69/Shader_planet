#iChannel0 "file://bufferA.glsl"

// BufferB
void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    fragColor = vec4(0, 0, 0, 0);

    for(int i = -1; i <= 1; i++)
    {
        for(int j = -1; j <= 1; j++)
        {
            vec2 coord = fragCoord + vec2(float(i), float(j));

            coord = mod(coord, iResolution.xy);

            vec2 uv = coord / iResolution.xy;

            vec4 data = texture(iChannel0, uv);

            vec2 pos = data.xy;

            if(data.x > 0.001 && abs(pos.x - fragCoord.x) < 0.5 && abs(pos.y - fragCoord.y) < 0.5)
            {
                fragColor = data;
            }
        }
    }
}