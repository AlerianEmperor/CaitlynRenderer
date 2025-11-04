#version 460

out vec4 color;
in vec2 tex;

uniform sampler2D pathTraceTexture;
uniform float invSampleCounter;

vec4 ToneMap(in vec4 c, float limit)
{
	float luminance = 0.3*c.x + 0.6*c.y + 0.1*c.z;

	return c * 1.0/(1.0 + luminance/limit);
}

void main()
{
	color = texture(pathTraceTexture, tex) * invSampleCounter;
	color = pow(ToneMap(color, 2.0), vec4(1.0 / 2.2));

	//color = vec4(0.0f, 1.0f, 0.0f, 1.0f);

	//color = texture(pathTraceTexture, tex);
}