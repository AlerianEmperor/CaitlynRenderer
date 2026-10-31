#version 460

out vec4 color;
in vec2 tex;

//sampler chi la ham de sample chu ko bind voi 1 texture cu the
//va no cung khong chua data cua texture

//uniform sampler2D pathTraceTexture2;
//uniform sampler2D pathTraceTexture4;
//uniform sampler2D pathTraceTexture16;

uniform sampler2D pathTraceTexture;

void main()
{
	color = texture(pathTraceTexture, tex);
}