#version 460

layout (location = 0) in vec2 position;
layout (location = 1) in vec2 texCoords;

out vec2 tex;

void main()
{
    gl_Position = vec4(position.x,position.y,0.0,1.0);
	tex = texCoords;
}