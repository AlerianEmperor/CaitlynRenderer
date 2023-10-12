#ifndef SHADER_H
#define SHADER_H

#include "glad.h"//<glad\glad.h>
#include <glm/glm.hpp>

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
//#include <gl\freeglut.h>

using namespace std;
using namespace glm;

class Shader
{
public:
	unsigned int id;
	// constructor generates the shader on the fly
	// ------------------------------------------------------------------------
	Shader() {}
	//Shader(const char* vertexPath, const char* fragmentPath, const char* geometryPath = nullptr)
	Shader(string vs_name, string fs_name, string gs_name = "")
	{
		GLuint vs = CreateShader(GL_VERTEX_SHADER, vs_name);
		GLuint fs = CreateShader(GL_FRAGMENT_SHADER, fs_name);

		id = glCreateProgram();

		glAttachShader(id, vs);
		glAttachShader(id, fs);

		if (gs_name != "")
		{
			GLuint gs = CreateShader(GL_GEOMETRY_SHADER, gs_name);

			glAttachShader(id, gs);
		}

		id = LinkingProgram(id);
	}
	Shader(string cs_name)//compute shader
	{
		GLuint cs = CreateShader(GL_COMPUTE_SHADER, cs_name);

		id = glCreateProgram();

		glAttachShader(id, cs);

		id = LinkingProgram(id);
	}
	string Read_Shader_Source(const string& filename)
	{
		ifstream file(filename, ios::in);
		string line = "";
		string result = "";

		while (!file.eof())
		{
			getline(file, line);
			result += line + "\n";
		}
		//cout << result;
		return result;
	}

	GLuint CreateShader(int shaderType, const string& ShaderPath)
	{
		GLint ShaderCompiled;
		string ShaderSource = Read_Shader_Source(ShaderPath);
		const char* source = ShaderSource.c_str();
		//cout << shaderType << "\n";

		GLuint ShaderReference = glCreateShader(shaderType);
		glShaderSource(ShaderReference, 1, &source, NULL);

		glCompileShader(ShaderReference);
		CheckOpenGLError();
		glGetShaderiv(ShaderReference, GL_COMPILE_STATUS, &ShaderCompiled);

		if (ShaderCompiled != 1)
		{
			cout << ShaderPath << " Error!\n";
			if (shaderType == 35633) cout << "Vertex ";
			if (shaderType == 36488) cout << "Tess Control ";
			if (shaderType == 36487) cout << "Tess Eval ";
			if (shaderType == 36313) cout << "Geometry ";
			if (shaderType == 35632) cout << "Fragment ";
			cout << "shader compilation error." << endl;
			PrintShaderLog(ShaderReference);
		}
		return ShaderReference;
	}

	//Check Errors
	bool CheckOpenGLError()
	{
		bool foundError = false;
		int glErr = glGetError();
		if (glErr != GL_NO_ERROR)
		{
			//cout << "glErr: " << glErr << "\n";
			foundError = true;
			glErr = glGetError();
		}
		return foundError;
	}

	void PrintShaderLog(GLuint shader)
	{
		int len = 0;
		int char_Written = 0;
		char* log;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
		if (len > 0)
		{
			log = (char*)malloc(len);
			glGetShaderInfoLog(shader, len, &char_Written, log);
			cout << "Shader Info Log " << log << "\n";
			free(log);
		}
	}

	void PrintProgramLog(GLuint program)
	{
		int len = 0;
		int char_Written = 0;
		char* log;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
		if (len > 0)
		{
			log = (char*)malloc(len);
			glGetProgramInfoLog(program, len, &char_Written, log);
			cout << "Program Info Log " << log << "\n";
			free(log);
		}
	}

	GLuint LinkingProgram(GLuint program)
	{
		GLint linked;
		glLinkProgram(program);
		CheckOpenGLError();
		glGetProgramiv(program, GL_LINK_STATUS, &linked);
		if (!linked)
		{
			cout << "linking failed!" << "\n";
			PrintProgramLog(program);
		}
		return program;
	}
	// activate the shader
	// ------------------------------------------------------------------------
	void use()
	{
		glUseProgram(id);
	}
	// utility uniform functions
	// ------------------------------------------------------------------------
	void setBool(const std::string& name, bool value) const
	{
		glUniform1i(glGetUniformLocation(id, name.c_str()), (int)value);
	}
	// ------------------------------------------------------------------------
	void setInt(const std::string& name, int value) const
	{
		glUniform1i(glGetUniformLocation(id, name.c_str()), value);
	}
	// ------------------------------------------------------------------------
	void setFloat(const std::string& name, float value) const
	{
		glUniform1f(glGetUniformLocation(id, name.c_str()), value);
	}
	// ------------------------------------------------------------------------
	void setVec2(const std::string& name, const glm::vec2& value) const
	{
		glUniform2fv(glGetUniformLocation(id, name.c_str()), 1, &value[0]);
	}
	void setVec2(const std::string& name, float x, float y) const
	{
		glUniform2f(glGetUniformLocation(id, name.c_str()), x, y);
	}
	// ------------------------------------------------------------------------
	void setVec3(const std::string& name, const glm::vec3& value) const
	{
		glUniform3fv(glGetUniformLocation(id, name.c_str()), 1, &value[0]);
	}
	void setVec3(const std::string& name, float x, float y, float z) const
	{
		glUniform3f(glGetUniformLocation(id, name.c_str()), x, y, z);
	}
	// ------------------------------------------------------------------------
	void setVec4(const std::string& name, const glm::vec4& value) const
	{
		glUniform4fv(glGetUniformLocation(id, name.c_str()), 1, &value[0]);
	}
	void setVec4(const std::string& name, float x, float y, float z, float w)
	{
		glUniform4f(glGetUniformLocation(id, name.c_str()), x, y, z, w);
	}
	// ------------------------------------------------------------------------
	void setMat2(const std::string& name, const glm::mat2& mat) const
	{
		glUniformMatrix2fv(glGetUniformLocation(id, name.c_str()), 1, GL_FALSE, &mat[0][0]);
	}
	// ------------------------------------------------------------------------
	void setMat3(const std::string& name, const glm::mat3& mat) const
	{
		glUniformMatrix3fv(glGetUniformLocation(id, name.c_str()), 1, GL_FALSE, &mat[0][0]);
	}
	// ------------------------------------------------------------------------
	void setMat4(const std::string& name, const glm::mat4& mat) const
	{
		glUniformMatrix4fv(glGetUniformLocation(id, name.c_str()), 1, GL_FALSE, &mat[0][0]);
	}
	void setIndex(const string& name, int index, int value)
	{
		string s = name + "[" + to_string(index) + "]";
		setInt(s, value);
	}
	void dispatch(u32 sx, u32 sy, u32 sz)
	{
		glDispatchCompute(sx, sy, sz);
	}
	void reset()
	{
		glUseProgram(0);
	}
private:
	// utility function for checking shader compilation/linking errors.
	// ------------------------------------------------------------------------
	void checkCompileErrors(GLuint shader, std::string type, const char *path)
	{
		GLint success;
		int infoLogLength;
		printf("checking %s ...\n", path);
		if (type != "PROGRAM")
		{
			glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
			if (infoLogLength > 0) {
				std::vector<char> errorMessage(infoLogLength + 1);
				glGetShaderInfoLog(shader, infoLogLength, NULL, &errorMessage[0]);
				printf("%s\n", &errorMessage[0]);
			}
		}
		else
		{
			glGetProgramiv(shader, GL_LINK_STATUS, &success);
			glGetProgramiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
			if (infoLogLength > 0) {
				std::vector<char> ProgramErrorMessage(infoLogLength + 1);
				glGetProgramInfoLog(shader, infoLogLength, NULL, &ProgramErrorMessage[0]);
				printf("%s\n", &ProgramErrorMessage[0]);
			}
		}
	}
};
#endif