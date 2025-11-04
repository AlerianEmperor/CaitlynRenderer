//#include "Scene_SSBO.h"
//#include "Scene_SSBO_SOA.h"
#include "Scene_SSBO_extend_center.h"
#include <gl\glfw3.h>
#include <glm\glm.hpp>
#include <time.h>
#include <iostream>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_internal.h"
 

float moveSpeed = 2.5f;
float mouseSensitivity = 0.05f;
bool keyPressed = false;

//ivec2 screen_size(500, 500);
ivec2 screen_size(screen_width, screen_height);

//RenderOptions renderOptions;

//Quad* quad;
//Camera* camera;

//FBO
//GLuint path_trace_fbo;

//texture of FBO
//GLuint path_trace_texture;

Scene scn;

//string shader_direction = "D:/a_gpu_path_tracing/4.GLSL_Final/GLSL_Final_Version/";
//string shader_direction = "D:\\a_a_opengl\\ky thuat quan trong\\Caitlyn Renderer\\test\\GLSL - Path - Tracer - main\\GLSL - Path - Tracer - main\\6b_SBVH_1_Leaf_Node\\build\\glsl_sbvh\\";
//string shader_direction = "D:/a_a_opengl/ky thuat quan trong/Caitlyn Renderer/test/GLSL-Path-Tracer-main/GLSL-Path-Tracer-main/6b_SBVH_1_Leaf_Node/build/glsl_sbvh/";//"D:/a_gpu_path_tracing_glsl/4.GLSL_Final/6.GLSL_SBVH/SBVH/SBVH/";
//string scene_direction = "D:/a_gpu_path_tracing/2._Robert_beckanb/OpenGL-PathTracer-master/assets/";//"D:/a_gpu_path_tracing/assets/";

string shader_direction = "";

string scene_direction = "D:\\a_a_opengl\\model\\cornell-box.obj";
//string scene_direction = "E:\\Models\\Cornell_Box\\cornell-box.obj";  //cornell_box_correct.obj";
//string scene_direction = "E:\\Models\\fireplace_room\\FireRoom\\fireplace_room.obj";
//string scene_direction = "E:\\Models\\Bath_Room\\contemporary_bathroom_lux_bikini.obj";
//string scene_direction = "E:\\Models\\Lamp\\blender\\Little_Lamp3.obj";

void init_scene(int index)
{
	string sceneFilenames[] =
	{
		"cornell.scene",
		"ajax.scene",
		"bathroom.scene",
		"boy.scene",
		"coffee.scene",
		"diningroom.scene",
		"glassBoy.scene",
		"hyperion.scene",
		"rank3police.scene",
		"spaceship.scene",
		"staircase.scene"
	};

	scn = Scene(scene_direction, shader_direction);

	//getchar();
	/*cout << scn.screen_width << " " << scn.screen_height << "\n";
	cout << scn.camera.position.x << " " << scn.camera.position.y << " " << scn.camera.position.z << "\n";
	cout << scn.camera.forward.x << " " << scn.camera.forward.y << " " << scn.camera.forward.z << "\n";
	cout << "Mat Size:\n";
	cout << scn.mats.size() << "\n";
	cout << "light:\n";
	cout << scn.lights.size() << "\n";
	cout << scn.lights[0].p.x << "\n";
	cout << scn.lights[0].e.x << "\n";
	cout << "Mesh:\n";
	cout << scn.vertices.size() << "\n";
	cout << scn.normal_texcoords.size() << "\n";
	cout << scn.triangles.size() << "\n";

	cout << "Test Mat:\n";

	for (int i = 0; i < scn.mats.size(); ++i)
	cout << scn.mats[i].color.x << " " << scn.mats[i].color.y << " " << scn.mats[i].color.z << "\n";*/

	//scn.Read_Scene("E:\Models\Cornell_Box\cornell-box.obj");

	//scn = Scene(scene_direction, shader_direction, sceneFilenames[index]);

	/*scn = Scene(scene_direction, sceneFilenames[index]);
	//Test
	cout << scn.vertices.size() << "\n";
	cout << scn.triangle_indices.size() << "\n";
	cout << scn.normals_texcoords.size() << "\n";
	cout << scn.mats.size() << "\n";
	cout << scn.lights.size() << "\n";*/
}

void render(GLFWwindow *window)
{
	scn.Render();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, screen_size.x, screen_size.y);
	//scn.display();

	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

	glfwSwapBuffers(window);
	glClear(GL_COLOR_BUFFER_BIT);
}

ImGuiIO mouse_io;


void update(float second, GLFWwindow* window)
{
	scn.update(second);

	//keyPressed = false;
	
	scn.camera.isMoving = false;

	if (glfwGetKey(window, 'W'))
	{
		scn.camera.offsetPosition(second * moveSpeed * scn.camera.forward);
		//keyPressed = true;
		scn.camera.isMoving = true;
		//glClear(GL_COLOR_BUFFER_BIT);
	}
	else if (glfwGetKey(window, 'S'))
	{
		scn.camera.offsetPosition(second * moveSpeed * -scn.camera.forward);
		//keyPressed = true;
		scn.camera.isMoving = true;
		//glClear(GL_COLOR_BUFFER_BIT);
	}
	if (glfwGetKey(window, 'A'))
	{
		scn.camera.offsetPosition(second * moveSpeed * -scn.camera.right);
		//keyPressed = true;
		scn.camera.isMoving = true;
		//glClear(GL_COLOR_BUFFER_BIT);
	}
	else if (glfwGetKey(window, 'D'))
	{
		scn.camera.offsetPosition(second * moveSpeed * scn.camera.right);
		//keyPressed = true;
		scn.camera.isMoving = true;
		//glClear(GL_COLOR_BUFFER_BIT);
	}

	

	//if (!ImGui::IsWindowHovered() && ImGui::IsMouseDown(1))
	//bool b = ImGui::IsWindowHovered();
	
	//ImVec2 min_window = ImGui::GetWindowContentRegionMin();
	//ImVec2 max_window = ImGui::GetWindowContentRegionMax();
	//scn.camera.isMoving = false;
	//https://github.com/ocornut/imgui/issues/1382
	if(!ImGui::IsWindowHovered(4) && ImGui::IsMouseDown(0))
	{
		ImVec2 mouse_delta = ImGui::GetMouseDragDelta();
		if (mouse_delta.x != 0 || mouse_delta.y != 0)
		{
			scn.camera.offsetOrientation(mouseSensitivity * mouse_delta.x, mouseSensitivity * mouse_delta.y);
			scn.camera.isMoving = true;
			ImGui::ResetMouseDragDelta();
		}
	}
}


void main()
{
	srand(time(0));

	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window;
	window = glfwCreateWindow(screen_size.x, screen_size.y, "Caitlyn", NULL, NULL);

	if (window == NULL)
	{
		fprintf(stderr, "Failed to Create OpenGL Context");
		getchar();
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		getchar();
		return;
	}

	GLint max;
	glGetIntegerv(GL_MAX_SHADER_STORAGE_BUFFER_BINDINGS, &max);
	printf("Max SSBO bindings = %d\n", max);

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();
	(void)io;

	ImGui::StyleColorsDark();

#if __APPLE__
	// GL 3.2 + GLSL 150
	const char* glsl_version = "#version 150";
	/*glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
	*/
#else
	// GL 3.0 + GLSL 130
	const char* glsl_version = "#version 130";
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	//glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
	//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif

	ImGui_ImplGlfw_InitForOpenGL(window, true);
	ImGui_ImplOpenGL3_Init(glsl_version);



	//everything must happend below this line!

	//set up
	//string shader_source = "D:/a_gpu_path_tracing/1_Sang_GLSL_Final/GLSL_Final/";
	//output_shader = new Shader(shader_source + "output.vs", shader_source + "output.fs");
	//path_trace_shader = new Shader(shader_source + "path_trace.vs", shader_source + "path_trace.fs");

	int currentSceneIndex = 0;

	init_scene(currentSceneIndex);

	//quad = new Quad();




	//set_up_camera(scn.path_trace_shader);

	//init_scene();

	//end set up

	double lastTime = glfwGetTime();
	double startTime = glfwGetTime();

	cout << "Reach This Line\n";

	while (!glfwWindowShouldClose(window))
	{
		glfwPollEvents();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();
		{
			ImGui::Begin("Caitlyn PathTracer");

			ImGui::Text("Aplication average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
			ImGui::Text("Total Frame %d", scn.frame_count);
			ImGui::Text("Render Time %f", glfwGetTime() - startTime);

			if (ImGui::Combo("Scene", &currentSceneIndex, "FireRoom\0ajax\0bathroom\0boy\0coffee\0diningroom\0glassBoy\0hyperion\0rank3police\0spaceship\0staircase\0"))
			{

			}
			bool render_options_changed = false;
			render_options_changed |= ImGui::Combo("Render Type", &renderOptions.rendererType, "Progressive\0Tiled\0");
			render_options_changed |= ImGui::InputInt2("Resolution", &renderOptions.resolution.x);
			render_options_changed |= ImGui::InputInt("Max Samples", &renderOptions.maxSamples);
			render_options_changed |= ImGui::InputInt("Max Depth", &renderOptions.maxDepth);
			render_options_changed |= ImGui::InputInt("Tiles X", &renderOptions.numTilesX);
			render_options_changed |= ImGui::InputInt("Tiles Y", &renderOptions.numTilesY);
			render_options_changed |= ImGui::Checkbox("Use envmap", &renderOptions.useEnvMap);
			render_options_changed |= ImGui::InputFloat("HDR multiplier", &renderOptions.hdrMultiplier);

			ImGui::End();
		}

		if (glfwGetKey(window, GLFW_KEY_ESCAPE))
			glfwSetWindowShouldClose(window, GL_TRUE);
		double presentTime = glfwGetTime();

		update((float)(presentTime - lastTime), window);
		lastTime = presentTime;
		render(window);
	}

	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGlfw_Shutdown();
	ImGui::DestroyContext();

	scn.delete_tex_data();
	glfwTerminate();

}

//Test Functions
//resize image to the same size, so they can be put in texture array

void main_test_resize_image()
{
	int w, h, n;

	string s = "E:\\Models\\fireplace_room\\FireRoom\\wood.png";

	unsigned char* image = stbi_load(s.c_str(), &w, &h, &n, 0);

	/*ofstream ofs("original_image.ppm");

	ofs << "P3\n" << w << " " << h << "\n255\n";

	vector<vec3> original;

	float inv_255 = 1.0f / 255.0f;

	for (int i = 0; i < w * h * 3; i += 3)
	{
		float x = 255 * (image[i] * inv_255);
		float y = 255 * (image[i + 1] * inv_255);
		float z = 255 * (image[i + 2] * inv_255);

		original.emplace_back(vec3(x, y, z));

		ofs << x << " " << y << " " << z << "\n";
	}*/

	int require_width = 628, require_height = 500;

	vector<vec3> resize;// (require_width * require_height);
	
	resize_image(image, w, h, resize, require_width, require_height);

	ofstream ofs2("resize_image_628_501.ppm");

	ofs2 << "P3\n" << require_width << " " << require_height << "\n255\n";

	for (int i = 0; i < require_width * require_height; ++i)
	{
		vec3 v = resize[i];
		
		ofs2 << v.x << " " << v.y << " " << v.z << "\n";
	}
}

