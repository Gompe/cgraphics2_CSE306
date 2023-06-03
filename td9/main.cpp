#define _CRT_SECURE_NO_WARNINGS 1

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <atomic>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../stb/stb_image.h"

#ifndef N_BOUNCES
	#define N_BOUNCES 5
#endif

#include "gx_camera.h"
#include "gx_scene.h"
#include "tutte.h"

struct RefScene {
    Camera *camera_ptr;
    Scene *scene_ptr;
    int W;
    int H;

    RefScene(Camera *camera_ptr, Scene *scene_ptr, int W, int H)
    : camera_ptr(camera_ptr), scene_ptr(scene_ptr), W(W), H(H)
    {}

};

std::vector<unsigned char> renderImage(RefScene refScene)
{
	int W = refScene.W;
	int H = refScene.H;
	Camera *camera_ptr = refScene.camera_ptr;
	Scene *scene_ptr = refScene.scene_ptr;

	std::vector<unsigned char> image(W * H * 3, 0);

	int NUMBER_OF_RAYS;
	std::cout << "Select the desired number of rays per pixel: ";
	std::cin >> NUMBER_OF_RAYS;

	auto start = std::chrono::high_resolution_clock::now();

	#pragma omp parallel for
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			Vector color = Vector(0,0,0);

			for(int _counter=0; _counter<NUMBER_OF_RAYS; _counter++){
				Ray r = camera_ptr->generate_ray(i, j);

				Vector vec = scene_ptr->getColor(r, N_BOUNCES);
				color+= vec;
			}

			color = color / NUMBER_OF_RAYS;
			color = color.clip(0., 255.);
			color /= 255.;
			
			color[0] = std::pow(color[0], 0.45);
			color[1] = std::pow(color[1], 0.45);
			color[2] = std::pow(color[2], 0.45);

			color *= 255.;

			image[(i * W + j) * 3 + 0] = std::min(color[0], 255.);
			image[(i * W + j) * 3 + 1] = std::min(color[1], 255.);
			image[(i * W + j) * 3 + 2] = std::min(color[2], 255.);
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto time_elapsed = std::chrono::duration<double>(end - start).count();
	printf("---Rendering complete. Waited %.1fs.\n", time_elapsed);

	return image;
}

int main(int argc, char **argv) {

	int W = 512;
	int H = 512;

	TriangleMesh *mmesh = new TriangleMesh("./goethe.obj", Material(Vector(0.2, 0.8, 0.8)));
	Object *mesh = mmesh;
	mesh->transformScale(0.15);

	TutteEmbedding(*mmesh, "goethe_tutte.obj");
	
	return 0;
}