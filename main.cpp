/*
Author: Nicholas Layden, December 2022
Compilation command:

clang++ raylibtest2.cpp `pkg-config --libs --cflags raylib` -std=c++20 -o test

-------------------------------------------------------------------------------
Testing framework for real time simulation of particles.
Particles contain multiple interaction settings (forces),
different boundary conditions (reflective, toroidal, spherical, etc)




-------------------------------------------------------------------------------*/
// Standard libraries
#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <format>

// Graphics libraries
#include "raylib.h"

// Custom libraries

double unirandomval(double min, double max)
{
  /// This function uses a Mersenne Twister to generate
  /// a uniform real distribution between min and max
  std::random_device rand_device;
  std::mt19937 gen(rand_device());
  std::uniform_real_distribution<double> distr(min, max);
  return distr(gen);
}

int main(){

	// Setup window context
	const int screenWidth = 800;
	const int screenHeight = 800;
	const int frameRate = 60;
	double frameTime = 1/((float)frameRate);
	uint FrameCounter = 0;
	SetTargetFPS(frameRate);
	InitWindow(screenWidth,screenHeight,"Particle Simulator");

	// Particle data
	int Nparticles = 2000;
	int ParticleSize = 4;
	float vmin = -2;
	float vmax = 2;
	float pmin_x = screenWidth*0.45;
	float pmax_x = screenWidth*0.55;
	float pmin_y = screenHeight*0.45;
	float pmax_y = screenHeight*0.55;
	
	// Particle containers
	std::vector<Vector2> particles(Nparticles);
	std::vector<Vector2> particles_v(Nparticles);
	std::vector<Vector2> particles_f(Nparticles);
	std::vector<Vector2> particles_m(Nparticles);
	Vector2 particle;
	std::vector<Vector2> particle0_position;

	// numerical properties
	float distance,dx,dy;
	float softener = 10;

	// Initial data
	for (int i = 0; i < Nparticles; ++i) {
		particles[i].x = unirandomval(pmin_x,pmax_x);
		particles[i].y = unirandomval(pmin_y,pmax_y);
		particles_v[i].x = unirandomval(vmin,vmax);
		particles_v[i].y = 0;//unirandomval(vmin,vmax);
		particles_f[i].x = 0;
		particles_f[i].y = 0;
		particles_m[i].x = 1;
		particles_m[i].y = 0;
	}
	particles_m[0].x = 50;
	// Tracking particle 0 with a line
	particle0_position.push_back(particles[0]);

	// Main window loop
	while(!WindowShouldClose()) {
		// Main update loop
		FrameCounter++;
		// Calculate any particle-particle interactions
		for (int i = 0; i < Nparticles; ++i) {
			for (int j = 0; j < Nparticles; ++j) {
				if (i != j) {
					// distance between particles i and j
					dx = particles[i].x - particles[j].x;
					dy = particles[i].y - particles[j].y;
					distance = sqrt(dx*dx + dy*dy) + softener	;
					// Compute forces, no self force!!
					particles_f[i].x -= 100.*particles_m[i].x*particles_m[j].x/(distance*distance*distance)*dx;
					particles_f[i].y -= 100.*particles_m[i].x*particles_m[j].x/(distance*distance*distance)*dy;
				} // end if
			} // end j loop
		} // end i loop



		// Update particle positions from velocities/forces
		for (int i = 0; i < Nparticles; ++i) {
			// Update velocities
			// particles_v[i].y += 20*frameTime; 
			particles_v[i].x += frameTime*particles_f[i].x;
			particles_v[i].y += frameTime*particles_f[i].y;
			// Update positions
			particles[i].x += particles_v[i].x;
			particles[i].y += particles_v[i].y;
			// Check if particles 'hit' boundary
			// x boundaries
			if ((particles[i].x >= (GetScreenWidth() - ParticleSize) )|| ( particles[i].x <= ParticleSize )) {
				particles_v[i].x *= -1.0f;
			}
			// y boundaries
			if ((particles[i].y >= (GetScreenWidth() - ParticleSize) )|| ( particles[i].y <= ParticleSize )) {
				particles_v[i].y *= -1.0f;
			}
			// Reset forces after velocities changed
			particles_f[i].x = 0.;
			particles_f[i].y = 0.;
		}
		particle0_position.push_back(particles[0]);

		// Main drawing loop
		BeginDrawing();
		ClearBackground(BLACK);
		DrawFPS(10,10);
		// Line for the 'heavy' particle
		if (particle0_position.size() >= 2) {
			for (int i = 0; i < particle0_position.size()-1; ++i) {
				float xs,ys,xe,ye;
				xs = particle0_position[i].x;
				ys = particle0_position[i].y;
				xe = particle0_position[i+1].x;
				ye = particle0_position[i+1].y;
				DrawLine(xs,ys,xe,ye,RAYWHITE);
			} // END I
		} // END IF

		// Draw the particles' current positions 
		for (int i = 0; i < Nparticles; ++i) {
			DrawCircle(particles[i].x,particles[i].y,ParticleSize,MAROON);
		} // end particle loop

		EndDrawing();

	} // end window loop

	// Free variables and allocated memory

	CloseWindow();

	return 0;
}

