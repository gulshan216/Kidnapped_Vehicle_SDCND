/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 5000;

	default_random_engine gen;

	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for(int i=0; i<num_particles;i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;

		particles.push_back(particle);
		weights.push_back(1.0);
	}
	is_initialized = true;	

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	for(int i=0; i<num_particles; i++)
	{
		double new_x,new_y,new_theta;
		if(yaw_rate == 0)
		{
			new_x = particles[i].x + velocity*cos(particles[i].theta)*delta_t;
			new_x = particles[i].y + velocity*sin(particles[i].theta)*delta_t;
			new_theta = particles[i].theta;
		}
		else
		{
			new_x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			new_y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			new_theta = particles[i].theta + yaw_rate*delta_t; 
		}
		normal_distribution<double> dist_x(new_x, std_pos[0]);
		normal_distribution<double> dist_y(new_y, std_pos[1]);
		normal_distribution<double> dist_theta(new_theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	 //  double highest_weight = -1.0;
	 //  Particle best_particle;
	 //  double weight_sum = 0.0;
	 //  for (int i = 0; i < num_particles; ++i) {
		// if (particles[i].weight > highest_weight) {
		// 	highest_weight = particles[i].weight;
		// 	best_particle = particles[i];
		// }
		// weight_sum += particles[i].weight;
	 //  }
	 //  cout << "highest w before update weights " << highest_weight << endl;
	 //  cout << "average w before update weights" << weight_sum/num_particles << endl;

	for(int i=0; i<particles.size();i++)
	{
		vector<LandmarkObs> trans_particle_obs;
		for(int j=0; j<observations.size(); j++)
		{	
			LandmarkObs transformed_map_obs;

			transformed_map_obs.x = particles[i].x + observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta);
			transformed_map_obs.y = particles[i].y + observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta);
			trans_particle_obs.push_back(transformed_map_obs);
		}
		if(i == 0)
		{
			cout << "observations x: " << observations[0].x << " y:" <<  observations[0].y << endl;
			cout << "transformed observations x: " << trans_particle_obs[0].x << " y:" <<  trans_particle_obs[0].y << endl;
		}

		vector<LandmarkObs> landmarks_in_range;
		for(int k=0; k<map_landmarks.landmark_list.size(); k++)
		{
			if(dist(particles[i].x,particles[i].y, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f) < sensor_range)
			{
				LandmarkObs landmark_obs;
				landmark_obs.id = map_landmarks.landmark_list[k].id_i;
				landmark_obs.x = map_landmarks.landmark_list[k].x_f;
				landmark_obs.y = map_landmarks.landmark_list[k].y_f;

				landmarks_in_range.push_back(landmark_obs);
			}
		}

		if(i == 0)
		{
			cout << "particle position x: " << particles[i].x << " y:" <<  particles[i].y << endl;
			cout << "Landmark in range x: " << landmarks_in_range[0].x << " y:" <<  landmarks_in_range[0].y << endl;
		}

		for(int l=0; l<trans_particle_obs.size();l++)
		{
			LandmarkObs obs;
			obs = trans_particle_obs[l];
			double min_dist = dist(obs.x,obs.y,landmarks_in_range[0].x,landmarks_in_range[0].y);
			int assoc_landmark_id = 0;
			for(int m=0; m<landmarks_in_range.size();m++)
			{
				double euclidean_dist = dist(obs.x,obs.y,landmarks_in_range[m].x,landmarks_in_range[m].y);
				if(euclidean_dist < min_dist)
				{
					min_dist = euclidean_dist;
					assoc_landmark_id = m;
				}
			}
			trans_particle_obs[l].id = assoc_landmark_id;
		}

		particles[i].weight = 1.0;
		for(int n=0;n<trans_particle_obs.size();n++)
		{

			double mu_x = landmarks_in_range[trans_particle_obs[n].id].x;
			double mu_y = landmarks_in_range[trans_particle_obs[n].id].y;

			double exp_a = pow((trans_particle_obs[n].x - mu_x),2)/(2*std_landmark[0]*std_landmark[0]);
			double exp_b = pow((trans_particle_obs[n].y - mu_y),2)/(2*std_landmark[1]*std_landmark[1]);
			particles[i].weight = particles[i].weight*((exp(-exp_a-exp_b))/(2*M_PI*std_landmark[0]*std_landmark[1]));
		}
		weights[i] = particles[i].weight;
	}

	 //  highest_weight = -1.0;
	 //  weight_sum = 0.0;
	 //  for (int i = 0; i < num_particles; ++i) {
		// if (particles[i].weight > highest_weight) {
		// 	highest_weight = particles[i].weight;
		// 	best_particle = particles[i];
		// }
		// weight_sum += particles[i].weight;
	 //  }
	 //  cout << "highest w after update weights " << highest_weight << endl;
	 //  cout << "average w after update weights" << weight_sum/num_particles << endl;

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// double highest_weight = -1.0;
	//   Particle best_particle;
	//   double weight_sum = 0.0;
	//   for (int i = 0; i < num_particles; ++i) {
	// 	if (particles[i].weight > highest_weight) {
	// 		highest_weight = particles[i].weight;
	// 		best_particle = particles[i];
	// 	}
	// 	weight_sum += particles[i].weight;
	//   }
	//   cout << "highest w before resample " << highest_weight << endl;
	//   cout << "average w before resample" << weight_sum/num_particles << endl;

	default_random_engine gen;
	discrete_distribution<int> dist(weights.begin(), weights.end());

	vector<Particle> resample_particles;
	for(int i=0; i<num_particles; i++)
	{
		resample_particles.push_back(particles[dist(gen)]);
	}

	particles = resample_particles;

	// highest_weight = -1.0;
	// weight_sum = 0.0;
	//   for (int i = 0; i < num_particles; ++i) {
	// 	if (particles[i].weight > highest_weight) {
	// 		highest_weight = particles[i].weight;
	// 		best_particle = particles[i];
	// 	}
	// 	weight_sum += particles[i].weight;
	//   }
	//   cout << "highest w after resample " << highest_weight << endl;
	//   cout << "average w after resample" << weight_sum/num_particles << endl;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
