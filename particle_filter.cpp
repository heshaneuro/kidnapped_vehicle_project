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
	// TODO: Set the number of particles. 
	// Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	default_random_engine gen;

	// Set standard deviations for x, y, and theta
	 std_x = std[0];
	 std_y = std[1];
	 std_theta = std[2];

	 num_particles = 10;
	 
	particles = std::vector<Particle>(num_particles);

	 for (int i = 0; i < num_particles; ++i) {
		double sample_x, sample_y, sample_theta;
		// create a normal (Gaussian) distribution for x, y and Theta
		normal_distribution<double> dist_x(x, std_x);
		normal_distribution<double> dist_y(y, std_y);
		normal_distribution<double> dist_theta(theta, std_theta);
		
		// Sample  and from these normal distrubtions like this: 
		//	 where "gen" is the random engine initialized earlier.
		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);

		particles[i].id = i;
		particles[i].x = sample_x;
		particles[i].y = sample_y;
		particles[i].theta = sample_theta;
		particles[i].weight = 1;

		particles[i].sense_x = std::vector<double>(11);
		particles[i].sense_y = std::vector<double>(11);
		particles[i].associations = std::vector<int>(11);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	double std_x, std_y, std_theta; // Standard deviations for x, y, and theta
	default_random_engine gen;
	double vbythetadot;
	double tempangle; 

	// Set standard deviations for x, y, and theta
	 std_x = std_pos[0];
	 std_y = std_pos[1];
	 std_theta = std_pos[2];


	 // Calculate velocity divided by Yaw rate
	 vbythetadot = velocity / yaw_rate;

	// Re-initialize for all particles before you start a new cycle
	for (int i = 0; i < num_particles; ++i) {
		particles[i].weight = 1;
		}


	// calculate predicted particles position and orientation
	for (int i = 0; i < num_particles; ++i) {
		
		double prev_x, prev_y, prev_theta, temp_angle;
		prev_x = particles[i].x;
		prev_y = particles[i].y;
		prev_theta = particles[i].theta;

		temp_angle = prev_theta + (yaw_rate * delta_t);

		if (fabs(yaw_rate) < 0.0001){
			particles[i].x = prev_x + (vbythetadot * (sin(temp_angle) - sin(prev_theta)));
		particles[i].y = prev_y + (vbythetadot * (cos(prev_theta) - cos(temp_angle)));
		//particles[i].theta = prev_theta + (0.0001*delta_t);


		}
		else{
		particles[i].x = prev_x + (vbythetadot * (sin(temp_angle) - sin(prev_theta)));
		particles[i].y = prev_y + (vbythetadot * (cos(prev_theta) - cos(temp_angle)));
		particles[i].theta = temp_angle;


		}		
		// generate gaussian distributions for each dimension
		std::normal_distribution<double> dist_x(particles[i].x, std_x);
		std::normal_distribution<double> dist_y(particles[i].y, std_y);
		std::normal_distribution<double> dist_theta(particles[i].theta, std_theta);

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

	// Initialize
	std::vector<LandmarkObs> tmpobservations = observations;
	int numobservations = tmpobservations.size();


	Map tmpmap = map_landmarks;
	int numlandmarks = tmpmap.landmark_list.size();

	// First Step a
	//get number of observations
	// loop on the observations per each particle 
	// transform each observation coordinate to the global coordinates
	// output a vector in each particle with landmark observations according to the global coordinates assigned to each particle

	for (int i = 0; i < num_particles; ++i) {
		for(int j = 0; j < numobservations; ++j) {
			particles[i].sense_x[j] = particles[i].x + (cos(particles[i].theta) * tmpobservations[j].x) - (sin(particles[i].theta) * tmpobservations[j].y);
			particles[i].sense_y[j] = particles[i].y + (sin(particles[i].theta) * tmpobservations[j].x) + (cos(particles[i].theta) * tmpobservations[j].y);
		}
	}

	// First Step b 
	// Filter out the landmarks in sensor range	
	/*Map filtered_out_Map = tmpmap;
	for (int i = 0; i < num_particles; ++i) {
		for(int j = 0; j < numobservations; ++j) {
			double obs_x, obs_y;
			obs_x = particles[i].sense_x[j];
			obs_y = particles[i].sense_y[j];
			difference_x = obs_x - sensor_range;
			difference_y = obs_y - sensor_range;
			 

		}
	} 
	*/

	
	// Second Step
	// for each particle measure the distance between the observation measurements and the landmarks positions from the map
	// find the lowest distance and assign each observation measurement with the associated landmark id

	double distance [numlandmarks] = {};
	

	for (int i = 0; i < num_particles; ++i) {
		for(int j = 0; j < numobservations; ++j) {
			double obs_x, obs_y;
			obs_x = particles[i].sense_x[j];
			obs_y = particles[i].sense_y[j];

			for(int k = 0; k < numlandmarks; ++k) {	
				double landmark_x, landmark_y;
				landmark_x = (double) tmpmap.landmark_list[k].x_f;
				landmark_y = (double) tmpmap.landmark_list[k].y_f;
				
				distance[k] = dist(obs_x, obs_y, landmark_x, landmark_y);
			}
			// get the smallest distance
			double smallest_distance = distance[0];
			int index = 0;
			for (int k = 1; k < numlandmarks; ++k) {
				if (distance[k] < smallest_distance) {
					smallest_distance = distance[k];
					index = k;
				}
			}
			// associate the measurement with the nearest neighbor index
			particles[i].associations[j] = index;
		}
	
	}

	
	// Third Step  Update weights
	// For each particle loop on the observations and their associated landmarks 
	// Calculate for each observation the multivariate gaussian probability 
	// the final weight is the multiplication of all the observation outputs
	double sigma_x = std_landmark[0];
	double sigma_y = std_landmark[1];
	
	double gaussian_norm = 1 / (2 * M_PI * sigma_x * sigma_y);
	
	for (int i = 0; i < num_particles; ++i) {
		for(int j = 0; j < numobservations; ++j) {
			double mu_x, mu_y;
			double exponent; 
			double obs_x, obs_y;
			double expo1, expo2; 
			int index;
			
			obs_x = particles[i].sense_x[j];
			obs_y = particles[i].sense_y[j];
			
			mu_x = tmpmap.landmark_list[particles[i].associations[j]].x_f;
			mu_y = tmpmap.landmark_list[particles[i].associations[j]].y_f;

			index = particles[i].associations[j];
			particles[i].associations[j] = tmpmap.landmark_list[index].id_i;
			
			expo1 =  (pow((obs_x - mu_x) / (sigma_x), 2)) / 2;
			expo2 = (pow((obs_y - mu_y) / (sigma_y), 2)) / 2;
			
			exponent = expo1 + expo2;
			particles[i].weight *= gaussian_norm * exp(-1 * exponent);
		}	
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> new_particles;
	new_particles.resize(num_particles);
	std::default_random_engine gen;

	vector<double> current_weights;
	
	for (int j = 0; j < num_particles; j++)
	{
		current_weights.push_back(particles[j].weight);
	}
	

	discrete_distribution<> sample(current_weights.begin(), current_weights.end());

	for (int i = 0; i < num_particles; i++)
	{
		//new_particles.push_back(particles[sample(random_generator)]);
		new_particles[i] = particles[sample(gen)];
	}
	
	particles = new_particles;
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
