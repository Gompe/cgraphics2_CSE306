#include "gx_random.h"

double random_uniform(){
    // Idea from
    // https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers
    
    static thread_local std::default_random_engine engine(RANDOM_SEED);
    std::uniform_real_distribution<double> uniform(0,1);
    return uniform(engine);
}

void boxMuller(double stdev, double &x, double &y) {
	double r1 = random_uniform();
	double r2 = random_uniform();

	x = sqrt(-2 * log(r1))*cos(2 * M_PI *r2)*stdev;
	y = sqrt(-2 * log(r1))*sin(2 * M_PI *r2)*stdev;
}