#include <pdb.hpp>	
#include <grid.hpp>
#include <string>

using namespace std;

void Grid::exhaustive_search(double z_min, double z_max, int t){  

	for(double x = 0; x < 2*M_PI; x+= 0.0174532925){
		for(double y = 0; y < 2*M_PI; y+= 0.0174532925){
			for(double z = z_min; z < z_max; z+= 0.25){
				double lowest_e = protein->orientate(x,y,z);
				boost::mutex::scoped_lock lock(result_mutex);	
				ncalls++;
				if(lowest_e < prevbest){
					prevbest = lowest_e;
				    protein->set_optparams(x,y,z);
				    if(threads > 1){
				    	printf("Lowest energy %f after %d function evaluations (thread %d).\n", lowest_e, ncalls, t+1);				    	
				    }else{	
						printf("Lowest energy %f after %d function evaluations.\n", lowest_e, ncalls);
					}
				}	
			}
		}
	}
} 

int Grid::run(){
	// Use boost threads
	if(threads){		
		for(int i = 0; i < threads; i++){
			double z_min = (i*z_range)-15.0;
			double z_max = ((i+1)*z_range)-15.0;
			// thread_group should take care of destroying these			
			boost::thread *t1 = new boost::thread(boost::bind(&Grid::exhaustive_search,this,z_min,z_max,i));
			g.add_thread(t1);
		}
		g.join_all();

	}else{
		exhaustive_search(2*M_PI,2*M_PI,1);		
	}
	return(ncalls);
}