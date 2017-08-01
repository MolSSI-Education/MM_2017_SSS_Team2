#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <iostream>
#include <pybind11/numpy.h>
#include <cmath>

namespace py = pybind11;

double one_particle (py::array_t <double> x_coords,
				    py::array_t <double> y_coords,
				    py::array_t <double> z_coords, 
				    double box_length, double cutoff2,
				    int i) 
{
  py::buffer_info x_coords_info = x_coords.request();
  py::buffer_info y_coords_info = y_coords.request();
  py::buffer_info z_coords_info = z_coords.request();

  int num_particles = x_coords_info.shape[0];

  // Maybe add checks to make sure arrays are same size

//  double system_energy(const std::vector<double> & xy,
   //                const std::vector<double> & v2)
//		   const std::vector<double> & v2)

 // {
//    if(v1.size() != v2.size())
  //      throw std::runtime_error("Vectors are different lengths");

    //if(v1.size() == 0)
      //  throw std::runtime_error("Zero-length vector");

    //double dot = 0.0;
    //for(size_t i = 0; i < v1.size(); i++)
      //  dot += v1[i] * v2[i];

   // return dot;
 // }

  const double * x_data = static_cast<double *>(x_coords_info.ptr);
  const double * y_data = static_cast<double *>(y_coords_info.ptr);
  const double * z_data = static_cast<double *>(z_coords_info.ptr);

  double i_x, i_y, i_z;
  double j_x, j_y, j_z; 
  double rijx, rijy, rijz, rij2;			
  double sig_by_r6;
  double sig_by_r12;
  double lj_pot, e_total; 

  e_total = 0.0;
  for (int i = 0; i < num_particles; i++)
  {
	i_x = x_data[3*i];
	i_y = y_data[3*i];
	i_z = z_data[3*i];

	for (int j = i + 1; j < num_particles
; j++)
	{
		j_x = x_data[3*j];
		j_y = y_data[3*j];
		j_z = z_data[3*j];

		rijx = j_x - i_x;
		rijy = j_y - i_y;
		rijz = j_z - i_z;

		rijx = rijx - box_length * round(rijx / box_length);
		rijy = rijy - box_length * round(rijy / box_length);
		rijz = rijz - box_length * round(rijz / box_length);

		rij2 = pow(rijx,2) + pow(rijy,2) + pow(rijz,2);

		if (rij2 < cutoff2)
		{
			sig_by_r6 = pow((1 / rij2),3);
    			sig_by_r12 = pow(sig_by_r6,2);
    			lj_pot = 4.0 * (sig_by_r12 - sig_by_r6);
			e_total += lj_pot;	
		}


	}
  }

  return e_total;
}

	int e_total = 0.0;
        for (int j = 1; j < num_particles; j++)
        {
                j_x = x_data[3*j];
                j_y = y_data[3*j];
                j_z = z_data[3*j];

                rijx = j_x - i_x;
                rijy = j_y - i_y;
                rijz = j_z - i_z;

                rijx = rijx - box_length * round(rijx / box_length);
                rijy = rijy - box_length * round(rijy / box_length);
                rijz = rijz - box_length * round(rijz / box_length);

                rij2 = pow(rijx,2) + pow(rijy,2) + pow(rijz,2);

                if (rij2 < cutoff2)
                {
                        sig_by_r6 = pow((1 / rij2),3);
                        sig_by_r12 = pow(sig_by_r6,2);
                        lj_pot = 4.0 * (sig_by_r12 - sig_by_r6);
                        e_total += lj_pot;
                }


        }
  


PYBIND11_PLUGIN(mc_demo)
{
   py::module m("mc_demo", "Dibyendu's basic module");
 
   m.def("one_particle", &one_particle);
   return m.ptr();

}
