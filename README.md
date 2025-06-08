This is a simple example of the typical N-Body simulation using C. This is meant to show how such a simulation could be paralelized, using OpenMP, OpenACC, CUDA... How to get the compiler to use SIMD...  
This is a personal project to learn and apply HPC stuff. My plan is to first make a simple and direct version (though that does not mean bad, i'll still try to get a good speedup) using tools like the ones mentioned before and then move onto the complex Barnes-hut implementation.  
To keep track of what is done and what is to be implemented here is a TODO list:  

# Direct version:  
    - [x] Vectorised
    - [x] OpenMP
    - [ ] OpenMPI  
    - [x] OpenACC  
    - [ ] CUDA  
  
# Barnes-Hut:  
    - [x] Vectorised  
    - [x] OpenMP  
    - [ ] OpenMPI  
    - [ ] OpenACC  
    - [ ] CUDA  
  
If for some reason you stubled onto this repo and want to learn how it works there is a presentation (in spanish). Though is consists basically of images and colors, so it should be easy enough to understand  
