#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#include "user-defined.hpp"



void makeColdUniformSphere(const PS::F64 mass_glb,
                           const PS::S64 n_glb,
                           const PS::S64 n_loc,
                           PS::F64 *& mass,
                           PS::F64vec *& pos,
                           PS::F64vec *& vel,
                           const PS::F64 eng = -0.25,
                           const PS::S32 seed = 0) {

  assert(eng < 0.0);
  {
    PS::MTTS mt;
    mt.init_genrand(0);
    for(PS::S32 i = 0; i < n_loc; i++){
      mass[i] = mass_glb / n_glb;
      const PS::F64 radius = 3.0;
      do {
	pos[i][0] = (2. * mt.genrand_res53() - 1.) * radius;
	pos[i][1] = (2. * mt.genrand_res53() - 1.) * radius;
	pos[i][2] = (2. * mt.genrand_res53() - 1.) * radius;
      }while(pos[i] * pos[i] >= radius * radius);
      vel[i][0] = 0.0;
      vel[i][1] = 0.0;
      vel[i][2] = 0.0;
    }
  }

  PS::F64vec cm_pos  = 0.0;
  PS::F64vec cm_vel  = 0.0;
  PS::F64    cm_mass = 0.0;
  for(PS::S32 i = 0; i < n_loc; i++){
    cm_pos  += mass[i] * pos[i];
    cm_vel  += mass[i] * vel[i];
    cm_mass += mass[i];
  }
  cm_pos /= cm_mass;
  cm_vel /= cm_mass;
  for(PS::S32 i = 0; i < n_loc; i++){
    pos[i] -= cm_pos;
    vel[i] -= cm_vel;
  }
}

template<class Tpsys>
void setParticlesColdUniformSphere(Tpsys & psys,
                                   const PS::S32 n_glb,
                                   PS::S32 & n_loc) {

  n_loc = n_glb;
  psys.setNumberOfParticleLocal(n_loc);

  PS::F64    * mass = new PS::F64[n_loc];
  PS::F64vec * pos  = new PS::F64vec[n_loc];
  PS::F64vec * vel  = new PS::F64vec[n_loc];
  const PS::F64 m_tot = 1.0;
  const PS::F64 eng   = -0.25;
  makeColdUniformSphere(m_tot, n_glb, n_loc, mass, pos, vel, eng);
  for(PS::S32 i = 0; i < n_loc; i++){
    psys[i].mass = mass[i];
    psys[i].pos  = pos[i];
    psys[i].vel  = vel[i];
    psys[i].id   = i;
  }
  delete [] mass;
  delete [] pos;
  delete [] vel;
}

template<class Tpsys>
void kick(Tpsys & system,
	  const PS::F64 dt) {
  PS::S32 n = system.getNumberOfParticleLocal();
  for(PS::S32 i = 0; i < n; i++) {
    //    system[i].vel += 0.50 * (system[i].force + system[i].prevforce) * dt / system[i].mass;
    system[i].vel  += system[i].force * dt;
  }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
  PS::S32 n = system.getNumberOfParticleLocal();
  for(PS::S32 i = 0; i < n; i++) {
    //    system[i].pos += system[i].vel * dt + 0.50 * dt * dt* system[i].prevforce / system[i].mass;
    system[i].pos  += system[i].vel * dt;
  }
}


template<class Tpsys,class Tdinfo>
void PeriodicBoundaryCondition(Tpsys & system,const Tdinfo &dinfo){
  PS::S32 n = system.getNumberOfParticleLocal();
  const PS::F64ort domain = dinfo.getPosRootDomain();
  const PS::F64vec length = domain.getFullLength();

  for(int i=0; i<n; i++){
    for(int k=0;k<3;k++){
      if(system[i].pos.x <  domain.low_.x){
        system[i].pos.x  += length.x;
      }
      if(system[i].pos.x >= domain.high_.x){
        system[i].pos.x  -= length.x;
      }
      if(system[i].pos.y <  domain.low_.y){
        system[i].pos.y  += length.y;
      }
      if(system[i].pos.y >= domain.high_.y){
        system[i].pos.y  -= length.y;
      }
      if(system[i].pos.z <  domain.low_.z){
        system[i].pos.z  += length.z;
      }
      if(system[i].pos.z >= domain.high_.z){
        system[i].pos.z  -= length.z;
      }
    }
  }
}



//PS::F64 FPLj::rcut  = 1.0;
PS::F64 FPLj::eps   = 1.0;
PS::F64 FPLj::sigma = 1.0;


int main(int argc, char *argv[]) {


  PS::Initialize(argc, argv);
  PS::F32 theta         = 0.5;
  PS::S32 n_leaf_limit  = 8;
  PS::S32 n_group_limit = 64;
  PS::F32 time_end      = 10.0;
  PS::F32 dt            = 1.0 / 128.0;
  PS::F32 dt_diag       = 1.0 / 8.0;
  PS::F32 dt_snap       = 1.0;
  char    dir_name[1024];
  PS::S64 n_tot         = 1024;
  PS::S32 c;

  PS::ParticleSystem<FPLj> system_lj;
  system_lj.initialize();
  PS::S32 n_loc    = 0;
  PS::F32 time_sys = 0.0;
  sprintf(dir_name,"./result");

  if(PS::Comm::getRank() == 0) {
    setParticlesColdUniformSphere(system_lj, n_tot, n_loc);
  } else {
    system_lj.setNumberOfParticleLocal(n_loc);
  }


  PS::F64vec box_size;
  box_size = 4 * powf(4.0/3.0,1.0/3.0);
  const PS::F32 coef_ema = 0.3;
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
                         PS::F64vec(box_size.x,box_size.y,box_size.z));
  PeriodicBoundaryCondition(system_lj,dinfo);
  dinfo.decomposeDomainAll(system_lj);
  system_lj.exchangeParticle(dinfo);
  n_loc = system_lj.getNumberOfParticleLocal();

  PS::TreeForForceShort<FPLj, FPLj, FPLj>::Gather tree_lj;
  for(int i=0;i<system_lj.getNumberOfParticleLocal();i++)
    system_lj[i].search_radius = CUTOFF_LENGTH;
  tree_lj.initialize(n_tot, theta, n_leaf_limit, n_group_limit);

  tree_lj.calcForceAllAndWriteBack(CalcLj<FPLj>,
				   system_lj,
				   dinfo);

  //output file                                                                                           
  std::ofstream fpo1("placs.xyz");

  PS::S64 n_loop = 0;

  while(time_sys < time_end){



    if(n_loop % 100 == 0 ){

      fpo1<<system_lj.getNumberOfParticleGlobal()<<std::endl;
      fpo1<<n_loop<<std::endl;
      for(PS::S32 i=0;i<system_lj.getNumberOfParticleGlobal();i++){
      fpo1<<"H"<<" "<<system_lj[i].pos.x<<" "<<system_lj[i].pos.y<<" "<<system_lj[i].pos.z<<std::endl;
      }
      
    }




    kick(system_lj, dt * 0.5);//why "*0.5" ?

    time_sys += dt;
    drift(system_lj, dt);
    PeriodicBoundaryCondition(system_lj,dinfo);
    system_lj.adjustPositionIntoRootDomain(dinfo);

    if(n_loop % 4 == 0){
      dinfo.decomposeDomainAll(system_lj);
    }

    system_lj.exchangeParticle(dinfo);

    tree_lj.calcForceAllAndWriteBack(CalcLj<FPLj>,
				     system_lj,
				     dinfo);


    kick(system_lj, dt * 0.5);
    
    n_loop++;
  }

  



  PS::Finalize();
  return 0;

}
