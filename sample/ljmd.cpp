#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#include "user-defined.hpp"

#include <cassert>



template<class System>
void initialPlace(System &sys, const PS::S32 num_par_axis, const PS::F64 initDistance){
  PS::S32 nx,ny,nz,jx,jy,jz,n=0;
  PS::MTTS mt;
  mt.init_genrand(0);
  nx=num_par_axis; 
  ny=num_par_axis; 
  nz=num_par_axis; 
  for (jx=1;jx<=nx;jx++){
    for (jy=1;jy<=ny;jy++){
      for (jz=1;jz<=nz;jz++){
	sys[n].pos.x=0.0+(jx-1)*initDistance;
	sys[n].pos.y=0.0+(jy-1)*initDistance;
	sys[n].pos.z=0.0+(jz-1)*initDistance;
	sys[n].vel.x=mt.genrand_res53() * 2.0 - 1.0;
	sys[n].vel.y=mt.genrand_res53() * 2.0 - 1.0;
	sys[n].vel.z=mt.genrand_res53() * 2.0 - 1.0;
	n+=1;
	sys[n].pos.x=0.0+(jx-1)*initDistance;
	sys[n].pos.y=initDistance/2.0+(jy-1)*initDistance;
	sys[n].pos.z=initDistance/2.0+(jz-1)*initDistance;
	sys[n].vel.x=mt.genrand_res53() * 2.0 - 1.0;
	sys[n].vel.y=mt.genrand_res53() * 2.0 - 1.0;
	sys[n].vel.z=mt.genrand_res53() * 2.0 - 1.0;
	n+=1;
	sys[n].pos.x=initDistance/2.0+(jx-1)*initDistance;
	sys[n].pos.y=initDistance/2.0+(jy-1)*initDistance;
	sys[n].pos.z=0.0+(jz-1)*initDistance;
	sys[n].vel.x=mt.genrand_res53() * 2.0 - 1.0;
	sys[n].vel.y=mt.genrand_res53() * 2.0 - 1.0;
	sys[n].vel.z=mt.genrand_res53() * 2.0 - 1.0;
	n+=1;
	sys[n].pos.x=initDistance/2.0+(jx-1)*initDistance;
	sys[n].pos.y=0.0+(jy-1)*initDistance;
	sys[n].pos.z=initDistance/2.0+(jz-1)*initDistance;
	sys[n].vel.x=mt.genrand_res53() * 2.0 - 1.0;
	sys[n].vel.y=mt.genrand_res53() * 2.0 - 1.0;
	sys[n].vel.z=mt.genrand_res53() * 2.0 - 1.0;
	n+=1;
      }
    }
  }
}



template<class System>
void
make_conf(System &sys, const PS::F64 L) {
  PS::MTTS mt;
  mt.init_genrand(0);
  //std::mt19937 mt(1);
  //std::uniform_real_distribution<double> ud(0.0, 1.0);
  const double V0 = 1.0;
  const int il = static_cast<int>(L);
  PS::S32 n_total = L * L * L;
  sys.setNumberOfParticleLocal(n_total);
  int n = 0;
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      for (int iz = 0; iz < L; iz++) {
        double z = mt.genrand_res53() * 2.0 - 1.0;
        double phi = mt.genrand_res53() * M_PI;
        sys[n].vel.x = V0 * sqrt(1 - z * z) * cos(phi);
        sys[n].vel.y = V0 * sqrt(1 - z * z) * sin(phi);
        sys[n].vel.z = V0 * z;
        sys[n].pos.x = ix + 0.5;
        sys[n].pos.y = iy + 0.5;
        sys[n].pos.z = iz + 0.5;
        n++;
      }
    }
  }
}


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
  PS::F64 dt            = 0.0050;
  PS::F32 dt_diag       = 1.0 / 8.0;
  PS::F32 dt_snap       = 1.0;
  char    dir_name[1024];
  PS::S32 c;
  PS::S32 num_par_axis  =5;
  PS::S64 n_tot         =4*num_par_axis*num_par_axis*num_par_axis;
  PS::F64 rho           =0.3;
  PS::F64 initDistance  =pow(4.0/rho,1.0/3.0); 
  PS::S32 nstep         =1000;

  PS::ParticleSystem<FPLj> system_lj;
  system_lj.initialize();
  PS::S32 n_loc    = 0;
  PS::F32 time_sys = 0.0;
  sprintf(dir_name,"./result");

  if(PS::Comm::getRank() == 0) {
    //    setParticlesColdUniformSphere(system_lj, n_tot, n_loc);
    system_lj.setNumberOfParticleLocal(n_tot);
    initialPlace< PS::ParticleSystem<FPLj> >(system_lj, num_par_axis, initDistance);

  } else {
    system_lj.setNumberOfParticleLocal(n_loc);
  }
  //  return 0;
  PS::F64vec box_size;

  box_size = num_par_axis * initDistance;
  const PS::F32 coef_ema = 0.3;
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);

  dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
                         PS::F64vec(box_size.x,box_size.y,box_size.z));
  //  system_lj.adjustPositionIntoRootDomain(dinfo);
  //PeriodicBoundaryCondition(system_lj,dinfo);
  dinfo.decomposeDomainAll(system_lj);
  system_lj.exchangeParticle(dinfo);
  n_loc = system_lj.getNumberOfParticleLocal();

  PS::TreeForForceShort<ForceLj, EPLj, EPLj>::Gather tree_lj;
  for(int i=0;i<system_lj.getNumberOfParticleLocal();i++)
    system_lj[i].search_radius = CUTOFF_LENGTH;
  tree_lj.initialize(n_tot, theta, n_leaf_limit, n_group_limit);

  tree_lj.calcForceAllAndWriteBack(CalcLj<EPLj>,
				   system_lj,
				   dinfo);

  //output file                                                                                           
  std::ofstream fpo1("placs.xyz");

  PS::S64 n_loop = 0;
  PS::S32 tt;
  for(tt=0;tt<nstep;tt++){
    //  while(time_sys < time_end){
    if(n_loop % 1 == 0 ){
    //    if(n_loop ==  1 ){
      fpo1<<system_lj.getNumberOfParticleGlobal()<<std::endl;
      fpo1<<n_loop<<std::endl;
      for(PS::S32 i=0;i<system_lj.getNumberOfParticleGlobal();i++){
      fpo1<<"H"<<" "<<system_lj[i].pos.x<<" "<<system_lj[i].pos.y<<" "<<system_lj[i].pos.z<<std::endl;
      }
    }
    kick(system_lj, dt * 0.5);//why "*0.5" ?
    time_sys += dt;
    drift(system_lj, dt);
    //PeriodicBoundaryCondition(system_lj,dinfo);
    system_lj.adjustPositionIntoRootDomain(dinfo);
    if(n_loop % 4 == 0){
      dinfo.decomposeDomainAll(system_lj);
    }
    system_lj.exchangeParticle(dinfo);

    PS::F64ort boundary = dinfo.getPosRootDomain();
    for(int i=0;i<system_lj.getNumberOfParticleLocal();i++){
      assert(boundary.low_.x <= system_lj[i].pos.x);
      assert(boundary.low_.y <= system_lj[i].pos.y);
      assert(boundary.low_.z <= system_lj[i].pos.z);

      assert(boundary.high_.x > system_lj[i].pos.x);
      assert(boundary.high_.y > system_lj[i].pos.y);
      assert(boundary.high_.z > system_lj[i].pos.z);
    }
    tree_lj.calcForceAllAndWriteBack(CalcLj<EPLj>,
				     system_lj,
				     dinfo);

    kick(system_lj, dt * 0.5);
    std::cout<<"time="<<tt<<std::endl;
    n_loop++;
  }

  



  PS::Finalize();
  return 0;

}
