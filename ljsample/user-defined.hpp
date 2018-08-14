#pragma once

class FileHeader{
public:
  PS::S64 n_body;
  PS::F64 time;
  PS::S32 readAscii(FILE * fp) {
    fscanf(fp, "%lf\n", &time);
    fscanf(fp, "%lld\n", &n_body);
    return n_body;
  }
  void writeAscii(FILE* fp) const {
    // fprintf(fp, "%e\n", time);
    fprintf(fp, "%lld\n", n_body);
    fprintf(fp, "%e\n", time);
  }
};

class ForceLj {
public:
  PS::F64vec force;
  PS::F64    pot;
  void clear(){
    force = 0.0;
    pot = 0.0;
  }
};

const PS::F64 CUTOFF_LENGTH = 3.0;
class FPLj{
public:
  PS::S64    id;
  static  PS::F64    mass;//一時的
  PS::F64vec pos;
  PS::F64vec vel;
  PS::F64vec force;
  PS::F64    pot;
  PS::F64    search_radius;

  // static PS::F64 rcut;
  static PS::F64 eps;
  static PS::F64 sigma;

  PS::F64 getRSearch() const {
    return CUTOFF_LENGTH;
  }

  PS::F64vec getPos() const {
    return pos;
  }

  PS::F64 getCharge() const {
    return mass;
  }
  
  void copyFromFP(const FPLj & fp){
    mass = fp.mass;
    pos  = fp.pos;
    search_radius = fp.search_radius;
  }
  
  void copyFromForce(const ForceLj & f) {
    force = f.force;
    pot = f.pot;
  }

  //"https://qiita.com/kaityo256/items/36d267be911ac45e848e"を参考
  void setPos(PS::F64vec npos){
    pos = npos;
  }

  void clear() {
    force = 0.0;
    pot = 0.0;
  }

  void writeAscii(FILE* fp) const {
    fprintf(fp, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	    this->id, this->mass,
	    this->pos.x, this->pos.y, this->pos.z,
	    this->vel.x, this->vel.y, this->vel.z);
  }

  void readAscii(FILE* fp) {
    fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
	   &this->id, &this->mass,
	   &this->pos.x, &this->pos.y, &this->pos.z,
	   &this->vel.x, &this->vel.y, &this->vel.z);
  }

  
};




class EPLj {
public:
  PS::S32 id;
  PS::F64vec pos;
  PS::F64 eps;
  PS::F64 sigma;
  PS::F64 search_radius;

  PS::F64vec getPos() const {
    return pos;
  }

  void setPos(const PS::F64vec & npos){
    pos = npos;
  }

  PS::F64 getRSearch() const {
    return CUTOFF_LENGTH;
  }

  void copyFromFP(const FPLj &fp){
    id    = fp.id;
    pos   = fp.pos;
    sigma = fp.sigma;
    eps   = fp.eps;
  }
};



template<class TParticleJ>
void CalcLj(const TParticleJ *ep_i,
	    const PS::S32 n_ip,
	    const TParticleJ *ep_j,
	    const PS::S32 n_jp,
	    ForceLj * force){
  PS::F64 ce12 = 4 * FPLj::eps * pow(FPLj::sigma,12);
  PS::F64 ce06 = 4 * FPLj::eps * pow(FPLj::sigma,6);
  PS::F64 cf12 = ce12 * 12;
  PS::F64 cf06 = ce06 * 6;
  PS::F64 rcut_sq = pow(CUTOFF_LENGTH,2.0);
  for(PS::S32 i = 0; i < n_ip; i++){
    PS::F64vec xi    = ep_i[i].pos;
    PS::F64vec forcei    = 0.0;
    PS::F64    poti  = 0.0;
    for(PS::S32 j = 0; j < n_jp; j++){
      PS::F64vec rij  = xi - ep_j[j].pos;
      PS::F64    r2   = rij * rij;
      if(0.0 < r2 && r2<=rcut_sq){
	PS::F64    r2i  = 1 / r2;
	PS::F64    r06i = r2i * r2i * r2i;
	PS::F64    r12i = r06i * r06i;
	PS::F64    ep   = ce12 * r12i - ce06 * r06i;
	//	PS::F64    fc   = - (cf12 * r12i - cf06 * r06i)*r2i;
	PS::F64    fc   =  (cf12 * r12i - cf06 * r06i)*r2i;
	forcei += fc * rij; 
	poti += ep * 0.5; 
      }
    }
    //force[i].prevforce = force[i].force;
    force[i].force = forcei;
    force[i].pot = poti;
  }

};


