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



class FPLj{
public:
  PS::S64    id;
  PS::F64    mass;
  PS::F64vec pos;
  PS::F64vec vel;
  //PS::F64vec prevforce;
  PS::F64vec acc;
  PS::F64    pot;

  static PS::F64 eps;
  static PS::F64 sigma;

  PS::F64vec getPos() const {
    return pos;
  }

  PS::F64 getCharge() const {
    return mass;
  }
  void copyFromFP(const FPLj & fp){
    mass = fp.mass;
    pos  = fp.pos;
  }
  
  void copyFromForce(const FPLj & force) {
    acc = force.acc;
    pot = force.pot;
  }
  
  void clear() {
    acc = 0.0;
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




template<class TParticleJ>
void CalcLj(const FPLj *ep_i,
	    const PS::S32 n_ip,
	    const TParticleJ *ep_j,
	    const PS::S32 n_jp,
	    FPLj * force){
  PS::F64 ce12 = 4 * FPLj::eps * pow(FPLj::sigma,12);
  PS::F64 ce06 = 4 * FPLj::eps * pow(FPLj::sigma,6);
  PS::F64 cf12 = ce12 * 12;
  PS::F64 cf06 = ce06 * 6;
  for(PS::S32 i = 0; i < n_ip; i++){
    PS::F64vec xi    = ep_i[i].getPos();
    PS::F64vec ai    = 0.0;
    PS::F64    poti  = 0.0;
    for(PS::S32 j = 0; j < n_jp; j++){
      PS::F64vec rij  = xi - ep_j[j].getPos();
      PS::F64    r2   = rij * rij;
      PS::F64    r2i  = 1 / r2;
      PS::F64    r06i = r2i * r2i * r2i;
      PS::F64    r12i = r06i * r06i;
      PS::F64    ep   = ce12 * r12i - ce06 * r06i;
      PS::F64    fc   = (cf12 * r12i - cf06 * r06i)*r2i;
      ai += fc * rij; 
      poti += ep; 
    }
    //force[i].prevforce = force[i].force;
    force[i].acc = ai;
    force[i].pot = poti;
  }

}
