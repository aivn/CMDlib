#include <sstream>
#include <aiwlib/binhead>
#include "full_model.hpp"

using namespace aiw;

//------------------------------------------------------------------------------
void Model::init(std::string path){
	for(int st=0; st<4; st++) data[st].resize(mesh_sz.prod());
	for(cell_t& c: data[0]){  c.m = M0; c.eta = eta0; }
	_dr = vecf(1)/dr;
	offsets.clear(); offsets.reserve(data[0].size()*6);
	for(int k=0; k<mesh_sz[2]; k++)
		for(int j=0; j<mesh_sz[1]; j++)
			for(int i=0; i<mesh_sz[0]; i++){
				offsets.push_back(i==0? (periodic&1? mesh_sz[0]-1: 0): -1);
				offsets.push_back(i==mesh_sz[0]-1? (periodic&1? -(mesh_sz[0]-1): 0): 1);
				offsets.push_back(j==0? (periodic&2? mesh_sz[0]*(mesh_sz[1]-1): 0): -mesh_sz[0]);
				offsets.push_back(j==mesh_sz[1]-1? (periodic&2? -mesh_sz[0]*(mesh_sz[1]-1):0): mesh_sz[0]);
				offsets.push_back(k==0? (periodic&4? mesh_sz[0]*mesh_sz[1]*(mesh_sz[2]-1): 0): -mesh_sz[0]*mesh_sz[1]);
				offsets.push_back(k==mesh_sz[2]-1? (periodic&4? -mesh_sz[0]*mesh_sz[1]*(mesh_sz[2]-1):0): mesh_sz[0]*mesh_sz[1]);

				// std::cout<<i<<' '<<j<<' '<<k<<':'; for(int l=0; l<6; l++){ std::cout<<offsets[6*(i+(j+k*mesh_sz[1])*mesh_sz[0])+l]<<','; } std::cout<<'\n';
			}
	M = M0; eta = eta0; t = 0; W[1] = -J*eta*n_b/2; W[2] = -Hext*M; W[3] = -K*(1-2*CMD::calc_Mu(nK*M)); W[0] = W[1]+W[2]+W[3];

	ftvals.open(path+"/tvals.dat"); ftvals<<"#:t M Mx My Mz eta W Wexch Wext Waniso\n"<<t<<' '<<M.abs()<<' '<<M<<' '<<eta<<' '<<W<<std::endl;
	fdata.open(path+"/data.msh");

	coeff_cmd.n_b = n_b; coeff_cmd.eta_c = eta_c;
	for(int i=0; i<3; i++) coeff_llb.nK[i] = coeff_cmd.nK[i] = nK[i];
}
//------------------------------------------------------------------------------
Model::cell_t Model::calc_dcdt(int st, int i){
	const int* offs = &(offsets[i*6]); const cell_t &c0 = data[st][i];
	Vecf<3> H = Hext;
	if(offs[0] && offs[1]) H += (data[st][i+offs[0]].m + data[st][i+offs[1]].m - 2*c0.m)*_dr[0]; 
	if(offs[2] && offs[3]) H += (data[st][i+offs[2]].m + data[st][i+offs[3]].m - 2*c0.m)*_dr[1]; 
	if(offs[4] && offs[5]) H += (data[st][i+offs[4]].m + data[st][i+offs[5]].m - 2*c0.m)*_dr[2];

	cell_t res;  res.m = -gamma*(c0.m%H + 2*alpha*T*c0.m);  res.eta = 0;
	if(use_cmd){
		coeff_cmd.eta = c0.eta; for(int k=0; k<3; k++){ coeff_cmd.M[k] = c0.m[k]; coeff_cmd.H[k] = H[k]; }
		coeff_cmd.calc();
		// M   += 2*adt*(K*coeff.Theta + .5*coeff.Xi*H + (n_b*coeff.U - T)*coeff.M);
		for(int k=0; k<3; k++) res.m[k] += gamma*( 2*K*(coeff_cmd.PHI[k]+alpha*coeff_cmd.THETA[k]) + alpha*coeff_cmd.XiH[k] + 2*n_b*alpha*J*coeff_cmd.U*c0.m[k] );
		// eta += 4*adt*(coeff.U*H*coeff.M*1.333 + K*coeff.Psi + coeff.Q - T*coeff.eta);
		res.eta = 4*gamma*alpha*(H*c0.m*coeff_cmd.U + K*coeff_cmd.Psi + J*coeff_cmd.Q - T*c0.eta);
	} else {
		for(int k=0; k<3; k++){ coeff_llb.M[k] = c0.m[k]; coeff_llb.H[k] = H[k] + n_b*eG*J*c0.m[k]; }
		coeff_llb.calc();
		// coeff.M += 2*adt*(K*coeff.Theta + .5*coeff.Xi*(H + n_b*eG*M) - T*M);
		for(int k=0; k<3; k++) res.m[k] += gamma*( 2*K*(coeff_llb.PHI[k]+alpha*coeff_llb.THETA[k]) + alpha*coeff_llb.XiH[k] );
		// for(int k=0; k<3; k++) res.m[k] += gamma*( 2*K*(coeff_llb.PHI[k]+alpha*coeff_llb.THETA[k]) + alpha*coeff_llb.XiH[k] );
	}
	return res;
}
//------------------------------------------------------------------------------
const float _3 = 1./3;
void Model::calc(int steps){
	float sghT = sqrt(2*dt*alpha*T/cell_sz);

	for(int it=0; it<steps; it++){
		M = vecf(0.); W = vec(0.); eta = 0;
		//----------------------------------------------------------------------
		for(int i=0, sz=data[0].size(); i<sz; ++i){
			cell_t &c0 = data[0][i], &c1 = data[1][i], &dc = data[3][i];
			cell_t dcdt = calc_dcdt(0, i);
			c1 = c0 + dcdt*(.5f*dt);
			dc = dcdt*.5f;
		}
		//----------------------------------------------------------------------
		for(int i=0, sz=data[0].size(); i<sz; ++i){
			cell_t &c0 = data[0][i], &c2 = data[2][i], &dc = data[3][i];
			cell_t dcdt = calc_dcdt(1, i);
			c2 = c0 + dcdt*(.5f*dt);
			dc += dcdt;
		}
		//----------------------------------------------------------------------
		for(int i=0, sz=data[0].size(); i<sz; ++i){
			cell_t &c0 = data[0][i], &c1 = data[1][i], &dc = data[3][i];
			cell_t dcdt = calc_dcdt(2, i);
			c1 = c0 + dcdt*dt;
			dc += dcdt;
		}
		//----------------------------------------------------------------------
		for(int i=0, sz=data[0].size(); i<sz; ++i){
			cell_t &c0 = data[0][i], &dc = data[3][i];
			cell_t dcdt = calc_dcdt(1, i);
			c0 += (dc+dcdt*.5f)*(dt*_3);
			c0.m = rotate(c0.m, randN01.V<3>()*sghT);
			if(!use_cmd) c0.eta = c0.m*c0.m;
			M += c0.m; W[3] += -K*(1-2*CMD::calc_Mu(nK*c0.m));  eta += data[0][i].eta;
		}
		//----------------------------------------------------------------------
		t += dt;
		M /= data[0].size();
		eta /= data[0].size();
		W[3] /= data[0].size();
		W[1] = -J*eta*n_b/2; W[2] = -Hext*M; W[0] = W[1]+W[2]+W[3];
		ftvals<<t<<' '<<M.abs()<<' '<<M<<' '<<eta<<' '<<W<<std::endl;
	}
}
//------------------------------------------------------------------------------
void Model::dump(){
	BinaryHead bh; bh.dim = 3; bh.szT = 16; bh.type = BinaryHead::mesh;
	for(int i=0; i<3; i++){ bh.bmin[i] = 0; bh.bmax[i] = mesh_sz[i]; bh.bbox[i] = mesh_sz[i]; }
	std::stringstream S; S<<"t="<<t; bh.head = S.str();
	bh.dump(fdata); fdata.write((const char*)(data[0].data()), data[0].size()*16);
}
//------------------------------------------------------------------------------


