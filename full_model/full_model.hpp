#include <fstream>
#include <vector>

#include <aiwlib/vec>
#include <aiwlib/gauss>

#include "../include/cmd_coeff.hpp"

class Model{
	struct cell_t{
		aiw::Vecf<3> m;
		float eta;
		
		void operator += (const cell_t& other){ m += other.m; eta += other.eta; }
		cell_t operator + (const cell_t& other) const { cell_t res; res.m = m+other.m; res.eta = eta+other.eta; return res; }
		cell_t operator * (float x) const { cell_t res; res.m = m*x; res.eta = eta*x; return res; }
	};
	cell_t calc_dcdt(int st, int i);
	std::vector<cell_t> data[4]; std::vector<int> offsets;
	aiw::Vecf<3> _dr;
	
	std::ofstream ftvals, fdata;
	aiw::RandN01<float> randN01;

	CMD::CoeffCMD<float> coeff_cmd;
	CMD::CoeffLLB<float> coeff_llb;
public:
	float J = 1.f;                          ///< обменный интеграл (межатомный!)
	int n_b;                                ///< число ближайших соседей
	float eG;                               ///< фактор Гаранина
	float eta_c;                            ///< значение eta в точке фазового перехода
	float gamma = 1.f;                      ///< гиромагнитное отношение
	float alpha = .1f;                      ///< коэффициент диссипации

	float dt = .01f;                        ///< шаг по времени
	double T = 1.;                          ///< температура
	double K = 0.;                          ///< анизотропия
	aiw::Vecf<3> Hext;                      ///< внешнее поле
	aiw::Vecf<3> nK = aiw::vecf(0.f, 0.f, 1.f);  ///< направление анизотропии
	aiw::Vecf<3> dr;                         ///< шаг сетки, в расстояниях между атомами
	aiw::Ind<3> mesh_sz;                    ///< размер сетки в ячейках
	int periodic = 7;                       ///< периодические граничные условия (битовая маска)
	int cell_sz = 0;                        ///< число атомов в ячейке сетки
	bool use_cmd = false;                   ///< CMD или LLB

	aiw::Vecf<3> M0;                        ///< начальная намагниченность
	double eta0;                            ///< начальная eta
	
	void init(std::string path);
	void calc(int steps);
	void dump();

	float t = 0;      ///< текущеее время
	aiw::Vecf<3> M;   ///< средняя намагниченность
	double eta;       ///< среднее eta
	aiw::Vec<4> W;    ///< энергия: W (полная), Wexch (обменная), Wext (вклад внешнего поля), Wanis (анизотропии)
};
