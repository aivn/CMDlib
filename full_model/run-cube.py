#!/usr/bin/python2
# -*- coding: utf-8 -*-

import os, math
from full_model import *
from aiwlib.vec import *
from aiwlib.racs import *

calc = Calc(t_max=(50., '#максимальное время расчета'), data_rank=(6, '# размер задачи в периодах кристаллической решетки'), sz=(4, '# размер расчетной сетки'),
            steps=(100, '#число шагов между сбросом данных'), _repo='repo', Hz=(0., '#внешнее поле (для задания из командной строки)'),
            M0z=(.99, '#начальная намагниченность (для задания из командной строки)'), mode=('SC', '#тип кристаллической решетки (SC|BCC|FCC)'))
calc.tags.add(calc.mode)

model = calc.wrap(Model(), '', lambda: model.t/calc.t_max)
model.M0[2] = calc.M0z
model.eta0 = calc.M0z**2
model.Hext[2] = calc.Hz

model.n_b, model.eG, model.eta_c, Tc, cell_sz = {'SC':(6, .725, .33,  1.445, 1), 'BCC':(8, .775, .27,  2.05, 2), 'FCC':(12, .7975, .25,  3.17, 4)}[calc.mode]

model.cell_sz = cell_sz*2**(3*calc.data_rank)/calc.sz
model.mesh_sz = ind(*[calc.sz]*3)
model.dr = vecf(*[2**calc.data_rank/calc.sz]*3)

if Tc-.8<=model.T and model.T<=Tc+1.4: calc.t_max *= 2
if Tc-.4<=model.T and model.T<=Tc+.7: calc.t_max *= 2
if Tc-.2<=model.T and model.T<=Tc+.3: calc.t_max *= 4

model.init(calc.path)

while model.t<calc.t_max:
    model.calc(calc.steps)    
    model.dump()
