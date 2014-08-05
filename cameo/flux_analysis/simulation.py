# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import time
from functools import partial
import sympy
from cameo.util import TimeMachine
from cameo.exceptions import SolveError


def fba(model, objective=None, *args, **kwargs):
    """Perform flux balance analysis."""
    tm = TimeMachine()
    if objective is not None:
        tm(do=partial(setattr, model, 'objective', objective),
           undo=partial(setattr, model, 'objective', model.objective))
    try:
        solution = model.solve()
        tm.reset()
        return solution
    except SolveError as e:
        tm.reset()
        raise e


def pfba(model, objective=None, *args, **kwargs):
    tm = TimeMachine()
    tm(do=partial(setattr, model, 'reversible_encoding', 'split'),
       undo=partial(setattr, model, 'reversible_encoding', model.reversible_encoding))
    try:
        if objective is not None:
            tm(do=partial(setattr, model, 'objective', objective),
               undo=partial(setattr, model, 'objective', model.objective))
        try:
            obj_val = model.solve().f
        except SolveError as e:
            print "pfba could not determine maximum objective value for\n%s." % model.objective
            raise e
        if model.objective.direction == 'max':
            fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression, lb=obj_val)
        else:
            fix_obj_constraint = model.solver.interface.Constraint(model.objective.expression, ub=obj_val)
        tm(do=partial(model.solver._add_constraint, fix_obj_constraint),
           undo=partial(model.solver._remove_constraint, fix_obj_constraint))

        pfba_obj = model.solver.interface.Objective(sympy.Add._from_args(
            [sympy.Mul._from_args((sympy.singleton.S.One, variable)) for variable in model.solver.variables.values()]),
                                                    direction='min', sloppy=True)
        # tic = time.time()
        tm(do=partial(setattr, model, 'objective', pfba_obj),
           undo=partial(setattr, model, 'objective', model.objective))
        # print "obj: ", time.time() - tic
        try:
            solution = model.solve()
            tm.reset()
            return solution
        except SolveError as e:
            tm.reset()
            print "pfba could not determine an optimal solution for objective %s" % model.objective
            raise e
    except Exception as e:
        tm.reset()
        raise e


def moma(model, objective=None, *args, **kwargs):
    pass


def lmoma(model, wt_reference=None):
    tm = TimeMachine()

    obj_terms = list()
    for rid, flux_value in wt_reference.iteritems():
        of_term = model.solver.interface.Variable("u_%s" % rid)
        tm(do=partial(model.solver._add_variable, of_term),
           undo=partial(model.solver._remove_variable, of_term))
        obj_terms.append(of_term)
        constraint1 = (model.solver.interface.Constraint(model.solver.variables[of_term.name] - model.reactions.get_by_id(rid).variable,
                                                         lb=-flux_value))
        constraint2 = (model.solver.interface.Constraint(model.solver.variables[of_term.name] - model.reactions.get_by_id(rid).variable,
                                                         lb=flux_value))

        tm(do=partial(model.solver._add_constraint, constraint1),
           undo=partial(model.solver._remove_constraint, constraint1))

        tm(do=partial(model.solver._add_constraint, constraint2),
           undo=partial(model.solver._remove_constraint, constraint2))

    lmoma_obj = model.solver.interface.Objective(sympy.Add._from_args(obj_terms), direction='min')

    tm(do=partial(setattr, model, 'objective', lmoma_obj),
       undo=partial(setattr, model, 'objective', model.objective))

    try:
        solution = model.solve()

        tm.reset()
        return solution
    except SolveError as e:
        tm.reset()
        print "lmoma could not determine an optimal solution for objective %s" % model.objective
        print model.solver
        raise e


def _cycle_free_flux(model, fluxes, fix=[]):
    """Remove cycles from a flux-distribution (http://cran.r-project.org/web/packages/sybilcycleFreeFlux/index.html)."""
    tm = TimeMachine()
    exchange_reactions = model.exchanges
    exchange_ids = [exchange.id for exchange in exchange_reactions]
    internal_reactions = [reaction for reaction in model.reactions if reaction.id not in exchange_ids]
    for exchange in exchange_reactions:
        exchange_flux = fluxes[exchange.id]
        tm(do=partial(setattr, exchange, 'lower_bound', exchange_flux),
           undo=partial(setattr, exchange, 'lower_bound', exchange.lower_bound))
        tm(do=partial(setattr, exchange, 'upper_bound', exchange_flux),
           undo=partial(setattr, exchange, 'upper_bound', exchange.upper_bound))
    obj_terms = list()
    for internal_reaction in internal_reactions:
        internal_flux = fluxes[internal_reaction.id]
        if internal_flux >= 0:
            obj_terms.append(sympy.Mul._from_args([sympy.S.One, internal_reaction.variable]))
            tm(do=partial(setattr, internal_reaction, 'lower_bound', 0),
               undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
            tm(do=partial(setattr, internal_reaction, 'upper_bound', internal_flux),
               undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
        elif internal_flux < 0:
            obj_terms.append(sympy.Mul._from_args([sympy.S.NegativeOne, internal_reaction.variable]))
            tm(do=partial(setattr, internal_reaction, 'lower_bound', internal_flux),
               undo=partial(setattr, internal_reaction, 'lower_bound', internal_reaction.lower_bound))
            tm(do=partial(setattr, internal_reaction, 'upper_bound', 0),
               undo=partial(setattr, internal_reaction, 'upper_bound', internal_reaction.upper_bound))
        else:
            pass
    for reaction_id in fix:
        reaction_to_fix = model.reactions.get_by_id(reaction_id)
        tm(do=partial(setattr, reaction_to_fix, 'lower_bound', fluxes[reaction_id]),
           undo=partial(setattr, reaction_to_fix, 'lower_bound', reaction_to_fix.lower_bound))
        tm(do=partial(setattr, reaction_to_fix, 'upper_bound', fluxes[reaction_id]),
           undo=partial(setattr, reaction_to_fix, 'upper_bound', reaction_to_fix.upper_bound))
    tm(do=partial(setattr, model, 'objective',
                  model.solver.interface.Objective(sympy.Add._from_args(obj_terms), name='Flux minimization',
                                                   direction='min', sloppy=True)),
       undo=partial(setattr, model, 'objective', model.objective))
    solution = model.optimize()
    tm.reset()
    return solution.x_dict


if __name__ == '__main__':
    import time
    from cobra.io import read_sbml_model
    from cobra.flux_analysis.parsimonious import optimize_minimal_flux
    from cameo import load_model

    sbml_path = '../../tests/data/EcoliCore.xml'
    #sbml_path = '../../tests/data/iJO1366.xml'

    cb_model = read_sbml_model(sbml_path)
    model = load_model(sbml_path)

    # model.solver = 'glpk'

    print "cobra fba"
    tic = time.time()
    cb_model.optimize(solver='cglpk')
    print "flux sum:", sum([abs(val) for val in cb_model.solution.x_dict.values()])
    print "cobra fba runtime:", time.time() - tic

    print "cobra pfba"
    tic = time.time()
    optimize_minimal_flux(cb_model, solver='cglpk')
    print "flux sum:", sum([abs(val) for val in cb_model.solution.x_dict.values()])
    print "cobra pfba runtime:", time.time() - tic

    print "pfba"
    tic = time.time()
    solution = pfba(model)
    print "flux sum:",
    print sum([abs(val) for val in solution.x_dict.values()])
    print "cameo pfba runtime:", time.time() - tic

    print "lmoma"
    tic = time.time()
    solution = pfba(model)
    solution = lmoma(model, wt_reference=solution.x_dict)
    print "flux sum:",
    print sum([abs(val) for val in solution.x_dict.values()])
    print "cameo lmoma runtime:", time.time() - tic

    print "lmoma w/ ko"
    tic = time.time()
    solution = pfba(model)
    model.reactions.FBA.lower_bound = 0
    model.reactions.FBA.upper_bound = 0
    solution = lmoma(model, wt_reference=solution.x_dict)
    print "flux sum:",
    print sum([abs(val) for val in solution.x_dict.values()])
    print "cameo lmoma runtime:", time.time() - tic

    # print model.solver