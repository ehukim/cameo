# Copyright 2016 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
from functools import partial

import cobra
import optlang
import six
from cobra import DictList, Object
from sympy.core.singleton import S

from cameo.core import SolverBasedModel, Gene, Reaction
from cameo.core.metabolite import Metabolite
from cameo.core.solution import LazySolution
from cameo.util import AutoVivification, flatten

NON_ENZYMATIC = "_non_enzymatic_"

logger = logging.getLogger(__name__)

__all__ = ["GPRBasedModel"]


class GPRReaction(Reaction):
    @classmethod
    def clone(cls, reaction, model=None):
        """Clone a reaction.

        Parameters
        ----------
        reaction : Reaction, cobra.core.Reaction.Reaction
        model : model, optional

        Returns
        -------
        Reaction

        """
        new_reaction = super(GPRReaction, cls).clone(reaction, model=model)
        new_reaction._forward_variables = {}
        new_reaction._reverse_variables = {}
        new_reaction._forward_isozyme_constraint = None
        new_reaction._reverse_isozyme_constraint = None
        return new_reaction

    def _clone_genes(self, model):
        cloned_genes = []
        for gene in self._genes:
            cloned_gene = GPRGene.clone(gene)
            cloned_genes.append(cloned_gene)
            if model is not None:
                model.genes._replace_on_id(cloned_gene)
        self._genes = set(cloned_genes)

    @property
    def forward_variables(self):
        if len(self._forward_variables) < len(self.isozymes):
            self._split_isozyme_variables()
        return self._forward_variables

    @property
    def reverse_variables(self):
        if len(self._reverse_variables) < len(self.isozymes):
            self._split_isozyme_variables()
        return self._reverse_variables

    @property
    def isozymes(self):
        if len(self.gene_reaction_rule) > 0:
            isozymes = self.gene_reaction_rule.split(" or ")
            isozymes = [enzyme.replace("(", "").replace(")", "").strip() for enzyme in isozymes]
            if "s0001" in isozymes:
                isozymes.remove("s0001")
                isozymes.append(NON_ENZYMATIC)
            return isozymes
        else:
            return [NON_ENZYMATIC]

    @property
    def isozyme_ids(self):
        return ["+".join(g.strip() for g in sorted(isozyme.split(" and "))) for isozyme in self.isozymes]

    @property
    def _isozyme_constraints(self):
        if self._model:
            solver = self._model.solver
            if self._forward_isozyme_constraint is None:
                fwd_id = "iso_%s_fwd" % self.id
                if fwd_id not in solver.constraints:
                    fwd_vars = self._forward_variables.values()
                    expression = sum(fwd_vars) - self.forward_variable
                    self._forward_isozyme_constraint = solver.interface.Constraint(expression, lb=0, ub=0, name=fwd_id)
                    solver.add(self._forward_isozyme_constraint, sloppy=True)
                else:
                    self._forward_isozyme_constraint = solver.constraints[fwd_id]

            if self._reverse_isozyme_constraint is None:
                rev_id = "iso_%s_rev" % self.id
                if rev_id not in solver.constraints:
                    rev_vars = self._reverse_variables.values()
                    expression = sum(rev_vars) - self.reverse_variable
                    self._reverse_isozyme_constraint = solver.interface.Constraint(expression, lb=0, ub=0, name=rev_id)
                    solver.add(self._reverse_isozyme_constraint, sloppy=True)

            return self._forward_isozyme_constraint, self._reverse_isozyme_constraint
        else:
            return None, None

    def _split_isozyme_variables(self):
        for isozyme_id in self.isozyme_ids:
            if isozyme_id == "s0001":
                continue
            if isozyme_id not in self._forward_variables:
                fwd_var_id = "%s_%s_fwd" % (self.id, isozyme_id)
                rev_var_id = "%s_%s_rev" % (self.id, isozyme_id)
                fwd_ = self._forward_variables[isozyme_id] = self._model.solver.interface.Variable(fwd_var_id, lb=0)
                rev_ = self._reverse_variables[isozyme_id] = self._model.solver.interface.Variable(rev_var_id, lb=0)
                self._model.solver.add([fwd_, rev_])

    @cobra.core.Reaction.gene_reaction_rule.setter
    def gene_reaction_rule(self, rule):
        old_isozyme_ids = self.isozyme_ids
        cobra.core.Reaction.gene_reaction_rule.fset(self, rule)
        self._clone_genes(self.model)
        self._split_isozyme_variables()
        for old_id in old_isozyme_ids:
            if old_id not in self._forward_variables:
                self._remove_isozyme(old_id)
        self._isozyme_constraints

    def _remove_isozyme(self, isozyme_id):
        self.model.solver.remove([self._forward_variables[isozyme_id], self._reverse_variables[isozyme_id]])

    @property
    def reversibility(self):
        return self._lower_bound < 0 < self._upper_bound

    @property
    def reactions(self):
        subreactions = []
        for isozyme_id in self.isozyme_ids:
            if isozyme_id == NON_ENZYMATIC:
                subreactions.append("%s : %s" % (self.id, self.reaction))
            else:
                reactants = []
                for metabolite in self.reactants:
                    if self.metabolites[metabolite] != -1:
                        reactants.append("%i %s" % (-self.metabolites[metabolite], metabolite.id))
                    else:
                        reactants.append(metabolite.id)

                products = []
                for metabolite in self.products:
                    if self.metabolites[metabolite] != 1:
                        products.append("%i %s" % (self.metabolites[metabolite], metabolite.id))
                    else:
                        products.append(metabolite.id)

                genes = isozyme_id.split("+")

                subreactions.append("fwd_%s : %s + (%s) --> %s" %
                                    (isozyme_id, " + ".join(reactants), " + ".join(genes), " + ".join(products)))
                subreactions.append("rev_%s : %s + (%s) --> %s" %
                                    (isozyme_id, " + ".join(products), " + ".join(genes), " + ".join(reactants)))

        return subreactions

    def _get_reverse_id(self):
        """Generate the id of reverse_variable from the reaction's id."""
        return '_'.join((self.id, 'reverse'))

    def _get_forward_id(self):
        """Generate the id of forward_variable from the reaction's id."""
        return self.id
        # return '_'.join((self.id, 'forward', hashlib.md5(self.id.encode('utf-8')).hexdigest()[0:5]))

    @property
    def flux_expression(self):
        """An optlang variable representing the forward flux (if associated with model), otherwise None.
        Representing the net flux if model.reversible_encoding == 'unsplit'"""
        model = self.model
        if model is not None:
            return 1. * self.forward_variable - 1. * self.reverse_variable
        else:
            return None

    @property
    def forward_variable(self):
        """An optlang variable representing the forward flux (if associated with model), otherwise None."""
        model = self.model
        if model is not None:
            if self._forward_variable is None:
                self._forward_variable = model.solver.variables[self._get_forward_id()]
            assert self._forward_variable.problem is self.model.solver

            return self._forward_variable
            # return model.solver.variables[self._get_forward_id()]
        else:
            return None

    @property
    def reverse_variable(self):
        """An optlang variable representing the reverse flux (if associated with model), otherwise None."""
        model = self.model
        if model is not None:
            if self._reverse_variable is None:
                self._reverse_variable = model.solver.variables[self._get_reverse_id()]
            assert self._reverse_variable.problem is self.model.solver
            return self._reverse_variable
            # return model.solver.variables[self._get_reverse_id()]
        else:
            return None

    def _reset_var_cache(self):
        self._forward_variables.clear()
        self._reverse_variables.clear()
        self._forward_variable = None
        self._reverse_variable = None
        self._forward_isozyme_constraint = None
        self._reverse_isozyme_constraint = None

    def add_metabolites(self, metabolites, combine=True, **kwargs):
        if combine:
            old_coefficients = self.metabolites
        super(Reaction, self).add_metabolites(metabolites, combine=combine, **kwargs)
        model = self.model
        if model is not None:
            for metabolite, coefficient in six.iteritems(metabolites):

                if isinstance(metabolite, six.string_types):  # support metabolites added as strings.
                    metabolite = model.metabolites.get_by_id(metabolite)
                if combine:
                    try:
                        old_coefficient = old_coefficients[metabolite]
                    except KeyError:
                        pass
                    else:
                        coefficient = coefficient + old_coefficient
                for isozyme_id in self.isozymes:
                    model.solver.constraints[metabolite.id].set_linear_coefficients({
                        self.forward_variables[isozyme_id]: coefficient,
                        self.reverse_variables[isozyme_id]: -coefficient
                    })

    @property
    def id(self):
        return getattr(self, "_id", None)  # Returns None if _id is not set

    @id.setter
    def id(self, value):
        if value == self.id:
            pass
        elif not isinstance(value, six.string_types):
            raise TypeError("ID must be a string")
        elif getattr(self, "_model", None) is not None:  # (= if hasattr(self, "_model") and self._model is not None)
            if value in self.model.reactions:
                raise ValueError("The model already contains a reaction with the id:", value)
            forward_variable = self.forward_variable
            reverse_variable = self.reverse_variable

            self._id = value
            self.model.reactions._generate_index()

            forward_variable.name = self._get_forward_id()
            reverse_variable.name = self._get_reverse_id()
        else:
            self._id = value

    def pop(self, metabolite_id):
        """Removes a given metabolite from the reaction stoichiometry, and returns the coefficient.
        """
        if self._model is None:
            return super(Reaction, self).pop(metabolite_id)
        else:
            if isinstance(metabolite_id, six.string_types):
                met = self.model.metabolites.get_by_id(metabolite_id)
            else:
                met = metabolite_id
            coef = self.metabolites[met]
            self.add_metabolites({met: -coef}, combine=True)
            return coef

    def remove_from_model(self, model=None, remove_orphans=False):
        reaction_model = self.model
        forward_variables = [self.forward_variable] + list(self._forward_variables.values())
        reverse_variables = [self.reverse_variable] + list(self._reverse_variables.values())
        super(Reaction, self).remove_from_model(model, remove_orphans)
        reaction_model.solver.remove(forward_variables + reverse_variables + list(self._isozyme_constraints))
        self.model = None  # Trigger model setter, since cobrapy only sets _model

    def delete(self, remove_orphans=False):
        reaction_model = self.model
        forward_variables = [self.forward_variable] + list(self._forward_variables.values())
        reverse_variables = [self.reverse_variable] + list(self._reverse_variables.values())
        super(Reaction, self).delete(remove_orphans)
        reaction_model.solver.remove(forward_variables + reverse_variables + list(self._isozyme_constraints))
        self.model = None  # Trigger model setter, since cobrapy only sets _model
        # if remove_orphans:
        #     model.solver.remove([metabolite.model.solver for metabolite in self.metabolites.keys()])

    def change_bounds(self, lb=None, ub=None, time_machine=None):
        """Changes one or both of the reaction bounds and allows the changes to be reversed with a TimeMachine"""
        if time_machine is not None:
            time_machine(do=int,
                         undo=partial(setattr, self, "lower_bound", self.lower_bound))
            time_machine(do=int,
                         undo=partial(setattr, self, "upper_bound", self.upper_bound))
        if lb is not None:
            self.lower_bound = lb
        if ub is not None:
            self.upper_bound = ub

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Id</strong></td><td>%s</td>
            </tr>
            <tr>
                <td><strong>Name</strong></td><td>%s</td>
            </tr>
            <tr>
                <td><strong>Stoichiometry</strong></td><td>%s</td>
            </tr>
            <tr>
                <td><strong>GPR</strong></td><td>%s</td>
            </tr>
            <tr>
                <td><strong>Lower bound</strong></td><td>%f</td>
            </tr>
            <tr>
                <td><strong>Upper bound</strong></td><td>%f</td>
            </tr>
        </table>
        """ % (self.id, self.name, self.reaction, self.gene_reaction_rule, self.lower_bound, self.upper_bound)


class GPRGene(Gene):

    @property
    def usage(self):
        if self.model is not None:
            return self.variable.primal
        else:
            return None

    @property
    def variable(self):
        if self.model:
            var_id = "u_%s" % self.id
            if var_id not in self.model.solver.variables:
                var = self.model.solver.interface.Variable(var_id, lb=0)
                self.model.solver.add(var, sloppy=True)
            else:
                var = self.model.solver.variables[var_id]
            return var
        else:
            raise AssertionError("Cannot retrieve a variable from a Gene without a model")

    @property
    def constraint(self):
        constraint_id = "u_%s" % self.id
        if constraint_id not in self.model.solver.constraints:
            variables = flatten(list(r.forward_variables.values()) + list(r.reverse_variables.values())
                                for r in self.reactions)
            exp = sum(-var for var in variables) + self.variable
            constraint = self.model.solver.interface.Constraint(exp, name=constraint_id, lb=0)
            self.model.solver.add(constraint, sloppy=True)
        else:
            constraint = self.model.solver.constraints[constraint_id]
        return constraint


class Protein(Object):
    def __init__(self, id="", name=None, genes=None, model=None):
        self._id = id
        super(Protein, self).__init__(id=id, name=name)
        self.model = model
        self.genes = DictList(genes or [])
        self._constraints = {}
        if genes:
            self._constraints.update({gene.id: gene.constraint for gene in genes})

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, id):
        self._id = id

    @property
    def constraints(self):
        if self.model:
            for gene in self.genes:
                if gene.id not in self._constraints:
                    self._constraints[gene.id] = gene.constraint
            return self._constraints
        else:
            raise AssertionError("Cannot retrieve a variable from a Protein without a model")


class GPRBasedModel(SolverBasedModel):
    def __init__(self, description=None, solver_interface=optlang, **kwargs):
        super(SolverBasedModel, self).__init__(description, **kwargs)
        self._solver = solver_interface.Model()
        self._solver.objective = solver_interface.Objective(S.Zero)

        self.proteins = DictList()
        cleaned_reactions = cobra.core.DictList()
        for reaction in self.reactions:
            reaction = GPRReaction.clone(reaction, model=self)
            cleaned_reactions.append(reaction)
        self.reactions = cleaned_reactions

        cleaned_genes = cobra.core.DictList()
        for gene in self.genes:
            cleaned_genes.append(GPRGene.clone(gene, model=self))
        self.genes = cleaned_genes

        cleaned_metabolites = cobra.core.DictList()
        for metabolite in self.metabolites:
            cleaned_metabolites.append(Metabolite.clone(metabolite, model=self))
        self.metabolites = cleaned_metabolites

        for metabolite in self.metabolites:
            metabolite._model = self
            metabolite._reaction = {self.reactions.get_by_id(re.id) for re in metabolite.reactions}

        for gene in self.genes:
            gene._model = self
            gene._reaction = {self.reactions.get_by_id(re.id) for re in gene.reactions}

        for reaction in self.reactions:
            reaction._genes = {self.genes.get_by_id(gene.id) for gene in reaction.genes}
            reaction._metabolites = {self.metabolites.get_by_id(met.id): coeff for
                                     met, coeff in six.iteritems(reaction.metabolites)}

        self._populate_solver(self.reactions, self.metabolites, self.genes)
        self._timestamp_last_optimization = None
        self.solution = LazySolution(self)
        for reaction in self.reactions:
            self._add_enzymes_for_reaction(reaction)

    def _add_enzymes_for_reaction(self, reaction):
        if len(reaction.gene_reaction_rule) > 0:
            for enzyme_id in reaction.isozyme_ids:
                gene_complex = set([g.strip() for g in enzyme_id.split("+")])
                try:
                    self.proteins.get_by_id(enzyme_id)
                except KeyError:
                    genes = [self.genes.get_by_id(gene) for gene in gene_complex if gene != NON_ENZYMATIC]
                    protein = Protein(id=enzyme_id, genes=genes, model=self)
                    self.proteins.append(protein)

    def _populate_solver(self, reaction_list, metabolite_list=None, genes_list=None):
        constraint_terms = AutoVivification()
        metabolite_constraints = AutoVivification()
        gene_constraints = AutoVivification()

        if metabolite_list is not None:
            for met in metabolite_list:
                constraint = self.solver.interface.Constraint(S.Zero, name=met.id, lb=0, ub=0)
                metabolite_constraints[met.id] = constraint
                self.solver.add(constraint, sloppy=True)

        if genes_list is not None:
            for gene in genes_list:
                constraint_id = var_id = "u_%s" % gene.id
                var = self.solver.interface.Variable(var_id, lb=0)
                constraint = self.solver.interface.Constraint(var, name=constraint_id, lb=0)
                gene_constraints[gene.id] = constraint
                self.solver.add([var, constraint], sloppy=True)

        for r in reaction_list:

            if r.reversibility:
                forward_variable = self.solver.interface.Variable(r._get_forward_id(), lb=0, ub=r._upper_bound)
                reverse_variable = self.solver.interface.Variable(r._get_reverse_id(), lb=0, ub=-1 * r._lower_bound)
            elif 0 == r.lower_bound and r.upper_bound == 0:
                forward_variable = self.solver.interface.Variable(r._get_forward_id(), lb=0, ub=0)
                reverse_variable = self.solver.interface.Variable(r._get_reverse_id(), lb=0, ub=0)
            elif r.lower_bound >= 0:
                forward_variable = self.solver.interface.Variable(r._get_forward_id(), lb=r._lower_bound,
                                                                  ub=r._upper_bound)
                reverse_variable = self.solver.interface.Variable(r._get_reverse_id(), lb=0, ub=0)
            elif r.upper_bound <= 0:
                forward_variable = self.solver.interface.Variable(r._get_forward_id(), lb=0, ub=0)
                reverse_variable = self.solver.interface.Variable(r._get_reverse_id(),
                                                                  lb=-1 * r._upper_bound,
                                                                  ub=-1 * r._lower_bound)

            setattr(r, "_forward_variable", forward_variable)
            setattr(r, "_reverse_variable", reverse_variable)

            self.solver.add(forward_variable)
            self.solver.add(reverse_variable)
            self.solver.update()

            if r._forward_isozyme_constraint:
                assert r._reverse_isozyme_constraint is None
                isozyme_constraint_fwd = r._forward_isozyme_constraint
                isozyme_constraint_rev = r._reverse_isozyme_constraint
            else:
                iso_fwd_id, iso_rev_id = "iso_%s_fwd" % r.id, "iso_%s_rev" % r.id
                isozyme_constraint_fwd = self.solver.interface.Constraint(S.Zero, name=iso_fwd_id, lb=0, ub=0)
                isozyme_constraint_rev = self.solver.interface.Constraint(S.Zero, name=iso_rev_id, lb=0, ub=0)
                constraint_terms[isozyme_constraint_fwd][forward_variable] = 1
                constraint_terms[isozyme_constraint_rev][reverse_variable] = 1
                self.solver.add([isozyme_constraint_fwd, isozyme_constraint_rev], sloppy=True)
                r._forward_isozyme_constraint = isozyme_constraint_fwd
                r._reverse_isozyme_constraint = isozyme_constraint_rev

            for isozyme_id in r.isozyme_ids:
                if isozyme_id not in r._forward_variables:
                    assert isozyme_id not in r._reverse_variables
                    fwd_var_id, rev_var_id = "%s_%s_fwd" % (r.id, isozyme_id), "%s_%s_rev" % (r.id, isozyme_id)

                    r._forward_variables[isozyme_id] = fwd_var = self.solver.interface.Variable(fwd_var_id, lb=0)
                    r._reverse_variables[isozyme_id] = rev_var = self.solver.interface.Variable(rev_var_id, lb=0)

                    constraint_terms[isozyme_constraint_fwd][fwd_var] = -1
                    constraint_terms[isozyme_constraint_rev][rev_var] = -1

                    self.solver.add([fwd_var, rev_var])
                else:
                    fwd_var = r._forward_variables[isozyme_id]
                    rev_var = r._reverse_variables[isozyme_id]

                for metabolite, coeff in six.iteritems(r.metabolites):
                    if metabolite.id in metabolite_constraints:
                        constraint = metabolite_constraints[metabolite.id]
                    else:
                        constraint = self.solver.interface.Constraint(S.Zero, name=metabolite.id, lb=0, ub=0)
                        self.solver.add(constraint, sloppy=True)

                    constraint_terms[constraint][fwd_var] = coeff
                    constraint_terms[constraint][rev_var] = -coeff

                for gene in r.genes:
                    if gene.id == "s0001":
                        continue
                    if gene.id in gene_constraints:
                        constraint = gene_constraints[gene.id]
                    else:
                        constraint = gene.constraint

                    constraint_terms[constraint][fwd_var] = -1
                    constraint_terms[constraint][rev_var] = -1

            objective_coeff = r._objective_coefficient
            if objective_coeff != 0.:
                if self.solver.objective is None:
                    self.solver.objective = self.solver.interface.Objective(0, direction='max')
                if self.solver.objective.direction == 'min':
                    self.solver.objective.direction = 'max'
                self.solver.objective.set_linear_coefficients({forward_variable: objective_coeff,
                                                               reverse_variable: -objective_coeff})

        self.solver.update()
        for constraint, terms in six.iteritems(constraint_terms):
            constraint.set_linear_coefficients(terms)

    def add_reactions(self, reaction_list):
        clean_reaction_list = []
        for reaction in reaction_list:
            clean_reaction_list.append(GPRReaction.clone(reaction, model=self))
        super(GPRBasedModel, self).add_reactions(clean_reaction_list)

        for reaction in clean_reaction_list:
            self._add_enzymes_for_reaction(reaction)
