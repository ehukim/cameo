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
from warnings import warn

import optlang
from cobra import DictList, Object

from cameo.core.solver_based_model import SolverBasedModel


class Enzyme(Object):
    def __init__(self, id=None, name=None, genes=None, reactions=None, model=None):
        super(Enzyme, self).__init__(id, name)
        self._id = id
        self.genes = DictList(genes or [])
        self.model = model
        self._reactions = DictList()
        for reaction in reactions:
            self.add_reaction(reaction)

    @property
    def reactions(self):
        return self._reactions

    def add_reaction(self, reaction):
        if reaction not in self.reactions:
            self.constraint.set_linear_coefficients({
                reaction.forward_variable: -1,
                reaction.reverse_variable: -1
            })
            self._reactions.append(reaction)

    def remove_reaction(self, reaction):
        if reaction in self.reactions:
            self.constraint.set_linear_coefficients({
                reaction.forward_variable: -0,
                reaction.reverse_variable: -0
            })
        else:
            warn("Reaction %s is not present in enzyme" % reaction.id)

    @property
    def variable(self):
        if self.model:
            if self.id not in self.model.solver.variables:
                var = self.model.solver.interface.Variable(self.id)
                self.model.solver.add(var, sloppy=True)
            else:
                var = self.model.solver.variables[self.id]
            return var

        raise AssertionError("Cannot retrieve a variable from a Enzyme without a model")

    @property
    def constraint(self):
        if self.model:
            if self.id not in self.model.solver.constraints:
                expression = sum(-r.forward_variable - r.reverse_variable for r in self.reactions) + self.variable
                constraint = self.model.solver.interface.Constraint(expression, name=self.id, lb=0, ub=0)
                self.model.solver.add(constraint, sloppy=True)

            return self.model.solver.constraints[self.id]

        raise AssertionError("Cannot retrieve a constraint from a Enzyme without a model")

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, id):
        self._id = id


class GPRBasedModel(SolverBasedModel):
    def __init__(self, description=None, solver_interface=optlang, **kwargs):
        self.enzymes = DictList()
        super(GPRBasedModel, self).__init__(description, **kwargs)

        for reaction in self.reactions:
            self._add_enzyme(reaction)

    def _add_enzyme(self, reaction):
        if len(reaction.gene_reaction_rule) == 0:
            return
        isozymes = reaction.gene_reaction_rule.split(" or ")
        for enzyme in isozymes:
            enzyme = enzyme.replace("(", "").replace(")", "")
            genes = [self.genes.get_by_id(gene) for gene in enzyme.split(" and ")]
            enzyme = "_".join(g.id for g in sorted(genes, key=lambda g: g.id))
            try:
                self.enzymes.get_by_id(enzyme).add_reaction(reaction)
            except KeyError:
                self.enzymes.append(Enzyme(id=enzyme, genes=genes, reactions=[reaction], model=self))

    def add_reactions(self, reaction_list):
        super(GPRBasedModel, self).add_reactions(reaction_list)

        for reaction in reaction_list:
            self._add_enzyme(reaction)

    def remove_reactions(self, the_reactions, delete=True, remove_orphans=False):
        super(GPRBasedModel, self).remove_reactions(the_reactions, delete=delete, remove_orphans=remove_orphans)

        for reaction in the_reactions:
            for enzyme in self.enzymes:
                if reaction in enzyme.reactions:
                    enzyme.remove_reaction(reaction)
