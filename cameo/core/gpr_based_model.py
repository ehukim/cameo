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

import optlang
from cobra import DictList

from cameo.core import SolverBasedModel, Gene, Reaction
from cameo.util import flatten

__all__ = ["GPRBasedModel"]


class GPRGene(Gene):

    @property
    def variable(self):
        if self.model:
            if self.id not in self.model.solver.variables:
                var = self.model.solver.interface.Variable(self.id)
                self.model.solver.add(var, sloppy=True)
            else:
                var = self.model.solver.variables[self.id]
            return var
        else:
            raise AssertionError("Cannot retrieve a variable from a Gene without a model")

    @property
    def constraint(self):
        if self.model:
            if self.id not in self.model.solver.constraints:
                constraint = self._build_constraint()
                self.model.solver.add(constraint, sloppy=True)

            return self.model.solver.constraints[self.id]
        raise AssertionError("Cannot retrieve a constraint from a Gene without a model")

    def _build_constraint(self):
        expression = sum(-r.forward_variable - r.reverse_variable for r in self.reactions) + self.variable
        return self.model.solver.interface.Constraint(expression, name=self.id, lb=0, ub=0)


class GPRBasedModel(SolverBasedModel):
    def __init__(self, description=None, solver_interface=optlang, **kwargs):
        super(GPRBasedModel, self).__init__(description, **kwargs)
        old_reactions = self.reactions
        self.reactions = DictList()
        for reaction in old_reactions:
            self.reactions.append(Reaction.clone(reaction, model=self, gene_class=GPRGene))

        self._populate_solver_genes(self.genes)

    def _populate_solver_genes(self, genes):
        for gene in genes:
            self.solver.add(gene._build_constraint(), sloppy=True)

    def add_reactions(self, reaction_list):
        for index, reaction in enumerate(reaction_list):
            reaction_list[index] = Reaction.clone(reaction, model=self, gene_class=GPRBasedModel)
        super(GPRBasedModel, self).add_reactions(reaction_list)

        genes = flatten([reaction.genes for reaction in reaction_list])
        for gene in genes:
            if gene.id not in self.solver.variables:
                self.solver.variables.add(gene.variable)
