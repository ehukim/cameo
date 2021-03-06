# -*- coding: utf-8 -*-
# Copyright 2014 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import absolute_import, print_function

import numpy as np
from cobra import Reaction
from inspyred.ec.emo import Pareto

from cameo import config, flux_variability_analysis

__all__ = ['biomass_product_coupled_yield', 'product_yield', 'number_of_knockouts']


class ObjectiveFunction(object):
    """
    Blueprint for objective function.

    All objective functions must extend this class and override the __call__ method.
    Because they implement the __call__ method, they will behave like functions.

    Methods
    -------
    __call__(model, solution, decoded_representation)
        Calculates the fitness of the solution

    """

    def __init__(self, *args, **kwargs):
        super(ObjectiveFunction, self).__init__(*args, **kwargs)

    def __call__(self, model, solution, decoded_representation):
        raise NotImplementedError

    @property
    def reactions(self):
        """
        Returns the reactions ids that require flux data for this evaluation function.
        """
        return []

    def _repr_latex_(self):
        return self.name

    @property
    def name(self):
        return self.__class__.__name__

    def worst_fitness(self, maximize=True):
        """

        Parameters
        ----------
        maximize: bool
            The sense of the optimization problem.

        Returns
        -------
        float
            The worst fitness value that a given objective can have based on the sense of the optimization.
        """
        raise NotImplementedError

    def __len__(self):
        return 1


class MultiObjectiveFunction(ObjectiveFunction):
    def __init__(self, objectives, *args, **kwargs):
        super(MultiObjectiveFunction, self).__init__(*args, **kwargs)
        failed = []
        for i, objective in enumerate(objectives):
            if not isinstance(objective, ObjectiveFunction):
                failed.append(i)

        if len(failed) > 0:
            raise ValueError("Objectives %s are not instance of ObjectiveFunction" %
                             ", ".join(str(objectives[i]) for i in failed))

        self.objectives = objectives

    def __getitem__(self, item):
        return self.objectives[item]

    def __call__(self, model, solution, targets):
        return Pareto(values=[o(model, solution, targets) for o in self.objectives])

    def worst_fitness(self, maximize=True):
        return Pareto(values=[o.worst_fitness(maximize) for o in self.objectives], maximize=maximize)

    @property
    def reactions(self):
        return set(sum([o.reactions for o in self.objectives], []))

    def __len__(self):
        return len(self.objectives)

    @property
    def name(self):
        return "MO: " + "|".join([of.name for of in self.objectives])

    def _repr_latex_(self):
        return "\\begin{gather}\n" + "\\\\".join(o._repr_latex_() for o in self.objectives) + "\n\\end{gather}\n"


class YieldFunction(ObjectiveFunction):
    def __init__(self, product, substrate, *args, **kwargs):
        super(YieldFunction, self).__init__(*args, **kwargs)
        self.n_c_substrate = 1
        self.n_c_product = 1
        if isinstance(product, Reaction):
            self.n_c_product = product.reactants[0].elements.get('C', 1)
            product = product.id
        self.product = product
        if isinstance(substrate, Reaction):
            self.n_c_substrate = substrate.reactants[0].elements.get('C', 1)
            substrate = substrate.id
        self.substrate = substrate
        self.__name__ = self.__class__.__name__

    def __call__(self, model, solution, targets):
        raise NotImplementedError

    def worst_fitness(self, maximize=True):
        if maximize:
            return 0
        else:
            return -np.inf


class biomass_product_coupled_yield(YieldFunction):
    """
    Biomass-Product Coupled Yield: (v[biomass] * v[product]) / v[substrate] [1]

    Parameters
    ----------
    biomass: str or Reaction
        biomass reaction identifier
    product: str or Reaction
        product reaction identifier
    substrate: str or Reaction
        substrate reaction identifier

    Returns
    -------
    float
        fitness value

    References
    ----------
    [1] Patil, K. R., Rocha, I., Förster, J., & Nielsen, J. (2005). "Evolutionary programming as a
    platform for in silico metabolic engineering". BMC Bioinformatics, 6, 308.
    doi:10.1186/1471-2105-6-308
    """

    def __init__(self, biomass, product, substrate, *args, **kwargs):
        super(biomass_product_coupled_yield, self).__init__(product, substrate, *args, **kwargs)
        if isinstance(biomass, Reaction):
            biomass = biomass.id
        self.biomass = biomass

    def __call__(self, model, solution, targets):
        try:
            biomass_flux = round(solution.fluxes[self.biomass], config.ndecimals)
            product_flux = round(solution.fluxes[self.product], config.ndecimals)
            substrate_flux = round(abs(solution.fluxes[self.substrate]), config.ndecimals)
            return round((biomass_flux * product_flux) / substrate_flux, config.ndecimals)

        except ZeroDivisionError:
            return 0.0

    def _repr_latex_(self):
        return "$$bpcy = \\frac{(%s * %s)}{%s}$$" % (
            self.biomass.replace("_", "\\_"), self.product.replace("_", "\\_"), self.substrate.replace("_", "\\_"))

    @property
    def name(self):
        return "bpcy = (%s * %s) / %s" % (self.biomass, self.product, self.substrate)

    @property
    def reactions(self):
        return [self.biomass, self.product, self.substrate]


class biomass_product_coupled_min_yield(YieldFunction):
    """
    Biomass-Product Coupled Minimum Yield: (v[biomass] * min(v[product])) / v[substrate] [1]

    Parameters
    ----------
    biomass: str or Reaction
        biomass reaction identifier
    product: str or Reaction
        product reaction identifier
    substrate: str or Reaction
        substrate reaction identifier

    Returns
    -------
    float
        fitness value

    """

    def __init__(self, biomass, product, substrate, *args, **kwargs):
        super(biomass_product_coupled_min_yield, self).__init__(product, substrate, *args, **kwargs)
        if isinstance(biomass, Reaction):
            biomass = biomass.id
        self.biomass = biomass

    def __call__(self, model, solution, targets):
        try:
            biomass_flux = round(solution.fluxes[self.biomass], config.ndecimals)
            fva_res = flux_variability_analysis(model, reactions=[self.product], fraction_of_optimum=1)
            min_product_flux = fva_res["lower_bound"][self.product]

            substrate_flux = round(abs(solution.fluxes[self.substrate]), config.ndecimals)
            return round((biomass_flux * min_product_flux) / substrate_flux, config.ndecimals)

        except ZeroDivisionError:
            return 0.0

    def _repr_latex_(self):
        return "$$bpcy = \\frac{(%s * min(%s))}{%s}$$" % (
            self.biomass.replace("_", "\\_"), self.product.replace("_", "\\_"), self.substrate.replace("_", "\\_"))

    @property
    def name(self):
        return "bpcy = (%s * min(%s)) / %s" % (self.biomass, self.product, self.substrate)

    @property
    def reactions(self):
        return [self.biomass, self.product, self.substrate]


class product_yield(YieldFunction):
    """
    Product Yield Objective function: v[product]/v[substrate]

    scaled to carbon content if the product and substrate are given as `Reaction` objects

    Parameters
    ----------
    product: str or Reaction
        product reaction identifier
    substrate: str or Reaction
        substrate reaction identifier

    Returns
    -------
    float
        fitness value
    """

    def __init__(self, product, substrate, *args, **kwargs):
        super(product_yield, self).__init__(product, substrate, *args, **kwargs)

    def __call__(self, model, solution, targets):
        try:
            product_flux = round(solution.fluxes[self.product], config.ndecimals) * self.n_c_product
            substrate_flux = round(abs(solution.fluxes[self.substrate]), config.ndecimals) * self.n_c_substrate
            return round(product_flux / substrate_flux, config.ndecimals)
        except ZeroDivisionError:
            return 0.0

    def _repr_latex_(self):
        return "$$yield = \\frac{%s}{%s}$$" % (self.product.replace('_', '\\_'), self.substrate.replace('_', '\\_'))

    @property
    def name(self):
        return "yield = (%s / %s)" % (self.product, self.substrate)

    @property
    def reactions(self):
        return [self.product, self.substrate]


class number_of_knockouts(ObjectiveFunction):
    """
    Number of Knockouts objective function.
    If sense is maximize then fitness is the number of knockouts, otherwise 1/#knockouts

    Parameters
    ----------
    sense: str
        'max' or 'min'

    Returns
    -------
    float
        fitness value
    """

    def __init__(self, sense='min', *args, **kwargs):
        super(number_of_knockouts, self).__init__(*args, **kwargs)
        self.sense = sense

    def __call__(self, model, solution, targets):
        if self.sense == 'max':
            return len(targets)
        else:
            return round(1.0 / len(targets), config.ndecimals)

    def _repr_latex_(self):
        return "$$ %s\\:\\#knockouts $$" % self.sense

    @property
    def name(self):
        if self.sense == 'max':
            return "max knockouts"
        else:
            return "min knockouts"

    def worst_fitness(self, maximize=True):
        if maximize:
            return 0
        else:
            return np.inf
