import numba
import numpy as np
import pandas as pd
import tskit

from .base import _check_instance, _check_dataframe, _check_non_decreasing  # noreorder


@numba.njit
def _traversal_genotype(
    nodes_individual,
    left_child_array,
    right_sib_array,
    stack,
    has_mutation,
    num_individuals,
    num_nodes,
):  # pragma: no cover
    """
    Numba to speed up the tree traversal algorithm to determine the genotype of
    individuals.
    Stack has to be Typed List in numba to use numba.
    """

    genotype = np.zeros(num_individuals)
    while len(stack) > 0:
        parent_node_id = stack.pop()
        if parent_node_id == num_nodes:
            individual_id = -1
        else:
            individual_id = nodes_individual[parent_node_id]
        if individual_id > -1:
            genotype[individual_id] += 1
        child_node_id = left_child_array[parent_node_id]
        while child_node_id != -1:
            if not has_mutation[child_node_id]:
                stack.append(child_node_id)
            child_node_id = right_sib_array[child_node_id]

    return genotype


class _GeneticValue:
    """GeneticValue class to compute genetic values of individuals.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Tree sequence data with mutation
    trait_df : pandas.DataFrame
        Dataframe that includes causal site ID, causal state, simulated effect
        size, and trait ID.
    alpha : float
        Parameter that determines the relative weight on rarer variants.
    """

    def __init__(self, ts, trait_df):
        self.trait_df = trait_df[["site_id", "effect_size", "trait_id", "causal_state"]]
        self.ts = ts

    def _individual_genotype(self, tree, site, causal_state, num_nodes):
        """
        Returns a numpy array that describes the number of causal mutation in an
        individual.
        """
        has_mutation = np.zeros(num_nodes + 1, dtype=bool)
        state_transitions = {tree.virtual_root: site.ancestral_state}
        for m in site.mutations:
            state_transitions[m.node] = m.derived_state
            has_mutation[m.node] = True
        stack = numba.typed.List()
        for node, state in state_transitions.items():
            if state == causal_state:
                stack.append(node)

        if len(stack) == 0:  # pragma: no cover
            genotype = np.zeros(self.ts.num_individuals)
        else:
            genotype = _traversal_genotype(
                nodes_individual=self.ts.nodes_individual,
                left_child_array=tree.left_child_array,
                right_sib_array=tree.right_sib_array,
                stack=stack,
                has_mutation=has_mutation,
                num_individuals=self.ts.num_individuals,
                num_nodes=num_nodes,
            )

        return genotype

    def _run(self):
        """Computes genetic values of individuals.

        Returns
        -------
        pandas.DataFrame
            Dataframe with simulated genetic value, individual ID, and trait ID.
        """

        num_ind = self.ts.num_individuals
        num_trait = np.max(self.trait_df.trait_id) + 1
        genetic_val_array = np.zeros((num_trait, num_ind))
        num_nodes = self.ts.num_nodes
        tree = tskit.Tree(self.ts)

        for data in self.trait_df.itertuples():
            site = self.ts.site(data.site_id)
            tree.seek(site.position)
            individual_genotype = self._individual_genotype(
                tree=tree,
                site=site,
                causal_state=data.causal_state,
                num_nodes=num_nodes,
            )
            genetic_val_array[data.trait_id, :] += (
                individual_genotype * data.effect_size
            )

        df = pd.DataFrame(
            {
                "trait_id": np.repeat(np.arange(num_trait), num_ind),
                "individual_id": np.tile(np.arange(num_ind), num_trait),
                "genetic_value": genetic_val_array.flatten(),
            }
        )

        return df


def genetic_value(ts, trait_df):
    """
    Obtains genetic value from a trait dataframe.

    Parameters
    ----------
    ts : tskit.TreeSequence
        The tree sequence data that will be used in the quantitative trait
        simulation.
    trait_df : pandas.DataFrame
        Trait dataframe.

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe that includes genetic value of individuals in the
        tree sequence data.

    See Also
    --------
    trait_model : Return a trait model, which can be used as `model` input.
    sim_trait : Return a trait dataframe, whch can be used as a `trait_df` input.
    GeneticResult : Dataclass object that will be used as an output.
    sim_env : Genetic value dataframe output can be used as an input to simulate
        environmental noise.

    Notes
    -----
    The `trait_df` input has some requirements that will be noted below.

    1. Columns

    The following columns must be included in `trait_df`:

        * **site_id**: Site IDs that have causal mutation.
        * **effect_size**: Simulated effect size of causal mutation.
        * **causal_state**: Derived state of causal mutation.
        * **trait_id**: Trait ID.

    2. Data requirements

        * Site IDs in **site_id** column must be sorted in an ascending order. Please
          refer to :py:meth:`pandas.DataFrame.sort_values` for details on sorting
          values in a :class:`pandas.DataFrame`.

        * Trait IDs in **trait_id** column must start from zero and be consecutive.

    The phenotype output of the simulation is given as a :py:class:`pandas.DataFrame`.

    The genetic value dataframe contains the following columns:

        * **trait_id**: Trait ID.
        * **individual_id**: Individual ID inside the tree sequence input.
        * **genetic_value**: Simulated genetic values.

    Examples
    --------
    See :ref:`genetic_value` for worked examples.
    """

    ts = _check_instance(ts, "ts", tskit.TreeSequence)
    if ts.num_individuals == 0:
        raise ValueError("No individuals in the provided tree sequence dataset")
    trait_df = _check_dataframe(
        trait_df, ["site_id", "effect_size", "trait_id", "causal_state"], "trait_df"
    )
    _check_non_decreasing(trait_df["site_id"], "site_id")

    trait_id = trait_df["trait_id"].unique()

    if np.min(trait_id) != 0 or np.max(trait_id) != len(trait_id) - 1:
        raise ValueError("trait_id must be consecutive and start from 0")

    genetic = _GeneticValue(ts=ts, trait_df=trait_df)

    genetic_result = genetic._run()

    return genetic_result