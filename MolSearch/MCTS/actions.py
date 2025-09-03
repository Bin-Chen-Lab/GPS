# coding=utf-8
# Copyright 2020 The Google Research Authors.
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

# Lint as: python2, python3
"""Defines the Markov decision process of generating a molecule.

The problem of molecule generation as a Markov decision process, the
state space, action space, and reward function are defined.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import copy
import itertools

import utils

from scores import *
import bottleneck
import numpy as np

from design_moves import DesignMove
replaceRule = DesignMove("chemblDB3.sqlitdb")


def top_k(valid_actions, score_func, small=False, k=3):
    if len(valid_actions) <= k:
        return valid_actions
    scores = []
    for s in valid_actions:
        scores.append(score_func(Chem.MolFromSmiles(s)))
    scores = np.array(scores)
    assert len(scores)==len(valid_actions)

    if len(scores) <= 5:
        if small:
            top5_tuple = [(scores[idx], valid_actions[idx]) for idx in range(len(scores))]
        else:
            top5_tuple = [(-scores[idx], valid_actions[idx]) for idx in range(len(scores))]
    else:
        if small:
            indexes = bottleneck.argpartition(scores, 5)[:5]
            top5_tuple = [(scores[idx], valid_actions[idx]) for idx in indexes] # top k method cannot used for tuples
        else:
            indexes = bottleneck.argpartition(-scores, 5)[:5]
            top5_tuple = [(-scores[idx], valid_actions[idx]) for idx in indexes]

    topk_tuple = sorted(top5_tuple)[:k]
    topk_actions = [t[1] for t in topk_tuple]
    return topk_actions

# change by Han Meng@2023Mar24
# def get_constraint_actions(actions, constraint_func, threshold, t=1.0):
#     valid_actions = []
#     for s in actions:
#         mol = Chem.MolFromSmiles(s)
#         score = constraint_func(mol)
#         if score > t * threshold:
#             valid_actions.append(s)
#     print("valid actions after constraint {:d}".format(len(valid_actions)))
#     return valid_actions


def get_constraint_actions(actions, functions, thresholds, t=1.0, l=3, L=5):
    valid_actions = []
    n_f = len(functions)
    for s in actions:
        mol = Chem.MolFromSmiles(s)
        scores = np.zeros(n_f)
        for i in range(n_f):
            scores[i] = functions[i](mol)
        if np.all(scores >= (t**(L-l)) * thresholds):
            valid_actions.append(s)
    print("t:%s,l:%s,L:%s, (t**(L-l)):%s"%(str(t),str(l),str(L),str(t**(L-l))))
    print("valid actions after constraint {:d}".format(len(valid_actions)))
    return valid_actions
# end of change

def get_mo_actions(actions, functions, thresholds, t=1.0):
    valid_actions = []
    n_f = len(functions)
    for s in actions:
        mol = Chem.MolFromSmiles(s)
        scores = np.zeros(n_f)
        for i in range(n_f):
            scores[i] = functions[i](mol)
        if np.all(scores >= t * thresholds):
            valid_actions.append(s)
    print("valid actions after constraint {:d}".format(len(valid_actions)))
    return valid_actions


def get_mo_stage2_actions(current_smiles, actions, functions, thresholds, t=1.0):
    ref_size = Chem.MolFromSmiles(current_smiles).GetNumAtoms()
    valid_actions = []
    n_f = len(functions)
    for s in actions:
        mol = Chem.MolFromSmiles(s)
        mol_size = mol.GetNumAtoms()
        scores = np.zeros(n_f)
        for i in range(n_f):
            scores[i] = functions[i](mol)
        if np.all(scores >= t * thresholds) and mol_size < ref_size:
            valid_actions.append(s)
    print("valid actions after constraint {:d}".format(len(valid_actions)))
    return valid_actions


def constraint_top_k(valid_actions, score_func, k=10):
    if len(valid_actions) <= k:
        return valid_actions
    scores = []
    for s in valid_actions:
        scores.append(score_func(Chem.MolFromSmiles(s)))
    scores = np.array(scores)
    assert len(scores)==len(valid_actions)

    all_tuple = [(-scores[idx], valid_actions[idx]) for idx in range(len(scores))]
    topk_tuple = sorted(all_tuple)[:k]
    topk_actions = [t[1] for t in topk_tuple]
    return topk_actions


def get_actions(
    state,
    atom_types=["C", "O", "N"],
    allow_removal=False,
    allow_no_modification=False,    #avoid same state as go deeper
    allowed_ring_sizes=[5,6,7],
    allow_bonds_between_rings=False,
    allow_ring_addition=False,
    allow_substitution=False,
    allow_atom_addition=True,
    allow_bond_addition=True,
):
    """Computes the set of valid actions for a given state.

  Args:
    state: String SMILES; the current state. If None or the empty string, we
      assume an "empty" state with no atoms or bonds.
    atom_types: Set of string atom types, e.g. {'C', 'O'}.
    allow_removal: Boolean whether to allow actions that remove atoms and bonds.
    allow_no_modification: Boolean whether to include a "no-op" action.
    allowed_ring_sizes: Set of integer allowed ring sizes; used to remove some
      actions that would create rings with disallowed sizes.
    allow_bonds_between_rings: Boolean whether to allow actions that add bonds
      between atoms that are both in rings.

  Returns:
    Set of string SMILES containing the valid actions (technically, the set of
    all states that are acceptable from the given state).

  Raises:
    ValueError: If state does not represent a valid molecule.
  """
    if not state:
        # Available actions are adding a node of each type.
        return copy.deepcopy(atom_types)
    mol = Chem.MolFromSmiles(state)
    if mol is None:
        raise ValueError("Received invalid state: %s" % state)
    atom_valences = {
        atom_type: utils.atom_valences([atom_type])[0] for atom_type in atom_types
    }
    atoms_with_free_valence = {}
    for i in range(1, max(atom_valences.values())):
        # Only atoms that allow us to replace at least one H with a new bond are
        # enumerated here.
        atoms_with_free_valence[i] = [
            atom.GetIdx() for atom in mol.GetAtoms() if atom.GetNumImplicitHs() >= i
        ]
    valid_actions = set()
    
    if allow_atom_addition:
        valid_actions.update(
            _atom_addition(
                mol,
                atom_types=atom_types,
                atom_valences=atom_valences,
                atoms_with_free_valence=atoms_with_free_valence,
            )
        )

    if allow_bond_addition:
        valid_actions.update(
            _bond_addition(
                mol,
                atoms_with_free_valence=atoms_with_free_valence,
                allowed_ring_sizes=allowed_ring_sizes,
                allow_bonds_between_rings=allow_bonds_between_rings,
            )
        )
    if allow_removal:
        valid_actions.update(_bond_removal(mol))
    if allow_no_modification:
        valid_actions.add(Chem.MolToSmiles(mol))
    if allow_ring_addition:                             # add rings option
        print(_ring_addition(mol))
        valid_actions.update(_ring_addition(mol))
    if allow_substitution:
        try:
            print("replace action calculation")
            valid_tmp = _frag_substitution(state, replaceRule)
            print("possible actions: {:d}".format(len(valid_tmp)))
            valid_actions.update(valid_tmp)
        except:
            pass

    return list(valid_actions)


def _frag_substitution(smi, rule, min_pairs=1):
    substitution_actions = rule.one_step_move(query_smi=smi, min_pairs=min_pairs)
    return set(substitution_actions)


def _atom_addition(state, atom_types, atom_valences, atoms_with_free_valence):
    """Computes valid actions that involve adding atoms to the graph.

  Actions:
    * Add atom (with a bond connecting it to the existing graph)

  Each added atom is connected to the graph by a bond. There is a separate
  action for connecting to (a) each existing atom with (b) each valence-allowed
  bond type. Note that the connecting bond is only of type single, double, or
  triple (no aromatic bonds are added).

  For example, if an existing carbon atom has two empty valence positions and
  the available atom types are {'C', 'O'}, this section will produce new states
  where the existing carbon is connected to (1) another carbon by a double bond,
  (2) another carbon by a single bond, (3) an oxygen by a double bond, and
  (4) an oxygen by a single bond.

  Args:
    state: RDKit Mol.
    atom_types: Set of string atom types.
    atom_valences: Dict mapping string atom types to integer valences.
    atoms_with_free_valence: Dict mapping integer minimum available valence
      values to lists of integer atom indices. For instance, all atom indices in
      atoms_with_free_valence[2] have at least two available valence positions.

  Returns:
    Set of string SMILES; the available actions.
  """
    bond_order = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    atom_addition = set()
    for i in bond_order:
        for atom in atoms_with_free_valence[i]:
            for element in atom_types:
                if atom_valences[element] >= i:
                    new_state = Chem.RWMol(state)
                    idx = new_state.AddAtom(Chem.Atom(element))
                    new_state.AddBond(atom, idx, bond_order[i])
                    sanitization_result = Chem.SanitizeMol(new_state, catchErrors=True)
                    # When sanitization fails
                    if sanitization_result:
                        continue
                    atom_addition.add(Chem.MolToSmiles(new_state))
    return atom_addition

from rdkit.Chem import AllChem
import json
def _ring_addition(state):
    rules = json.load(open('rings/add_rings_simple.json'))
    ring_addition = set()
    for rule in rules:
        rxn = AllChem.ReactionFromSmarts(rule['smarts'])
        products = rxn.RunReactants((state,))
        for p in products:
            sanitization_result = Chem.SanitizeMol(p[0], catchErrors=True)
            if sanitization_result:
                continue
            ring_addition.add(Chem.MolToSmiles(p[0]))
            #print(Chem.MolToSmiles(p[0]))
    return ring_addition


def _bond_addition(
    state, atoms_with_free_valence, allowed_ring_sizes, allow_bonds_between_rings
):
    """Computes valid actions that involve adding bonds to the graph.

  Actions (where allowed):
    * None->{single,double,triple}
    * single->{double,triple}
    * double->{triple}

  Note that aromatic bonds are not modified.

  Args:
    state: RDKit Mol.
    atoms_with_free_valence: Dict mapping integer minimum available valence
      values to lists of integer atom indices. For instance, all atom indices in
      atoms_with_free_valence[2] have at least two available valence positions.
    allowed_ring_sizes: Set of integer allowed ring sizes; used to remove some
      actions that would create rings with disallowed sizes.
    allow_bonds_between_rings: Boolean whether to allow actions that add bonds
      between atoms that are both in rings.

  Returns:
    Set of string SMILES; the available actions.
  """
    bond_orders = [
        None,
        Chem.BondType.SINGLE,
        Chem.BondType.DOUBLE,
        Chem.BondType.TRIPLE,
    ]
    bond_addition = set()
    for valence, atoms in atoms_with_free_valence.items():
        for atom1, atom2 in itertools.combinations(atoms, 2):
            # Get the bond from a copy of the molecule so that SetBondType() doesn't
            # modify the original state.
            bond = Chem.Mol(state).GetBondBetweenAtoms(atom1, atom2)
            new_state = Chem.RWMol(state)
            # Kekulize the new state to avoid sanitization errors; note that bonds
            # that are aromatic in the original state are not modified (this is
            # enforced by getting the bond from the original state with
            # GetBondBetweenAtoms()).
            Chem.Kekulize(new_state, clearAromaticFlags=True)
            if bond is not None:
                if bond.GetBondType() not in bond_orders:
                    continue  # Skip aromatic bonds.
                idx = bond.GetIdx()
                # Compute the new bond order as an offset from the current bond order.
                bond_order = bond_orders.index(bond.GetBondType())
                bond_order += valence
                if bond_order < len(bond_orders):
                    idx = bond.GetIdx()
                    bond.SetBondType(bond_orders[bond_order])
                    new_state.ReplaceBond(idx, bond)
                else:
                    continue
            # If do not allow new bonds between atoms already in rings.
            elif not allow_bonds_between_rings and (
                state.GetAtomWithIdx(atom1).IsInRing()
                and state.GetAtomWithIdx(atom2).IsInRing()
            ):
                continue
            # If the distance between the current two atoms is not in the
            # allowed ring sizes
            elif (
                allowed_ring_sizes is not None
                and len(Chem.rdmolops.GetShortestPath(state, atom1, atom2))
                not in allowed_ring_sizes
            ):
                continue
            else:
                new_state.AddBond(atom1, atom2, bond_orders[valence])
            sanitization_result = Chem.SanitizeMol(new_state, catchErrors=True)
            # When sanitization fails
            if sanitization_result:
                continue
            bond_addition.add(Chem.MolToSmiles(new_state))
    return bond_addition


def _bond_removal(state):
    """Computes valid actions that involve removing bonds from the graph.

  Actions (where allowed):
    * triple->{double,single,None}
    * double->{single,None}
    * single->{None}

  Bonds are only removed (single->None) if the resulting graph has zero or one
  disconnected atom(s); the creation of multi-atom disconnected fragments is not
  allowed. Note that aromatic bonds are not modified.

  Args:
    state: RDKit Mol.

  Returns:
    Set of string SMILES; the available actions.
  """
    bond_orders = [
        None,
        Chem.BondType.SINGLE,
        Chem.BondType.DOUBLE,
        Chem.BondType.TRIPLE,
    ]
    bond_removal = set()
    for valence in [1, 2, 3]:
        for bond in state.GetBonds():
            # Get the bond from a copy of the molecule so that SetBondType() doesn't
            # modify the original state.
            bond = Chem.Mol(state).GetBondBetweenAtoms(
                bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            )
            if bond.GetBondType() not in bond_orders:
                continue  # Skip aromatic bonds.
            new_state = Chem.RWMol(state)
            # Kekulize the new state to avoid sanitization errors; note that bonds
            # that are aromatic in the original state are not modified (this is
            # enforced by getting the bond from the original state with
            # GetBondBetweenAtoms()).
            Chem.Kekulize(new_state, clearAromaticFlags=True)
            # Compute the new bond order as an offset from the current bond order.
            bond_order = bond_orders.index(bond.GetBondType())
            bond_order -= valence
            if bond_order > 0:  # Downgrade this bond.
                idx = bond.GetIdx()
                bond.SetBondType(bond_orders[bond_order])
                new_state.ReplaceBond(idx, bond)
                sanitization_result = Chem.SanitizeMol(new_state, catchErrors=True)
                # When sanitization fails
                if sanitization_result:
                    continue
                bond_removal.add(Chem.MolToSmiles(new_state))
            elif bond_order == 0:  # Remove this bond entirely.
                atom1 = bond.GetBeginAtom().GetIdx()
                atom2 = bond.GetEndAtom().GetIdx()
                new_state.RemoveBond(atom1, atom2)
                sanitization_result = Chem.SanitizeMol(new_state, catchErrors=True)
                # When sanitization fails
                if sanitization_result:
                    continue
                smiles = Chem.MolToSmiles(new_state)
                parts = sorted(smiles.split("."), key=len)
                # We define the valid bond removing action set as the actions
                # that remove an existing bond, generating only one independent
                # molecule, or a molecule and an atom.
                if len(parts) == 1 or len(parts[0]) == 1:
                    bond_removal.add(parts[-1])
    return bond_removal


if __name__ == '__main__':
    # s = 'CCOC'
    # mol = Chem.MolFromSmiles(s)
    # _ring_addition(mol)
    a = np.array([2, 5, 3, 10, 4, 14, 7, 9])
    print(bottleneck.argpartition(-a, 3)[:3])
    print(np.argsort(a))

    a = ['C=C(C(C)=C(CCCCCCCC)CCCCCCCC)C(CCC)=C(C)CCCCCCCCC', 'C=C(C(C)=C(CCCCCCCC)CCCCCCCC)C(CC)=C(C)CCCCCCCCCC', 'C=C(C(C)=C(CCCCCCCC)CCCCCCCCC)C(CC)=C(C)CCCCCCCCC']
    for s in a:
        print(plogp(Chem.MolFromSmiles(s)))
    b = ['C=C(C(C)=C(CCCCCCC)CCCCCCCCC)C(CC)=C(C)CCCCCCCCCC', 'C=C(C(C)=C(CCCCCCC)CCCCCCCCCC)C(CC)=C(C)CCCCCCCCC', 'C=C(C(C)=C(CCCCCCCC)CCCCCCCCC)C(CC)=C(C)CCCCCCCCC']
    for s in b:
        print(plogp(Chem.MolFromSmiles(s)))

    s = 'C=C(C(C)=C(CCCCCCC)CCCCCCCC)C(CC)=C(C)CCCCCCCCC'
    valid_actions = get_valid_actions(s)
    topk_actions = top_k(valid_actions, plogp)
    # scores = []
    # for a in valid_actions:
    #     scores.append(plogp(Chem.MolFromSmiles(a)))
    # print(sorted(scores)[-10:])
    for a in topk_actions:
        print(a, plogp(Chem.MolFromSmiles(a)))

