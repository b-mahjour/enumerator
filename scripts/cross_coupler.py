import rdkit
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem
import numpy as np
from scripts.util import (
    copy_edit_mol,
    get_reactant_map,
    map_atoms,
    get_reacting_atoms,
    adjust_atom_maps,
)
import copy
import itertools


def enumerate_bonds_between_templates(mol1, mol2):
    map1indices = []
    new_mol_template = Chem.RWMol(Chem.MolFromSmiles(""))
    new_mol1 = Chem.RWMol(Chem.MolFromSmiles(""))
    map1 = 1
    for atom in mol1.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetAtomMapNum(map1)
        if atom.GetIsAromatic() and atom.GetSymbol() == "N":
            new_atom.SetNumExplicitHs(atom.GetTotalNumHs())

        new_mol1.AddAtom(new_atom)
        new_mol_template.AddAtom(new_atom)
        map1indices.append(map1)
        map1 = map1 + 1
    map2indices = []
    new_mol2 = Chem.RWMol(Chem.MolFromSmiles(""))
    map2 = map1
    for atom in mol2.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetAtomMapNum(map2)
        if atom.GetIsAromatic() and atom.GetSymbol() == "N":
            new_atom.SetNumExplicitHs(atom.GetTotalNumHs())

        new_mol2.AddAtom(new_atom)
        new_mol_template.AddAtom(new_atom)
        map2indices.append(map2)
        map2 = map2 + 1

    for bond in mol1.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol1.AddBond(a1, a2, bt)
        new_mol_template.AddBond(a1, a2, bt)
    for bond in mol2.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol2.AddBond(a1, a2, bt)
    for bond in new_mol2.GetBonds():
        a1 = bond.GetBeginAtom().GetAtomMapNum() - 1
        a2 = bond.GetEndAtom().GetAtomMapNum() - 1
        bt = bond.GetBondType()
        new_mol_template.AddBond(a1, a2, bt)

    mesh = np.array(np.meshgrid(map1indices, map2indices))
    combinations = mesh.T.reshape(-1, 2)
    products = []
    new_bonds = []
    for k in combinations:
        mol = copy_edit_mol(new_mol_template)
        bt = rdkit.Chem.rdchem.BondType.SINGLE
        mol.AddBond(int(k[0]) - 1, int(k[1]) - 1, bt)
        new_bonds.append([mol.GetBondBetweenAtoms(int(k[0]) - 1, int(k[1]) - 1).GetIdx()])
        products.append(mol)
    return new_mol1, new_mol2, products, new_bonds


def combinate(set1, set2):
    set1_mols = [Chem.MolFromSmiles(m) for m in set1]
    set2_mols = [Chem.MolFromSmiles(m) for m in set2]

    if len(set1_mols) > 1:
        mcs1 = rdFMCS.FindMCS(set1_mols)
        mol1 = Chem.MolFromSmarts(mcs1.smartsString)
    else:
        mcs1 = Chem.MolToSmarts(set1_mols[0])
        mol1 = set1_mols[0]
    if len(set2_mols) > 1:
        mcs2 = rdFMCS.FindMCS(set2_mols)
        mol2 = Chem.MolFromSmarts(mcs2.smartsString)
    else:
        mcs2 = Chem.MolToSmarts(set2_mols[0])
        mol2 = set2_mols[0]

    mol1, mol2, products, hit_bonds = enumerate_bonds_between_templates(mol1, mol2)

    mol1sm = Chem.MolToSmarts(mol1)
    mol2sm = Chem.MolToSmarts(mol2)
    rxns = []
    for k in products:
        prod_sm = Chem.MolToSmarts(k)
        rxns.append(AllChem.ReactionFromSmarts(mol1sm + "." + mol2sm + ">>" + prod_sm))

    all_reaction_inputs = list(itertools.product(set1_mols, set2_mols))
    prods = []
    bond_highlight = []
    for ko in all_reaction_inputs:
        map_atoms(ko)
        rmap = get_reactant_map(ko)
        for r in rxns:
            k = copy.deepcopy(ko)
            p = r.RunReactants((k[0], k[1]))
            for pp in p:
                for ppp in pp:
                    adjust_atom_maps(ppp, r, k)
                    Chem.SanitizeMol(ppp, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
                    ratoms = get_reacting_atoms(ppp, rmap)
                    prods.append(ppp)
                    atom1 = 0
                    atom2 = 0
                    for atom in ppp.GetAtoms():
                        if atom.GetAtomMapNum() == ratoms[0]:
                            atom1 = atom.GetIdx()
                        if atom.GetAtomMapNum() == ratoms[1]:
                            atom2 = atom.GetIdx()
                    bond_idx = ppp.GetBondBetweenAtoms(atom1, atom2).GetIdx()
                    bond_highlight.append([bond_idx])

    return mol1, mol2, prods, bond_highlight
