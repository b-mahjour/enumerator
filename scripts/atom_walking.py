import rdkit
from rdkit import Chem
import numpy as np
from scripts.util import copy_edit_mol
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def atom_walk(mol, atom):
    new_mol = Chem.RWMol(Chem.MolFromSmiles(""))
    added_atom = Chem.Atom(atom)
    cnt = 1
    og_atms = []
    for atom in mol.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetAtomMapNum(cnt)
        og_atms.append(cnt)
        cnt = cnt + 1
        new_mol.AddAtom(new_atom)

    new_atom = Chem.Atom(added_atom.GetSymbol())
    new_atom.SetFormalCharge(added_atom.GetFormalCharge())
    new_atom.SetAtomMapNum(cnt)

    new_mol.AddAtom(new_atom)

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)

    walked, new_bonds = [], []
    for ind in og_atms:
        dup = copy_edit_mol(new_mol, True)
        bt = rdkit.Chem.rdchem.BondType.SINGLE
        dup.AddBond(cnt-1, ind-1, bt)
        try:
            Chem.SanitizeMol(dup)
            new_bonds.append([dup.GetBondBetweenAtoms(cnt-1, ind-1).GetIdx()])
        except:
            continue
        walked.append(dup)

    return walked, new_bonds


def atom_swap(mol, atom):
    added_atom = Chem.Atom(atom)
    cnt = 1
    og_atms = []
    for atom in mol.GetAtoms():
        og_atms.append(cnt)
        cnt = cnt + 1

    walked, changed_atom = [], []
    changed_atom_saved = []
    for k in og_atms:
        dup = copy_edit_mol(mol)
        new_mol = Chem.RWMol(Chem.MolFromSmiles(""))
        cnt = 1
        for atom in dup.GetAtoms():
            if cnt == k:
                new_atom = Chem.Atom(added_atom.GetSymbol())
                new_atom.SetFormalCharge(added_atom.GetFormalCharge())
                changed_atom.append([k-1])
            else:
                new_atom = Chem.Atom(atom.GetSymbol())
                new_atom.SetFormalCharge(atom.GetFormalCharge())
            new_atom.SetAtomMapNum(cnt)
            new_mol.AddAtom(new_atom)
            cnt = cnt + 1

        for bond in dup.GetBonds():
            a1 = bond.GetBeginAtom().GetIdx()
            a2 = bond.GetEndAtom().GetIdx()
            bt = bond.GetBondType()
            new_mol.AddBond(a1, a2, bt)

        try:
            Chem.SanitizeMol(new_mol)
        except:
            continue
        walked.append(new_mol)
        changed_atom_saved.append(changed_atom[-1])
    return walked, changed_atom_saved


def ch_coupling(mol1, mol2):
    new_mol = Chem.RWMol(Chem.MolFromSmiles(""))
    cnt = 0
    c1s = []
    for atom in mol1.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetAtomMapNum(cnt)

        if atom.GetSymbol() == "C" and atom.GetNumImplicitHs() > 0:
            c1s.append(cnt)
        cnt = cnt + 1
        new_mol.AddAtom(new_atom)

    c2s = []
    cnt2 = 0
    mol2connect = {}
    for atom in mol2.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        mol2connect[cnt2] = cnt
        new_atom.SetAtomMapNum(cnt)
        if atom.GetSymbol() == "C" and atom.GetNumImplicitHs() > 0:
            c2s.append(cnt)
        cnt = cnt + 1
        cnt2 = cnt2 + 1
        new_mol.AddAtom(new_atom)

    for bond in mol1.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)

    for bond in mol2.GetBonds():
        a1 = mol2connect[bond.GetBeginAtom().GetIdx()]
        a2 = mol2connect[bond.GetEndAtom().GetIdx()]
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)

    mesh = np.array(np.meshgrid(c1s, c2s))
    combinations = mesh.T.reshape(-1, 2)

    new_mols = []
    for k in combinations:
        edit = copy_edit_mol(new_mol)
        bt = rdkit.Chem.rdchem.BondType.SINGLE
        edit.AddBond(int(k[0]), int(k[1]), bt)
        try:
            Chem.SanitizeMol(edit)
        except:
            continue
        new_mols.append(edit)

    return new_mols, []


bond_order_map = {
    "single": rdkit.Chem.rdchem.BondType.SINGLE,
    "double": rdkit.Chem.rdchem.BondType.DOUBLE,
    "triple": rdkit.Chem.rdchem.BondType.TRIPLE,
}

def bond_order_enumerator(mol, initial_bond_order="single", target_bond_order="double"):

    initial_bond = bond_order_map[initial_bond_order]
    target_bond = bond_order_map[target_bond_order]
    new_mol = Chem.RWMol(Chem.MolFromSmiles(""))
    cnt = 1
    bonds_to_change = []
    for atom in mol.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetAtomMapNum(cnt)
        cnt = cnt + 1
        new_mol.AddAtom(new_atom)

    bnd_count = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)
        if bt == initial_bond:
            bonds_to_change.append(bnd_count)
        bnd_count += 1

    new_mols = []
    changed_bond = []
    for b in bonds_to_change:
        edit = copy_edit_mol(new_mol, True)
        bond = edit.GetBondWithIdx(b)
        if bond.GetIdx() == b:
            edit.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            edit.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), target_bond)
        try:
            Chem.SanitizeMol(edit)
        except:
            continue
        changed_bond.append([bond.GetIdx()])
        new_mols.append(edit)

    return new_mols, changed_bond