from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np
from scripts.util import copy_edit_mol, map_atoms
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def enumerate_fusions_between_templates(mol1, mol2):
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
        a1 = bond.GetBeginAtom().GetAtomMapNum()-1
        a2 = bond.GetEndAtom().GetAtomMapNum()-1
        bt = bond.GetBondType()
        new_mol_template.AddBond(a1, a2, bt)

    mesh = np.array(np.meshgrid(map1indices, map2indices, map2indices))
    combinations = mesh.T.reshape(-1, 3)

    products = []
    new_bonds = []
    indices = []
    for k in combinations:
        if k[1] == k[2]:
            continue
        mol = copy_edit_mol(new_mol_template)
        mol.AddBond(int(k[0])-1, int(k[1])-1, bt)
        mol.AddBond(int(k[0])-1, int(k[2])-1, bt)
        try:
            Chem.SanitizeMol(mol)
        except:
            continue
        indices.append(k)
        new_bonds.append(
            [
                mol.GetBondBetweenAtoms(int(k[0])-1, int(k[1])-1).GetIdx(),
                mol.GetBondBetweenAtoms(int(k[0])-1, int(k[2])-1).GetIdx(),
            ]
        )
        products.append(mol)

    mesh = np.array(np.meshgrid(map2indices, map1indices, map1indices))
    combinations = mesh.T.reshape(-1, 3)
    for k in combinations:
        if k[1] == k[2]:
            continue
        mol = copy_edit_mol(new_mol_template)
        mol.AddBond(int(k[0])-1, int(k[1])-1, bt)
        mol.AddBond(int(k[0])-1, int(k[2])-1, bt)
        try:
            Chem.SanitizeMol(mol)
        except:
            continue
        indices.append(k)
        new_bonds.append(
            [
                mol.GetBondBetweenAtoms(int(k[0])-1, int(k[1])-1).GetIdx(),
                mol.GetBondBetweenAtoms(int(k[0])-1, int(k[2])-1).GetIdx(),
            ]
        )
        products.append(mol)

    return new_mol1, new_mol2, products, new_bonds, indices


def generateReactionTemplates(mol1, mol2, products):
    sm1 = Chem.MolToSmarts(mol1)
    sm2 = Chem.MolToSmarts(mol2)
    reactionTemplates = []
    for k in products:
        pr = Chem.MolToSmarts(k)
        st = sm1 + "." + sm2 + ">>" + pr
        reactionTemplates.append(Chem.rdChemReactions.ReactionFromSmarts(st))
    return reactionTemplates


def generateProducts(mol1, mol2, reactionTemplates, indices):
    products, invalids, prodInd, invInd, new_bonds = [], [], [], [], []
    for k, i in zip(reactionTemplates, indices):
        prod = k.RunReactants((mol1, mol2))[0][0]
        try:
            Chem.SanitizeMol(prod)
            products.append(prod)
            for k in prod.GetAtoms():
                props = k.GetPropsAsDict()
                if "old_mapno" in props:
                    k.SetAtomMapNum(k.GetPropsAsDict()["old_mapno"])

            prodInd.append(i)
        except:
            invalids.append(prod)
            invInd.append(i)
            continue
    return products, invalids, prodInd, invInd, new_bonds


# def mapAtoms(set1_mols, set2_mols, mcs_mol1, mcs_mol2):
#     for m1 in set1_mols:
#         ss = m1.GetSubstructMatches(mcs_mol1)
#         for MCS_atom, ss_ind in zip(mcs_mol1.GetAtoms(), ss[0]):
#             at = m1.GetAtomWithIdx(int(ss_ind))
#             at.SetAtomMapNum(MCS_atom.GetAtomMapNum())

#     for m2 in set2_mols:
#         ss = m2.GetSubstructMatches(mcs_mol2)
#         for MCS_atom, ss_ind in zip(mcs_mol2.GetAtoms(), ss[0]):
#             at = m2.GetAtomWithIdx(int(ss_ind))
#             at.SetAtomMapNum(MCS_atom.GetAtomMapNum())


# def mapProductAtoms(prod, mcs_mol1, mcs_mol2):
#     ss = prod.GetSubstructMatches(mcs_mol1)
#     for MCS_atom, ss_ind in zip(mcs_mol1.GetAtoms(), ss[0]):
#         at = prod.GetAtomWithIdx(int(ss_ind))
#         at.SetAtomMapNum(MCS_atom.GetAtomMapNum())
#     ss = prod.GetSubstructMatches(mcs_mol2)
#     for MCS_atom, ss_ind in zip(mcs_mol2.GetAtoms(), ss[0]):
#         at = prod.GetAtomWithIdx(int(ss_ind))
#         at.SetAtomMapNum(MCS_atom.GetAtomMapNum())


def fusionCoupling(set1, set2):

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

    map_atoms([mol1, mol2])

    enum_products1, ats1 = fusionCouplingSingle(set1_mols, set2_mols, mol1, mol2)

    return enum_products1, ats1


def fusionCouplingSingle(set1_mols, set2_mols, mol1, mol2):

    mol1, mol2, products, new_bonds, indices = enumerate_fusions_between_templates(
        mol1, mol2
    )

    reactionTemplates = generateReactionTemplates(mol1, mol2, products)
    enum_products, invalids, enum_ind, inv_ind, new_bonds = generateProducts(
        set1_mols[0], set2_mols[0], reactionTemplates, indices
    )

    ats = []
    sms = []
    out = []
    for k, p in zip(enum_ind, enum_products):
        inds = []
        for atm in p.GetAtoms():
            props = atm.GetPropsAsDict()
            if "old_mapno" in props:
                if props["old_mapno"] in k:
                    inds.append(atm.GetIdx())
            if atm.GetIdx() in k and atm.GetIdx() == 0:
                inds.append(0)
        sm = Chem.MolToSmiles(p)
        if sm not in sms:
            out.append(p)
            sms.append(sm)
            ats.append(inds)

    return out, ats
