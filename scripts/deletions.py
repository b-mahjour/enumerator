from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def deleteAtom(mol):
    new_mol = Chem.RWMol(Chem.MolFromSmiles(""))
    cnt = 0
    bond_counter = {}
    for atom in mol.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetAtomMapNum(cnt)
        if atom.GetIsAromatic() and atom.GetSymbol() == "N":
            new_atom.SetNumExplicitHs(atom.GetTotalNumHs())
        bond_counter[cnt] = []
        cnt = cnt + 1
        new_mol.AddAtom(new_atom)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bond_counter[a1].append(a2)
        bond_counter[a2].append(a1)
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)

    Chem.Kekulize(new_mol, clearAromaticFlags=True)

    new_molecules = []
    removed = []
    for k in bond_counter:
        edit_mol = Chem.RWMol(Chem.MolFromSmiles(""))

        if len(bond_counter[k]) != 2:
            continue

        for atom in new_mol.GetAtoms():
            new_atom = Chem.Atom(atom.GetSymbol())
            new_atom.SetFormalCharge(atom.GetFormalCharge())
            new_atom.SetAtomMapNum(atom.GetAtomMapNum()+1)
            edit_mol.AddAtom(new_atom)

        for bond in new_mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomMapNum()-1 == int(k) or a2.GetAtomMapNum()-1 == int(k):
                continue
            bt = bond.GetBondType()
            edit_mol.AddBond(a1.GetIdx(), a2.GetIdx(), bt)

        edit_mol.AddBond(bond_counter[k][0], bond_counter[k][1], bt)

        edit_mol.RemoveAtom(int(k))

        try:
            Chem.SanitizeMol(edit_mol)
        except:
            continue
        new_molecules.append(edit_mol)
        removed.append([list(edit_mol.GetBonds())[-1].GetIdx()])

    return new_molecules, new_mol, bond_counter, removed
