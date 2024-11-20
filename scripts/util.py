from rdkit import Chem
from rdkit.Chem import Draw

def render_enumeration(products, hit_bonds=[], hit_atoms=[], numMolsPerRow=5, titles=[]):

    if len(titles) == 0:
        titles = ['' for _ in products]
    else:
        if len(titles) != len(products):
            raise ValueError('Number of titles must match number of products')

    
    img = Draw.MolsToGridImage(
        products,
        molsPerRow=numMolsPerRow,
        subImgSize=(200, 200),
        legends=titles,
        highlightAtomLists=hit_atoms,
        highlightBondLists=hit_bonds,
    )

    return img

def copy_edit_mol(mol, san=False):
    new_mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for atom in mol.GetAtoms():
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetAtomMapNum(atom.GetAtomMapNum())
        if atom.GetIsAromatic() and atom.GetSymbol() == 'N':
            if san:
                Chem.SanitizeMol(mol)
            new_atom.SetNumExplicitHs(atom.GetTotalNumHs())

        new_mol.AddAtom(new_atom)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)
    return new_mol 


def map_atoms(reactant_objs):
    n = 1
    for robj in reactant_objs:
        for a in robj.GetAtoms():
            a.SetAtomMapNum(n)
            n += 1


def get_reactant_map(reactant_objs):
    atom_map_to_reactant_atom = {}
    for robj in reactant_objs:
        Chem.SanitizeMol(robj, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
        for a in robj.GetAtoms():
            atom_map_to_reactant_atom[a.GetAtomMapNum()] = a

    return atom_map_to_reactant_atom



def adjust_atom_maps(
    i, rxn, rxt_in
):
    for atom in i.GetAtoms():
        if (
            "molAtomMapNumber" not in atom.GetPropsAsDict()
            and "old_mapno" in atom.GetPropsAsDict()
        ):
            old_atom_idx = atom.GetPropsAsDict()["old_mapno"]
            react_atom_idx = atom.GetPropsAsDict()["react_atom_idx"]
            hit = 0
            for i2, rx in enumerate(rxn.GetReactants()):
                found = False
                for i3, j in enumerate(rx.GetAtoms()):
                    if j.GetAtomMapNum() == old_atom_idx:
                        found = True
                        break
                if found:
                    hit = i2
                    break
            atom_map = rxt_in[hit].GetAtomWithIdx(react_atom_idx).GetAtomMapNum()
            atom.SetAtomMapNum(atom_map)

            # atom_map_to_product_atom[atom_map] = atom


def get_reacting_atoms(i, atom_map_to_reactant_atom):
    rxt_atoms = []
    for atom in i.GetAtoms():
        old_atom = atom_map_to_reactant_atom[atom.GetAtomMapNum()]
        old_atom_info = [
            old_atom.GetFormalCharge(),
            old_atom.GetNumExplicitHs(),
            old_atom.GetNumImplicitHs(),
            old_atom.GetDegree(),
            old_atom.GetTotalDegree(),
            old_atom.GetTotalValence(),
            old_atom.GetTotalNumHs(),
            old_atom.GetExplicitValence(),
            old_atom.GetImplicitValence(),
        ]
        new_atom_info = [
            atom.GetFormalCharge(),
            atom.GetNumExplicitHs(),
            atom.GetNumImplicitHs(),
            atom.GetDegree(),
            atom.GetTotalDegree(),
            atom.GetTotalValence(),
            atom.GetTotalNumHs(),
            atom.GetExplicitValence(),
            atom.GetImplicitValence(),
        ]
        # print(old_atom.GetAtomMapNum(), old_atom.GetSmarts(), len(old_atom.GetNeighbors()))
        # print(atom.GetAtomMapNum(), old_atom.GetSmarts(), len(old_atom.GetNeighbors()))
        # print()
        if (
            new_atom_info[0] == old_atom_info[0]
            and new_atom_info[4] != old_atom_info[4]
            and new_atom_info[5] != old_atom_info[5]
            and new_atom_info[6] != old_atom_info[6]
        ):
            old_bonds = [x.GetBondType() for x in old_atom.GetBonds()]
            new_bonds = [x.GetBondType() for x in atom.GetBonds()]

            if old_bonds == new_bonds:
                # basically in situations where ([N;H1,H2;+0:1]-[c;H0;+0:2].[O;H0;+0:3]=[C;!$([C](=[O])[O]);H0;+0:4])>>([N;+1;H1,H2:1](-[c;H0;+0:2])-[C;!$([C](=[O])[O]);H0;+0:4]-[O;-1;H0:3])
                # the intramolecular reaction forms a bond that's already there.
                # rdkit will run the reaction anyways, and mess up the protonation state,
                # creating a molecule with improper valence
                return None

        if (
            new_atom_info[3] == old_atom_info[3]
            and new_atom_info[4] == old_atom_info[4]
            and new_atom_info[5] == old_atom_info[5]
            and new_atom_info[6] == old_atom_info[6]
            and new_atom_info[0] != old_atom_info[0]
        ):
            old_bonds = [x.GetBondType() for x in old_atom.GetBonds()]
            new_bonds = [x.GetBondType() for x in atom.GetBonds()]
            if old_bonds == new_bonds:
                # basically in situations where ([N;H1,H2;+0:1]-[c;H0;+0:2].[O;H0;+0:3]=[C;!$([C](=[O])[O]);H0;+0:4])>>([N;+1;H1,H2:1](-[c;H0;+0:2])-[C;!$([C](=[O])[O]);H0;+0:4]-[O;-1;H0:3])
                # the intramolecular reaction forms a bond that's already there.
                # rdkit will run the reaction anyways, and mess up the protonation state,
                # creating a molecule with improper valence
                return None

        if old_atom_info != new_atom_info:
            rxt_atoms.append(atom.GetAtomMapNum())

    return rxt_atoms
