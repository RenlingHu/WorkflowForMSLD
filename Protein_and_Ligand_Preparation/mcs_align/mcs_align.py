import copy
import rdkit
import warnings
import numpy as np
from rdkit import Chem

from . import rmsd
from . import zmatrix

__version__ = "0.1"
__license__ = "WTFPL Version 2"

class MCSALIGN(object):
    r"""
    Align small molecules according to their maximum common substructure (MCS).

    Parameters
    ----------
    mol_list : list(rdkit.Chem.rdchem.Mol)
        The molecules to align.
    """

    def __init__(self, mol_list):
        r"""
        Initialize some internal data and try find the MCS.
        """
        self._mol = copy.deepcopy(mol_list)
        self._mol_aligned = None

        ps = Chem.rdFMCS.MCSParameters()
        results = Chem.rdFMCS.FindMCS(self._mol, ps)
        self._mcs = Chem.MolFromSmarts(results.smartsString)
        self._mcs_indices = [x.GetSubstructMatch(self._mcs) for x in self._mol]

        if (self._mcs.GetNumAtoms() == 0):
            raise RuntimeError('Can not find MCS for given molecules!')

    def _substructure_kabsch(self, mol_1, mol_2, indices_1, \
            indices_2) -> rdkit.Chem.rdchem.Mol:
        r"""
        The kabsch method based on partial atoms of molecules.

        Parameters
        ----------
        mol_1 : rdkit.Chem.rdchem.Mol
            The reference molecule.
        mol_2 : rdkit.Chem.rdchem.Mol
            The aligned molecule.
        indices_1 : list(int)
            The indices of the atoms that belong to the substructure of mol_1.
        indices_2 : list(int)
            The indices of the atoms that belong to the substructure of mol_2.

        Notes
        -----
        The first conformation of both molecules will be used. And the
        alignment will NOT happend in-placed.
        """
        if (len(indices_1) != len(indices_2)):
            raise RuntimeError('Can not align two substructures with ' + \
                'different atom numbers!')

        P = mol_1.GetConformer().GetPositions()
        Q = mol_2.GetConformer().GetPositions()
        P_c = rmsd.centroid(P[indices_1])
        Q_c = rmsd.centroid(Q[indices_2])
        U = rmsd.kabsch((Q[indices_2] - Q_c), (P[indices_1] - P_c))
        Q_1 = np.dot((Q - Q_c), U) + P_c

        result = copy.deepcopy(mol_2)
        conf = result.GetConformer()
        for i in range(0, conf.GetNumAtoms()):
            conf.SetAtomPosition(i, Q_1[i])

        return result

    def _mol_to_zmatrix(self, mol_id, ref_indices) -> (list, list):
        r"""
        Convert a molecule to a zmatrix using three given atoms as the
        references.

        Parameters
        ----------
        mol_id : int
            The index of converted molecule.
        ref_indices : list(int)
            The indices of the reference atoms.
        """
        atom_map = []
        cartesian = []
        mol = self._mol[mol_id]

        P = mol.GetConformer().GetPositions()
        for i in ref_indices:
            atom = mol.GetAtomWithIdx(i)
            cartesian.append([atom.GetSymbol(), P[i], atom.GetMass()])
            atom_map.append(i)

        count = 3
        while (count != mol.GetNumAtoms()):
            for i in range(0, mol.GetNumAtoms()):
                if (i not in atom_map):
                    atom = mol.GetAtomWithIdx(i)
                    for n in atom.GetNeighbors():
                        if (n.GetIdx() in atom_map):
                            line = [atom.GetSymbol(), P[i], atom.GetMass()]
                            cartesian.append(line)
                            atom_map.append(i)
                            count += 1
                            break

        c = zmatrix.Converter()
        c.cartesian = cartesian
        c.add_first_three_to_zmatrix()
        for i, atom in enumerate(c.cartesian[3:], start=3):
            neighbors = []
            if (atom_map[i] in self.mcs_indices[mol_id]):
                neighbors = [atom_map.index(x) for x in ref_indices]
            else:
                for n in mol.GetAtomWithIdx(atom_map[i]).GetNeighbors():
                    if (n.GetIdx() in atom_map[:i]):
                        neighbors.append(atom_map[:i].index(n.GetIdx()))
                        break
                for n in \
                    mol.GetAtomWithIdx(atom_map[neighbors[0]]).GetNeighbors():
                    if (n.GetIdx() in atom_map[:i]):
                        neighbors.append(atom_map[:i].index(n.GetIdx()))
                        break
                for n in \
                    mol.GetAtomWithIdx(atom_map[neighbors[1]]).GetNeighbors():
                    if ((n.GetIdx() in atom_map[:i]) and \
                        (n.GetIdx() != atom_map[neighbors[0]])):
                        neighbors.append(atom_map[:i].index(n.GetIdx()))
                        break
                if (len(neighbors) == 2):
                    for j in range(0, i):
                        if (j not in neighbors):
                            neighbors.append(j)
                            break
            c.add_atom_to_zmatrix(i, atom, neighbors)

        for i in range(2, mol.GetNumAtoms()):
            c.zmatrix[i][1][1][1] = np.radians(c.zmatrix[i][1][1][1])
        for i in range(3, mol.GetNumAtoms()):
            c.zmatrix[i][1][2][1] = np.radians(c.zmatrix[i][1][2][1])

        atom_map = np.array(atom_map).argsort().tolist()

        return (c.zmatrix, atom_map)

    def _align_mol_pair(self, mol_id_1, mol_id_2, strict_align=True):
        r"""
        Align a pair of input molecules.

        Parameters
        ----------
        mol_id_1 : int
            The index of the reference molecule.
        mol_id_2 : int
            The index of the aligned molecule.
        strict_align : bool
            If use the zmatrix method to perform strict structure
            transformation.
        """
        result = self._substructure_kabsch(self._mol[mol_id_1], \
            self._mol[mol_id_2], list(self.mcs_indices[mol_id_1]), \
            list(self.mcs_indices[mol_id_2]))

        if (strict_align):
            zmat_1, atom_map_1 = self._mol_to_zmatrix(mol_id_1, \
                self.mcs_indices[mol_id_1][0:3])
            zmat_2, atom_map_2 = self._mol_to_zmatrix(mol_id_2, \
                self.mcs_indices[mol_id_2][0:3])
            for i in range(0, self.mcs.GetNumAtoms()):
                index_1 = atom_map_1[self.mcs_indices[mol_id_1][i]]
                index_2 = atom_map_2[self.mcs_indices[mol_id_2][i]]
                for j in range(0, 3):
                    try:
                        zmat_2[index_2][1][j][1] = zmat_1[index_1][1][j][1]
                    except:
                        pass

            c = zmatrix.Converter()
            c.zmatrix = zmat_2
            c.zmatrix_to_cartesian()

            conf = result.GetConformer()
            for i in range(0, conf.GetNumAtoms()):
                conf.SetAtomPosition(i, c.cartesian[atom_map_2[i]][1])

            result = self._substructure_kabsch(self._mol[mol_id_1], result, \
                list(self.mcs_indices[mol_id_1]), \
                list(self.mcs_indices[mol_id_2]))

        return result

    def align(self, ref_mol_index=0, strict_align=True):
        r"""
        Align input molecules to a given molecule.

        Parameters
        ----------
        ref_mol_index : int
            The index of the reference molecule.
        strict_align : bool
            If use the zmatrix method to perform strict structure
            transformation.
        """
        if (self.mcs.GetNumAtoms() < 3):
            warnings.warn('Too few atoms in the MCS! ' + \
                'Strict alignment will be disabled!')
            strict_align = False

        self._mol_aligned = []
        for i in range(len(self.molecules_original)):
            if (i != ref_mol_index):
                m = self._align_mol_pair(ref_mol_index, i, strict_align)
            else:
                m = copy.deepcopy(self._mol[ref_mol_index])
            self._mol_aligned.append(m)
        return self._mol_aligned

    @property
    def mcs(self) -> rdkit.Chem.rdchem.Mol:
        r"""
        The maximum common substructure of the given molecules.
        """
        return self._mcs

    @property
    def mcs_indices(self) -> list:
        r"""
        Indices of the maximum common substructure of the given molecules.
        """
        return self._mcs_indices

    @property
    def molecules_original(self) -> list:
        r"""
        The original molecules.
        """
        return self._mol

    @property
    def molecules_aligned(self) -> list:
        r"""
        The aligned molecules.
        """
        return self._mol_aligned
